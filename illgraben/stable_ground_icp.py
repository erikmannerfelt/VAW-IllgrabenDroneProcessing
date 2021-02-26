"""Stable ground ICP coalignment analysis."""

from __future__ import annotations

import os
from collections import namedtuple

import Metashape as ms

from illgraben import main, processing_tools
from illgraben.files import INPUT_FILEPATHS, PROCESSING_FOLDER

# What command to run pdal as. TODO: Remove?
PDAL_PATH = "pdal"


def get_stable_ground_locations() -> list[tuple[float, float, float]]:
    """
    Read the hand-picked stable ground locations file.

    return: coords: A list of X/Y/Z coordinates for the stable ground locations.
    """

    with open(INPUT_FILEPATHS["stable_ground_points"]) as infile:
        lines = infile.read().splitlines()

    header = lines.pop(0)
    assert header == "X,Y,Z", "stable_ground_points.xyz has unexpected first row. Expected 'X,Y,Z', got {}".format(
        header)

    stable_points = []
    for line in lines:
        stable_points.append(tuple([float(string) for string in line.split(",")]))

    assert len(stable_points) > 0, "No points read from stable_ground_points.xyz!"
    return stable_points  # type: ignore


# Create a helper namedtuple for bounds checking
Bounds = namedtuple("Bounds", ["x_min", "x_max", "y_min", "y_max", "z_min", "z_max"])


def get_bounding_boxes(points: list[tuple[float, float, float]], radius: float = 15) -> list[Bounds]:
    """
    Convert point coordinates to bounding boxes for point extraction.

    param: points: The points to create bounding boxes around.
    param: radius: The radii of the output bounding boxes.

    return: bounds: A list of bounding boxes
    """
    bounds: list[Bounds] = []
    for x_coord, y_coord, z_coord in points:
        bound = Bounds(
            x_min=x_coord - radius,
            x_max=x_coord + radius,
            y_min=y_coord - radius,
            y_max=y_coord + radius,
            z_min=z_coord - radius,
            z_max=z_coord + radius
        )
        bounds.append(bound)

    return bounds


def extract_features(chunk: ms.Chunk, bounds: list[Bounds]):
    """
    Extract subsets of a chunk's dense point cloud using a given list of bounding boxes.

    param: chunk: The input chunk.
    type: chunk. Metashape.Chunk
    param: bounds: A list of bounds
    type: bounds: List[Bounds]
    """
    # Make the directory in which to save the extracted features
    features_dir = os.path.join(main.PROCESSING_FOLDER, chunk.label, "features")
    os.makedirs(features_dir, exist_ok=True)

    # Pipeline to provide PDAL with
    extraction_pipeline = '''
    [
        "INPUT_FILENAME",
        {
            "type": "filters.crop",
            "bounds": [
                        BOUNDS
                      ]
        },
        "OUTPUT_FILENAME_TEMPLATE"
    ]'''

    # Loop through each bounds and format them to PDAL's standard
    bounds_string = ''
    for bound in bounds:
        bounds_string += ' "([{0}, {1}], [{2}, {3}], [{4}, {5}])",\n\t\t\t'.format(
            bound.x_min,
            bound.x_max,
            bound.y_min,
            bound.y_max,
            bound.z_min,
            bound.z_max)

    # Keep everything until the last comma
    bounds_string = bounds_string[:bounds_string.rindex(",")]

    # Run with point streaming turned off
    # Why? Because otherwise it doesn't work, that's why!
    processing_tools.run_pdal_pipeline(
        pipeline=extraction_pipeline,
        stream=False,
        parameters={
            "INPUT_FILENAME": os.path.join(PROCESSING_FOLDER, chunk.label, "dense_cloud_for_ICP.ply"),
            "OUTPUT_FILENAME_TEMPLATE": os.path.join(features_dir, "feature_#.las"),
            "BOUNDS": bounds_string
        }
    )


def compare_features(reference_chunk: ms.Chunk, aligned_chunk: ms.Chunk, bounds: list[Bounds]):
    """
    Coregister features from two chunks. 

    param: reference_chunk: The chunk to act as reference.
    param: aligned_chunk: The chunk to be aligned.
    param: bounds: Bounding boxes of the features.

    """

    # Make appropriate folder names for the features.
    reference_feature_dir = os.path.join(main.PROCESSING_FOLDER, reference_chunk.label, "features")
    aligned_feature_dir = os.path.join(main.PROCESSING_FOLDER, aligned_chunk.label, "features")

    # Loop through all the features that exist in the reference folder
    features = [feature for feature in os.listdir(reference_feature_dir) if feature.endswith(".las")]

    # Make a helper namedtuple for point start and destination coordinates
    Point = namedtuple("Point", ["start", "destination"])

    points: list[Point] = []
    print("Running feature-wise ICP for chunk: {}".format(aligned_chunk.label))
    for feature in features:
        icp_output_meta_file = os.path.join(aligned_feature_dir, feature.replace(".las", "_icp_meta.json"))
        # Use a preexisting transform if it exists
        if os.path.isfile(icp_output_meta_file):
            print("Using cached ICP")
        else:
            processing_tools.run_icp(
                os.path.join(reference_feature_dir, feature),
                os.path.join(aligned_feature_dir, feature),
            )

        # Use the first point in the reference vs. aligned point cloud as a source-destination pair
        starting_point_0 = processing_tools.get_first_point_info(os.path.join(aligned_feature_dir, feature))
        destination_point_0 = processing_tools.get_first_point_info(os.path.join(
            aligned_feature_dir, feature.replace(".las", "_ICP_aligned.las")))

        # Append the pair to the points list
        points.append(
            Point(
                start=[starting_point_0["X"], starting_point_0["Y"], starting_point_0["Z"]],
                destination=[destination_point_0["X"], destination_point_0["Y"], destination_point_0["Z"]]
            )
        )

    return points
