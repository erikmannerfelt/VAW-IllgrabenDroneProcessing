#!/usr/bin/env python

import json
import os
from collections import namedtuple

import Metashape as ms

from . import main, pdal

# What command to run pdal as
PDAL_PATH = "pdal"


def get_stable_ground_locations(chunk: ms.Chunk):
    group_label = "stable_ground"
    coords = []
    for shape in chunk.shapes:
        if not shape.group.label == group_label:
            continue
        if not shape.type == ms.Shape.Type.Point:
            continue
        # Extract the first (and only) vertex of the point
        x_coord, y_coord, z_coord = shape.vertices[0]
        coords.append([x_coord, y_coord, z_coord])

    assert len(coords) > 0, "No points in the group: {}".format(group_label)

    return coords


def get_stable_ground_locations():
    with open("input/stable_ground_points.xyz") as infile:
        lines = infile.read().splitlines()

    header = lines.pop(0)
    assert header == "X,Y,Z", "stable_ground_points.xyz has unexpected first row. Expected 'X,Y,Z', got {}".format(
        header)

    coords = []
    for line in lines:
        coords.append([float(string) for string in line.split(",")])

    assert len(coords) > 0, "No points read from stable_ground_points.xyz!"
    return coords


Bounds = namedtuple("Bounds", ["x_min", "x_max", "y_min", "y_max", "z_min", "z_max"])


def get_bounding_boxes(coords, radius=15):
    bounds = []
    for x_coord, y_coord, z_coord in coords:
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


def extract_features(chunk, bounds):
    """
    Extract subsets of a chunk's dense point cloud using a given list of bounding boxes.

    param: chunk: The input chunk.
    type: chunk. Metashape.Chunk
    param: bounds: A list of bounds
    type: bounds: List[Bounds]
    """
    # Make the directory in which to save the extracted features
    features_dir = os.path.join(main.TEMP_FOLDER, chunk.label, "features")
    if not os.path.isdir(features_dir):
        os.makedirs(features_dir)

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
    ]'''\
        .replace("INPUT_FILENAME", os.path.join(main.TEMP_FOLDER, chunk.label, "dense_cloud_lowres.ply"))\
        .replace("OUTPUT_FILENAME_TEMPLATE", os.path.join(features_dir, "feature_#.las"))

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
    # Add the bounds_string to the pipeline
    extraction_pipeline = extraction_pipeline.replace("BOUNDS", bounds_string[:bounds_string.rindex(",")])

    # Run with point streaming turned off
    # Why? Because otherwise it doesn't work, that's why!
    pdal.run_pipeline(extraction_pipeline, stream=False)


def compare_features(reference_chunk, aligned_chunk, bounds):

    reference_feature_dir = os.path.join(main.TEMP_FOLDER, reference_chunk.label, "features")
    aligned_feature_dir = os.path.join(main.TEMP_FOLDER, aligned_chunk.label, "features")
    # assert os.listdir(reference_feature_dir) == os.listdir(aligned_feature_dir),\
    #    "{} and {} have different feature counts".format(reference_feature_dir, aligned_feature_dir)

    features = []
    for filename in os.listdir(reference_feature_dir):
        if filename.endswith(".las"):
            features.append(filename)
    Point = namedtuple("Point", ["start", "destination"])

    points = []
    print("Running feature-wise ICP for chunk: {}".format(aligned_chunk.label))
    for i, feature in enumerate(features):
        icp_file = os.path.join(aligned_feature_dir, feature.replace(".las", "_icp_meta.json"))
        print(icp_file)
        if os.path.isfile(icp_file):
            with open(icp_file) as infile:
                print("Loading cached ICP")
                transform = json.load(infile)["stages"]["filters.icp"]["composed"].replace("\n", " ")
        else:
            transform = pdal.icp_coregistration(
                os.path.join(reference_feature_dir, feature),
                os.path.join(aligned_feature_dir, feature),
                composed=True)

        # bound = bounds[i]
        # starting_point = [(bound.x_max + bound.x_min) / 2, (bound.y_max +
        #                                                    bound.y_min) / 2, (bound.z_max + bound.z_min) / 2]
        # transformed_point = pdal.transform_point(starting_point, transform, invert=True)
        # t_float = [float(string) for string in transform.split(" ") if string]
        # transformed_point = [starting_point[0] + t_float[3],
        #                     starting_point[1] + t_float[7],
        #                     starting_point[2] + t_float[11]]

        starting_point_0 = pdal.get_first_point_info(os.path.join(aligned_feature_dir, feature))
        destination_point_0 = pdal.get_first_point_info(os.path.join(
            aligned_feature_dir, feature.replace(".las", "_ICP_aligned.las")))

        points.append(Point([starting_point_0["X"], starting_point_0["Y"], starting_point_0["Z"]],
                            [destination_point_0["X"], destination_point_0["Y"], destination_point_0["Z"]]))

    return points


def add_correction_marker(chunk, initial_position, corrected_position, name):
    local_initial_position = chunk.transform.matrix.inv().mulp(chunk.crs.unproject(initial_position))

    for previous_marker in chunk.markers:
        if previous_marker.label == name:
            chunk.remove(previous_marker)
    marker = chunk.addMarker(local_initial_position)
    marker.label = name

    marker.reference.location = corrected_position
    marker.reference.enabled = True


def align_stable_ground_locations(reference_chunk, aligned_chunk):
    bounds = get_bounding_boxes(get_stable_ground_locations())

    if not os.path.isdir(os.path.join(main.TEMP_FOLDER, reference_chunk.label, "features")):
        print("Extracting reference chunk features")
        extract_features(reference_chunk, bounds)

    # TODO: Dangerous assumption that they exist
    if not os.path.isdir(os.path.join(main.TEMP_FOLDER, aligned_chunk.label, "features")):
        print("Extracting features from chunk: {}".format(aligned_chunk.label))
        extract_features(aligned_chunk, bounds)

    points = compare_features(reference_chunk, aligned_chunk, bounds)

    for i, point in enumerate(points):
        add_correction_marker(aligned_chunk, point.start, point.destination, "auto_ICP_{}".format(str(i).zfill(3)))
