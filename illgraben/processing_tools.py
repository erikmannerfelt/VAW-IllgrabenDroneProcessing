"""Processing functions using tools unrelated to Metashape. Mostly PDAL wrappers."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
import warnings
from typing import Any, Optional

from illgraben.constants import CONSTANTS


def run_pdal_pipeline(pipeline: str, output_metadata_file: Optional[str] = None,
                      parameters: Optional[dict[str, str]] = None, stream: bool = True) -> dict[str, Any]:
    """
    Run a PDAL pipeline.

    param: pipeline: The pipeline to run.
    param: output_metadata_file: Optional. The filepath for the pipeline metadata.
    param: parameters: Optional. Parameters to fill the pipeline with, e.g. {"FILEPATH": "/path/to/file"}.
    param: stream: Whether to use stream mode. Will not work with certain tools.

    return: output_meta: The metadata produced by the output.
    """
    # Create a temporary directory to save the output metadata in
    temp_dir = tempfile.TemporaryDirectory()
    # Fill the pipeline with
    if parameters is not None:
        for key in parameters:
            # Warn if the key cannot be found in the pipeline
            if key not in pipeline:
                warnings.warn(
                    f"{key}:{parameters[key]} given to the PDAL pipeline but the key was not found", RuntimeWarning)
            # Replace every occurrence of the key inside the pipeline with its corresponding value
            pipeline = pipeline.replace(key, str(parameters[key]))

    try:
        json.loads(pipeline)  # Throws an error if the pipeline is poorly formatted
    except json.decoder.JSONDecodeError as exception:
        raise ValueError("Pipeline was poorly formatted: \n" + pipeline + "\n" + str(exception))

    # Run PDAL with the pipeline as the stdin
    commands = ["pdal", "pipeline", "--stdin", "--metadata", os.path.join(temp_dir.name, "meta.json")]
    if not stream:
        commands.append("--nostream")
    subprocess.run(commands, input=pipeline, check=True, encoding="utf-8")

    # Load the temporary metadata file
    with open(os.path.join(temp_dir.name, "meta.json")) as infile:
        output_meta = json.load(infile)

    # Save it with a different name if one was provided
    if output_metadata_file is not None:
        with open(output_metadata_file, "w") as outfile:
            json.dump(output_meta, outfile)

    return output_meta


def run_icp(reference_cloud_filepath: str, aligned_cloud_filepath: str) -> str:
    """
    Run ICP registration between a reference and (to-be-) aligned point cloud.

    param: reference_cloud: The path to the point cloud acting reference.
    param: aligned_cloud: The path to the point cloud to be aligned.

    return: transform: The calculated transform between the two point clouds.
    """
    # Name for the ICP registration metadata file.
    icp_meta_filename = os.path.join(os.path.dirname(aligned_cloud_filepath), "icp_meta.json")
    transformed_cloud_output_filepath = os.path.join(
        os.path.dirname(aligned_cloud_filepath), "point_cloud_registered.las")

    # The pipeline to feed PDAL
    icp_pipeline = '''
    [
        "CLOUD1",
        "CLOUD2",
        {
           "type":"filters.icp"
        },
        "OUTPUT"
    ]
    '''

    # Run PDAL with the given pipeline and save the metadata to a file.
    metadata = run_pdal_pipeline(
        pipeline=icp_pipeline,
        output_metadata_file=icp_meta_filename,
        parameters={
            "CLOUD1": reference_cloud_filepath,
            "CLOUD2": aligned_cloud_filepath,
            "OUTPUT": transformed_cloud_output_filepath
        }
    )

    # Extract the resultant transform from the metadata
    transform = metadata["stages"]["filters.icp"]["composed"].replace("\n", " ")

    return transform


def transform_camera_locations(locations_filepath: str, output_filepath: str, transform: str):
    """
    Take a list of camera locations and translate them using a transformation matrix.

    param: cameras_filename: The filename of the camera location file.
    param: cameras_output_filename: The filename for the transformed camera locations.
    param: transform: The 4x4 transformation matrix in string format.
    """
    temp_dir = tempfile.TemporaryDirectory()
    pdal_infile_path = os.path.join(temp_dir.name, "pdal_cameras_in.txt")
    pdal_outfile_path = os.path.join(temp_dir.name, "pdal_cameras_out.txt")

    # Open the camera location file
    header = "label,easting,northing,altitude"
    with open(locations_filepath) as infile:
        camera_locations = infile.read().splitlines()

    # Make sure that the cameras have a correct column order
    assert camera_locations[0] == header, "{} has wrong header type".format(locations_filepath)

    # Remove the header from the data
    camera_locations.pop(0)

    # Remove the label column for PDAL so only the X,Y,Z columns are left.
    pdal_file = [",".join(line.split(",")[1:]) for line in camera_locations]

    # Write the PDAL-formatted locations file
    with open(pdal_infile_path, "w") as outfile:
        outfile.write("X,Y,Z\n")
        for line in pdal_file:
            outfile.write(line + "\n")

    # The pipeline for PDAL to process
    translation_pipeline = '''
    [
        {
            "type":"readers.text",
            "filename":"PDAL_INFILE"
        },
        {
            "type":"filters.transformation",
            "matrix": "ICP_MATRIX"
        },
        {
            "type":"writers.text",
            "filename":"PDAL_OUTFILE"
        }
    ]
    '''

    run_pdal_pipeline(
        pipeline=translation_pipeline,
        parameters={
            "PDAL_INFILE": pdal_infile_path,
            "PDAL_OUTFILE": pdal_outfile_path,
            "ICP_MATRIX": transform
        }
    )

    # Read the result and remove the header
    with open(pdal_outfile_path) as infile:
        pdal_result = infile.read().splitlines()
        pdal_result.pop(0)

    # Write the transformed camera locations file.
    with open(output_filepath, "w") as outfile:
        outfile.write(header + "\n")
        for i, camera_location_row in enumerate(camera_locations):
            label = camera_location_row.split(",")[0]
            outfile.write(label + "," + pdal_result[i] + "\n")


def calculate_dem_extent(chunk):
    """
    Calculate a reasonable extent for a DEM based on the chunk's camera positions.

    param: chunk: The chunk on which to calculate.

    return: extent: The reasonable extent as (xmin, xmax, ymin, ymax).
    """
    x_positions = []
    y_positions = []
    for camera in chunk.cameras:
        x_position, y_position, _ = camera.reference.location
        x_positions.append(x_position)
        y_positions.append(y_position)

    min_x = min(x_positions)
    min_y = min(y_positions)
    max_x = max(x_positions)
    max_y = max(y_positions)

    extent = [
        min_x - (min_x % CONSTANTS.dem_gridsize) - CONSTANTS.dem_gridsize,
        max_x - (max_x % CONSTANTS.dem_gridsize) + CONSTANTS.dem_gridsize,
        min_y - (min_y % CONSTANTS.dem_gridsize) - CONSTANTS.dem_gridsize,
        max_y - (max_y % CONSTANTS.dem_gridsize) + CONSTANTS.dem_gridsize
    ]

    return extent


def get_first_point_info(point_cloud_filepath: str) -> dict[str, Any]:
    """
    Get the data of the first point (index=0) of a point cloud.

    param: point_cloud_filepath: The path to the point cloud.
    return: point_data: The data of the first point.
    """
    result = subprocess.run(["pdal", "info", point_cloud_filepath, "-p", "0"],
                            stdout=subprocess.PIPE, check=True, encoding="utf-8")
    point_data = json.loads(result.stdout)["points"]["point"]

    return point_data
