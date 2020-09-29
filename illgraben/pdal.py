"""PDAL wrappers."""
import subprocess
import os
import json

# If an absolute path is needed, this should be changed to it.
PDAL_COMMAND = "pdal"


def run_pipeline(pipeline, metadata_output_file=None, stream=True):
    """
    Run a PDAL pipeline.

    param: pipeline: The input pipeline string.
    type: pipeline: str
    param: metadata_output_file: Optional. Output filename of the metadata of the pipeline.
    type: Optional[str]
    param: stream: Whether to use stream mode. Will not work with certain tools.
    type: bool

    return: result: The stdout of PDAL, usually just a progress update.
    type: result: str
    """

    commands = [PDAL_COMMAND, "pipeline", "--stdin"]
    if metadata_output_file is not None:
        commands.append("--metadata")
        commands.append(metadata_output_file)

    if not stream:
        commands.append("--nostream")

    result = subprocess.run(commands, check=True, stdout=subprocess.PIPE,
                            input=pipeline, encoding="utf-8").stdout
    return result


def icp_coregistration(reference_filename, aligned_filename, composed=True):

    aligned_basename, aligned_extension = os.path.splitext(aligned_filename)
    icp_pipeline = '''
    [
        "REFERENCE_FILENAME",
        "ALIGNED_FILENAME",
        {
            "type": "filters.icp"
        },
        "OUTPUT_FILENAME"
    ]'''\
            .replace("REFERENCE_FILENAME", reference_filename)\
            .replace("ALIGNED_FILENAME", aligned_filename)\
            .replace("OUTPUT_FILENAME", aligned_basename + "_ICP_aligned" + aligned_extension)

    output_metadata = os.path.splitext(aligned_filename)[0] + "_icp_meta.json"

    run_pipeline(icp_pipeline, metadata_output_file=output_metadata)

    with open(output_metadata) as infile:
        metadata = json.load(infile)

    if composed:
        transform = metadata["stages"]["filters.icp"]["composed"].replace("\n", " ")
    else:
        transform = metadata["stages"]["filters.icp"]["transform"].replace("\n", " ")

    return transform


def transform_point(point, transform, invert=False):

    with open("temp/temp_point.xyz", "w") as outfile:
        outfile.write("X,Y,Z\n")
        outfile.write("{},{},{}\n".format(*point))

    transform_pipeline = '''
    [
        {
            "type": "readers.text",
            "filename": "temp/temp_point.xyz"
        },
        {
            "type": "filters.transformation",
            "matrix": "TRANSFORM",
            "invert": "INVERT"
        },
        "temp/temp_point_transformed.xyz"
    ]'''.replace("TRANSFORM", transform).replace("INVERT", "true" if invert else "false")

    run_pipeline(transform_pipeline)

    with open("temp/temp_point_transformed.xyz") as infile:
        transformed_point_str = infile.read().splitlines()[1]
        transformed_point = [float(coord) for coord in transformed_point_str.split(",")]

    return transformed_point


def get_first_point_info(point_cloud_filename):

    result = subprocess.run(["pdal", "info", point_cloud_filename, "-p", "0"],
                            stdout=subprocess.PIPE, check=True, encoding="utf-8")
    point_data = json.loads(result.stdout)["points"]["point"]

    return point_data
