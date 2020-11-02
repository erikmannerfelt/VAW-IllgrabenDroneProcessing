"""Align images of the Illigraben dataset autonomously."""

import json
import os
import shutil
import subprocess
import time

import Metashape as ms

from . import stable_ground
from .utilities import no_stdout

TEMP_FOLDER = "temp"
DEM_GRIDSIZE = 1
DEPTH_MAP_DOWNSCALING = 4
TEMPORARY_DEPTH_MAP_DOWNSCALING = 8
ORTHOMOSAIC_RESOLUTION = 0.25

if not os.path.isdir(TEMP_FOLDER):
    os.mkdir(TEMP_FOLDER)


def big_print(string):
    """
    Print a string with space above and below it.

    param: string: The string to print.
    type: string: str
    """
    separator = "============================="

    current_time = time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(time.time()))

    print("\n", separator, "\n")
    print(current_time.replace("T", "  "), "\t", string)

    log_filepath = os.path.join(TEMP_FOLDER, "process.log")
    if os.path.isfile(log_filepath):
        with open(log_filepath, "a+") as outfile:
            outfile.write("{}\t{}\n".format(current_time, string))
    print("\n", separator, "\n")


def run_icp(reference_cloud, aligned_cloud):
    """
    Run ICP registration between a reference and (to-be-) aligned point cloud.

    param: reference_cloud: The path to the point cloud acting reference.
    type: reference_cloud: str
    param: aligned_cloud: The path to the point cloud to be aligned.
    type: reference_cloud: str

    return: transform: The calculated transform between the two point clouds.
    """
    # Name for the ICP registration metadata file.
    icp_meta_filename = os.path.join(os.path.dirname(aligned_cloud), "icp_meta.json")
    transformed_cloud_output_filename = os.path.join(os.path.dirname(aligned_cloud), "point_cloud_registered.las")

    # Remove any previous metadata file
    # If PDAL errors out, it won't stop the program. If no metadata file exists, however, the program stops.
    if os.path.isfile(icp_meta_filename):
        os.remove(icp_meta_filename)

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
    '''.replace("CLOUD1", reference_cloud).replace("CLOUD2", aligned_cloud).replace("OUTPUT", transformed_cloud_output_filename)

    # Run PDAL with the given pipeline and save the metadata to a file.
    subprocess.run("echo '{}' | pdal pipeline --stdin --metadata {}".format(icp_pipeline,
                                                                            icp_meta_filename), shell=True, check=True)

    # Open the resultant metadata file and parse it
    with open(icp_meta_filename) as infile:
        metadata = json.load(infile)

    # Extract the resultant transform from the metadata
    transform = metadata["stages"]["filters.icp"]["composed"].replace("\n", " ")

    return transform


def transform_cameras(cameras_filename, cameras_output_filename, transform):
    """
    Take a list of camera locations and translate them using a transformation matrix.

    param: cameras_filename: The filename of the camera location file.
    type: cameras_filename: str
    param: cameras_output_filename: The filename for the transformed camera locations.
    type: cameras_output_filename: str
    param: transform: The 4x4 transformation matrix in string format.
    type: transform: str
    """
    temp_pdal_infile = os.path.join(os.path.dirname(cameras_filename),  "pdal_cameras_in.txt")
    temp_pdal_outfile = os.path.join(os.path.dirname(cameras_filename), "pdal_cameras_out.txt")

    # Remove previous results if they exist.
    # PDAL does not return a proper error if it fails, so the next step (read results) is its error-trigger.
    for filepath in [temp_pdal_infile, temp_pdal_outfile]:
        if os.path.isfile(filepath):
            os.remove(filepath)

    # Open the camera location file
    header = "label,easting,northing,altitude"
    with open(cameras_filename) as infile:
        camera_locations = infile.read().splitlines()

    # Make sure that the cameras have a correct column order
    assert camera_locations[0] == header, "{} has wrong header type".format(cameras_filename)

    # Remove the header from the data
    camera_locations.pop(0)

    # Remove the label column for PDAL so only the X,Y,Z columns are left.
    pdal_file = [",".join(line.split(",")[1:]) for line in camera_locations]

    # Write the PDAL-formatted locations file
    with open(temp_pdal_infile, "w") as outfile:
        outfile.write("X,Y,Z\n")
        for line in pdal_file:
            outfile.write(line + "\n")

    # The pipeline for PDAL to process
    translate_pipeline = '''
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
    '''.replace("PDAL_INFILE", temp_pdal_infile).replace("PDAL_OUTFILE", temp_pdal_outfile).replace("ICP_MATRIX", transform)

    # Run the PDAL pipeline
    subprocess.run("echo '{}' | pdal pipeline --stdin".format(translate_pipeline), shell=True, check=True)

    # Read the result and remove the header
    with open(temp_pdal_outfile) as infile:
        pdal_result = infile.read().splitlines()
        pdal_result.pop(0)

    # Write the transformed camera locations file.
    with open(cameras_output_filename, "w") as outfile:
        outfile.write(header + "\n")
        for i in range(len(camera_locations)):
            label = camera_locations[i].split(",")[0]
            outfile.write(label + "," + pdal_result[i] + "\n")


def export_camera_reference(chunk, filename):
    """
    Export the camera reference information from a chunk in CSV format.

    param: chunk: The Metashape chunk to export from.
    type: chunk: Metashape Chunk
    param: filename: The output filename.
    type: filename: str
    """
    header = "label,easting,northing,altitude"

    with open(filename, "w") as outfile:
        # Write the header
        outfile.write(header + "\n")
        # Loop through each camera
        for cam in chunk.cameras:
            # Get its label from its filepath
            # Importing needs the extension, e.g. .JPG, which is not always in the chunk label
            label = os.path.basename(cam.photo.path)
            easting, northing, altitude = cam.reference.location
            outfile.write("{},{},{},{}".format(label, easting, northing, altitude) + "\n")


def import_camera_reference(chunk, filename):
    """
    Import camera reference data from a CSV.

    It is assumed the the entries are only for cameras and that they share the chunk crs.

    param: chunk: The chunk to import the reference to.
    type: chunk: Metashape Chunk
    param: filename: The input filename.
    type: filename: str
    """
    with no_stdout():
        chunk.importReference(path=filename, delimiter=",", columns="nxyz", create_markers=False,
                              crs=chunk.crs, items=ms.ReferenceItemsCameras)
        chunk.updateTransform()


def align_chunk(reference_chunk, aligned_chunk):
    """
    Run all functions to align one chunk to a reference chunk.

    param: reference_chunk: The chunk to act as reference.
    type: reference_chunk: Metashape Chunk.
    param: aligned_chunk: The chunk to be aligned.
    type: aligned_chunk: Metashape Chunk.

    """
    # Check that they share the same CRS.
    assert reference_chunk.crs.name == aligned_chunk.crs.name, "Chunk CRS's differ!"
    # Export the sparse point cloud for both chunks.
    save_point_cloud(aligned_chunk)

    # Estimate the transform (offset) using ICP registration.
    transform = run_icp(os.path.join(TEMP_FOLDER, reference_chunk.label, "point_cloud.las"),
                        os.path.join(TEMP_FOLDER, aligned_chunk.label, "point_cloud.las"))

    # Export the camera locations (before registration)
    cameras_filename = os.path.join(TEMP_FOLDER, aligned_chunk.label, "cameras.csv")
    export_camera_reference(aligned_chunk, filename=cameras_filename)

    # Transform the camera locations and save the output file.
    translated_cameras_filename = os.path.join(TEMP_FOLDER, aligned_chunk.label, "translated_cameras.csv")
    transform_cameras(cameras_filename=cameras_filename,
                      cameras_output_filename=translated_cameras_filename, transform=transform)

    # Import the tranformed camera reference information and update the chunk transform.
    import_camera_reference(aligned_chunk, translated_cameras_filename)


def save_point_cloud(chunk):
    """
    Save a sparse point cloud in the chunk's temp folder.

    param: chunk: What chunk to export the point cloud from.
    type: chunk: Metashape Chunk
    """
    with no_stdout():
        chunk.exportPoints(os.path.join(TEMP_FOLDER, chunk.label, "point_cloud.las"),
                           source_data=ms.DataSource.PointCloudData,
                           save_normals=False,
                           save_colors=False,
                           save_classes=False,
                           crs=chunk.crs)


def init_chunk(doc, label):
    """
    Initialise a chunk.

    Sets the appropriate settings and aligns the images.

    param: doc: The Metashape document instance.
    type: doc: Metashape Document
    param: label: What label to assign the chunk.
    type: label: str
    """
    big_print("Importing survey: {}".format(label))

    # Create the chunk's temporary folder
    if not os.path.isdir(os.path.join("temp", label)):
        os.makedirs(os.path.join("temp", label))

    chunk = doc.addChunk()
    chunk.label = label

    # Add the chunk's images
    photo_dir = os.path.join("input", label, "images")
    photos = [os.path.join(photo_dir, photo) for photo in os.listdir(photo_dir)]
    with no_stdout():
        chunk.addPhotos(photos)

    # Set the x/y/z location accuracy
    chunk.camera_location_accuracy = [2.5] * 3

    # Convert camera coordinates to CH1903 / LV03
    out_crs = ms.CoordinateSystem("EPSG::21781")
    for cam in chunk.cameras:
        if cam.reference.location is None:
            continue
        cam.reference.location = ms.CoordinateSystem.transform(
            cam.reference.location, chunk.crs, out_crs)
    chunk.crs = out_crs

    # Remove cameras at a low altitude (for the initial ascent and final descent)
    min_altitude = min(
        [camera.reference.location.z for camera in chunk.cameras])
    for camera in chunk.cameras:
        if camera.reference.location.z < (min_altitude + 80):
            chunk.remove(camera)

    # Remove the rotation information as it seems to be erroneous
    for cam in chunk.cameras:
        cam.reference.rotation = None

    # Enable rolling shutter compensation on all cameras.
    for sensor in chunk.sensors:
        sensor.rolling_shutter = True

    # Align the cameras
    with no_stdout():
        chunk.matchPhotos()
        chunk.alignCameras()


def initialise_chunks(doc, reference_chunk_label, redo):

    # Run alignment if it should be redone, or if it has never been done
    align = redo or len(doc.chunks) == 0

    if not align:
        for chunk in doc.chunks:
            # Check if any chunk does not have a sparse point cloud
            align = chunk.point_cloud is None

            if align:
                big_print("Chunks exist but some or all are not aligned. Running the alignment again.")
                break

    if align:
        big_print("Aligning all chunks")
        # List all of the surveys in the input folder
        surveys = []
        for item in os.listdir("input/"):
            if os.path.isdir(os.path.join("input/", item)):
                surveys.append(item)

        # Check that at least two surveys exist
        assert len(surveys) > 1, "At least two surveys are needed to process. Found {}".format(len(surveys))

        # Initialise a chunk (import, set params, align) for each survey.
        for survey in surveys:
            init_chunk(doc, label=survey)
    else:
        big_print("Chunks already exist. Skipping to next step")

    # Separate the chunks to be aligned with the reference chunk.
    chunks_to_be_aligned = []
    reference_chunk = None
    for chunk in doc.chunks:
        if reference_chunk_label == chunk.label:
            reference_chunk = chunk
        else:
            chunks_to_be_aligned.append(chunk)

    if align:
        # Save the reference chunk point cloud, so it can be used with the chunk alignment
        save_point_cloud(reference_chunk)

        # Align all chunks using ICP
        for chunk_to_be_aligned in chunks_to_be_aligned:
            align_chunk(reference_chunk, chunk_to_be_aligned)

    return reference_chunk, chunks_to_be_aligned


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
        min_x - (min_x % DEM_GRIDSIZE) - DEM_GRIDSIZE,
        max_x - (max_x % DEM_GRIDSIZE) + DEM_GRIDSIZE,
        min_y - (min_y % DEM_GRIDSIZE) - DEM_GRIDSIZE,
        max_y - (max_y % DEM_GRIDSIZE) + DEM_GRIDSIZE
    ]

    return extent


def generate_dem(chunk, redo):
    """
    Generate a DEM using PDAL.

    Generating a DEM in PDAL is better than in Metashape since the grid size can be specified and interpolation is off.

    param: chunk: The input chunk.
    param: redo: Whether to redo the analysis even if it is partially completed
    """
    extent = calculate_dem_extent(chunk)

    dense_cloud_path = os.path.join(TEMP_FOLDER, chunk.label, "dense_cloud.ply")

    if not os.path.isfile(dense_cloud_path) or redo:
        with no_stdout():
            chunk.exportPoints(
                dense_cloud_path,
                source_data=ms.DataSource.DenseCloudData,
                crs=chunk.crs,
                save_confidence=True)

    output_raster_path = os.path.join(os.path.dirname(dense_cloud_path), "dem.tif")

    variables = {
        "DENSE_CLOUD_PATH": dense_cloud_path,
        "GRID_SIZE": DEM_GRIDSIZE,
        "MIN_X": extent[0],
        "MAX_X": extent[1],
        "MIN_Y": extent[2],
        "MAX_Y": extent[3],
        "OUTPUT_RASTER_PATH": output_raster_path
    }
    dem_pipeline = '''
    [
        "DENSE_CLOUD_PATH",
        {
            "type": "filters.range",
            "limits": "confidence[2:]"
        },
        {
            "resolution": GRID_SIZE,
            "bounds": "([MIN_X, MAX_X],[MIN_Y, MAX_Y])",
            "output_type": "mean",
            "filename": "OUTPUT_RASTER_PATH"
        }
    ]
    '''
    # Fill the data into the template
    for key, value in variables.items():
        dem_pipeline = dem_pipeline.replace(key, str(value))

    # Run PDAL
    print("Running PDAL")
    if not os.path.isfile(output_raster_path) or redo:
        subprocess.run("echo '{}' | pdal pipeline --stdin".format(dem_pipeline), shell=True, check=True)

    # Import the raster
    chunk.importRaster(path=output_raster_path, crs=chunk.crs)


def save(doc, filename=None):
    """Save the document"""
    if filename is None:
        with no_stdout():
            doc.save()
        print("Saved project")
    else:
        doc.save(filename)


def main(redo=False):
    """
    Run the entire processing pipeline.

    1. Import the images in separate chunks and set appropriate settings
    2. Align the images
    3. Align the chunks to a reference chunk using ICP
    4. Generate dense clouds, DEMs and orthomosaics
    """
    if redo:
        big_print("Redo flag set to True. Redoing steps that already seem to exist")
    else:
        big_print("Redo flag set to False. Not redoing steps that already seem to exist")
    # Remove all temporary results if the analysis should be redone.
    if redo:
        shutil.rmtree("temp")
        os.mkdir("temp")

    # Instantiate a Metashape document
    doc = ms.Document()
    # Set a fitting name (with its full path)
    document_name = os.path.join(os.getcwd(), "Illgraben.psx")
    # Check that an input folder is present
    if not os.path.isdir("input"):
        raise ValueError("No input/ folder in the working directory!")

    # Load an already existing document if redo is not False
    if os.path.isfile(document_name) and not redo:
        doc.open(document_name)
    # Otherwise, make a new one
    else:
        save(doc, document_name)
        with open(os.path.join(TEMP_FOLDER, "process.log"), "w") as outfile:
            outfile.write("New process started.")

    # Check that the document is not in readonly mode
    assert not doc.read_only, "Document is in read-only mode."

    # Load the reference chunk label
    reference_chunk_label = open("input/reference.txt").read().strip()

    reference_chunk, chunks_to_be_aligned = initialise_chunks(doc, reference_chunk_label, redo)

    # Check that a reference chunk exists
    assert reference_chunk is not None, "Reference chunk {} could not be found".format(reference_chunk_label)

    save(doc)

    # Generate low-resolution dense point clouds for fine-grained ICP
    for chunk in doc.chunks:
        if chunk.dense_cloud is not None and not redo and\
                chunk.dense_cloud.meta["BuildDepthMaps/downscale"] == str(TEMPORARY_DEPTH_MAP_DOWNSCALING):
            continue
        big_print("Generating small dense cloud for {}".format(chunk.label))
        with no_stdout():
            chunk.buildDepthMaps(downscale=TEMPORARY_DEPTH_MAP_DOWNSCALING,
                                 filter_mode=ms.FilterMode.AggressiveFiltering)
            chunk.buildDenseCloud()
            chunk.exportPoints(
                os.path.join(TEMP_FOLDER, chunk.label, "dense_cloud_lowres.ply"),
                source_data=ms.DataSource.DenseCloudData,
                crs=chunk.crs)
        save(doc)

    big_print("Running fine grained ICP on stable ground features")
    for chunk in chunks_to_be_aligned:
        # Check if an automatic ICP exists.
        if "auto_ICP_000" in (marker.label for marker in chunk.markers) and not redo:
            print("auto_ICPs already seem to exist. Skipping {}".format(chunk.label))
        stable_ground.align_stable_ground_locations(reference_chunk, chunk)
        with no_stdout():
            chunk.optimizeCameras()
        save(doc)

    # Generate dense point clouds
    if not redo:
        big_print("Checking for dense clouds")
    for chunk in doc.chunks:
        # Check if dense cloud exists, if it should be redone, and whether the dense cloud has the right resolution
        # If any of those criteria are false, it rebuilds the dense cloud
        if chunk.dense_cloud is not None and not redo and\
                chunk.dense_cloud.meta["BuildDepthMaps/downscale"] == str(DEPTH_MAP_DOWNSCALING):
            continue
        big_print("Generating dense cloud for {}".format(chunk.label))
        with no_stdout():
            chunk.buildDepthMaps(downscale=DEPTH_MAP_DOWNSCALING, filter_mode=ms.FilterMode.AggressiveFiltering)
            chunk.buildDenseCloud(point_confidence=True)
        save(doc)

    # Generate DEMs
    if not redo:
        big_print("Checking for DEMs")
    for chunk in doc.chunks:
        if chunk.elevation is not None and not redo:
            continue
        big_print("Generating DEM for {}".format(chunk.label))
        generate_dem(chunk, redo)
        save(doc)

    # Generate orthomosaics
    if not redo:
        big_print("Checking for orthomosaics")
    for chunk in doc.chunks:
        if chunk.orthomosaic is not None and not redo:
            continue
        big_print("Generating orthomosaic for {}".format(chunk.label))
        with no_stdout():
            chunk.buildOrthomosaic(surface_data=ms.DataSource.ElevationData, resolution=ORTHOMOSAIC_RESOLUTION)
        save(doc)
