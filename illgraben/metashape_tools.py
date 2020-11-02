"""Tools to interact with the Metashape API."""
import os
import subprocess

import Metashape as ms
import statictypes

from illgraben import processing_tools, stable_ground_icp
from illgraben.constants import CONSTANTS
from illgraben.files import PROCESSING_FOLDER, log
from illgraben.utilities import big_print, no_stdout


@statictypes.enforce
def import_camera_reference(chunk: ms.Chunk, filepath: str) -> None:
    """
    Import camera reference data from a CSV.

    It is assumed the the entries are only for cameras and that they share the chunk crs.

    param: chunk: The chunk to import the reference to.
    param: filepath: The input filename.
    """
    with no_stdout():
        chunk.importReference(path=filepath, delimiter=",", columns="nxyz", create_markers=False,
                              crs=chunk.crs, items=ms.ReferenceItemsCameras)
        chunk.updateTransform()


# @statictypes.enforce
def align_chunk(reference_chunk: ms.Chunk, aligned_chunk: ms.Chunk) -> None:
    """
    Run all functions to align one chunk to a reference chunk.

    param: reference_chunk: The chunk to act as reference.
    param: aligned_chunk: The chunk to be aligned.
    """
    # Check that they share the same CRS.
    assert reference_chunk.crs.name == aligned_chunk.crs.name, "Chunk CRS's differ!"
    # Export the sparse point cloud for both chunks.
    save_point_cloud(aligned_chunk)

    # Estimate the transform (offset) using ICP registration.
    transform = processing_tools.run_icp(os.path.join(PROCESSING_FOLDER, reference_chunk.label, "point_cloud.las"),
                                         os.path.join(PROCESSING_FOLDER, aligned_chunk.label, "point_cloud.las"))

    # Export the camera locations (before registration)
    cameras_filename = os.path.join(PROCESSING_FOLDER, aligned_chunk.label, "cameras.csv")
    export_camera_reference(aligned_chunk, filepath=cameras_filename)

    # Transform the camera locations and save the output file.
    translated_cameras_filename = os.path.join(PROCESSING_FOLDER, aligned_chunk.label, "translated_cameras.csv")
    processing_tools.transform_camera_locations(locations_filepath=cameras_filename,
                                                output_filepath=translated_cameras_filename, transform=transform)

    # Import the tranformed camera reference information and update the chunk transform.
    import_camera_reference(aligned_chunk, translated_cameras_filename)


def save_point_cloud(chunk: ms.Chunk) -> None:
    """
    Save a sparse point cloud in the chunk's temp folder.

    param: chunk: What chunk to export the point cloud from.
    """
    with no_stdout():
        chunk.exportPoints(os.path.join(PROCESSING_FOLDER, chunk.label, "point_cloud.las"),
                           source_data=ms.DataSource.PointCloudData,
                           save_normals=False,
                           save_colors=False,
                           save_classes=False,
                           crs=chunk.crs)


def init_chunk(doc: ms.Document, label: str):
    """
    Initialise a chunk.

    Sets the appropriate settings and aligns the images.

    param: doc: The Metashape document instance.
    param: label: What label to assign the chunk.
    """
    big_print("Importing survey: {}".format(label))

    # Create the chunk's temporary folder
    os.makedirs(os.path.join(PROCESSING_FOLDER, label), exist_ok=True)

    chunk = doc.addChunk()
    chunk.label = label

    # Add the chunk's images
    photo_dir = os.path.join("input/surveys", label, "images")
    photos = [os.path.join(photo_dir, photo) for photo in os.listdir(photo_dir)]
    with no_stdout():
        chunk.addPhotos(photos)

    # Set the x/y/z location accuracy
    chunk.camera_location_accuracy = [2.5] * 3

    # Convert camera coordinates a projected coordinate system
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
        if camera.reference.location.z < (min_altitude + CONSTANTS.low_height_threshold):
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

    print(f"Finished aligning: {label}")
    log(f"Finished aligning: {label}")


def initialise_chunks(doc: ms.Document, reference_chunk_label: ms.Chunk) -> tuple[ms.Chunk, list[ms.Chunk]]:
    """
    Initialise all chunks in a document.

    param: doc: The Metashape document.
    param: reference_chunk_label: The label of the chunk acting reference.
    param: redo: Whether to redo the analysis.
    """
    # Run alignment if it has never been done
    align = len(doc.chunks) == 0

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
        for item in os.listdir("input/surveys"):
            if os.path.isdir(os.path.join("input/surveys", item)):
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


def export_camera_reference(chunk: ms.Chunk, filepath: str) -> None:
    """
    Export the camera reference information from a chunk in CSV format.

    param: chunk: The Metashape chunk to export from.
    param: filepath: The output filename.
    """
    header = "label,easting,northing,altitude"

    with open(filepath, "w") as outfile:
        # Write the header
        outfile.write(header + "\n")
        # Loop through each camera
        for cam in chunk.cameras:
            # Get its label from its filepath
            # Importing needs the extension, e.g. .JPG, which is not always in the chunk label
            label = os.path.basename(cam.photo.path)
            easting, northing, altitude = cam.reference.location
            outfile.write(",".join((label, easting, northing, altitude)) + "\n")


def generate_dem(chunk: ms.Chunk, redo: bool = False):
    """
    Generate a DEM using PDAL.

    Generating a DEM in PDAL is better than in Metashape since the grid size can be specified and interpolation is off.

    param: chunk: The input chunk.
    param: redo: Whether to redo the analysis even if it is partially completed
    """
    extent = processing_tools.calculate_dem_extent(chunk)

    dense_cloud_path = os.path.join(PROCESSING_FOLDER, chunk.label, "dense_cloud.ply")

    if not os.path.isfile(dense_cloud_path) or redo:
        with no_stdout():
            chunk.exportPoints(
                dense_cloud_path,
                source_data=ms.DataSource.DenseCloudData,
                crs=chunk.crs,
                save_confidence=True)

    output_raster_path = os.path.join(os.path.dirname(dense_cloud_path), "dem.tif")

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

    # Run PDAL
    print("Running PDAL")
    processing_tools.run_pdal_pipeline(
        pipeline=dem_pipeline,
        parameters={
            "DENSE_CLOUD_PATH": dense_cloud_path,
            "GRID_SIZE": str(CONSTANTS.dem_gridsize),
            "MIN_X": extent[0],
            "MAX_X": extent[1],
            "MIN_Y": extent[2],
            "MAX_Y": extent[3],
            "OUTPUT_RASTER_PATH": output_raster_path
        }
    )

    # Import the raster
    chunk.importRaster(path=output_raster_path, crs=chunk.crs)


def save(doc, filename=None):
    """Save the document."""
    if filename is None:
        with no_stdout():
            doc.save()
        print("Saved project")
    else:
        doc.save(filename)


def build_dense_cloud(chunk: ms.Chunk, point_confidence: bool = False):
    """
    Build a dense cloud for the selected chunk.

    param: chunk: The chunk to be processed.
    param: point_confidence: Whether to calculate point confidences.
    """
    with no_stdout():
        chunk.buildDepthMaps(downscale=CONSTANTS.depth_map_downscaling,
                             filter_mode=ms.FilterMode.AggressiveFiltering)
        chunk.buildDenseCloud(point_confidence=point_confidence)


def export_dense_cloud(chunk: ms.Chunk, filename: str) -> None:
    """
    Export a chunk's dense cloud.

    param: chunk: The chunk to be processed.
    param: name_template: The name to give the point cloud in its appropriate chunk processing directory.

    """
    with no_stdout():
        chunk.exportPoints(
            os.path.join(PROCESSING_FOLDER, chunk.label, filename),
            source_data=ms.DataSource.DenseCloudData,
            crs=chunk.crs)


def build_orthomosaic(chunk: ms.Chunk) -> None:
    """Build an orthomosaic."""
    with no_stdout():
        chunk.buildOrthomosaic(surface_data=ms.DataSource.ElevationData,
                               resolution=CONSTANTS.orthomosaic_resolution)


def add_correction_marker(chunk: ms.Chunk, initial_position: tuple[float, float, float],
                          corrected_position: tuple[float, float, float], name: str) -> None:
    """
    Add a marker from its assumed position to its actual position.

    param: chunk: The chunk to analyse.
    param: initial_position: The initial X/Y/Z position of the point.
    param: corrected_position: The corrected X/Y/Z position of the point.
    param: name: The name to label the marker.
    """
    # Transform the initial position to local coordinates
    local_initial_position = chunk.transform.matrix.inv().mulp(chunk.crs.unproject(initial_position))

    # Remove any marker with the same name (if the analysis is run on a preexisting chunk)
    for previous_marker in chunk.markers:
        if previous_marker.label == name:
            chunk.remove(previous_marker)

    # Add the marker with the local initial (assumed) position to estimate projections from
    marker = chunk.addMarker(local_initial_position)

    # Pin the projections (to consistently keep them where they are)
    for camera in marker.projections.keys():
        marker.projections[camera].pinned = True

    marker.label = name
    marker.reference.location = corrected_position
    marker.reference.enabled = True


def align_stable_ground_locations(reference_chunk: ms.Chunk, aligned_chunk: ms.Chunk):
    """
    Use feature-wise ICP to align one chunk to another.

    param: reference_chunk: The chunk acting reference.
    param: aligned_chunk: The chunk to be aligned.

    """
    bounds = stable_ground_icp.get_bounding_boxes(stable_ground_icp.get_stable_ground_locations())

    if not os.path.isdir(os.path.join(PROCESSING_FOLDER, reference_chunk.label, "features")):
        print("Extracting reference chunk features")
        stable_ground_icp.extract_features(reference_chunk, bounds)

    # TODO: Fix dangerous assumption that they exist
    if not os.path.isdir(os.path.join(PROCESSING_FOLDER, aligned_chunk.label, "features")):
        print("Extracting features from chunk: {}".format(aligned_chunk.label))
        stable_ground_icp.extract_features(aligned_chunk, bounds)

    points = stable_ground_icp.compare_features(reference_chunk, aligned_chunk, bounds)

    for i, point in enumerate(points):
        add_correction_marker(aligned_chunk, point.start, point.destination, "auto_ICP_{}".format(str(i).zfill(3)))
