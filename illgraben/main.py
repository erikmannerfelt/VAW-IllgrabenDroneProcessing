"""Align images of the Illigraben dataset autonomously."""

import os
import shutil

import Metashape as ms

from illgraben import metashape_tools
from illgraben.constants import CONSTANTS
from illgraben.files import PROCESSING_FOLDER, log
from illgraben.utilities import big_print, no_stdout, notify


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
        shutil.rmtree(PROCESSING_FOLDER)
        os.mkdir(PROCESSING_FOLDER)

    # Instantiate a Metashape document
    doc = ms.Document()
    # Set a fitting name (with its full path)
    document_name = os.path.join(os.getcwd(), PROCESSING_FOLDER, "Illgraben.psx")
    # Check that an input folder is present
    if not os.path.isdir("input"):
        raise ValueError("No input/ folder in the working directory!")

    # Load an already existing document if it shouldn't be redone
    if os.path.isfile(document_name) and not redo:
        doc.open(document_name)
        log("Existing document loaded")
    # Otherwise, make a new one
    else:
        metashape_tools.save(doc, document_name)
        log("New document created")

    # Check that the document is not in readonly mode
    assert not doc.read_only, "Document is in read-only mode."

    # Load the reference chunk label
    reference_chunk_label = open("input/reference.txt").read().strip()

    # Load/create the reference chunk and a list of chunks to be aligned
    reference_chunk, chunks_to_be_aligned = metashape_tools.initialise_chunks(doc, reference_chunk_label)

    # Check that a reference chunk exists
    assert reference_chunk is not None, "Reference chunk {} could not be found".format(reference_chunk_label)

    metashape_tools.save(doc)

    # Generate low-resolution dense point clouds for fine-grained ICP
    # TODO: Dense clouds are right now created in the same resolution
    for chunk in doc.chunks:
        # Check if the dense cloud step should be skipped (if a dense cloud with the proper resolution exists)
        if chunk.dense_cloud is not None and\
                chunk.dense_cloud.meta["BuildDepthMaps/downscale"] == str(CONSTANTS.depth_map_downscaling):
            continue
        big_print("Generating small dense cloud for {}".format(chunk.label))
        # TODO: Change to build smaller dense clouds?
        try:
            metashape_tools.build_dense_cloud(chunk, point_confidence=True)
        except Exception as exception:
            if "Zero resolution" in str(exception):
                continue
            if "Assertion 23910910127 failed" in str(exception):
                continue

            raise exception
        metashape_tools.export_dense_cloud(chunk, "dense_cloud_for_ICP.ply")
        metashape_tools.save(doc)

    big_print("Running fine grained ICP on stable ground features")
    for chunk in chunks_to_be_aligned:
        # Check if an automatic ICP exists.
        if "auto_ICP_000" in (marker.label for marker in chunk.markers):
            print("auto_ICPs already seem to exist. Skipping {}".format(chunk.label))
            continue
        # Create automatic ICP tie points
        metashape_tools.align_stable_ground_locations(reference_chunk, chunk)
        with no_stdout():
            chunk.optimizeCameras()
        metashape_tools.save(doc)

    # Generate dense point clouds
    if not redo:
        big_print("Checking for dense clouds")
    for chunk in doc.chunks:
        # Check if dense cloud exists, and whether the dense cloud has the right resolution
        # If any of those criteria are false, it rebuilds the dense cloud
        if chunk.dense_cloud is not None and\
                chunk.dense_cloud.meta["BuildDepthMaps/downscale"] == str(CONSTANTS.depth_map_downscaling):
            continue
        big_print("Generating dense cloud for {}".format(chunk.label))

        try:
            metashape_tools.build_dense_cloud(chunk, point_confidence=True)
        except Exception as exception:
            if "Zero resolution" not in str(exception):
                raise exception
        metashape_tools.save(doc)

    # Generate DEMs
    if not redo:
        big_print("Checking for DEMs")
    for chunk in doc.chunks:
        if chunk.elevation is not None:
            continue
        big_print("Generating DEM for {}".format(chunk.label))
        metashape_tools.generate_dem(chunk, redo=redo)
        metashape_tools.save(doc)

    # Generate orthomosaics
    if not redo:
        big_print("Checking for orthomosaics")
    for chunk in doc.chunks:
        if chunk.orthomosaic is not None:
            continue
        big_print("Generating orthomosaic for {}".format(chunk.label))
        metashape_tools.build_orthomosaic(chunk)
        metashape_tools.save(doc)
    notify("Finished processing")
