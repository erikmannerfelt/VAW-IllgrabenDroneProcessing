"""Functions and constants to keep track of input and output files."""
import os
import time

PROCESSING_FOLDER = "processing/"
OUTPUT_FOLDER = "output/"


INPUT_FILEPATHS = {
    "stable_ground_points": "input/stable_ground_points.xyz",
    "reference_chunk_label": "input/reference.txt"
}


if not os.path.isdir(PROCESSING_FOLDER):
    os.mkdir(PROCESSING_FOLDER)


def log(text: str) -> None:
    """Write a line of text to the log."""
    # Get the current date and time in the local timezone
    now: str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

    with open(os.path.join(PROCESSING_FOLDER, "processing.log"), "a+") as outfile:
        outfile.write(f"{now} {text}")
