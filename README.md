# illgraben: Automated drone survey processing for Illgraben in Switzerland.

This script utilises the Agisoft Metashape API and feature-wise ICP registration to coalign and process an arbitrary amount of repeat datasets at Illgraben in Switzerland.
The feature-wise ICP registration works by extracting manually picked areas of ground assumed to be stable, and using the post-ICP transform to generate automatic Metashape markers which aligns one dataset to a reference.

## Requirements

* Conda
* Agisoft Metashape 1.6.5 (for another version, change `environment.yml` accordingly.


## Installation
First, setup the conda environment:
```bash
$ conda env create -f environment.yml  # this will create an environment called 'illgraben'

$ conda activate illgraben
```
Then, make sure that the environment variable `agisoft_LICENSE` points toward a valid `metashape.lic` file.

#### Data structure
The input folder structure is fixed:

```bash
input
|    reference.txt  #  A text file with the name of the reference survey (the folder name)
|    stable_ground_points.xyz  # A csv of X/Y/Z coordinates for stable points in the terrain.
└────surveys  # A folder of folders for the surveys
|    └────'survey1'  # Replace this with the name of the actual survey
|    |    └────images
|    |    |    |    'img1.jpg'
|    |    |    |    'img2.jpg'
|    |    |    |     ...
|    └────'survey2'  # Replace this with the name of the actual survey
|    |    └────images
|    |    |    |    'img1.jpg'
|    |    |    |    'img2.jpg'
|    |    |    |     ...
|    ...
```

## Usage
First, write the name of a reference survey in the `input/reference.txt` file.
This will be used to coalign the other surveys.
Then, define some stable points in the terrain and write their coordinates to `input/stable_ground_points.xyz` (in the same CRS as the project).

Make sure you are in the working directory with the `illgraben` environment activated, then run the `illgraben` module directly:
```bash
$ python -m illgraben
```

# Todo list
* Generate good error measures for the surveys
* Align the "reference survey" to the swissAlti3D "higher reference" dataset
* Add a filtering method for poor ICP matches.
* Additionally, skip extracted features that have a low or no point count.
* Make sure that the above fixes improve the alignment quality of all the surveys, especially in the top
* Remove the forest from the subsequent DEMs (GNDVI maybe?)

