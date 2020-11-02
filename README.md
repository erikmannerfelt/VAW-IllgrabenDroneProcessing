# illgraben: Automated drone survey processing for Illgraben in Switzerland.

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
