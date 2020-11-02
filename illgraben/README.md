# File descriptions

### main.py
The main script where the datasets are aligned.

## constants.py
Readonly constants (e.g. the CRS, DEM resolution etc.) that are read throughout the module.

### files.py
Handling of input and output files, as well as hardcoded filepaths.
Functions and variables are used throughout the module.

### metashape\_tools.py
Wrappers and tools around the Metashape API.
Functions are called from here to `main.py`.

### processing\_tools.py
Wrappers and tools around other CLI interfaces (mostly just PDAL)
Functions are called from here to `metashape_tools.py`

### stable\_ground\_icp.py
Functions for the feature-wise stable ground ICP registration (later run in `metashape_tools.py`).

### utilities.py
Python and C helper functions, used throughout the module.
