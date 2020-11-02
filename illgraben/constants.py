"""Constant parameters that are used throughout the module."""

from illgraben.utilities import ConstantType


class Constants(ConstantType):  # pylint: disable=R0903
    """Readonly constants."""

    low_height_threshold: float = 80.0  # Remove cameras in the lower X metres (for ascent and descent images)
    dem_gridsize: float = 1.0  # The output DEM resolution in metres
    orthomosaic_resolution: float = 0.25  # The output orthomosaic resolution in metres
    depth_map_downscaling = 4  # The downscaling of depth maps to create DEMs (4 == Medium)
    # temporary_depth_map_downscaling = 8  # The downscaling of depth maps to coalign surveys (8 == Low) CURRENTLY UNUSED
    coordinate_system: str = "EPSG::2056"  # CH1903/LV95


CONSTANTS = Constants()
