# classes.py
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table
from .functions import get_stars_in_region




class RegionStars:
    def __init__(self, coord: SkyCoord, radius: Angle) -> None:
        stars = get_stars_in_region(coord, radius)
        return stars

    def star(self, gaiaID):
        stars = self.stars
        star = stars[stars['Source'] == gaiaID]
        return star
