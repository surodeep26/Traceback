# classes.py
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table
from .functions import get_stars_in_region
import contextlib
import os
from astropy import units as u


class Star:
    """
    Example usage:
    >>> mystar = Star(3326026144355977472)
    """
    def __init__(self, sourceID) -> None:
        self.name = f"Gaia DR3 {sourceID}"
        self.source = sourceID
        self.coords = SkyCoord.from_name(self.name)
        # supress prints in the with block
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                _stars = get_stars_in_region(skycoordinate=self.coords, radius=0.5*u.arcmin)
        self.info = _stars[_stars['Source'] == self.source]
        self.SkyCoord = self.info['SkyCoord']
        self.ra = self.info['RA_ICRS_1'][0]*u.deg  # Right Ascension (deg)
        self.dec = self.info['DE_ICRS_1'][0]*u.deg  # Declination (deg)
        self.HIP = self.info['HIP']  # Hipparcos ID
        self.TYC2 = self.info['TYC2']  # Tycho-2 ID
        self.rgeo = self.info['rgeo'][0]*u.pc  # Geometric distance (pc)
        self.b_rgeo = self.info['b_rgeo'][0]*u.pc  # Lower bound geometric distance (pc)
        self.B_rgeo = self.info['B_rgeo'][0]*u.pc  # Upper bound geometric distance (pc)
        self.plx = self.info['Plx'][0]*u.mas  # Parallax (mas)
        self.e_plx = self.info['e_Plx'][0]*u.mas  # Parallax error (mas)
        self.pmRA = self.info['pmRA'][0]*u.mas/u.yr  # Proper motion in RA (mas/yr)
        self.e_pmRA = self.info['e_pmRA'][0]*u.mas/u.yr  # Proper motion error in RA (mas/yr)
        self.pmDE = self.info['pmDE'][0]*u.mas/u.yr  # Proper motion in Declination (mas/yr)
        self.e_pmDE = self.info['e_pmDE'][0]*u.mas/u.yr  # Proper motion error in Declination (mas/yr)
        self.RUWE = self.info['RUWE']  # Renormalized unit weight error
        self.Teff = self.info['Teff'][0]*u.K  # Effective temperature (K)
        self.logg = self.info['logg']  # Surface gravity (log(cm/s^2))
        self.Gmag = self.info['Gmag'][0]*u.mag  # G-band magnitude
        self.e_Gmag = self.info['e_Gmag'][0]*u.mag  # Error in G magnitude
        self.BP_RP = self.info['BP-RP'][0]*u.mag  # BP-RP color index
        self.e_BP_RP = self.info['e_BP-RP'][0]*u.mag  # Error in BP-RP color index
        self.BPmag = self.info['BPmag'][0]*u.mag  # BP magnitude
        self.e_BPmag = self.info['e_BPmag'][0]*u.mag  # Error in BP magnitude
        self.RPmag = self.info['RPmag'][0]*u.mag  # RP magnitude
        self.e_RPmag = self.info['e_RPmag'][0]*u.mag  # Error in RP magnitude
        self.RV = self.info['RV'][0]*u.km/u.s  # Radial velocity (km/s)
        self.e_RV = self.info['e_RV'][0]*u.km/u.s  # Radial velocity error (km/s)
        self.FG = self.info['FG'][0]  # Flux in G band
        self.e_FG = self.info['e_FG'][0]  # Error in G band flux
        self.FBP = self.info['FBP'][0]  # Flux in BP band
        self.e_FBP = self.info['e_FBP'][0]  # Error in BP band flux
        self.FRP = self.info['FRP'][0]  # Flux in RP band
        self.e_FRP = self.info['e_FRP'][0]  # Error in RP band flux


class RegionStars:
    def __init__(self, coord: SkyCoord, radius: Angle):
        self.stars = get_stars_in_region(coord, radius)
        return None

    def star(self, gaiaID):
        stars = self.stars
        star = stars[stars['Source'] == gaiaID]
        return star
