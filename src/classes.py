# classes.py
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table, vstack
from .functions import get_stars_in_region, fs4giesler
import contextlib
import os
from astropy import units as u
import configparser
import time

class Traceback:
    """
    Example usage
    >>> trace = Trace(195633325090480896, 195633320791780608)
    >>> trace.create_input_file()
    """
    def __init__(self, sourceID1, sourceID2) -> None:
        self.sourceID1 = sourceID1
        self.sourceID2 = sourceID2
        self.stars = (Star(sourceID1), Star(sourceID2))
        self.path = "giessler_traceback/"
        pass
    def create_input_file(self, path='two_trace'):
        two_stars = vstack([self.stars[0].info,
                            self.stars[1].info])
        g = fs4giesler(two_stars)
        output_file = os.path.join(self.path, path,
                                   "input.tsv"
                                   )
        
        g.write(output_file, format='csv', delimiter='\t', overwrite=True)
        print(f"Input file created at {output_file}")

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
class SimulationConfig:
    def __init__(self, config_file):
        self.config_file = config_file
        self.config = configparser.ConfigParser(allow_no_value=True, delimiters=("=",))
        self.config.optionxform = str  # Preserve case sensitivity of keys
        self.load_config()

    def load_config(self):
        """Loads the configuration file into the config parser."""
        with open(self.config_file, 'r') as file:
            self.config.read_file(file)

    def save_config(self, output_file=None):
        """Saves the modified configuration back to the file."""
        if output_file is None:
            output_file = self.config_file
        with open(output_file, 'w') as file:
            self.config.write(file, space_around_delimiters=False)

    def __getattr__(self, name):
        """Dynamically get attribute values from the [Simulation] section."""
        section = 'Simulation'
        if section in self.config:
            # Compare in uppercase to handle case-insensitive matching
            for key in self.config[section]:
                if key.upper() == name.upper():
                    return self.config[section][key]
        raise AttributeError(f"'SimulationConfig' object has no attribute '{name}'")

    def __setattr__(self, name, value):
        """Dynamically set attribute values in the [Simulation] section."""
        if name in ['config_file', 'config']:
            super().__setattr__(name, value)  # Allow normal setting of non-config attributes
        else:
            section = 'Simulation'
            if section in self.config:
                for key in self.config[section]:
                    if key.upper() == name.upper():
                        self.config[section][key] = str(value)
                        return
            raise AttributeError(f"'SimulationConfig' object has no attribute '{name}'")