from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table, join, vstack
from astropy.time import Time
from astropy import units as u
from astropy.utils.metadata import MergeConflictWarning
from astroquery.vizier import Vizier
import numpy as np
import time
import warnings

def get_catalog(name):
    '''
    ### example:
    green = get_catalog("VII/272/snrs")
    '''
    return Vizier(columns=["*", "+_r"], row_limit=-1).get_catalogs(name)[0]

def searchDR3(skycoordinate: SkyCoord, radius: Angle) -> Table:
    start_time = time.time()
    print(f"Searching in Gaia DR3 I/355/gaiadr3 {radius} around {skycoordinate.ra, skycoordinate.dec, skycoordinate.distance}")
    filters = {}
    stars_fromDR3 = Vizier(columns=["*", "+_r"], row_limit=-1).query_region(
        skycoordinate, 
        radius=radius, 
        catalog="I/355/gaiadr3",
        column_filters=filters
    )[0]  
    stars_fromDR3['SkyCoord1'] = SkyCoord(
        ra=stars_fromDR3['RA_ICRS'],
        dec=stars_fromDR3['DE_ICRS'],
        distance=(stars_fromDR3['Plx']).to(u.pc, equivalencies=u.parallax()),
        pm_ra_cosdec=stars_fromDR3['pmRA'],
        pm_dec=stars_fromDR3['pmDE'],
        obstime=(Time('J2000') + 1 * u.Myr)
    )
    end_time = time.time()
    print(f"found {len(stars_fromDR3):,} sources in {end_time - start_time:.2f} seconds")
    return stars_fromDR3

def searchDR3_dist(skycoordinate: SkyCoord, radius: Angle) -> Table:
    start_time = time.time()
    print(f"Searching in Gaia DR3 distances I/352/gedr3dis {radius} around {skycoordinate.ra, skycoordinate.dec, skycoordinate.distance}")
    stars_fromDR3_dist = Vizier(columns=["*", "+_r"], row_limit=-1).query_region(
        skycoordinate, 
        radius=radius, 
        catalog="I/352/gedr3dis"
    )[0]
    stars_fromDR3_dist['SkyCoord2'] = SkyCoord(
        ra=stars_fromDR3_dist['RA_ICRS'],
        dec=stars_fromDR3_dist['DE_ICRS'],
        distance=(stars_fromDR3_dist['rgeo']),
        obstime=(Time('J2000') + 1 * u.Myr)
    )
    end_time = time.time()
    print(f"found {len(stars_fromDR3_dist):,} sources in {end_time - start_time:.2f} seconds")
    return stars_fromDR3_dist

def merge_gaia_tables(stars_fromDR3: Table, stars_fromDR3_dist: Table) -> Table:
    start_time = time.time()
    print("Starting merge of DR3 and distance catalog data")

    warnings.simplefilter('ignore', MergeConflictWarning)
    merged = join(stars_fromDR3, stars_fromDR3_dist, keys='Source', join_type='inner')

    # Order the columns
    merged = merged['RA_ICRS_1', 'DE_ICRS_1', 'e_RA_ICRS', 'e_DE_ICRS', '_r_1',
                                       'HIP', 'TYC2', 'Source', 'rgeo', 'Plx', 'e_Plx',
                                       'pmRA', 'pmDE', 'e_pmRA', 'e_pmDE',
                                       'RUWE', 'Teff', 'logg', 'Gmag', 'BP-RP', 'BPmag', 'RPmag', 'RV', 'e_RV',
                                       'b_rgeo', 'B_rgeo', 'FG', 'e_FG', 'FBP', 'e_FBP', 'FRP', 'e_FRP', 'RAVE5', 'RAVE6']

    sigmaG_0 = 0.0027553202
    sigmaGBP_0 = 0.0027901700
    sigmaGRP_0 = 0.0037793818

    merged['e_Gmag'] = np.sqrt((-2.5 / np.log(10) * merged['e_FG'] / merged['FG'])**2 + sigmaG_0**2)
    merged['e_BPmag'] = np.sqrt((-2.5 / np.log(10) * merged['e_FBP'] / merged['FBP'])**2 + sigmaGBP_0**2)
    merged['e_RPmag'] = np.sqrt((-2.5 / np.log(10) * merged['e_FRP'] / merged['FRP'])**2 + sigmaGRP_0**2)
    merged['e_BP-RP'] = merged['e_BPmag'] + merged['e_RPmag']

    merged['SkyCoord'] = SkyCoord(
        ra=merged['RA_ICRS_1'],
        dec=merged['DE_ICRS_1'],
        distance=(merged['rgeo']),
        pm_ra_cosdec=merged['pmRA'],
        pm_dec=merged['pmDE'],
        obstime=(Time('J2000') + 1 * u.Myr)
    )

    end_time = time.time()
    print(f"{len(merged):,} sources found by merging in {end_time - start_time:.2f} seconds")
    return merged

def get_stars_in_region(skycoordinate: SkyCoord, radius: Angle) -> Table:
    c = skycoordinate
    t1 = searchDR3(c, radius)
    t2 = searchDR3_dist(c, radius)
    t3 = merge_gaia_tables(t1, t2)
    t3.sort('Gmag')
    return t3

def fs4giesler(sir, save=None):
        table = sir
        g = Table()
        g['TypeInput'] = np.ones_like(table['e_Plx'].value).astype(int)
        g['RA'] = table['SkyCoord'].ra.to_string(unit='hourangle', sep=' ', precision=3, pad=True)
        g['DE'] = table['SkyCoord'].dec.to_string(unit='degree', sep=' ', precision=3, pad=True)
        g['Plx'] = table['rgeo'].to(u.mas, u.parallax())
        g['e_Plx'] = table['e_Plx']
        g['RV'] = np.zeros_like(table['e_Plx'].value).astype(int)
        g['e_RV'] = np.zeros_like(table['e_Plx'].value).astype(int)
        g['RVdist'] = np.zeros_like(table['e_Plx'].value).astype(int)
        g['pmRA'] = table['pmRA']
        g['e_pmRA'] = table['e_pmRA']
        g['pmDE'] = table['pmDE']
        g['e_pmDE'] = table['e_pmDE']
        g['Source'] = table['Source'].astype(str)

        new_table = Table(names=g.colnames, dtype=[col.dtype for col in g.columns.values()])
        g = vstack([new_table,g])
        if save:
             g.write("giesler_input.tsv", format='csv', delimiter='\t', overwrite=True)
        return g

def read_traceback_out(matches_file):
    column_names = [
        'min_sep',        # Column 1: minimum separation between star 1 and 2 [pc]
        'time',           # Column 2: time since minimum separation [yr]
        'col3',           # Column 3: always 0
        'col4',           # Column 4: always 0
        'plx_1',          # Column 5: π [mas] of star 1
        'RV_1',           # Column 6: RV or U [km/s] of star 1
        'mu_alpha_1',     # Column 7: µ_α* [mas/yr] or V [km/s] of star 1
        'mu_delta_1',     # Column 8: µ_δ [mas/yr] or W [km/s] of star 1
        'l_1',            # Column 9: Galactic longitude [degrees] of star 1
        'b_1',            # Column 10: Galactic latitude [degrees] of star 1
        'd_sun_1',        # Column 11: Distance from Sun [pc] of star 1
        'plx_2',          # Column 12: π [mas] of star 2
        'RV_2',           # Column 13: RV or U [km/s] of star 2
        'mu_alpha_2',     # Column 14: µ_α* [mas/yr] or V [km/s] of star 2
        'mu_delta_2',     # Column 15: µ_δ [mas/yr] or W [km/s] of star 2
        'l_2',            # Column 16: Galactic longitude [degrees] of star 2
        'b_2',            # Column 17: Galactic latitude [degrees] of star 2
        'd_sun_2'         # Column 18: Distance from Sun [pc] of star 2
    ]
    # Read the table from the file using the custom column names
    table = Table.read(matches_file, format="ascii", delimiter=r'\s', names=column_names, guess=False)
    # Add units to the columns
    table['min_sep'].unit = u.pc                       # Parsec
    table['time'].unit = u.yr                          # Years
    table['plx_1'].unit = u.mas                        # Milliarcseconds
    table['RV_1'].unit = u.km / u.s                    # Kilometers per second
    table['mu_alpha_1'].unit = u.mas / u.yr            # Milliarcseconds per year
    table['mu_delta_1'].unit = u.mas / u.yr            # Milliarcseconds per year
    table['l_1'].unit = u.deg                          # Degrees
    table['b_1'].unit = u.deg                          # Degrees
    table['d_sun_1'].unit = u.pc                       # Parsec
    table['plx_2'].unit = u.mas                        # Milliarcseconds
    table['RV_2'].unit = u.km / u.s                    # Kilometers per second
    table['mu_alpha_2'].unit = u.mas / u.yr            # Milliarcseconds per year
    table['mu_delta_2'].unit = u.mas / u.yr            # Milliarcseconds per year
    table['l_2'].unit = u.deg                          # Degrees
    table['b_2'].unit = u.deg                          # Degrees
    table['d_sun_2'].unit = u.pc                       # Parsec
    # Now you can access the table and the columns will have the correct units
    return table

def best_match(traceback_out):
    idx = np.argmin(traceback_out['min_sep'])
    best_case = traceback_out[idx]
    return best_case

def read_trajectory_out(matches_file):
    table = Table.read(matches_file, format="ascii", delimiter=r'\s')
    # Add units to the columns
    table.rename_column('col1','time')
    table.rename_column('col2','l_1')
    table.rename_column('col3','b_1')
    table.rename_column('col4','dist_1')
    table.rename_column('col5','l_2')
    table.rename_column('col6','b_2')
    table.rename_column('col7','dist_2')
    table.rename_column('col8','sep3d')
    table['time'].unit = u.yr
    table['l_1'].unit = u.deg
    table['b_1'].unit = u.deg
    table['dist_1'].unit = u.pc
    table['l_2'].unit = u.deg
    table['b_2'].unit = u.deg
    table['dist_2'].unit = u.pc
    table['sep3d'].unit = u.pc
    return table

def best_params(matches_file):
    traceback_out = read_traceback_out(matches_file)
    best_case = best_match(traceback_out)
    __skycoords = SkyCoord( l=[best_case['l_1'], best_case['l_2']]*u.deg,
                            b=[best_case['b_1'], best_case['b_2']]*u.deg, 
                            distance=[best_case['d_sun_1'], best_case['d_sun_2']]*u.pc, 
                            frame='galactic')
    _skycoords = __skycoords.transform_to('icrs')
    skycoords = SkyCoord(ra=_skycoords.ra,
                        dec=_skycoords.dec,
                        distance=_skycoords.distance,
                        pm_ra_cosdec=[best_case['mu_alpha_1'], best_case['mu_alpha_2']]*u.mas/u.yr,
                        pm_dec=[best_case['mu_delta_1'], best_case['mu_delta_2']]*u.mas/u.yr,
                        )
    table = Table()
    table['RA'] = skycoords.ra.to_string(unit='hourangle', sep=' ', precision=3, pad=True)
    table['DE'] = skycoords.dec.to_string(unit='degree', sep=' ', precision=3, pad=True)
    table['Plx'] = np.around(skycoords.distance.to(u.mas, u.parallax()).value, 5)
    table['e_Plx'] = np.zeros_like(table['Plx'])
    table.add_column([1,1], index=0, name="TypeInput")
    table['RV'] = np.zeros_like(table['Plx'])
    table['e_RV'] = np.zeros_like(table['Plx'])
    table['RVdist'] = np.zeros_like(table['Plx'])
    table['pmRA'] = np.around(skycoords.pm_ra_cosdec.value, 5)
    table['e_pmRA'] = np.zeros_like(table['Plx'])
    table['pmDE'] = np.around(skycoords.pm_dec.value, 5)
    table['e_pmDE'] = np.zeros_like(table['Plx'])
    table['Source'] = np.zeros_like(table['Plx'])
    table['SkyCoord'] = skycoords
    table.remove_column('SkyCoord')
    return table