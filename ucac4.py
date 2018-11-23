from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import os
from astropy import units as u
from astropy.table import Table, Column
from astropy.time import Time
import os

data_loc = os.environ.get('UCAC4_DATA_LOC', '/scratch/ph1jrm/UCAC4/UCAC4')
star_dtype = np.dtype([
    ('ra', 'i4'), ('dec', 'i4'),
    ('model_mag', 'i2'), ('aperture_mag', 'i2'),
    ('mean_mag_sigma', 'u1'), ('obj_type', 'u1'),
    ('double_star_flag', 'u1'), ('e_ra', 'i1'),
    ('e_dec', 'i1'), ('n_ucac_total', 'u1'),
    ('n_ucac_used', 'u1'), ('n_cats_used', 'u1'),
    ('epoch_ra', 'u2'), ('epoch_dec', 'u2'),
    ('pm_ra', 'i2'), ('pm_dec', 'i2'),
    ('e_pm_ra', 'i1'), ('e_pm_dec', 'i1'),
    ('twomass_id', 'u4'), ('mag_j', 'i2'),
    ('mag_h', 'i2'), ('mag_k', 'i2'),
    ('j_flag', 'u1'), ('h_flag', 'u1'),
    ('k_flag', 'u1'), ('e_mag_j', 'u1'),
    ('e_mag_h', 'u1'), ('e_mag_k', 'u1'),
    ('mag_B', 'i2'), ('mag_V', 'i2'), ('mag_g', 'i2'),
    ('mag_r', 'i2'), ('mag_i', 'i2'),
    ('e_mag_B', 'u1'), ('e_mag_V', 'u1'),
    ('e_mag_g', 'u1'), ('e_mag_r', 'u1'),
    ('e_mag_i', 'u1'), ('yale_flag', 'u1'),
    ('catalog_flags', 'u4'), ('leda_flag', 'u1'),
    ('twomass_extended_flag', 'u1'),
    ('mpos1', 'u4'), ('ucac2_zone', 'u2'),
    ('ucac2_number', 'u4')
])


def star_array_from_zone(zone):
    zone_file = os.path.join(data_loc, 'u4b/z{:03d}'.format(zone))
    sarr = np.fromfile(zone_file, dtype=star_dtype)
    return sarr


def extract_star_array(ra, dec, width, height):
    dec1 = dec - height / 2
    dec2 = dec + height / 2
    ra1 = ra - width / 2
    ra2 = ra + width / 2
    zone_height = .2  # zones are .2 degrees each

    start_zone = int((dec1 + 90) / zone_height) + 1
    end_zone = int((dec2 + 90) / zone_height) + 1

    arrays = []
    for zone in range(start_zone, end_zone + 1):

        sarr = star_array_from_zone(zone)

        # ra in microarcseconds
        sra = sarr['ra'] / 3600 / 1000
        mask = (sra >= ra1) & (sra <= ra2)

        ids = [
            '{:03d}-{:06d}'.format(zone, i) for i in np.where(mask)[0]
        ]
        sarr = np.lib.recfunctions.append_fields(
            sarr[mask], 'ucac4', ids, usemask=False
        )
        arrays.append(sarr)

    return convert_star_array(np.hstack(arrays))


def convert_star_array(star_array):
    """
    Rationalise star array into astropy table, with units
    """
    converters = dict(
        ucac4=lambda x: x,
        ra=lambda x: u.deg * x / 3600 / 1000,
        dec=lambda x: u.deg * (-90 + x / 3600 / 1000),
        e_ra=lambda x: u.deg * x / 3600 / 1000,
        e_dec=lambda x: u.deg * x / 3600 / 1000,
        model_mag=lambda x: u.mag * x * 0.001,
        aperture_mag=lambda x: u.mag * x * 0.001,
        mean_mag_sigma=lambda x: u.mag * x * 0.001,
        obj_type=lambda x: x,
        double_star_flag=lambda x: x,
        epoch_ra=lambda x: Time(1900 + x/100, format='byear'),
        epoch_dec=lambda x: Time(1900 + x/100, format='byear'),
        n_ucac_total=lambda x: x,
        n_ucac_used=lambda x: x,
        n_cats_used=lambda x: x,
        pm_ra=lambda x: u.mas * 0.1 * x / u.yr,
        pm_dec=lambda x: u.mas * 0.1 * x / u.yr,
        e_pm_ra=lambda x: u.mas * 0.1 * x / u.yr,
        e_pm_dec=lambda x: u.mas * 0.1 * x / u.yr,
        mag_j=lambda x: x * u.mag * 0.001,
        e_mag_j=lambda x: x * u.mag * 0.001,
        j_flag=lambda x: x,
        mag_h=lambda x: x * u.mag * 0.001,
        e_mag_h=lambda x: x * u.mag * 0.001,
        h_flag=lambda x: x,
        mag_k=lambda x: x * u.mag * 0.001,
        e_mag_k=lambda x: x * u.mag * 0.001,
        k_flag=lambda x: x,
        mag_B=lambda x: x * u.mag * 0.001,
        e_mag_B=lambda x: x * u.mag * 0.001,
        mag_V=lambda x: x * u.mag * 0.001,
        e_mag_V=lambda x: x * u.mag * 0.001,
        mag_g=lambda x: x * u.mag * 0.001,
        e_mag_g=lambda x: x * u.mag * 0.001,
        mag_r=lambda x: x * u.mag * 0.001,
        e_mag_r=lambda x: x * u.mag * 0.001,
        mag_i=lambda x: x * u.mag * 0.001,
        e_mag_i=lambda x: x * u.mag * 0.001,
        twomass_extended_flag=lambda x: x,
        leda_flag=lambda x: x
    )

    columns = []
    for k, f in converters.items():
        columns.append(Column(f(star_array[k]), k))
    return Table(columns)