import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy import wcs
from astropy import units as u
from astropy.io import fits
from astropy.io import ascii
from astropy.stats import SigmaClip
from astropy.coordinates import SkyCoord
from astropy.visualization import ZScaleInterval, ImageNormalize
from ucac4 import extract_star_array as get_local_ucac4

# method to convert pixel coordinates to RA and Dec
'''
    Parameters:
    - file_name: file containing table of sources detected on image
    - fits_image: fits image

    Returns:
    - sources: table updated with RA and Dec of sources
'''
def pix2wcs(file_name, fits_image):
    # read in list of image sources into astropy table
    sources = ascii.read(file_name)

    # get pixel coordinates from pt5m image
    x_pix = np.array(sources['xcenter'].data)
    y_pix = np.array(sources['ycenter'].data)
    pix_coords = np.vstack((x_pix, y_pix)).T

    hdulist = fits.open(fits_image)   # Load FITS hdulist
    w = wcs.WCS(hdulist[0].header)  # parse WCS keywords in primary HDU

    # Convert pixel coordinates to RA and Dec
    world = w.wcs_pix2world(pix_coords, 1)
    sky_coords = SkyCoord(world, frame='icrs', unit='deg')

    sources['ra'] = sky_coords.ra
    sources['dec'] = sky_coords.dec
    return sources

# method to get centre of image
'''
    Returns RA and Dec of image centre
'''
def get_cr(fits_image):
    hdulist = fits.open(fits_image)   # Load FITS hdulist
    w = wcs.WCS(hdulist[0].header)  # parse WCS keywords in primary HDU

    cr_ra = hdulist[0].header['CRVAL1']
    cr_dec = hdulist[0].header['CRVAL2']
    return cr_ra, cr_dec

# method to get dimensions of image
'''
    Returns width and height of image in degrees, and 1 pixel separation in degrees
'''
def get_dims(fits_image):
    hdulist = fits.open(fits_image)   # Load FITS hdulist
    w = wcs.WCS(hdulist[0].header)  # parse WCS keywords in primary HDU

    x_range = hdulist[0].header['NAXIS1']
    y_range = hdulist[0].header['NAXIS2']
    # convert pixels to RA and Dec
    lim_coords = w.wcs_pix2world([[0,0], [x_range,y_range]], 1)
    lim_coords = SkyCoord(lim_coords, frame='icrs', unit='deg')
    width = (lim_coords.ra.max() - lim_coords.ra.min())
    height = (lim_coords.dec.max() - lim_coords.dec.min())
    pix_sep = ((width/x_range) + (height/y_range)) / 2
    return width, height, pix_sep

# method to calculate error in the mean
def std_error_mean(data, err_data):
    n = len(data)
    sq_sum = 0

    for i in range(n):
        sq_sum += err_data[i] ** 2
    
    error = np.sqrt(sq_sum)/n
    return error

# method to determine zero point for magnitude calibration using UCAC4 catalog
'''
   Parameters:
   - sources : table file for sources detected in image
   - fits_image: fits image of science frame
   - filter_name : filter used in image
   - radius : separation threshold to match images stars to catalog stars (in pixels)
   - mag : column header for instrumental magnitudes

   Returns:
   - zp: median clipped zero point for image magnitude calibration
   - zp_err: uncertainty in zero point
   - n_total: number of matched sources used
   - n_rej: number of rejected sources
   - delta_m = array of difference between catalog and instrumental magnitudes
   - ucac4_mag = array of matched catalog star magnitudes
'''
def find_zp(sources, fits_image, filter_name, radius, mag):
    sources = pix2wcs(sources, fits_image)
    # get RA and Dec for sources
    src_coords = SkyCoord(sources['ra'], sources['dec'], unit=u.deg)

    ra, dec = get_cr(fits_image)    # get centre of image
    width, height, pix_sep = get_dims(fits_image)   # get dimensions and unit separation on image

    min_sep = radius * pix_sep

    ucac4 = get_local_ucac4(ra, dec, width.deg, height.deg)
    good_mags = (ucac4[filter_name] < 18)
    ucac4 = ucac4[good_mags]
    ucac4_coords = SkyCoord(ucac4['ra'], ucac4['dec'], unit='deg')  # gets RA and Dec of stars from catalog
    
    idx, sep, _ = src_coords.match_to_catalog_sky(ucac4_coords)
    mask = sep < min_sep
 
    delta_m = ucac4[filter_name][idx] - sources[mag]
    delta_m = delta_m[mask]
    ntotal = len(delta_m)   
    ucac4_mag = ucac4[filter_name][idx]
    ucac4_mag = ucac4_mag[mask]

    filter_name_error = 'e_' + filter_name  # get magnitude error of catalog stars
    e_mag = ucac4[filter_name_error][idx]
    e_mag = e_mag[mask]

    # sigma clip data
    sc = SigmaClip(sigma=3, iters=4)
    filtered_delta_m = sc(delta_m)
    mask = ~filtered_delta_m.mask

    nrej = ntotal - np.sum(mask)    # number of rejected sources
    zp = np.ma.median(filtered_delta_m)     # sigma-clipped zero point

    zp_err = std_error_mean(ucac4_mag, e_mag)   # zero point error

    #zp_err = filtered_delta_m.std()/np.sqrt(np.sum(mask))   # zero point error

    return zp, zp_err, ntotal, nrej, delta_m, ucac4_mag

img_files = glob.glob('final_science*.fits')
source_files = glob.glob('*_sources.txt')

for fitsname in img_files:
    header = fits.getheader(fitsname)
    filter_name = header['FILTER']
    filter_label = filter_name
    filename = filter_name + '_sources.txt'

    if filter_name.upper() in ('B', 'V'):
        filter_name = 'mag_' + filter_name.upper()
    elif filter_name.lower() in ('r', 'i'):
        filter_name = 'mag_' + filter_name.lower()
    else:
        raise ValueError('cannot calibrate filter {} against UCAC4'.format(filter_name))

    zp, zp_err, ntot, nrej, delta_m, ucac4_mag = find_zp(filename, fitsname, filter_name, 4, 'mag_inst')

    print('zero point in ' + filter_label.upper() + ' band is: ' + str(np.round(zp, 2)) + ' +/-' + str(np.round(zp_err, 2)))

    # calculate object magnitudes
    sources = ascii.read(filename)
    mag_obj = sources['mag_inst'] + zp
    sources['mag_obj'] = mag_obj

    # calculate error in object magnitudes
    mag_inst_err = sources['mag_inst_err']
    mag_obj_err = np.sqrt((mag_inst_err ** 2) + (zp_err ** 2))
    sources['mag_obj_err'] = mag_obj_err

    # write magnitudes to file
    file = open(filename, 'w')
    sources.write(file, format='csv')
    file.close()

    plt.scatter(mag_obj, mag_obj_err, marker='.', color='#002db3')
    plt.xlabel('$m_{obj}$')
    plt.ylabel('$\sigma_m$')
    plt.title(filter_label.upper() + ' band')
    plt.show()