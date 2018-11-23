import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval, ImageNormalize
from photutils import datasets
from photutils import DAOStarFinder
from photutils import aperture_photometry, CircularAperture, CircularAnnulus

fits_imgs = glob.glob('final_science_2x2*.fits')

for frame in fits_imgs:
    # load image data
    img, hdr = fits.getdata(frame, header=True)
    filt_name = hdr['FILTER']

    # estimate the background and background noise using sigma clipped statistics
    # pixels that are above or below a specified sigma level from the median are discarded 
    # and the statistics are recalculated
    mean, median, std = sigma_clipped_stats(img, sigma=3, iters=5)

    # subtract background and use DAOStarFinder to locate stars in the image that have FWHMs
    # of around 5 pixels and have peaks approximately 4-sigma above the background
    daofind = DAOStarFinder(fwhm=5, threshold=4.*std)
    sources = daofind(img - median)
    for col in sources.colnames:
        sources[col].info.format = '%.8g'   # for consistent table output

    positions = (sources['xcentroid'], sources['ycentroid'])    # list of source pixel positions

    # create circular aperture and annulus aperture objects, radius in pixel scale
    apertures = CircularAperture(positions, r=5)
    annulus_apertures = CircularAnnulus(positions, r_in=10, r_out=20)

    # perform photometry in both apertures
    apers = [apertures, annulus_apertures]
    phot_table = aperture_photometry(img, apers)
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'    # for consistent table output

    # to calculate the mean local background within the circular annulus aperture, we need to divide
    # its sum by its area, which can be calculated using the area() method
    bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area()

    # the background sum within the circular aperture is then the mean local background times the
    # circular aperture area
    bkg_sum = bkg_mean * apertures.area()
    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    phot_table['residual_aperture_sum'] = final_sum
    phot_table['residual_aperture_sum'].info.format = '%.8g'    # for consistent table output

    # removes sources with negative counts
    rm_idx = []
    for row in phot_table:
        if np.sign(row['residual_aperture_sum']) == -1:
            idx = int(row['id']) - 1
            rm_idx.append(idx)
        else:
            pass
    phot_table.remove_rows(rm_idx)

    # calculate uncertainty in counts
    n_pix = apertures.area()
    n_sky = annulus_apertures.area()
    readnoise = 7
    c0 = n_pix * (1 + (n_pix/n_sky))

    count_err = np.sqrt(phot_table['aperture_sum_0'] + c0 * (phot_table['aperture_sum_1'] + readnoise ** 2))
    phot_table['count_err'] = count_err

    # calculate error in magnitde
    mag_err = (1.0857 / phot_table['aperture_sum_0']) * count_err
    phot_table['mag_inst_err'] = mag_err

    # calculate instrumental magnitudes (using counts)
    counts = phot_table['residual_aperture_sum']
    mag_inst = -2.5 * np.log10(counts)
    phot_table['mag_inst'] = mag_inst

    #print(phot_table)

    # Write table of sources to file
    filename = filt_name + '_sources.txt'

    try:
        os.remove(filename)
    except FileNotFoundError:
        pass
    
    print('Printing sources to .txt file - ' + filename)
    file = open(filename, 'w')
    phot_table.write(file, format='csv')
    file.close()

'''
    # plot the image and mark the location of detected sources
    norm = ImageNormalize(img, interval=ZScaleInterval())
    plt.figure(figsize=(9,12))          # set figure size
    plt.imshow(img, cmap='Greys', origin='lower', norm=norm)
    apertures.plot(color='green', lw=1., alpha=0.5)
    annulus_apertures.plot(color='red', lw=1., alpha=0.5)
    #plt.savefig(path + fname + '_apphot.pdf', format='pdf', bbox_inches='tight', dpi=300)
    plt.show()
'''