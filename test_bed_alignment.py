import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from data_reduction import make_stack, median_comb
import astroalign as aa 
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval, ImageNormalize
from photutils import DAOStarFinder

# method to sort files by filter, returns a 2d array of files sorted by filter
def filter_sort(file_list, filter_list):
    file_filt_list = []

    for n in filter_list:
        same_filt_list = []
        for f in file_list:
            if fits.getheader(f)['FILTER'] == n:
                same_filt_list.append(f)
            else:
                pass

        if len(same_filt_list) == 0:
            pass
        else:
            file_filt_list.append(same_filt_list)
        same_filt_list = []
    
    return file_filt_list

sci_frames = glob.glob('science_2*.fits')   # gets list of all science FITS images
sci_frames.sort()   # sort by filename

filter_names = ['B', 'V', 'R', 'I']

sci_list = filter_sort(sci_frames, filter_names)

sources_list = []

# Find stars in image
for filt_list in sci_list:
    same_filt_tables = []
    for frame in filt_list:

        img, header = fits.getdata(frame, header=True)

        w = wcs.WCS(header)     # parse WCS keywords in primary HDU

        mean, median, std = sigma_clipped_stats(img, sigma=3, iters=4)
        daofind = DAOStarFinder(fwhm=5, threshold=5.*std)
        sources = daofind(img - median)

        print('\r' + 'Finding stars in ' + frame + ' ', end='')

        for col in sources.colnames:
            sources[col].info.format = '%.8g'   # for consistent table output
        #print(sources)

        x_pix = np.array(sources['xcentroid'].data)
        y_pix = np.array(sources['ycentroid'].data)
        pix_coords = np.vstack((x_pix, y_pix)).T

        # Convert pixel coordinates to RA and Dec
        world = w.wcs_pix2world(pix_coords, 1)
        sky_coords = SkyCoord(world, frame='icrs', unit='deg')
        # Add RA and Dec to table
        sources['ra'] = sky_coords.ra
        sources['dec'] = sky_coords.dec
        sources['filename'] = frame

        sources = sources.to_pandas()
        sources = sources.sort_values(by=['flux'], ascending=False)
        sources = sources[:6]  # grabs 10 brightest sources

        trim_srcs = sources[['filename', 'flux', 'xcentroid', 'ycentroid', 'ra', 'dec']]

        same_filt_tables.append(trim_srcs)

    sources_list.append(same_filt_tables)
    same_filt_tables = []
print()

# Check for any misalignment in images

misaligned_frames = []

for filter_stack in sources_list:
    f_name = filter_stack[0]['filename'].iloc[0]
    filter_n = fits.getheader(f_name)['FILTER']

    bad_frames = []

    for i in range(len(filter_stack)):
        #print('New loop')
        source_frame = filter_stack[i]
        source_filename = source_frame['filename'].iloc[0]
        source_x = source_frame['xcentroid'].values
        source_y = source_frame['ycentroid'].values

        check_idx = np.arange(0, len(filter_stack))
        check_idx = np.delete(check_idx, i)

        # Check for misalignment by comparing pixel positions between frames
        for j in check_idx:

            target_frame = filter_stack[j]
            filename = target_frame['filename'].iloc[0]

            #print(source_filename, filename)
            target_x = target_frame['xcentroid'].values
            target_y = target_frame['ycentroid'].values

            diff_x = np.asarray(np.abs(source_x - target_x))
            diff_y = np.asarray(np.abs(source_y - target_y))

            if np.all(diff_x > 1) or np.all(diff_y > 1):
                #print('MISALIGNED FRAME - ' + filename)
                bad_frames.append(filename)
    
    # Get the misaligned frame
    if len(bad_frames) != 0:
        print('Bad frame found in ' + filter_n + ' band')
        mode = Counter(bad_frames).most_common(1)
        bad_frame = mode[0][0]
        misaligned_frames.append(bad_frame)
        print('Misaligned frame - ' + bad_frame)

# Correct misaligned image
for filter_stack in sources_list:
    
    temp_stack = []    # create placeholder stack for good images
    bad_img = []
    good_sources_x = []    # store pixel locations of sources from good images
    good_sources_y = []
    good_sources_ra = []    # store ra and dec locations of sources from good images
    good_sources_dec = []
    
    final_img_stk = []
    
    f_name = filter_stack[0]['filename'].iloc[0]
    filter_n = fits.getheader(f_name)['FILTER']
    #print(filter_n + ' band')
    
    bad_count = 0
    for frame in filter_stack:
            
        filename = frame['filename'].iloc[0]
        #print(filename)
        src_x = frame['xcentroid'].values
        src_y = frame['ycentroid'].values
        
        src_ra = frame['ra'].values
        src_dec = frame['dec'].values
        
        if filename in misaligned_frames:
            bad_count += 1
            bad_img.append(frame)
        else:
            temp_stack.append(filename)
            
            #src_pos = np.vstack((src_x, src_y)).T
            #print(src_pos)
            good_sources_x.append(np.round(src_x,1))
            good_sources_y.append(np.round(src_y,1))
            good_sources_ra.append(src_ra)
            good_sources_dec.append(src_dec)
    
    #print(temp_stack)
    #print(good_sources_x)
    if bad_count != 0:
        mean_src_x = np.median(good_sources_x, axis=0)
        mean_src_y = np.median(good_sources_y, axis=0)
        
        mean_src_ra = np.median(good_sources_ra, axis=0)
        mean_src_dec = np.median(good_sources_dec, axis=0)
        
        # Get median of good source positions
        good_src_pos = np.vstack((mean_src_x, mean_src_y)).T
        
        # Median combine good frames to create a reference to align with
        good_frame_stk = make_stack(temp_stack)
        ref_img = median_comb(good_frame_stk)
        
        # align bad frames
        for frame in bad_img:
            src_x = np.round(frame['xcentroid'].values,1)
            src_y = np.round(frame['ycentroid'].values,1)
            
            src_ra = frame['ra'].values
            src_dec = frame['dec'].values
            
            filename = frame['filename'].iloc[0]
            
            # Get bad source positions
            bad_src_pos = np.vstack((src_x, src_y)).T
            
            # Get difference in WCS positions between sources
            # from two frames to determine which ones will be
            # used for transformation
            diff_ra = np.abs(mean_src_ra - src_ra)
            diff_dec = np.abs(mean_src_dec - src_dec)
            
            src_rm_idx = np.where(diff_ra > 0.001)
            
            # Remove sources that don't appear in both frames
            if len(src_rm_idx) > 0:
                bad_src_pos = np.delete(bad_src_pos, src_rm_idx, axis=0)
                good_src_pos = np.delete(good_src_pos, src_rm_idx, axis=0)
            
            # Estimate transformation required to align image
            print('Correcting misaligned image ' + filename)
            tform = aa.estimate_transform('affine', bad_src_pos, good_src_pos)
            
            # Get misaligned frame data and applies transformation
            bad_img = fits.getdata(filename)
            aligned_img = aa.apply_transform(tform, bad_img, ref_img)

            good_frame_stk.append(aligned_img)
            
        
        final_img_stk = good_frame_stk
    else:
        final_img_stk = make_stack(temp_stack)
        
    final_img = np.sum(final_img_stk, axis=0)
    final_hdr = fits.getheader(temp_stack[0])

    filename = 'final_' + temp_stack[0][:11] + '_' + filter_n + '.fits'

    try:
        os.remove(filename)
    except FileNotFoundError:
        pass

    fits.writeto(filename, final_img, final_hdr)
    print('Creating science image in ' + str(filter_n) + ' band: ' + filename)

    '''
    img_toplot = final_img
        
    # creates ZScale interval as used by DS9
    interval = ZScaleInterval()
    interval.get_limits(img_toplot)
    norm = ImageNormalize(img_toplot, interval=ZScaleInterval())
    plt.figure(figsize=(9,12))          # set figure size
    #plt.scatter(src_x, src_y)
    plt.imshow(img_toplot, cmap='Greys', origin='lower', norm=norm)
    plt.show()
    ''' 