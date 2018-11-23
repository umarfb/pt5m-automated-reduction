import glob
import os
from data_reduction import rebin, make_stack, median_comb, ms_sub, dark_subtract
import numpy as np
from astropy.io import fits

# method to sort files by pixel binning, returns a 2d array of files sorted by pixel binning
def bin_sort(file_list, bin_list):
    file_bin_list = []
    for n in bin_list:      # organises bias frames by pixel binning
        same_bin_stack = []
        for f in file_list:
            if fits.getheader(f)['XBINNING'] == n:
                same_bin_stack.append(f)
            else:
                pass

        if len(same_bin_stack) == 0:
            pass
        else:
            file_bin_list.append(same_bin_stack)
        same_bin_stack = []

    return file_bin_list

# method to median combine frames, writes median combined image to file
def make_masterframe(file_list):
    for bin_list in file_list:  
        header = fits.getheader(bin_list[0])
        binsize = header['XBINNING']
        f_type = str(header['TYPE']).lower()
        binsize = str(binsize) + 'x' + str(binsize)
    
        frame_bin_stack = make_stack(bin_list)
        master_frame = median_comb(frame_bin_stack)
        filename = 'master_' + f_type + '_' + binsize + '.fits'
        print('Created master ' + f_type + ' frame for ' + binsize + ' pixel bin: ' + filename)

        try:
            os.remove(filename)
        except FileNotFoundError:
            pass
    
        header['HISTORY'] = 'Median combined'
        fits.writeto(filename, master_frame, header)

# method to carry out bias subtraction routine, writes subtracted image to file
def bias_subtract_routine(file_list, bin_list, mast_frame):
    for n in bin_list:
        for binned_frames in file_list:
            for f in binned_frames:
                header = fits.getheader(f)
                binsize = header['XBINNING']
                f_type = str(header['TYPE']).lower()
                bin_dim = str(binsize) + 'x' + str(binsize)
                master_frame = mast_frame + bin_dim + '.fits'

                '''
                if f == 'r*.fits':
                    filename = f_type + '_' + bin_dim + '_' + f[5:13]
                    print(filename)
                else:
                    filename = f
                    #print(filename)
                '''

                if f[0] == 'r':
                    filename = f_type + '_' + bin_dim + '_' + f[5:13]
                else:
                    filename = f

                if binsize == n:
                    bias_sub = ms_sub(master_frame, f)
                    header['HISTORY'] = 'Bias subtracted'
                    print('\r' + 'Bias subtracting ' + f_type + ' frames ', end='')
                    try:
                        os.remove(filename)
                    except FileNotFoundError:
                        pass
                    fits.writeto(filename, bias_sub, header)
                else:
                    pass
    print()

# method to carry out dark subtraction, writes subtracted image to file
def dark_subtract_routine(file_list, bin_list, mast_frame):
    for n in bin_list:
        for binned_frames in file_list:
            for f in binned_frames:
                header = fits.getheader(f)
                binsize = header['XBINNING']
                f_type = str(header['TYPE']).lower()
                bin_dim = str(binsize) + 'x' + str(binsize)
                master_frame = mast_frame + bin_dim + '.fits'
                tf_exp = header['EXPTIME']       # get exposure time in seconds for target frame
                mf_exp = fits.getheader(master_frame)['EXPTIME']

                scale = mf_exp / tf_exp     # gets scaling factor of exposure time

                if binsize == n:
                    dark_sub = dark_subtract(master_frame, f, scale)
                    header['HISTORY'] = 'Dark subtracted'
                    print('\r' + 'Dark subtracting ' + f_type + ' frames ', end='')
                    try:
                        os.remove(f)
                    except FileNotFoundError:
                        pass
                    fits.writeto(f, dark_sub, header)
                else:
                    pass
    print()

files = glob.glob('r*.fits')    # get list of all raw FITS images
files = np.asarray(files)
files.sort()                    # sort by filename

type_name  = []
filter_name = []
pix_bin_size = []

bias_frames = []
dark_frames = []
flat_fields = []
science_frames = []

# get start date of observations
date_start = str(fits.getheader(files[0])['DATE-OBS'])
print('Observations for the night of ' + date_start[0:10])

# sort FITS images by TYPE (BIAS, DARK, SCIENCE, SKY)
for f in files:
    filter = fits.getheader(f)['FILTER']
    f_type = fits.getheader(f)['TYPE']
    pixel_bin = fits.getheader(f)['XBINNING']

    type_name.append(f_type)
    filter_name.append(filter)
    pix_bin_size.append(pixel_bin)

    if f_type == 'BIAS':
        bias_frames.append(f)
    elif f_type == 'DARK':
        dark_frames.append(f)
    elif f_type == 'SKY':
        flat_fields.append(f)
    elif f_type == 'SCIENCE':
        science_frames.append(f)

# get the type of image, filter used, and pixel bins
type_name = np.unique(type_name)
filter_name = np.unique(filter_name)
pix_bin_size = np.unique(pix_bin_size)

print(str(len(files)) + ' files found: ')
print('- ' + str(len(bias_frames)) + ' bias frames')
print('- ' + str(len(dark_frames)) + ' dark frames')
print('- ' + str(len(flat_fields)) + ' flat fields')
print('- ' + str(len(science_frames)) + ' science frames')

# determine binning of science frames
sci_bins = []
for f in science_frames:
    pixel_bin = fits.getheader(f)['XBINNING']
    sci_bins.append(pixel_bin)
sci_bins = np.unique(sci_bins)

# rebin flat fields to match binning in science frames
temp_flats = []

for n in sci_bins:
    for f in flat_fields:
        rebin(f, n)
        f_type = str(fits.getheader(f)['TYPE']).lower()
        bin_dim = str(n) + 'x' + str(n)
        filename = f_type + '_' + bin_dim + '_' + f[5:13]
        temp_flats.append(filename)
    print()

flat_fields = np.concatenate([flat_fields, temp_flats]) # update list of flat fields

# make master bias frames
bias_list = bin_sort(bias_frames, sci_bins)     # sort bias frames by pixel bin
make_masterframe(bias_list)

# do bias subtraction from dark, science, and flat fields
# -----------------------------------------------------------------------------
dark_list = bin_sort(dark_frames, sci_bins)     # sort dark frames by pixel bin
flats_list = bin_sort(flat_fields, sci_bins)    # sort flat fields by pixel bin
science_list = bin_sort(science_frames, sci_bins)   # sort science frames by pixel bin

master_bias = 'master_bias_'    # define master bias frame to be used

bias_subtract_routine(dark_list, sci_bins, master_bias)  # subtract bias from dark frames
bias_subtract_routine(flats_list, sci_bins, master_bias) # subtract bias from flat fields
bias_subtract_routine(science_list, sci_bins, master_bias)   # subtract bias from science frames

# -----------------------------------------------------------------------------

# clear and load bias subtracted files into buffer arrays
# -----------------------------------------------------------------------------
dark_frames = []
flat_fields = []
science_frames = []

dark_frames = np.asarray(glob.glob('dark*.fits'))
flat_fields = np.asarray(glob.glob('sky*.fits'))
science_frames = np.asarray(glob.glob('science*.fits'))

dark_list = bin_sort(dark_frames, sci_bins)     # sort dark frames by pixel bin
flats_list = bin_sort(flat_fields, sci_bins)    # sort flat fields by pixel bin
science_list = bin_sort(science_frames, sci_bins)   # sort science frames by pixel bin
# -----------------------------------------------------------------------------

# make master dark frames
make_masterframe(dark_list)

# do dark subtraction from science frames and flat fields
master_dark = 'master_dark_'
dark_subtract_routine(flats_list, sci_bins, master_dark)
dark_subtract_routine(science_list, sci_bins, master_dark)

# make master flat fields
for bin_list in flats_list:
    binsize = fits.getheader(bin_list[0])['XBINNING']
    bin_dim = str(binsize) + 'x' + str(binsize)
    for filt in filter_name:

        same_filt_stack = []
        for frame in bin_list:  # create stack of frames sorted by filter
            frame_filter = fits.getheader(frame)['FILTER']

            if frame_filter == filt:
                same_filt_stack.append(frame)
        
        header = fits.getheader(same_filt_stack[0])
        # scale each image by median, and median combine and normalise to make master frame
        frame_stack = make_stack(same_filt_stack, scale=True)
        master_frame = median_comb(frame_stack, normalise=True)

        filename = 'master_flat_' + filt + '_' + bin_dim + '.fits'

        try:
            os.remove(filename)
        except FileNotFoundError:
            pass

        print('Created master flat field in ' + filt + ' band for ' + bin_dim + ' pixel bin: ' + filename)
        fits.writeto(filename, master_frame, header)
        same_filt_stack = []

# divide science frames by flat fields, for each filter

for bin_list in science_list:
    binsize = fits.getheader(bin_list[0])['XBINNING']
    bin_dim = str(binsize) + 'x' + str(binsize)

    for frame in bin_list:
        data, header = fits.getdata(frame, header=True)
        frame_filter = header['FILTER']

        master_flat = 'master_flat_' + frame_filter + '_' + bin_dim + '.fits'
        mflat_data = fits.getdata(master_flat)

        data_out = data / mflat_data
        header['HISTORY'] = 'Flat fielded'

        try:
            os.remove(frame)
        except FileNotFoundError:
            pass
        print('\r' + 'Flat fielding science frames for ' + str(frame_filter) + ' band ', end='')
        fits.writeto(frame, data_out, header)
print()