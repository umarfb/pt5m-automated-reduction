from astropy.io import fits
import numpy as np
import os

'''
Method to rebin FITS image by chose scale

Returns:
    new_data : rebinned FITS image file
'''
def rebin(file, scale):
    data = fits.getdata(file)
    hdr = fits.getheader(file)

    f_type = str(hdr['TYPE']).lower()
    binsize = hdr['XBINNING']

    bin_dim = str(scale) + 'x' + str(scale)
    filename = f_type + '_' + bin_dim + '_' + file[5:13]

    try:
        tf_binsize = fits.getheader(filename)['XBINNING']
        if binsize == scale or tf_binsize == scale:
            print('\r' + file + ' already pixel binned to ' + str(scale) + 'x' + str(scale), end='')
        else:
            oldDim = [data.shape[0], data.shape[1]]                   # gets dimensions of image to be rebinned
            newDim = [data.shape[0]//scale, data.shape[1]//scale]     # sets dimensions of new rebinned image
    
            # rebins old image by summing underlying pixels
            new_data = data.reshape([newDim[0], oldDim[0]//newDim[0], newDim[1], oldDim[1]//newDim[1]]).sum(3).sum(1)

            hdr['NAXIS1'] = hdr['NAXIS1']//scale
            hdr['NAXIS2'] = hdr['NAXIS2']//scale
            hdr['XBINNING'] = scale
            hdr['YBINNING'] = scale
            hdr['XPIXSZ'] = hdr['XPIXSZ'] * scale
            hdr['YPIXSZ'] = hdr['YPIXSZ'] * scale
            hdr['CRPIX1'] = hdr['CRPIX1']//scale
            hdr['CRPIX2'] = hdr['CRPIX2']//scale
            hdr['HISTORY'] = 'Rebinned to ' + str(scale) + 'x' + str(scale)

            fits.writeto(filename, new_data, hdr)      # writes rebinned data to FITS file
            print('\r' + 'Creating rebinned images ', end='')
    except OSError:
        if binsize == scale:
            print('\r' + file + ' already pixel binned to ' + str(scale) + 'x' + str(scale), end='')
        else:
            oldDim = [data.shape[0], data.shape[1]]                   # gets dimensions of image to be rebinned
            newDim = [data.shape[0]//scale, data.shape[1]//scale]     # sets dimensions of new rebinned image
    
            # rebins old image by summing underlying pixels
            new_data = data.reshape([newDim[0], oldDim[0]//newDim[0], newDim[1], oldDim[1]//newDim[1]]).sum(3).sum(1)

            hdr['NAXIS1'] = hdr['NAXIS1']//scale
            hdr['NAXIS2'] = hdr['NAXIS2']//scale
            hdr['XBINNING'] = scale
            hdr['YBINNING'] = scale
            hdr['XPIXSZ'] = hdr['XPIXSZ'] * scale
            hdr['YPIXSZ'] = hdr['YPIXSZ'] * scale
            hdr['CRPIX1'] = hdr['CRPIX1']//scale
            hdr['CRPIX2'] = hdr['CRPIX2']//scale
            hdr['HISTORY'] = 'Rebinned to ' + str(scale) + 'x' + str(scale)

            fits.writeto(filename, new_data, hdr)      # writes rebinned data to FITS file
            print('\r' + 'Creating rebinned images ', end='')

'''
Method to create frame stack of images with same pixel binning, with option to scale by median

Returns:
    data_stack : stack of FITS image data
'''
def make_stack(file_list, **kwargs):
    data_stack = []
    scale_median = kwargs.get('scale')

    for f in file_list:
        data = fits.getdata(f)

        if scale_median == True:
            data = data / np.median(data)
        else:
            pass

        data_stack.append(data)
    return data_stack

'''
Method to median combine frame stack, with option to normalise by mean

Returns:
    median_data : median combined FITS image of a stack
'''
def median_comb(file_stack, **kwargs):
    normalise = kwargs.get('norm')

    if normalise == True:
        m = np.mean(file_stack)
        file_stack = file_stack / m
    else:
        pass

    median_data = np.median(file_stack, axis=0)
    return median_data

'''
Method to scale each image by the median, and combine into a stack

Returns:
    data_stack : stack of median scaled FITS images
'''
# method to scale image data by median and combine into stack
def scale_stack(file_list, pixelbin):
    data_stack = []
    for f in file_list:
        binpix = fits.getheader(f)['XBINNING']     # identifies pixel binning in the frame
        if binpix == pixelbin:
            data = fits.getdata(f)
            data = data / np.median(data)
            data_stack.append(data)
    return data_stack

'''
Method to normalise FITS image by the mean

Returns:
    norm_data : mean normalised FITS image
'''
def mean_norm(file):
    mean = np.mean(file)
    norm_data = file/mean
    return norm_data

'''
Method to subtract by master frame

Returns:
    data_out : subtracted frame
'''
def ms_sub(ms_frame, frame):
    ms_data = fits.getdata(ms_frame)
    frame_data = fits.getdata(frame)
    data_out = frame_data - ms_data                    
    return data_out

'''
Method to subtract master dark frame

Returns:
    data_out : subtracted frame

'''
def dark_subtract(dark_frame, frame, scale):
    dark_data = fits.getdata(dark_frame)
    f_data = fits.getdata(frame)
    data_out = f_data - (dark_data//scale)                           
    return data_out