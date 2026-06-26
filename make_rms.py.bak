import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
from astropy.io import fits
import os
import re
import sys
import glob
import shutil


def make_file_dict():
    dates_list = glob.glob('slice*/dates')
    dates_list = np.sort(dates_list)
    fins = {}
    for dates in dates_list:
        files = np.genfromtxt(dates, usecols=(0), dtype=str)
        try:
            for f in files:
                fin = int(f.split('-')[2]) 
                fins[fin] = dates.split('/')[0]
        except TypeError:
            for f in files.reshape(1):
                fin = int(f.split('-')[2]) 
                fins[fin] = dates.split('/')[0]
            
    return fins

        
def make_rms(inlist ):

    if '/o' in os.getcwd():
        lookup_dict = make_file_dict()
    else:
        lookup_dict = None
        
    #assumes files in inlist are in current working directory
    avg = 0
    square = 0
    N = 0
    for infile in inlist:
        #dirlook = 'cam{}_ccd{}'.format(cam,ccd)
        #infile_use = glob.glob( dirlook + '/o1a/slice*' + '/conv_tess2023132091625-{:08d}-{}-crm-ffi_ccd{}.cal.fits'.format(int(infile),cam,ccd))[0]
        #print(infile_use)
        #infile_use = glob.glob(dirlook + '/o1a/tess*' + str(infile) + '*')[0]
        #print(dirlook, str(infile), infile_use)

        if lookup_dict is not None:
            fin = int(infile.split('-')[2])
            slice_use = lookup_dict[fin]

            print(os.getcwd(), slice_use)
            #untar images
            #important to give the destination directory, or all the images land in the wrong place
            evt = shutil.unpack_archive(os.path.join(slice_use, 'images.tar'),
                                        slice_use,)
            d = fits.open( os.path.join(slice_use, 'conv_' + infile)  )[0].data
            #delete images
            for prefix in ['interp_', 'conv_', 'bkg_']:
                flist1 = glob.glob( os.path.join(slice_use, prefix + '*fits'))
                for fl in flist1:
                    os.remove(fl)

            
        else:
            d = fits.open( infile )[0].data
        avg += d
        square += d**2
        N += 1
    rms = np.sqrt( square/N -  (avg/N)**2   )
    fits.writeto('rms.fits', rms,overwrite=True)

if __name__ == "__main__":
    inlist = np.genfromtxt(sys.argv[1],dtype=str)
    make_rms(inlist)
