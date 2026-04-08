#!/usr/bin/env python

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
from scipy.ndimage import median_filter
import matplotlib.pyplot as plt
import sys
import re
import glob
sys.path.insert(0,os.getenv('PYTHONPATH'))
from catalog2tess_px.catalogs.AsciiCol import AsciiCol
import argparse
sys.path.insert(0,os.getenv('PIPELINE_DIR'))
from tess_time.cut_ffi.cut_data import cut_data, cut_multisector_data
from tess_time.btjd.btjd_correction import btjd_correction

from multiprocessing import Pool


def get_meta_data(ifile,metafile,decimal=False):
    wdir  =  os.path.abspath( os.path.dirname(ifile))
    print(wdir, ifile, metafile)

    #sector number
    try:
        sector_search = re.search('s(\d\d\d\d)',wdir)
        sector = sector_search.group(1)
    except AttributeError:
        sector_search = re.search('sector(\d+)',wdir)
        sector = sector_search.group(1)

    #cam, ccd number
    cam_search = re.search('cam(\d)',wdir)
    cam = cam_search.group(1)
    ccd_search = re.search('ccd(\d)',wdir)
    ccd = ccd_search.group(1)
    
    print('infile',ifile)
    print('sector,cam,ccd',sector,cam,ccd)

    if decimal:
        print('decimal flag passed')
        cat = AsciiCol(metafile, sector, int(cam),sexagesimal=False,ignore_image_buffer=True)
    else:
        print('decimal flag NOT passed')
        cat = AsciiCol(metafile, sector, int(cam),sexagesimal=True,ignore_image_buffer=True)

    rel_ifile = ifile.replace('\\', '/')
    if '/lcGRB/' in rel_ifile:
        rel_ifile = 'lcGRB/' + rel_ifile.split('/lcGRB/', 1)[1]
    else:
        rel_ifile = os.path.basename(rel_ifile)

    base_ifile = os.path.basename(rel_ifile)
    raw_base = base_ifile.replace('_cleaned', '')

    candidates = [
        rel_ifile,
        base_ifile,
        raw_base,
    ]

    print('lookup candidates', candidates)

    m = np.zeros(len(cat.obj_name), dtype=bool)
    for cand in candidates:
        m = np.in1d(cat.obj_name, cand)
        if any(m):
            print('matched', cand)
            break

    print(any(m))
    print(cat.obj_name)
    return {'RA':float(cat.ra[m][0]),
            'DEC':float(cat.dec[m][0])
    }
            

def get_fluxcal_faster( lc_names, fluxes, light_curve_name):
    print("------------->",fluxes, type(fluxes))
    if fluxes.size == 1:
        fluxes = np.array([fluxes])

    #obj_list = []
    #for lc_name in lc_names:
    #    obj_search =re.search('lc_(.*)',os.path.basename(ifile))
    #    obj = obj_search.group(1)
    #    #print('lookup obj',obj)
    #    if obj == 'detrended':
    #        obj = ifile.split('_')[-2]
    #    #print(obj)
    #    obj_list.append(obj)
    #    
    #m= np.in1d(obj_list, object_use)


    m = np.in1d(lc_names, light_curve_name.replace('//','/').split('/',1) )
    print(lc_names, light_curve_name.split('/',1),m)
    print(fluxes[m])
    assert len(fluxes[m]) == 1

    return fluxes[m]


def get_fluxcal(fluxcal_file, light_curve_name):
    lc_names = np.genfromtxt(fluxcal_file,usecols=(0),dtype=str)
    fluxes   = np.genfromtxt(fluxcal_file,usecols=(1) )
    if fluxes.size == 1:
        fluxes = np.array([fluxes])

    #obj_list = []
    #for lc_name in lc_names:
    #    obj_search =re.search('lc_(.*)',os.path.basename(ifile))
    #    obj = obj_search.group(1)
    #    #print('lookup obj',obj)
    #    if obj == 'detrended':
    #        obj = ifile.split('_')[-2]
    #    #print(obj)
    #    obj_list.append(obj)
    #    
    #m= np.in1d(obj_list, object_use)

    m = np.in1d(lc_names, light_curve_name.split('/',1) )
    assert len(fluxes[m]) == 1

    return fluxes[m]


#def clean_lc(ifile, metadata,
#             multisector=False,
#             reference_flux = None):

def clean_lc_parallel(intuple):
    ifile = intuple[0]
    metadata = intuple[1]
    multisector= intuple[2]
    reference_flux = intuple[3]
    print(ifile, metadata, multisector, reference_flux)

    x,y,z,bkg = np.genfromtxt(ifile,unpack=1,usecols=(0,1,2,6))
    

    #remove gaps
    wdir = os.path.abspath(os.path.dirname(ifile))

    if not multisector:
        sector_idx = wdir.find('sector')
        sector     = wdir[sector_idx : sector_idx+8]
        #sector_idx = wdir.find('sector')
        #sector_search = re.search('s(\d\d\d\d)',wdir)
        
        cam_idx = wdir.find('cam')
        cam = wdir[cam_idx : cam_idx+4]
        ccd_idx = wdir.find('ccd')
        ccd = wdir[ccd_idx : ccd_idx+4]
        print(sector, cam, ccd)
        x2,y2,z2, = cut_data(x,y,z, sector, cam, ccd)
        x2,bkg2,z2, = cut_data(x,bkg,z, sector, cam, ccd)

        #correct by 15 minutes to mid exposure
        #convert to BTJD
        #would like to save x2, for look up later
        #print(metadata['RA'],metadata['DEC'])
        if int(sector[-2:]) < 56:
            if int(sector[-2:]) < 27:
                exptime = 30.0/60./24.0
            else:
                exptime = 10.0/60./24.0
        else:
            exptime = 200.0/3600./24.0

    else:
        x2,y2,z2,   = cut_multisector_data(x,y,z, )
        x2,bkg2,z2, = cut_multisector_data(x,bkg,z)
        
        exptime = np.ones(len(x2))
        #roughly the start of s27
        exptime[ x2  < 2036.10 ] = 30.0/60./24.0
        exptime[ (x2 >= 2036.10) & (x2 <=2825.4) ] = 10.0/60./24.0
        exptime[ x2 > 2825.4 ] = 200.0/3600./24.
            
    #load up the filtered background estimate, if it exists
    dstem,dtarget = os.path.split(wdir)
    ifile2 = os.path.join(dstem,'bkg_phot',dtarget,os.path.basename(ifile))
    if os.path.isfile(ifile2):
        x_bkg,y_bkg,z_bkg,bkg_bkg = np.genfromtxt(ifile2, unpack=1,usecols=(0,1,2,6))

        #remove gaps
        x_bkg2,y_bkg2,z_bkg2   = cut_data(x_bkg, y_bkg,z_bkg, sector, cam, ccd)
        x_bkg2,bkg_bkg2,z_bkg2 = cut_data(x_bkg, bkg_bkg,z_bkg, sector, cam, ccd)
        bkg_model = median_filter(y_bkg2, size=48, mode='reflect')
    else:
        y_bkg2 =   np.array([np.nan]*len(x2))
        z_bkg2 =   np.array([np.nan]*len(x2))
        bkg_bkg2 = np.array([np.nan]*len(x2))
        bkg_model = np.array([np.nan]*len(x2))


    x_correct = btjd_correction(x2 + exptime/2.0, metadata['RA'], metadata['DEC'] )

    #adding this Nov 22, 2022. If fluxcal is set, we want magnitudes
    #and counts per second in the light curves
    if reference_flux is not None:
        print('saving {}'.format(ifile+'_cleaned'))

        if reference_flux > 0:
            y2 = reference_flux - y2
        else:
            #this means reference flux is zero or negative
            #essentially, quite small
            y2 = -y2


        inttime = exptime*0.8*0.99*86400
        y2 = y2/inttime
        z2 = z2/inttime
        bkg2 = -bkg2/inttime
        bkg_model = -bkg_model/inttime
        y_bkg2 = -y_bkg2/inttime
        z_bkg2 = z_bkg2/inttime

        mag = np.zeros(len(y2))
        emag = np.zeros(len(y2))        

        #set mag limits based on 3sigma relative to uncertainties
        #negative flux will be treated as an upper limit (3sigma)
        mask = y2 < z2*3
        mag[mask]  = -2.5*np.log10(z2[mask]*3) + 20.44
        emag[mask] = 99.9

        mag[~mask]  = -2.5*np.log10(y2[~mask]) + 20.44
        emag[~mask] = z2[~mask]/y2[~mask]*2.5/np.log(10)

        
        np.savetxt(ifile+'_cleaned',np.c_[x_correct, x2, y2, z2, mag, emag, bkg2,  bkg_model, y_bkg2, z_bkg2],
                   fmt='%15.5f %15.5f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.6f %15.6f',
                   header='reference_flux: {:15.4f}\n{:>13s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s}'.format(
                       reference_flux[0]/inttime,'BTJD','TJD','cts_per_s','e_cts_per_s','mag','e_mag',
                       'bkg','bkg_model', 'bkg2','e_bkg2'))

    else:
        print('saving {}'.format(ifile+'_cleaned'))
        #it seems as though the bkg_bkg is very similar to bkg (as it
        #should be!).  so we won't bother saving it
        np.savetxt(ifile+'_cleaned',np.c_[x_correct, x2, -y2, z2,-bkg2, -bkg_model, -y_bkg2, z_bkg2],
                   fmt='%15.5f %15.5f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f',
                   header='{:>13s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s}'.format('BTJD','TJD','cts','e_cts','bkg','bkg_model', 'bkg2','e_bkg2'))


def get_inputs(args):
    parser = argparse.ArgumentParser( description="Specify TNS light curves to plot---will look up the source in the catalogs and display metadata.")
    parser.add_argument('infiles', nargs='+')
    parser.add_argument('--metafile', help = 'file with ra dec, etc')
#     parser.add_argument('--decimal',action='store_true',help='If set, metafile is assumed to have decimal coordinates')
    parser.add_argument('--multisector',action='store_true',help='If set, assumes light curves over sevral sectors and removes all time ranges saturated in cam4')
    parser.add_argument('--fluxcal', default=None, help='Set to name of file with fluxcalibration.  Expects to find light curve names identical to what is in phot.data')
    return parser.parse_args()
    
def main():
    args = get_inputs(sys.argv[1:])
    tmpra = np.genfromtxt(args.metafile,usecols=(1),dtype=str)
    if tmpra.size==1:
        if ":" in str(tmpra):
            decimal = False
        else:
            decimal = True
    else:
        if ":" in str(tmpra[0]):
            decimal = False
        else:
            decimal = True

    if args.multisector:
        obj_name = np.genfromtxt(args.metafile,usecols=(0),dtype=str)
        ra,dec = np.genfromtxt(args.metafile,usecols=(1,2),unpack=1)


#    p = Pool(150)

    if args.fluxcal is not None:
        lc_names = np.genfromtxt(args.fluxcal,  usecols=(0),dtype=str)
        fluxes   = np.genfromtxt(args.fluxcal,  usecols=(1) )

    metadata = get_meta_data(args.infiles[0], args.metafile,decimal=decimal)

    for ifile in args.infiles:
        print(ifile)
        if '.png' in ifile:
            continue
        if '_cleaned' in ifile:
            continue

        if os.path.isfile(ifile + '_cleaned'):
            continue

        #if not args.multisector:
        #    metadata = get_meta_data(ifile, args.metafile,decimal=decimal)
        #else:
        #
        #
        #    obj_search =re.search('lc_(.*)', os.path.basename(ifile))
        #    obj = obj_search.group(1)
        #    print('lookup obj',obj)
        #    if obj == 'detrended':
        #        obj = ifile.split('_')[-2]
        #    print(obj)
        #    m = np.in1d(obj_name, obj)
        #    print(any(m))
        #    metadata = {'RA':float(ra[m][0]),
        #                'DEC':float(dec[m][0])}
            


        #if args.fluxcal is not None:
        #    reference_flux = get_fluxcal(args.fluxcal, ifile )
        if args.fluxcal is not None:
            reference_flux = get_fluxcal_faster( lc_names, fluxes, ifile)
        else:
            reference_flux = None

        #clean_lc(ifile,metadata,
        #        multisector=args.multisector,
        #        reference_flux = reference_flux)
        #p.apply_async(clean_lc_parallel, (ifile,metadata,
        #                                  args.multisector,
        #                                  reference_flux,), )
        clean_lc_parallel((ifile,metadata,
                           args.multisector,
                           reference_flux,), )

    #p.close()
    #p.join()

if __name__== '__main__':
    main()

#for i in 1 3; do for  j in 1 2 3 4; do cd cam${i}_ccd${j}; pwd; for d in discovery postdiscovery prediscovery; do cp lc_${d}/phot.data . ; mkdir lc; python ~/image_sub/ap_phot.py; mv lc/ref.phot lc_${d}/ref.phot; done ; cd ..; done; done
