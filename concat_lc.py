#!/usr/bin/env python

import numpy as np
import glob
import os
import sys
import argparse



def get_inputs(args):
    parser = argparse.ArgumentParser( description="Tool for copying slice light curves into one light curve.  Will assume you run this in the directory where you want the outputs.")

    parser.add_argument('--sector')
    parser.add_argument('--cam')
    parser.add_argument('--ccd')
    parser.add_argument('--lcdir')
    parser.add_argument('--override',action='store_true',help="By defaul, script exits if there is no data in the first"
                        "slice, because this looks like a pipeline error. But, over ride this behavior (run any way)"
                        " if there really isn't any data in the first slice. For example, form scattered light.")
    
    parser.add_argument('--delete', action='store_true', help="Delete light curves from the slice directories. Important for saving inodes.")

    return parser.parse_args()
    

def read_data(indir,lcstem):
    infiles = glob.glob( os.path.join(
        indir,
        lcstem,
        "lc_*"
    ))
    infiles.sort()
    output = {}
    for infile in infiles:
        #print(infile)
        d = np.genfromtxt(infile)
        output[ os.path.basename(infile) ] = d

    infiles = glob.glob( os.path.join(
        indir,
        "bkg_phot/",
        lcstem,
        "lc_*"
        ))
    infiles.sort()
    output_bkg = {}
    for infile in infiles:
        #print(infile)
        d = np.genfromtxt(infile)
        output_bkg[ os.path.basename(infile) ] = d

    return output, output_bkg

def delete_data(indir,lcstem):
    infiles = glob.glob( os.path.join(
        indir,
        lcstem,
        "lc_*"
        ))
    infiles.sort()
    for infile in infiles:

        os.remove(infile)

    infiles = glob.glob( os.path.join(
        indir,
        "bkg_phot",
        lcstem,
        "lc_*"
        ))
    infiles.sort()

    for infile in infiles:
        os.remove(infile)
    return 0


if __name__ == '__main__':
    duse = os.environ['DATA_DIR']
    dhome = os.getcwd()

    args = get_inputs(sys.argv[1:])
    
    sector = args.sector
    sectoruse = 's{:04d}'.format(int(sector))
    print(sectoruse)

    cam = args.cam
    ccd = args.ccd
    lcdir = args.lcdir

    dout = os.path.join(dhome,
                        "sector"+str(sector),
                        "cam"+str(cam)+"_ccd"+str(ccd)
                        )

    dtarget = os.path.join(duse,
                           sectoruse,
                           "cam"+str(cam)+"-ccd"+str(ccd),
                           )

    print(dout)
    if os.path.isdir(os.path.join(dout, lcdir)):
        pass
    else:
        os.makedirs(os.path.join(dout, lcdir))

    if os.path.isdir(os.path.join(dout, 'bkg_phot',lcdir)):
        pass
    else:
        os.makedirs(os.path.join(dout, 'bkg_phot', lcdir))
    

    
    indirs = glob.glob( os.path.join(
        dtarget, "o*", "slice*"
    ))
    indirs.sort()



    results = []
    results_bkg = []
    for indir in indirs:
        #print(indir)
        out1, out2 = read_data(indir,lcdir)
        results.append( out1  )
        results_bkg.append( out2 )

        
    if bool(results[0]) == False:
        #over ride option if there really shouldn't be data in the first slice,
        #for example scattered light.  This happened is s87, cam1, ccd3, o1a.
        if args.override==False:
            print('Error found---first director (o1a/slice0000) is empty')
            sys.exit()


    #what to do for override here?
    for lc in results[0].keys():
        print(lc)
        output = []
        output_bkg = []
    
        for ii,r in enumerate(results):
            #print(ii,r)
            try:
                print(indirs[ii], lc)
                output.append( r[lc] )
                output_bkg.append( results_bkg[ii][lc] )
            except KeyError:
                print('error! no key for',indirs[ii], lc)
                print('passing over this directory')
                #print(os.getcwd())
                #print(lc)
                #raise

        output = np.vstack(output)
        output_bkg = np.vstack(output_bkg)

        idx = np.argsort(output[:,0])
        output = output[idx]
        idx = np.argsort(output_bkg[:,0])
        output_bkg = output_bkg[idx]
        
        output = np.unique(output, axis=0)
        output_bkg = np.unique(output_bkg, axis=0)

    
        output_file = os.path.join(dhome,
                                   "sector"+str(sector),
                                   "cam"+str(cam)+"_ccd"+str(ccd),
                                   lcdir,
                                   os.path.basename(lc)
                                   )

        np.savetxt(output_file, output)

        output_file2 = os.path.join(dhome,
                                    "sector"+str(sector),
                                    "cam"+str(cam)+"_ccd"+str(ccd),
                                    "bkg_phot",
                                    lcdir,
                                    os.path.basename(lc)
                                    )
        np.savetxt(output_file2, output_bkg)
        #check that number of entries in output and output_bkg
        #match what is in the save file
        
        n_output = np.shape(output)[0]
        n_output_bkg = np.shape(output_bkg)[0]
        
        with open(output_file,'r') as fin:
            lines = fin.readlines()
            n_output_saved = len(lines)
        with open(output_file2,'r') as fin:
            lines = fin.readlines()
            n_output_bkg_saved = len(lines)
            
        print('n saved (lc, bkg): ', n_output_saved, n_output_bkg_saved)
        print('n in data/tica/ (lc, bkg): ', n_output, n_output_bkg)
        assert (n_output_saved == n_output)
        assert (n_output_bkg_saved == n_output_bkg)
        assert (n_output_saved == n_output_bkg_saved)
        

    
        if args.delete:

            for indir in indirs:
                print(indir, lcdir)
                delete_data(indir,lcdir)
        
