import matplotlib as mpl
from io import StringIO
mpl.use('Agg')
import numpy as np
import warnings
import pandas as pd
import traceback
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit
import argparse
from multiprocessing import Pool
from astropy.time import Time
import astropy.units as astrou
from scipy.stats import chi2, mode
from scipy.stats import norm
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import glob
from tqdm import tqdm
from tess_stars2px import tess_stars2px_reverse_function_entry
import pyvo as vo
from astroquery.mast import Catalogs
from astroquery.vizier import Vizier
import matplotlib
import warnings

font = {'weight' : 'bold',
        'size'   : 6}

matplotlib.rc('font', **font)
# service = vo.dal.TAPService("https://datalab.noirlab.edu/tap")
parser = argparse.ArgumentParser()
parser.add_argument("--photfile", required=True,type=str, help='phot.data file')
parser.add_argument("--singleproc",action="store_true")
parser.add_argument("--sector",type=int,required=True)
parser.add_argument("--N",type=int, default=4)
parser.add_argument("--bursttab", type=str)
args = parser.parse_args()

cam = int(args.photfile.replace("/phot.data","").split("/")[-1][3])
ccd = int(args.photfile.replace("/phot.data","").split("/")[-1][-1])
print(cam,ccd)
sourcedata = np.loadtxt(args.photfile, dtype=[('fcol','f8'),('frow','f8'),('icol','i8'),('irow','i8'),('fname',"<U128"),('flag','i4')])
if sourcedata.size==1:
    sourcedata = sourcedata[np.newaxis]
def ctstomag(cts):
    return -2.5*np.log10(cts) + 20.44
    # return cts
# Burst time in TJD
alldates = ""
for datesfile in glob.glob(f"/lustre/scratch/sarchast/s{args.sector:04}/cam{cam}-ccd{ccd}/o??/slice*/dates"):
    thisorbit = datesfile.split("/")[-3]
    with open(datesfile, "r") as f:
        for line in f:
            alldates  = alldates + f"{thisorbit} {line}"
         
dates = np.unique(np.genfromtxt(StringIO(alldates),dtype=[('detorbit','<U8'),('name','<U128'),('dates','f8')]))
dates.sort(order='dates')
orbitbounds = dates['dates'][:-1][(dates['dates'][1:] - dates['dates'][:-1]) > (1/24)]
def gaussian_nobackground(xy, amplitude, xo, yo, sigma_x, sigma_y, theta):
    x, y = xy
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                            + c*((y-yo)**2)))
    return g.ravel()

def plot(lcfile,dat):
    if "_cleaned" not in lcfile:
        raise Exception("Expected a _cleaned lcfile got instead:",lcfile)
    d = np.genfromtxt(lcfile)
    timesort = np.argsort(d[:,0])
    d = d[timesort]
    xvals = np.sort(d[:,0])
    yvals = (d[:,2] - d[:,6])#  + (d[:,8] - d[:,6])
    yerr = d[:,3]
    ymag = d[:,4]
    ymagerr = d[:,5]
    ybkg = d[:,6]
    tdiff = np.amin(np.round((xvals[1:] - xvals[:-1])*24*60*60).astype(int))
    srcname = lcfile.split('/')[-1]
    GRBname = srcname.split("_")[1]
    tburst = Time(dat['Trigger'][dat['GCN Name']==GRBname].to_numpy().astype(str)).jd-2457000
    tdiffdays = tdiff/60/60/24.0
    detorbit = dates['detorbit'][((tburst - dates['dates']) <= tdiffdays) & ((tburst - dates['dates']) > 0)][0]
    src = sourcedata[np.isin([s+"_cleaned" for s in sourcedata['fname'].astype(str)],f"lcGRB/{srcname}")]
    orbitfolders = [o.split("/")[-1] for o in glob.glob(f"/lustre/scratch/sarchast/s{args.sector:04}/cam{cam}-ccd{ccd}/o??")]
    orbitnames = np.array(orbitfolders)
    xlims = []
    ylims = []
    if src.size !=1:
        raise Exception(f"no data for {srcname}")
    frow = src['frow'][0] 
    fcol = src['fcol'][0] 
    icol = int(round(src['icol'][0]))
    irow = int(round(src['irow'][0]))
    ra, dec, _ = tess_stars2px_reverse_function_entry(args.sector, cam, ccd, fcol, frow)
    try: 
        flaggedlist = []
        rmslen = 20
        percentmax = 0.25
        runrms = []
        fig = plt.figure(figsize=(10,10)) 
        gs = GridSpec(6,5,figure=fig)
        for i in range(len(xvals)-rmslen):
            curybkg = ybkg[i:i+rmslen]
            curxvals = xvals[i:i+rmslen]
            if ((curxvals.max() - curxvals.min())>1):
                currms = np.nan
            else:
                currms = np.sqrt(np.average(curybkg**2))
            runrms.append((np.average(curxvals),currms))
            if currms > 100:
                for c in curxvals:
                    flaggedlist.append(c)
        for c in xvals[0:rmslen]:
            flaggedlist.append(c)
        for c in xvals[-rmslen:]:
            flaggedlist.append(c)
            

        flaggedtimes = np.array(flaggedlist)

       #  for i in range(len(xvals)-rmslen):
       #      curxvals = xvals[i:i+rmslen]
       #      flagpercent = np.sum(np.isin(curxvals,flaggedtimes))/20
       #      if flagpercent > percentmax:
       #          for c in curxvals:
       #              if c not in flaggedlist:
       #                  flaggedlist.append(c)
       #  flaggedtimes = np.array(flaggedlist)
        #################### Running RMS plot ############################3        
        runrms = np.array(runrms)
        ax = fig.add_subplot(gs[-4,:])
        ax.scatter(runrms[:,0],runrms[:,1])
        ax.axhline(10)
        ax.set_yscale('log')
        

        #################### CTS PER SEC PLOT ################### 
        ax = fig.add_subplot(gs[-3,:])
        ax.errorbar(xvals[~np.isin(xvals,flaggedtimes)],yvals[~np.isin(xvals,flaggedtimes)],yerr[~np.isin(xvals,flaggedtimes)],ls='none',color='black')
        ax.errorbar(xvals[np.isin(xvals,flaggedtimes)],yvals[np.isin(xvals,flaggedtimes)],yerr[np.isin(xvals,flaggedtimes)],ls='none',color='red',marker='X',label='flagged')
        negunflag = np.sum(yvals[~np.isin(xvals,flaggedtimes)] < 0 )
        try:
            xlims = (xvals.min()-(2/24), xvals.max()+(2/24))
            ylims = (yvals[~np.isin(xvals[(xvals>xlims[0]) & (xvals <xlims[1])],flaggedtimes)].min()-100,yvals[~np.isin(xvals[(xvals>xlims[0]) & (xvals <xlims[1])],flaggedtimes)].max()+100)
        except ValueError:
            nosetlims = True 
        # recalculate rms, exclude negatives
        inROI =  ((xvals-tburst) > (-0.5/24)) & ((xvals - tburst) < (1/24))  
        if np.sum(inROI)==0:
            raise Exception("No data near Tburst")
        outROI = ~inROI
        subxvals = xvals[inROI]
        subyvals = yvals[inROI]
        subymagvals = ymag[inROI]
        subxmag = subxvals[(subyvals >0)]
        subymag = ymag[np.isin(xvals,subxmag)]
        orbits = None
        otherflux = yvals[(~np.isin(yvals,subyvals)) & (~np.isin(xvals,flaggedtimes))]
        otherrms = np.sqrt(np.average(otherflux**2))
        rms = np.sqrt(np.average(subyvals[subyvals>0]**2))
        ctssensitivity = 3*np.sqrt(np.average(subyvals[(tburst-subxvals)>(200/60/60/24)]**2))
        magsensitivity = ctstomag(ctssensitivity)
        xmax = subxvals[np.amin(np.abs(subxvals-tburst))==np.abs(subxvals-tburst)]
        ymax = subyvals[subxvals==xmax]
        flagaroundmax = np.sum((flaggedtimes > (tburst-(2/24))) & (flaggedtimes < (tburst + (2/24))))
        if (rms<(3*(otherrms))) or (flagaroundmax > 0) or (otherrms > rms) or (np.abs(xmax-tburst)*24*60*60 > 400):
             prefix="flag_"
        else:
             prefix=""
        errmax = yerr[yvals==ymax]
        if errmax.size>1:
            errmax = errmax[0]
        ax.errorbar(xmax,ymax,errmax,ls='none',color='black',marker='*')
        ax.axvline(tburst,ls='--',label="T Burst")
        ax.legend()
        if ylims:
            ax.set_ylim(ylims[0],ylims[1])
            ax.set_xlim(xlims[0],xlims[1])
        # ax.set_yscale('symlog')
        ############### BKG PLOT ################
        ax = fig.add_subplot(gs[-2,:])
        ax.scatter(xvals[~np.isin(xvals,flaggedtimes)],ybkg[~np.isin(xvals,flaggedtimes)],color='black')
        ax.scatter(xvals[np.isin(xvals,flaggedtimes)],ybkg[np.isin(xvals,flaggedtimes)],color='red',marker='X',label='flagged')
        ax.axvline(tburst,ls='--',label="T Burst")
        ax.set_ylabel('Flux (counts/second)')
        ax.legend()
        if ylims:
            ax.set_ylim(ylims[0],ylims[1])
            ax.set_xlim(xlims[0],xlims[1])
        
        # ax.set_yscale('symlog')
        ############ magnitude plot ################3
        axmag = fig.add_subplot(gs[-1,:])
        axmag.invert_yaxis()
                    
        axmag.axvline(tburst,ls='--',label="T Burst")
               #  if orbits is not None:
               #      orbitGAPlabel = ["o1a->o1b","o1b->o2a","o2a->o2b"]
               #      for o,olabel in zip(orbits,orbitGAPlabel):
               #          ax.axvline(xvals[o],ls=':',label=olabel,color='black',alpha=0.5)
               #          axmag.axvline(xvals[o],ls=':',label=olabel,color='black',alpha=0.5)
        maglimitx = xvals[(yvals<ctssensitivity)]
        maglimity = np.array([magsensitivity for m in maglimitx])
        axmag.scatter(xvals[(yvals>0) & (~np.isin(xvals,flaggedtimes)) & (~np.isin(xvals,maglimitx))],ctstomag(yvals[(yvals>0) & (~np.isin(xvals,flaggedtimes)) & (~np.isin(xvals,maglimitx))]),color='black')
        # subxmag = subxvals[(subyvals >0) & (subxvals > (xmax - (1/24))) & (subxvals < (xmax + 1/24))]
        regxmag = subxmag[((tburst-subxmag) < (0.5/24)) & ((subxmag-tburst) < 1/24)]
        regymag = subymag[((tburst-subxmag) < (0.5/24)) & ((subxmag-tburst) < 1/24)]
        axmag.scatter(regxmag[~np.isin(regxmag,maglimitx)],regymag[~np.isin(regxmag,maglimitx)],color='blue',label=f"Near Tburst")
        detxmag = subxmag[(subymag < ctstomag(3*otherrms)) & ((tburst-subxmag) < (0.5/24)) & ((subxmag-tburst) < 1/24)]
        detymag = subymag[(subymag < ctstomag(3*otherrms)) & ((tburst-subxmag) < (0.5/24)) & ((subxmag-tburst) < 1/24)]
        axmag.scatter(detxmag,detymag,color='goldenrod',label=f"mag < {ctstomag(3*otherrms)}")
        axmag.scatter(maglimitx,maglimity,color='black',marker='v',label=f"mag > {magsensitivity}")
        if len(detxmag)==1:
            numpts=1
        else: 
            numpts=len(detymag)
        if numpts==0:
            prefix = "flag_"
           
        axmag.set_xlim(tburst-(1/24),tburst+(2/24))
        try:
            axmag.set_ylim(20,0.9*detymag.min())
        except Exception as e:
            axmag.set_ylim(20,10)
                
            
        axmag.legend() 
        trigcand = Time(tburst+2457000,format='jd').iso
       
        # plt.gca().set_xlabel('JD $-$ 2,457,000 (days)')
        axmag.set_xlabel('TJD (days)')
        title =  'TESS Light Curve of '+lcfile+f"|| RMS:{rms:9.2f} || "+str(numpts)+" points || negatives: "+str(negunflag)
        # plt.xlim([3797.7, 3798.2])
        # ax.set_xlim(0.1,0.25)
        # ax.set_xscale("log")

        mindiff = np.amin(np.abs(dates['dates'] - tburst))
        thisFILE = dates['name'][np.abs(dates['dates'] - tburst)==mindiff]
        thisINDEX = np.where(dates['name']==thisFILE)[0]
        IMSIZE=20
        axs = [fig.add_subplot(gs[0,-2]),fig.add_subplot(gs[0,-1]),fig.add_subplot(gs[1,-2]),fig.add_subplot(gs[1,-1])]
        ax = axs[0]
        ax.set_title("data")
        try:
            for curfile in glob.glob(f"/lustre/scratch/sarchast/s{args.sector:04}/cam{cam}-ccd{ccd}/{dates['detorbit'][thisINDEX][0]}/slice*/*{dates['name'][thisINDEX][0]}"):
                if "conv_" in curfile:
                    thisCONV = curfile
                if "interp_" in curfile:
                    thisINTERP = curfile
        except Exception as e:
            print(e)
            print(f"/lustre/scratch/sarchast/s{args.sector:04}/cam{cam}-ccd{ccd}/o??/slice*/",thisFILE)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with fits.open(thisINTERP) as hdul:
                header = hdul[0].header
            with fits.open(thisCONV) as hdul:
                data = -hdul[0].data/(tdiff*0.8*0.99)
                detim = np.copy(data)
                detax = ax
                # maxval = ymax[0]
                subdata = np.copy(data[(irow-IMSIZE):(irow+IMSIZE+1),(icol-IMSIZE):(icol+IMSIZE+1)])
                flatdata = np.ravel(subdata)
                # X2 = np.sort(flatdata)
                # F2 = np.arange(len(X2))/len(X2)
                # qscale = 0.9
                # minq = X2[F2<=(1-qscale)].max()
                # maxq = X2[F2>=qscale].min()

                # ax.imshow(data[(icol-50):(icol+51),(irow-50):(irow+51)], cmap='Greys',origin='lower',vmin=minq,vmax=maxq)
                # ax.imshow(subdata, cmap='Greys',origin='lower',vmin=minq,vmax=maxq)
                
             
                constraint=4
                imdata = np.copy(detim[(irow-constraint):(irow+constraint+1),(icol-constraint):(icol+constraint+1)])
                if imdata.shape[0]!=imdata.shape[1]:
                    raise Exception(lcfile,"is too close to the edge",icol,irow)
                x0 = (frow-irow)+constraint
                y0 = (fcol-icol)+constraint
                f0 = detim[int(x0),int(y0)]
                bsizepix = 1.2
                initial_guess = [f0,x0,y0,bsizepix,bsizepix,0]
                minf0 = max(0,f0*0.5)
                maxf0 = min(1e6,f0*1.5)
                if minf0==maxf0:
                   minf0 = 0
                   maxf0 = 1e9
                initbounds = [(minf0,maxf0),(x0-constraint,x0+constraint),(y0-constraint,y0+constraint),(bsizepix*0.5,bsizepix*5),(bsizepix*0.5,bsizepix*5),(0,2*np.pi)]
                bounds0 = tuple([b[0] for b in initbounds])
                bounds1 = tuple([b[1] for b in initbounds])
                bounds = [bounds0,bounds1]
                fitdata = np.transpose(imdata)
                circgauss = lambda xy, amplitude, xo, yo, sigma_x, sigma_y, theta: gaussian_nobackground(xy,amplitude,xo,yo,sigma_x, sigma_y, theta)
                xarr = np.arange(2*constraint+1)
                yarr = np.arange(2*constraint+1)
                subx,suby = np.meshgrid(xarr,yarr)
                preshape=fitdata.shape
                try:
                    popt, pcov = curve_fit(circgauss, (subx, suby), fitdata.ravel(), p0=initial_guess, bounds=bounds)
                    trueROW = popt[1] + (irow-constraint)
                    trueCOL = popt[2] + (icol-constraint)
                    wcs = WCS(header)
                    srcSC = wcs.pixel_to_world(trueCOL,trueROW)
                    ra = srcSC.ra.deg
                    dec = srcSC.dec.deg
                    srcGAL = srcSC.transform_to('galactic')
                    ax = axs[1] 
                    ax.set_title("model")
                    ax.axis("off")
                    model = circgauss((subx,suby),*popt).reshape(preshape)
                    vmin=0
                    vmax=model.max()
                    axs[0].imshow(subdata, cmap='Greys',origin='lower',vmin=0,vmax=vmax)
                    circle = plt.Circle(((fcol-icol)+IMSIZE,(frow-irow)+IMSIZE),2.5,fill=False)
                    ax.add_artist(circle)
                    ax.axis('off')
                    ax.imshow(model, cmap='Greys', origin='lower',vmin=vmin, vmax=vmax)
                    resid = np.abs(fitdata-model)
                    ax = axs[2]
                    ax.set_title("residuals")
                    ax.imshow(resid, cmap='Greys', origin='lower',vmin=vmin, vmax=vmax)
                    ax.axis("off")
                    ax = axs[3]
                    ax.text(0,0.1, f"{popt[0]:5.4e},{trueCOL:6.1f},{trueROW:6.1f},\n{popt[3]:5.4e},{popt[4]:5.4e},{popt[5]:5.4e}")
                    ax.text(0,0.7, f"Gaussian fit values")
                    ax.text(0,0.5, f"Amp, col., row, rad1, rad2, theta")
                    ax.axis("off")
                  #   if (np.abs(srcGAL.b.deg) < 10.1) and  (srcGAL.l.deg > 124.1) and (srcGAL.l.deg < 6.1):
                  #       decapsresult = service.search(f"SELECT TOP 10 * FROM decaps_dr2.object WHERE ra BETWEEN {ra-(3/60)} AND {ra+(3/60)} AND dec BETWEEN {dec-(3/60)} AND {dec+(3/60)} ")
                  #       ax.text(0,0.15,decapsresult)
                  #   panstarrs = Catalogs.query_region(coordinates=f"{ra} {dec}", radius=2/60/60, catalog="PANSTARRS", table="mean", data_release="dr2",columns=["objName","distance"])
                  #   if panstarrs:
                  #       panstarrs['distance']*=3600
                  #       ax.text(0,0.2,panstarrs)
                  #       print(panstarrs)
                  #   vizier = Vizier()
                  #   desdr2 = vizier.query_region(SkyCoord(ra, dec, unit='deg'), radius=4*astrou.arcsec, catalog="II/371/des_dr2")
                  #   if desdr2:
                  #       ax.text(0,0.3,desdr2)
                  #       print(desdr2)
                    
                except Exception as e:
                   print(e)
                   srcGAL = SkyCoord(ra, dec, unit='deg').transform_to('galactic')
            rootfolder = f"/lustre/research/mfausnau/data/tica/s{args.sector:04}/cam{cam}-ccd{ccd}"
            reffile = f"{rootfolder}/ref.fits" 
            axs = [fig.add_subplot(gs[0,2])]
            for ax, ffile in zip(axs,[reffile]):
                try:
                    with fits.open(ffile) as hdul:
                        header = hdul[0].header
                        data = hdul[0].data/(tdiff*0.8*0.99)
                        ax.imshow(data[(irow-IMSIZE):(irow+IMSIZE+1),(icol-IMSIZE):(icol+IMSIZE+1)], cmap='Greys',origin='lower',vmin=0,vmax=400)
                        circle = plt.Circle(((fcol-icol)+IMSIZE,(frow-irow)+IMSIZE),2.5,fill=False)
                        ax.add_artist(circle)
                        ax.axis('off')
                except Exception as e:
                    print(e)
            axs = []
            rmsind = 0
            allrms = np.unique(glob.glob(f"{rootfolder}/o??/rms.fits"))
            detidx = 0
            
            if allrms.size==0:
                print("could not find rms images in",f"{rootfolder}/o??/")
            if allrms.size>1:
                for idx,a in enumerate(allrms):
                    if detorbit in a:
                        orbitdetected=True 
                        detidx= idx
                if detidx<4:
                    halfrms = allrms[:4]
                else:
                    halfrms = allrms[4:]
                for i,rmsfile in enumerate(halfrms):
                        rowno = (i)%2
                        colno =  int((i/2)%2)
                        axs.append(fig.add_subplot(gs[colno,rowno]))
                datalist = []
                vmax = 5
                otherim = []
                founddetim = False
                for ffile in allrms:
                    if True:
                        oname = ffile.split("/")[-2]
                        with fits.open(ffile) as hdul:
                            header = hdul[0].header
                            filedata = hdul[0].data/(tdiff*0.8*0.99)
                            datalist.append(np.copy(filedata))
                            if oname==detorbit:
                                detim = filedata
                                founddetim=True
                            else:
                                otherim.append(filedata)
                                # print("Other im:",ffile)
                if not founddetim:
                    warnings.warn("Did not find detrms in rms images. Double length sectors are currently not working properly")
                
                medianOTHER = np.median(np.array(otherim).reshape(len(otherim),otherim[0].shape[0],otherim[0].shape[1]),axis=0)
                detim -= medianOTHER
                constraint=4
                imdata = detim[(irow-constraint):(irow+constraint+1),(icol-constraint):(icol+constraint+1)]
                x0 = (frow-irow)+constraint
                y0 = (fcol-icol)+constraint
                f0 = detim[int(x0),int(y0)]
                bsizepix = 1.2
                initial_guess = [f0,x0,y0,bsizepix,bsizepix,0]
                minf0 = max(0,f0*0.5)
                maxf0 = min(1e6,f0*1.5)
                if minf0==maxf0:
                   minf0 = 0
                   maxf0 = 1e9
                initbounds = [(minf0,maxf0),(x0-constraint,x0+constraint),(y0-constraint,y0+constraint),(bsizepix*0.5,bsizepix*5),(bsizepix*0.5,bsizepix*5),(0,2*np.pi)]
                bounds0 = tuple([b[0] for b in initbounds])
                bounds1 = tuple([b[1] for b in initbounds])
                bounds = [bounds0,bounds1]
                fitdata = np.transpose(imdata)
                circgauss = lambda xy, amplitude, xo, yo, sigma_x, sigma_y, theta: gaussian_nobackground(xy,amplitude,xo,yo,sigma_x, sigma_y, theta)
                xarr = np.arange(2*constraint+1)
                yarr = np.arange(2*constraint+1)
                subx,suby = np.meshgrid(xarr,yarr)
                popt, pcov = curve_fit(circgauss, (subx, suby), fitdata.ravel(), p0=initial_guess, bounds=bounds)
                ax = fig.add_subplot(gs[1,2])
                ax.imshow(detim[(irow-IMSIZE):(irow+IMSIZE+1),(icol-IMSIZE):(icol+IMSIZE+1)], cmap='Greys',origin='lower',vmax=vmax,vmin=0)
                circle = plt.Circle(((fcol-icol)+IMSIZE,(frow-irow)+IMSIZE),2.5,fill=False)
                ax.add_artist(circle)
                ax.axis('off')
                ax.set_title("RMS subtracted")
                                
                   #  except Exception as e:
                   #      print(e)
                for ax, data in zip(axs, datalist):
                    try:
                        ax.imshow(data[(irow-IMSIZE):(irow+IMSIZE+1),(icol-IMSIZE):(icol+IMSIZE+1)], cmap='Greys',origin='lower',vmax=vmax,vmin=0)
                        circle = plt.Circle(((fcol-icol)+IMSIZE,(frow-irow)+IMSIZE),2.5,fill=False)
                        ax.add_artist(circle)
                        ax.axis('off')

                    except Exception as e:
                        print(e)
                # axs[0].set_title(f"{popt}") 
                axs[1].set_title(f"ratio: {popt[-3]/popt[-2]}")
                ellipRAT = popt[-3]/popt[-2]
                if (ellipRAT > 1.5) or (ellipRAT < 0.5):
                    prefix="flag_"
           #  elif "flag" not in prefix:    
           #      axs.append(fig.add_subplot(gs[0:2,0:2]))
           #      avg = np.zeros((15,15))
           #      sq = np.zeros((15,15))
           #      N = 0
           #      for datenear in dates[((dates['dates'] - tburst)>(-0.5/24)) & ((dates['dates'] - tburst)< (1/24)) ]:
           #          for curfile in glob.glob(f"/lustre/scratch/sarchast/s{args.sector:04}/cam{cam}-ccd{ccd}/{datenear['detorbit']}/slice*/*{datenear['name']}"): 
           #              if datenear['dates'].astype(float) not in flaggedtimes:
           #                  with fits.open(curfile) as hdul:
           #                      avg += hdul[0].data[(irow-7):(irow+8),(icol-7):(icol+8)]
           #                      sq += hdul[0].data[(irow-7):(irow+8),(icol-7):(icol+8)]**2
           #                      N+=1
           #      rmswin = np.sqrt(sq/N - (avg/N)**2)
           #      flatdata = np.ravel(rmswin)
           #      X2 = np.sort(flatdata)
           #      F2 = np.arange(len(X2))/len(X2)
           #      qscale = 0.9
           #      minq = X2[F2<=(1-qscale)].max()
           #      maxq = X2[F2>=qscale].min()
           #      axs[-1].imshow(rmswin, cmap='Greys',origin='lower',vmax=maxq,vmin=minq)
                                   
                            
                        
            outfolder = lcfile.rsplit('/',2)[0]+"/plots"
            try:
                if not os.path.exists(outfolder):
                    os.mkdir(outfolder)
            except Exception as e:
                print(e)
            outfile = lcfile.split('/')[-1].replace(".txt","") + "cam" + str(cam) + "_" + "ccd" + str(ccd)+ ".jpg"
            if numpts==1:
                prefix = prefix +"single_"
            fig.suptitle(f"{title}\n CLASS: {prefix} || (RA,Dec):({ra:7.4f},{dec:7.4f}) || (col,row):({fcol:6.1f},{frow:6.1f})\nDatafile:{lcfile} || (l,b): ({srcGAL.l.deg:7.4f},{srcGAL.b.deg:7.4f})")
            plt.savefig(outfolder+"/"+prefix+outfile)
            plt.close()
            # with open(outfolder+"/"+prefix+outfile.replace(".png",".txt"),"w") as f:
            #     f.write(f"{outfile}\t{ra}\t{dec}\t{trigcand}\t{srcGAL.l.deg}\t{srcGAL.b.deg}")

    except Exception as e:
        print(traceback.format_exc())
def strip_spaces(a_str_with_spaces):
    return a_str_with_spaces.replace(' ', '')
if __name__=="__main__":
    fulldat = pd.read_csv(args.bursttab,converters={'GCN Name': strip_spaces})
    dat = fulldat[(fulldat["TESS Sector"]==int(args.sector)) & (fulldat["Enclosed Probability"] > 0)].copy(deep=True)
    dat['GCN Name'] = dat['GCN Name'].fillna('')
    dat.loc[dat["GCN Name"]=="","GCN Name"]=dat.loc[dat["GCN Name"]=="","Fermi Name"]
    for ind,row in dat.iterrows():
        lcfiles = glob.glob(f"{args.photfile.replace('phot.data','')}/lcGRB/lc_{row['GCN Name'].replace(' ','')}*_cleaned")
        tburst = Time(row['Trigger']).jd - 2457000
 
        if args.singleproc:
            for lc in tqdm(lcfiles):
                print(lc)
                plot(lc,dat)
        else:
            with Pool(args.N) as pool:
                pool.starmap(plot,[(lc,dat) for lc in lcfiles])

