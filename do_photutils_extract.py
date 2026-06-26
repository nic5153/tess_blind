import numpy as np 
from astropy.stats import sigma_clipped_stats, SigmaClip
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint
from photutils.background import MMMBackground, Background2D, MedianBackground, MeanBackground, ModeEstimatorBackground, SExtractorBackground, BiweightLocationBackground, LocalBackground
from photutils.detection import DAOStarFinder
from photutils.psf import IterativePSFPhotometry, CircularGaussianPRF, PSFPhotometry
from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel, detect_sources, deblend_sources, SourceFinder, SourceCatalog
from astropy.io import fits
from astropy.table import setdiff
import matplotlib.pyplot as plt
import datetime
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("file",type=str)
args = parser.parse_args()

with fits.open(args.file) as hdul:
    data = hdul[0].data
sigma_clip = SigmaClip(sigma=3.0, maxiters=100)
threshold = detect_threshold(data, nsigma = 2.0, sigma_clip=sigma_clip)
segment_img = detect_sources(data, threshold, npixels=10)
footprint = circular_footprint(radius=10)
mask = segment_img.make_source_mask(footprint=footprint)
mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
print(np.array((mean, median, std)))

fig, axs = plt.subplots(nrows=4,ncols=1, figsize=(25,100))
estimator_names  = [ "MMMBackground"]
bkg_estimator_list =  [MMMBackground]
bkg_est = MMMBackground
ax = axs
e_name = estimator_names[0]
t1 = datetime.datetime.now() 
print("Model Background")
flatdata = np.ravel(data)
X2 = np.sort(flatdata)
F2 = np.arange(len(X2))/len(X2)
qscale = 0.9
minq = X2[F2<=(1-qscale)].max()
maxq = X2[F2>=qscale].min()
ax[0].imshow(data, origin='lower', cmap='Greys_r', vmin=minq, vmax=maxq, interpolation='nearest')
if "rms" not in args.file:
    bkg = Background2D(data, (50,50), filter_size=(3, 3), 
                       sigma_clip=sigma_clip, 
                       bkg_estimator=bkg_est())
    data -= bkg.background 
    noise = np.amin(bkg.background_rms)
else:
    bkg = Background2D(data, (50,50), filter_size=(3, 3), 
                       sigma_clip=sigma_clip, 
                       bkg_estimator=bkg_est())
    noise = np.amin(bkg.background_rms)
ax[1].imshow(bkg.background, origin='lower', cmap='Greys_r', vmin=minq, vmax=maxq, interpolation='nearest')
ax[2].imshow(data, origin='lower', cmap='Greys_r', vmin=minq, vmax=maxq, interpolation='nearest')
ax[0].set_title(e_name, fontsize=48)
fits.writeto("bgsub_"+args.file, data-bkg.background, overwrite=True)
print("Run DAOFind")
if args.file[0:4]=="hlsp":
    fwhm=1.5
else:
    fwhm=2.4
if "rms" in args.file:
    daofind = DAOStarFinder(fwhm=fwhm, threshold=10*np.amin(bkg.background_rms),
                        exclude_border=True)
else:
    daofind = DAOStarFinder(fwhm=fwhm, threshold=3*np.amin(bkg.background_rms),
                        exclude_border=True)
phot = daofind(data)

phot['xcentroid']-=0.5
phot['ycentroid']-=0.5
phot = phot[(phot['xcentroid'] > 48) & (phot['xcentroid'] < 2090) & (phot['ycentroid'] > 15) & (phot['ycentroid'] < 2053)]
with open("outcat"+args.file.replace(".fits",".tsv"),"w") as outfile:
    outfile.write(phot.to_pandas().to_csv(sep='\t',index=False))
t2 = datetime.datetime.now()
print(e_name,"elapsed:",t2-t1)
ax[3].imshow(data, origin='lower', cmap='Greys_r', vmin=minq, vmax=maxq, interpolation='nearest')
ax[3].scatter(phot['xcentroid'],phot['ycentroid'],color='red',marker='x',alpha=0.4)
plt.savefig(args.file.replace(".fits","")+"bkg_est.png")
plt.close()

