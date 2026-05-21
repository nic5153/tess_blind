#!/usr/bin/env bash


duse="/data/tess/image_sub/"
dhome=$(pwd);

function do_phot(){

sector=$1
cam=$2
ccd=$3

cp "sector$sector""/cam$cam""_ccd$ccd""/phot.data"   "$duse""/sector$sector""/cam$cam""_ccd$ccd""/"
cd "$duse""/sector$sector""/cam$cam""_ccd$ccd"
mkdir lc
~/isis/phot2.csh #&
cd "$duse""/sector$sector""/cam$cam""_ccd$ccd""/bkg_phot"
cp ../phot.data .
pwd
mkdir lc
~/isis/phot2.csh 
cd $dhome
return 0
}

#usage
#do_phot "01"    4   2 

function copy_phot(){

mv  "/data/tess/image_sub/sector$1""/cam$2""_ccd$3""/lc" "sector$1""/cam$2""_ccd$3"
mkdir "sector$1""/cam$2""_ccd$3""/bkg_phot"
mv  "/data/tess/image_sub/sector$1""/cam$2""_ccd$3""/bkg_phot/lc" "sector$1""/cam$2""_ccd$3""/bkg_phot"

rm  "/data/tess/image_sub/sector$1""/cam$2""_ccd$3""/phot.data"
rm  "/data/tess/image_sub/sector$1""/cam$2""_ccd$3""/bkg_phot/phot.data"
}
