#!/usr/bin/env bash

set -euo

dhome=$(pwd);

function copy_phot(){

    [[ -d "$dhome""/sector$1""/cam$2""_ccd$3""/lcGRB" ]] || { mkdir -p "$dhome""/sector$1""/cam$2""_ccd$3""/lcGRB" ; }
    [[ -d "$dhome""/sector$1""/cam$2""_ccd$3""/bkg_phot" ]] || { mkdir -p "$dhome""/sector$1""/cam$2""_ccd$3""/bkg_phot" ; }
    [[ -d "$dhome""/sector$1""/cam$2""_ccd$3""/bkg_phot/lcGRB" ]] || { mkdir -p "$dhome""/sector$1""/cam$2""_ccd$3""/bkg_phot/lcGRB" ; }

    sectoruse=$(printf "s%04d" $1)
    duse=/lustre/scratch/${USER}
    dtarget=$duse"/$sectoruse""/cam$2""-ccd$3"
    #mkdir "$dtarget""/lcGRB"
    #mkdir "$dtarget""/bkg_phot/"
    #mkdir "$dtarget""/bkg_phot/lcGRB/"
    for o in ${dtarget}/o??; do
        echo $o
	for slice in $(ls -d "$o"/slice*); do
	    cd $slice"/lcGRB"
            pwd
	    for f in $(ls lc_*); do 
                echo $f
		cat $f >> "$dhome""/sector$1""/cam$2""_ccd$3""/lcGRB/$f" &
		#cat $f >> "$dtarget""/lcGRB/$f" &
	    done
	    wait

	    cd ../bkg_phot/lcGRB
	    for f in $(ls lc_*); do
		cat $f >> "$dhome""/sector$1""/cam$2""_ccd$3""/bkg_phot/lcGRB/$f" &
		#cat $f >> "$dtarget""/bkg_phot/lcGRB/$f" &
	    done
	    wait
	    cd ../../../
	done
    done


    cd  "$dhome""/sector$1""/cam$2""_ccd$3""/lcGRB/"
    #cd  "$dtarget""/lc/"
    for f in $(ls lc_* | grep -v cleaned | grep -v png); do
        sort -nuk1 $f > tmp
        mv tmp $f
    done
    cd ../bkg_phot/lcGRB
    for f in $(ls lc_* | grep -v cleaned | grep -v png); do
        sort -nuk1 $f > tmp
        mv tmp $f
    done
    cd "$dhome"

     #mv  "/data/tess/image_sub/sector$1""/cam$2""_ccd$3""/lc" "sector$1""/cam$2""_ccd$3"
	#mkdir "sector$1""/cam$2""_ccd$3""/bkg_phot"
	#mv  "/data/tess/image_sub/sector$1""/cam$2""_ccd$3""/bkg_phot/lc" "sector$1""/cam$2""_ccd$3""/bkg_phot"

#rm  "/data/tess/image_sub/sector$1""/cam$2""_ccd$3""/phot.data"
#rm  "/data/tess/image_sub/sector$1""/cam$2""_ccd$3""/bkg_phot/phot.data"

}

function clean_phot(){
    sectoruse=$(printf "s%04d" $1)
    dtarget=$DATA_DIR"/$sectoruse""/cam$2""-ccd$3"
    
    for f in $(ls "$dhome""/sector$1""/cam$2""_ccd$3""/lcGRB"); do

	echo $f
	ncheck=$(wc "$dhome""/sector$1""/cam$2""_ccd$3""/lcGRB"/$f | awk '{print $1}')
	ncheck2=$(wc "$dhome""/sector$1""/cam$2""_ccd$3""/bkg_phot/lcGRB/"$f | awk '{print $1}')
	echo "N_lc_home N_bkg_lc_home"
	echo $ncheck $ncheck2


	#make sure that all lines in /data/tess/image_sub have been concated
	nlines=0
	cd "$dtarget"
	for o in ${dtarget}/o??; do
	    cd $o
	    for slice in $(ls -d slice*); do
		cd $slice
		if [ -e "lcGRB/$f" ]; then
		    ii=$(wc "lcGRB/$f" | awk '{print $1}')
		    let nlines=$(($nlines + $ii))
		fi
		cd ..
	    done
	    cd ..
	done

	echo "N_lc_home N_bkg_lc_home N_lines_found_in_slices"
	echo "$ncheck $ncheck2 $nlines"
	if [[ $ncheck == $ncheck2 ]]; then
	    if [[ $ncheck == $nlines ]]; then
		echo "$f passed, deleted from slice directories"
		for o in ${dtarget}/o??; do 
		    cd $o
		    for slice in $(ls -d slice*); do
			cd $slice
			rm lcGRB/"$f"
			if [ -d bkg_phot ]; then
			    cd bkg_phot
			    rm lcGRB/"$f"
			    cd ../
			fi
			cd ../
		    done
		    cd ..
		done
		    
	    fi
	fi
	
    done
    
}
