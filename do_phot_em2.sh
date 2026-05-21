#!/usr/bin/env bash

set -euo

duse=$DATA_DIR
dhome=$(pwd);

function do_phot(){

    sector=$1
    sectoruse=$(printf "s%04d" $1)
    echo "sector use:",$sectoruse
    cam=$2
    ccd=$3

    for o in o1a o1b o2a o2b; do 
    #for o in o2b; do 

	cp "sector$sector""/cam$cam""_ccd$ccd""/phot.data"   "$duse""/$sectoruse""/cam$cam""-ccd$ccd""/$o"
    done
    for o in o1a o1b o2a o2b; do 
	cd "$duse""/$sectoruse""/cam$cam""-ccd$ccd""/$o"
	for slice in $(ls -d slice*); do
	    cd $slice
	    for p in $(ls ../psf_file*fits); do
                [[ -e ${p:3:30} ]] || { ln -s $p ; }
            done
	    [[ -e psf_table ]] || { ln -s ../psf_table ; }
	    
	    if [ ! -d lc ]; then
		mkdir lc
	    fi

	    [[ -e kernel_data.tgz ]] && { tar -xzf kernel_data.tgz; }
	    
	    srun --job-name="phot2_cam${cam}-ccd${ccd}" \
                 --output=${LOG_DIR}/%x.o%j --error=${LOG_DIR}/%x.e%j \
                 --partition=nocona \
                 --nodes 1 --cpus-per-task 1 \
                 --ntasks-per-node=1\
                 ${ISIS_DIR}/phot2.csh &
	    
	    [[ -e 'kernel_data.tgz' ]] || {
		srun --job-name="phot2_cam${cam}-ccd${ccd}" \
		     --output=/dev/null --error=/dev/null \
		     --partition=nocona \
		     --nodes 1 --cpus-per-task 1 \
		     --ntasks-per-node=1\
		     tar -czf kernel_data.tgz {kt_,kc_}*fits &
	    }

	    cd ..
	done
	wait

	#clean up kernel files
	for slice in $(ls -d slice*); do
	    cd "$slice"
	    rm kt_*fits kc_*fits
	    cd ..
	done

	for slice in $(ls -d slice*); do
	    cd "$slice"
	    if [ ! -d bkg_phot ]; then
                mkdir bkg_phot
                cd bkg_phot
                bash $PIPELINE_DIR/setup/make_bkg_phot_dir_em2
            else
                cd bkg_phot
            fi

	    if [ ! -d lc ]; then
		mkdir lc
	    fi
	    [[ -e kernel_data.tgz ]] && { tar -xzf kernel_data.tgz; }

	    srun --job-name="bkg_phot2_cam${cam}-ccd${ccd}" \
                 --output=${LOG_DIR}/%x.o%j --error=${LOG_DIR}/%x.e%j \
                 --partition=nocona \
                 --nodes 1 --cpus-per-task 1 \
                 --ntasks-per-node=1\
                 ${ISIS_DIR}/phot2.csh &
	    
	    [[ -e 'kernel_data.tgz' ]] || {
                srun --job-name="phot2_cam${cam}-ccd${ccd}" \
                     --output=/dev/null --error=/dev/null \
                     --partition=nocona \
                     --nodes 1 --cpus-per-task 1 \
                     --ntasks-per-node=1\
                     tar -czf kernel_data.tgz {kt_,kc_}*fits &
            }
	
	    cd ../../;
	done
	wait
	
	#clean up kernel files
	for slice in $(ls -d slice*); do
	    cd "$slice"/bkg_phot
	    rm kt_*fits kc_*fits
	    cd ../../
	done

	cd ../ ;
    done
    
    cd $dhome
    return 0
}

#usage
#do_phot "01"    4   2 

function copy_phot(){

    [[ -d "$dhome""/sector$1""/cam$2""_ccd$3""/lc" ]] || { mkdir -p "$dhome""/sector$1""/cam$2""_ccd$3""/lc" ; }
    [[ -d "$dhome""/sector$1""/cam$2""_ccd$3""/bkg_phot" ]] || { mkdir -p "$dhome""/sector$1""/cam$2""_ccd$3""/bkg_phot" ; }
    [[ -d "$dhome""/sector$1""/cam$2""_ccd$3""/bkg_phot/lc" ]] || { mkdir -p "$dhome""/sector$1""/cam$2""_ccd$3""/bkg_phot/lc" ; }

    sectoruse=$(printf "s%04d" $1)
    dtarget=$DATA_DIR"/$sectoruse""/cam$2""-ccd$3"
    #mkdir "$dtarget""/lc"
    #mkdir "$dtarget""/bkg_phot/"
    #mkdir "$dtarget""/bkg_phot/lc/"
    for o in o1a o1b o2a o2b; do
	for slice in $(ls -d "$dtarget""/$o"/slice*); do
	    cd $slice"/lc"
	    for f in $(ls lc_*); do
		cat $f >> "$dhome""/sector$1""/cam$2""_ccd$3""/lc/$f" &
		#cat $f >> "$dtarget""/lc/$f" &
	    done
	    wait

	    cd ../bkg_phot/lc
	    for f in $(ls lc_*); do
		cat $f >> "$dhome""/sector$1""/cam$2""_ccd$3""/bkg_phot/lc/$f" &
		#cat $f >> "$dtarget""/bkg_phot/lc/$f" &
	    done
	    wait
	    cd ../../../
	done
    done


    cd  "$dhome""/sector$1""/cam$2""_ccd$3""/lc/"
    #cd  "$dtarget""/lc/"
    for f in $(ls lc_* | grep -v cleaned | grep -v png); do
        sort -nuk1 $f > tmp
        mv tmp $f
    done
    cd ../bkg_phot/lc
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
    
    for f in $(ls "$dhome""/sector$1""/cam$2""_ccd$3""/lc"); do

	echo $f
	ncheck=$(wc "$dhome""/sector$1""/cam$2""_ccd$3""/lc"/$f | awk '{print $1}')
	ncheck2=$(wc "$dhome""/sector$1""/cam$2""_ccd$3""/bkg_phot/lc/"$f | awk '{print $1}')
	echo "N_lc_home N_bkg_lc_home"
	echo $ncheck $ncheck2


	#make sure that all lines in /data/tess/image_sub have been concated
	nlines=0
	cd "$dtarget"
	for o in o1a o1b o2a o2b; do
	    cd $o
	    for slice in $(ls -d slice*); do
		cd $slice
		if [ -e "lc/$f" ]; then
		    ii=$(wc "lc/$f" | awk '{print $1}')
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
		for o in o1a o1b o2a o2b; do 
		    cd $o
		    for slice in $(ls -d slice*); do
			cd $slice
			rm lc/"$f"
			if [ -d bkg_phot ]; then
			    cd bkg_phot
			    rm lc/"$f"
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
