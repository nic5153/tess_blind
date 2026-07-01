import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
from astropy.io import fits
import sys
import glob
import tarfile
import argparse

def get_frame_id(filename):
    """Extracts the integer frame ID from a TESS filename."""
    try:
        base = os.path.basename(filename)
        return int(base.split('-')[2])
    except (IndexError, ValueError):
        return None

def make_file_dict(lustre_ccd_dir, orbit=None):
    orbit_glob = orbit if orbit is not None else 'o*'
    dates_pattern = os.path.join(lustre_ccd_dir, orbit_glob, 'slice*/dates')
    dates_list = glob.glob(dates_pattern)
    
    file_info = {}
    for dates_path in sorted(dates_list):
        if "bkg_phot" in dates_path:
            continue
            
        path_parts = dates_path.split('/')
        # Identifies the orbit (e.g., 'o1a') from the path
        orbit_name = next((p for p in path_parts if p.startswith('o') and len(p) == 3), 'o1a')
        slice_dir = os.path.dirname(dates_path)
        
        try:
            with open(dates_path, 'r') as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 2:
                        f_id = get_frame_id(parts[0])
                        if f_id is not None:
                            file_info[f_id] = {
                                'full_name': os.path.basename(parts[0]),
                                'lustre_slice_path': slice_dir,
                                'orbit': orbit_name,
                                'tjd': float(parts[1])
                            }
        except Exception:
            pass
    return file_info

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--cam", type=int, required=True)
    parser.add_argument("--ccd", type=int, required=True)
    parser.add_argument("--orbit", type=str, default=None)
    args = parser.parse_args()

    work_base = "/home/nimcclur/work/TESS/photometry"
    lustre_base = "/lustre/research/mfausnau/data/tica"
    
    lustre_ccd_dir = f"{lustre_base}/s{args.sector:04}/cam{args.cam}-ccd{args.ccd}"
    work_ccd_dir = f"{work_base}/sector{args.sector}/cam{args.cam}_ccd{args.ccd}"
    
    orbit_glob = args.orbit if args.orbit is not None else "o*"
    tica_lists = glob.glob(f"{lustre_ccd_dir}/{orbit_glob}/rms_list_tica")
    lookup_dict = make_file_dict(lustre_ccd_dir, args.orbit)
    
    if not lookup_dict:
        return

    daily_bins = {}
    for t_list in tica_lists:
        try:
            with open(t_list, 'r') as f:
                for line in f:
                    item = line.strip().split()[0]
                    try:
                        f_id = int(item)
                    except ValueError:
                        f_id = get_frame_id(item)

                    if f_id in lookup_dict:
                        day = int(np.floor(lookup_dict[f_id]['tjd']))
                        if day not in daily_bins:
                            daily_bins[day] = []
                        daily_bins[day].append(f_id)
        except Exception:
            pass

    # Processing Loop
    for day, ids_in_day in sorted(daily_bins.items()):
        # Determine orbit folder from the first image of the day
        sample_info = lookup_dict[ids_in_day[0]]
        orbit_dir = os.path.join(work_ccd_dir, sample_info['orbit'])
        
        if not os.path.exists(orbit_dir):
            os.makedirs(orbit_dir, exist_ok=True)
        
        output_path = os.path.join(orbit_dir, f'rms_day_{day}.fits')
        if os.path.exists(output_path):
            print(f"Skipping existing: {output_path}")
            continue

        avg, square, N = 0, 0, 0
        
        for f_id in ids_in_day:
            info = lookup_dict[f_id]
            tar_path = os.path.join(info['lustre_slice_path'], 'images.tar')
            target_fits = 'conv_' + info['full_name']
            
            if os.path.exists(tar_path):
                try:
                    with tarfile.open(tar_path, "r") as tar:
                        f_obj = tar.extractfile(target_fits)
                        if f_obj:
                            with fits.open(f_obj) as hdul:
                                d = hdul[0].data
                                avg += d
                                square += d**2
                                N += 1
                except Exception:
                    pass

        if N > 0:
            rms_img = np.sqrt(square/N - (avg/N)**2)
            fits.writeto(output_path, rms_img, overwrite=True)
            print(f"Saved: {output_path}")

if __name__ == "__main__":
    main()
