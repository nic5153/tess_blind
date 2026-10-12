import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("photfile", type=str, nargs="*")
args = parser.parse_args()

print(os.getcwd())
# Drop candidates near the CCD edges from each phot.data, and delete the light
# curves of the dropped candidates only. (The old version globbed
# **/outcatrms_* for deletion, which removed the daily outcatrms_day_*.tsv
# source catalogs instead of light curves.)
for pf in args.photfile:
    camdir = os.path.dirname(pf)
    newpf = ""
    dropped = []
    with open(pf, "r") as pfh:
        for line in pfh:
            parts = line.split()
            if len(parts) != 6:
                print("skipping malformed line:", repr(line))
                continue
            col, row, icol, irow, lcfile, flag = parts
            if ((float(col) > 48) and (float(col) < 2080)) and ((float(row) > 30) and (float(row) < 2040)):
                newpf = newpf + line
            else:
                print(line, "removed")
                dropped.append(lcfile)
    with open(pf, "w") as pfh:
        pfh.write(newpf)
    for lcfile in dropped:
        for path in (os.path.join(camdir, lcfile), os.path.join(camdir, "bkg_phot", lcfile)):
            if os.path.exists(path):
                print("removing ", path)
                os.remove(path)
