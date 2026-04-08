import numpy as np
import glob
import argparse
import os
from tess_stars2px import tess_stars2px_reverse_function_entry
parser = argparse.ArgumentParser()
parser.add_argument("photfile",type=str,nargs="*")
args = parser.parse_args()

cwd = os.getcwd()
print(cwd)
sector = cwd.split("sector")[-1].replace("/","")
for pf in args.photfile:
    newpf = ""
    lcfiles = glob.glob(pf.replace("phot.data","")+"/"+"**/outcatrms_*",recursive=True)
    validlcf = []
    with open(pf,"r") as pfh:
        for line in pfh:
            parts = line.split()
            if len(parts) != 6:
                print("skipping malformed line:", repr(line))
                continue
            col,row,icol,irow,lcfile,flag = parts
            if ((float(col)>48) and (float(col)<2080)) and ((float(row)>30) and (float(row)<2040)):
                newpf= newpf+line
                validlcf.append(lcfile.split("/")[-1])
            else:
                print(line,"removed")
    with open(pf,"w") as pfh:
        pfh.write(newpf)
    for lcf in lcfiles:
        basename = lcf.split("/")[-1]
        if basename not in validlcf:
            print("removing ",lcf)
            os.remove(lcf)
                
            
        
    
           
       
