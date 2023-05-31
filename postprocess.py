# merge output files and remove old files

import glob
import os


files = glob.glob("out/*.dat")

records = []

for f in files:
    arr = f.split('-')

    if len(arr) != 2:
        continue
    
    p, r = arr
    iter_no = int(p.split('_')[1])

    rank_no = int(r.split('_')[1].split(".")[0])

    target = p+".dat"

    records.append((target, rank_no, f))

records = sorted(records)

# print(records)


for target, rank_no, f in records:

    # merge file
    os.system("type %s >> %s" %(f, target))

    # remove old file
    os.system("del %s" %f)
