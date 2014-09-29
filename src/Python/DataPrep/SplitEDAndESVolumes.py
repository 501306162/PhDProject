import sys
import os
import shutil

folder = sys.argv[1]
ed_folder = sys.argv[2]
es_folder = sys.argv[3]



files = [os.path.join(folder, f) for f in os.listdir(folder) if os.path.isfile(os.path.join(folder,f))]

ed_files = []
es_files = []

for f in files:
    # get file number
    parts = f.split('_')
    num = int(parts[1])

    if f.find("ED.nrrd") > -1:
        ed_files.append([num,f])
    else:
        es_files.append([num,f])



# sort the lists on the number
ed_files = sorted(ed_files, key=lambda x :x[0])
es_files = sorted(es_files, key=lambda x :x[0])


# create the listing files
lf = open(os.path.join(ed_folder, 'listing.txt'), 'w')

# copy the output
for fp in  ed_files:
    # make output filename
    fname = os.path.join(ed_folder, 'label_' + str(fp[0]) + '.nrrd')
    shutil.copy2(fp[1], fname)
    lf.write(fname + '\n')

lf.close()


# create the listing files
lf = open(os.path.join(es_folder, 'listing.txt'), 'w')

# copy the output
for fp in  es_files:
    # make output filename
    fname = os.path.join(es_folder, 'label_' + str(fp[0]) + '.nrrd')
    shutil.copy2(fp[1], fname)
    lf.write(fname + '\n')

lf.close()
