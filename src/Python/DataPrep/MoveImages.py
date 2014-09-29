import sys
import os
import os.path
import shutil

def do_originals(folder, output_folder):
    files = []
    fp = open(os.path.join(output_folder, '..', 'OriginalImageFilenames.txt'), 'w')
    for n in range(16):
        fin = os.path.join(folder, 'patient_' + str(n) + '_time_0.nrrd')
        
        # create the output path
        fout = os.path.join(output_folder, 'image_' + str(n) + '.nrrd')
        
        shutil.copy2(fin, fout)
        fp.write(fout + '\n')

    
def do_normalised(folder, output_folder):
    fp = open(os.path.join(output_folder, '..', 'NormalisedImageFilenames.txt'), 'w')
    for n in range(16):
        fin = os.path.join(folder, 'patient_' + str(n) + '_time_0.nrrd')
        
        # create the output path
        fout = os.path.join(output_folder, 'image_' + str(n) + '.nrrd')
        
        shutil.copy2(fin, fout)
        fp.write(fout + '\n')




original_folder = sys.argv[1]
output_originals = sys.argv[2]
normalised_folder = sys.argv[3]
output_normalised = sys.argv[4]

do_originals(original_folder, output_originals)
do_normalised(normalised_folder, output_normalised)

