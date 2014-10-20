import dicom
import sys
import os
import os.path

def main():
    # get all the dicom directories
    doc_root = sys.argv[1]

    # look into all the directories
    dirs  = [os.path.join(doc_root,x) 
            for x in os.listdir(doc_root) 
            if os.path.isdir(os.path.join(doc_root,x))]


    for d_dir in dirs:
        dicom_files = get_dicom_files(d_dir)
        descriptions = get_series_descriptions(dicom_files)
        study_ids = get_study_ids(dicom_files)

        print study_ids, descriptions
        
        



""" get the list of study ids and the list of series descriptions """
def get_series_descriptions(dicom_files):
    descriptions = []
    for dcm in dicom_files:
        if not dcm.SeriesDescription in descriptions:
            descriptions.append(dcm.SeriesDescription)

    return descriptions

def get_study_ids(dicom_files):
    ids = []
    for dcm in dicom_files:
        if not dcm.StudyInstanceUID in ids:
            ids.append(dcm.StudyInstanceUID)

    return ids



""" recursive function to get all the dicom files """
def get_dicom_files(root_dir):
    paths = [os.path.join(root_dir, x) for x in os.listdir(root_dir)]
    files = []
    for path in  paths:
        if os.path.isdir(path):
            files = files + get_dicom_files(path)
        else:
            dcm = load_dicom(path)
            if dcm: files.append(dcm)
    return files

""" function to load the dicom file """
def load_dicom(path):
    try:
        dcm = dicom.read_file(path)
        if check_dicom(dcm) : return dcm 
    except:
        return False

""" function to check the dicom file is an image file """
def check_dicom(dcm):
    if dcm.StudyInstanceUID and dcm.SeriesDescription:
        return True
    return False

if __name__ == '__main__':
    main()

