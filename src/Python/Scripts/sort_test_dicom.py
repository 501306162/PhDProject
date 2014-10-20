import os
import os.path
import sys
import dicom
import shutil


def main():
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    files = get_files(input_folder)

    for f in files: 
        # check that the file is a cine image
        try:
            dcm = dicom.read_file(f)
            basename = os.path.basename(f)
            new_path = os.path.join(output_folder,basename)
            shutil.copy2(f, new_path)

        except Exception, e:
            raise e




def get_files(doc_root):
    paths = [os.path.join(doc_root,x) for x in os.listdir(doc_root)]
    files = []
    for path in paths:
        if os.path.isdir(path):
            files = files + get_files(path)
        else:
            files.append(path)

    return files




if __name__ == '__main__':
    main()

