import h5py
import numpy as np
import sys
import sklearn.cross_validation
import sklearn.svm 
import sklearn.decomposition
import math

def extract_data(filename):
    f = h5py.File(filename, "r")
    X = f["X"].value
    y = f["labels"].value
    owners = f["owners"].value
    f.close()

    return X,y,owners

def get_positive(X,y):
    n_samples = X.shape[0]
    output_list = []
    for n in range(n_samples):
        if y[n,0] == 1:
            output_list.append(X[n,:])

    X_out = np.zeros((len(output_list), X.shape[1]))
    for n, s in enumerate(output_list):
        X_out[n,:] = s.ravel()

    return X_out



def main():
    # load the data set
    data_filename = sys.argv[1]
    X, y, owners = extract_data(data_filename)
    
    # get the positive data 
    pos = get_positive(X,y)
    
    


    pca = sklearn.decomposition.PCA(n_components=35, whiten=True);
    m = pca.fit(pos)

    out = pca.score_samples(X)
    for n in range(X.shape[0]):
        print math.exp(out[n]), y[n,:]
    



if __name__ == '__main__':
    main()
