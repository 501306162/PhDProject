import h5py
import sys
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import SpectralClustering
import numpy as np
import pylab as pl
import nrrd


f = h5py.File(sys.argv[1], 'r')
max_size = int(sys.argv[3])

'''
locations = f['locations'][:]
adjacency = f['adjacency'][:]
affinity = f['affinity'][:]



# apply kernel
similarity = np.exp(-affinity / affinity.std())

print similarity.shape


#af  = AffinityPropagation(damping=0.99,max_iter=400,affinity='precomputed', verbose=True)
#af.fit(similarity)
#cluster_centers_indices = af.cluster_centers_indices_
#labels = af.labels_
sc = SpectralClustering(n_clusters = 50,affinity='precomputed')
sc.fit(similarity)
labels = sc.labels_





counts = np.bincount(labels)
max_val = np.argmax(counts)


print max_val  
'''

im, tmp = nrrd.read(sys.argv[2])
implot = pl.imshow((im[:,:,7]).T)
implot.set_cmap('gray')
locations = f['locations'][:]
labels = f['index'][:]

count = np.bincount(labels.ravel());

n_clusters = 0
for n in count:
    if n >= max_size:
        n_clusters+=1

for n,l in enumerate(labels):
    if count[l] >= max_size:
        pl.plot(locations[n,0], locations[n,1], 'o', color=pl.cm.jet(float(l) / np.max(n_clusters)))
pl.show()


