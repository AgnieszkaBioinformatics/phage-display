# One hot encoding for clustering (k-means)

from Bio import SeqIO
import umap
import umap.plot
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from bokeh.plotting import save
import sklearn.cluster as cluster


### one hot encoding ###
nuc_codes ={'A': [1,0,0,0,0], 'C': [0,1,0,0,0], 'G': [0,0,1,0,0], 'T': [0,0,0,1,0], '-': [0,0,0,0,1]}
coded_matrix = []
seq_ids = []
fasta_file = 'clustal-bar2-750-1250-pel-myc.fasta'

with open(fasta_file, 'r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        # extracting the original seq
        org_seq = str(record.seq)

        # extracting the names of reads
        seq_id = record.description
       
        coded_nuc_seq = []

        # encoding for every nucleotide in the sequence
        for nuc in org_seq:
            coded_nuc_seq+=nuc_codes[nuc]

        # adding sequence to the list
        coded_matrix.append(coded_nuc_seq)
        seq_ids.append(seq_id)


### starting umap ###

# model learns about the data
mapper = umap.UMAP(n_neighbors=15, n_components=2, random_state=42).fit(coded_matrix)

# for interactive description plots
hover_data = pd.DataFrame({'id': seq_ids})
umap.plot.output_file('clustal-100.html')
p = umap.plot.interactive(mapper,
                          hover_data=hover_data,
                          point_size=5)
save(p)

# returning the transformed data as a numpy array
mapper=mapper.transform(coded_matrix)

# kmeans clustering
kmeans_labels = cluster.KMeans(n_clusters=10).fit_predict(coded_matrix)
plt.scatter(mapper[:, 0], mapper[:, 1], c=kmeans_labels, s=0.1, cmap='Spectral')
plt.savefig('kmeans-final.png')


### saving the result ###

# connecting the read names and cluster number to which the sequence was assigned
kmeans_labels_ready=[]
for n, c in zip(seq_ids, kmeans_labels):
    ready_k_tuple=()
    ready_k_tuple+=(n, c)
    kmeans_labels_ready.append(ready_k_tuple)

# storing the result in the dictionary in which the keys are the clusters, 
# the values are the reads that were classified to the particular cluster
un = np.unique(kmeans_labels)

ready_w_klastrach=dict()


for i in un:
    ready_w_klastrach[i]=[]

for k, v in kmeans_labels_ready:
    if v in ready_w_klastrach:
        ready_w_klastrach[v].append(k)

with open('kmeans-final.txt', 'w') as f:
    f.write(ready_w_klastrach)


# searching for the smallest cluster
for k, v in ready_w_klastrach.items():
    if v==min(ready_w_klastrach.values()):
        print(k)