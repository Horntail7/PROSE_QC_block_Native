import math
import numpy as np
import argparse
import os.path
import matplotlib.pyplot as plt
import pandas as pd
import warnings
import seaborn as sns

from sklearn.preprocessing import StandardScaler  
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN


#from dtw import *
#from dtaidistance import dtw as dtwd
#import statistics

parser = argparse.ArgumentParser()

parser.add_argument('--inp', '--input_file', type=str, default='data.csv',
                    help='Full path to the input file (data.csv)')
#parser.add_argument('--k', '--key_file', type=str, default='key.dat',
#                    help='Full path to the key file (key.dat)')
parser.add_argument('--out', '--out_name', type=str, default='TSNE.png',
                    help='Full path to the output file (TSNE.png)')

args = parser.parse_args()

#kf=open('args.k','r')
#kys=kf.readlines()

df = pd.read_csv(args.inp)

X_orig=df[df.columns[1:]]
print(X_orig.head())
y=X_orig.pop('class_id')
X=np.array(X_orig)

#scaler = StandardScaler()  
#scaler.fit(X)  
#X_trans = scaler.transform(X)  
#print(X_trans)

model = TSNE(n_components=2, perplexity=30, init="random")
Y = model.fit_transform(X) 
print(Y.shape)


df_tsne = pd.DataFrame({
    'tsne1': Y[:,0],
    'tsne2': Y[:,1],
    'label': y  # your actual labels
})


sns.scatterplot(data=df_tsne, x='tsne1', y='tsne2', hue='label')
plt.show()

#x_min=-40
#x_max=-15
#y_min=-10
#y_max=5
#
#
#filtered = df_tsne[
#    (df_tsne['tsne1'] > x_min) &
#    (df_tsne['tsne1'] < x_max) &
#    (df_tsne['tsne2'] > y_min) &
#    (df_tsne['tsne2'] < y_max)
#]
#
#filtered_indices = filtered.index.to_list()

clustering = DBSCAN(eps=5, min_samples=5).fit(Y)
df_tsne['tsne_cluster'] = clustering.labels_

sns.scatterplot(data=df_tsne, x='tsne1', y='tsne2', hue='tsne_cluster', palette='tab10')
plt.title("t-SNE Clusters")
plt.show()

separated_cluster = df_tsne[df_tsne['tsne_cluster'] == 0]

print(separated_cluster)
sns.scatterplot(data=separated_cluster, x='tsne1', y='tsne2', hue='label', palette='tab10')
plt.title("t-SNE Clusters")
plt.show()

print(separated_cluster['label'].value_counts())

indices = separated_cluster.index.to_list()


merged_df = pd.concat([X_orig, y], axis=1)

X_filtered=merged_df.loc[indices]



X_filtered.to_csv("cluster.csv", index=False)

#print(indices)

#plt.scatter(Y[:, 0], Y[:, 1], c=y, cmap='viridis')
#plt.scatter(Y[:, 0], Y[:, 1], c=y, cmap='Paired')
#plt.savefig(args.out)
#plt.show() 
