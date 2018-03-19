import sys
import json
import numpy as np
from sklearn.cluster import KMeans
from argparse import ArgumentParser

parser = ArgumentParser(description="Cluster correlation results.")
parser.add_argument("corr_file", type=str)
parser.add_argument("output_prefix", type=str)
parser.add_argument('--corr_type', nargs='?', const="P4", type=str, default="P4")
parser.add_argument('--xmin', nargs='?', const=3, type=int, default=3)
parser.add_argument('--xmax', nargs='?', const=150, type=int, default=150)
opts = parser.parse_args()
labels = []
X = []
corr_results = {}
with open(opts.corr_file) as reader:
    for line in reader:
        data = json.loads(line)
        labels.append(data['ID'])
        corr = []
        for res in data['Results']:
            if res['Type'] == 'P2' and int(res['Lag']) == 0 and int(res['N']) > 0:
                if len(corr) == 0:
                    corr.append(float(res['Mean']))
                else:
                    corr[0] = float(res['Mean'])
            elif res['Type'] == opts.corr_type and int(res['Lag']) > 0 and int(res['Lag']) < 150 and int(res['N']) > 0:
                idx = int(res['Lag']) / 3
                while len(corr) <= idx:
                    corr.append(0.0)
                corr[idx] = float(res['Mean'])
        while len(corr) < 50:
            corr.append(0.0)
        X.append(corr)
        corr_results[data['ID']] = data
X = np.array(X)
n_clusters = 2
kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(X)

for i in range(n_clusters):
    outfile = "%s_%d.json" % (opts.output_prefix, i)
    with open(outfile, 'w') as w:
        for idx, c in enumerate(kmeans.labels_):
            if c == i:
                w.write(json.dumps(corr_results[labels[idx]]) + "\n")
    outfile = "%s_%d.txt" % (opts.output_prefix, i)
    with open(outfile, 'w') as w:
        for idx, c in enumerate(kmeans.labels_):
            if c == i:
                w.write(labels[idx] + "\n")