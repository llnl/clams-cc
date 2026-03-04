# ClaMS: Clustering Comparison

## About

This is a stand-alone repository that is part of the [ClaMS](https://github.com/llnl/ClaMS) project.
It compares clusterings in which each point is assigned to a single cluster label.
It reads clustering files of the form `point_id cluster_id` in each line and calculates a variety
of metrics to compare two clusterings.
For example, it can compare clusters found using [ClaMS](https://github.com/llnl/ClaMS).

## Build

```shell
git clone https://github.com/llnl/clams-cc.git
cd clams-cc
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```

## Run

We use [YGM](https://github.com/llnl/ygm/), a C++ asynchronous communication library, to efficiently compare clusterings in distributed memory. To take advantage of this, run the program on an HPC system with multiple ranks. For example, on a slurm-based (srun and sbatch) HPC system, request more nodes and ranks to speed up the calculation. A sufficient amount of compute needs to be provided to avoid running out of memory (see [note below](#a-note-on-space-requirements)).

## Calculating clustering metrics

Clustering files must be of the form `point_id cluster_id` in each line. Clusterings can have noise points which will be ignored when calculating metrics. Points that are labeled as noise in either clustering are not included in the metrics calculated. By default, points with label `-1` are noise points. If the label type `cluster_id_type` is changed to an unsigned type, then the maximum value becomes the noise label.

Currently, the resulting clustering metrics are just printed to the terminal or output file.

The script `src/clustering_metrics.cpp` calculates all supported clustering metrics. 
For two clusterings $X$ and $Y$, the clustering metrics it calculates are:
- Information theoretic quantities without adjustment
    - Unnormalized mutual information $I(X,Y)$ 
    - Geometric normalized mutual information $I(X,Y) / \sqrt{H(X)H(Y)}$
    - Arithmetic normalized mutual information $I(X,Y) / 0.5(H(X) + H(Y))$
- Information theoretic quantities with adjustment
    - Geometric-normalized adjusted mutual information
    - Arithmetic-normalized adjusted mutual information (AMI) 
      (this is the default adjusted mutual information calculate by python scikit-learn)
- Clustering metrics calculated from the pair-confusion matrix
    - Adjusted Rand index (ARI)
    - Fowlkes Mallows index
    - Purity (if `-p` flag provided)
    - (Pair-confusion) Balanced accuracy
    - (Pair-confusion) Geometric mean

You can compare a "ground truth" clustering to one or more clusterings in a single run. Here is an example:

```shell
#Example for calculating clustering metrics for individual files to compare to the first clustering
src/clustering_metrics -v -g ../path/to/ground_truth_clusters_file.txt ../path/to/comparison_clusters_file.txt ../path/to/comparison_clusters_file2.txt 

#Example for calculating clustering metrics for a list of clustering files to compare to the first clustering
src/clustering_metrics -v -g ../path/to/ground_truth_clusters_file.txt -l ../path/to/list_of_clusters_files.txt 
```

In this example,
- The `-v` option indicates verbose printout, which will print more intermediate outputs.
- The `-g` option is required and indicates the file name of the "ground truth" clustering
- The `-l` option indicates the file name of a file containing a list of clustering files to compare to the ground truth. Each line in this list should be the file name of a comparison clustering. Each comparison clustering will be compared to the ground truth clustering.
- The `-p` option indicates that we should also calculate and print out the purity. Purity is calculated for the clustering(s) with respect to the ground truth clustering indicated by `-g`.

### Calculate pair-confusion metrics only

The script `src/clustering_metrics_no_mi.cpp` calculates clustering metrics that are only based
on the pair confusion matrix. These metrics are:
- Adjusted Rand index (ARI)
- Fowlkes Mallows index
- Purity (if `-p` flag provided)
- (Pair-confusion) Balanced accuracy
- (Pair-confusion) Geometric mean

`src/clustering_metrics_no_mi.cpp` runs with the same user options as `src/clustering_metrics.cpp`. 
However, for large clusterings, it is currently much faster and has smaller space requirements 
than `src/clustering_metrics.cpp` because it does not calculate adjusted mutual information. The slowest step in the full computation is calculating the expected mutual information.

### Example

The `example` directory contains two small clusterings `example_clusters1.txt` and `example_clusters2.txt`. From the main directory, running
```
build/src/clustering_metrics -g example/example_clusters1.txt example/example_clusters2.txt
```
will produce the output in `example/example_output.txt`.

### A note on space requirements

#### Pair-confusion based clustering metrics

To calculate the clustering comparison metrics that rely on the pair confusion matrix, there
are three YGM maps used. These are:
1. A map of point id -> a pair of cluster ids indicate which cluster in clusterings 1 and 2
  the point belongs to
2. A map of cluster id -> cluster size (number of points) for clustering 1
3. A map of cluster id -> cluster size (number of points) for clustering 2
4. A map of pair of cluster ids (i,j) -> the (nonzero) overlap size 
  (number of points in cluster i of clustering 1 and cluster j of clustering 2).
  This should (hopefully) be like a sparse contingency/confusion matrix
5. If purity is calculated, a map of cluster id in clustering 2 -> size of its largest overlap with a cluster in (ground truth) clustering 1

Enough ranks need to be provided to hold these maps.

#### Information theoretic clustering metrics

The most expensive (both in time and in space) calculation in the clustering metrics is the 
adjusted mutual information (specifically, calculating the expected mutual information to 
make the adjustment). In addition to the above YGM maps, calculating the adjusted mutual
information requires the following YGM maps:

6. A map of cluster size -> number of clusters with that size for clustering 1
7. A map of cluster size -> number of clusters with that size for clustering 2
8. A map of pair of cluster sizes -> pair of count and contribution to the expected mutual 
   information.
   Here, there is an entry for each distinct pair of cluster sizes where one size is for a 
   cluster in clustering 1 and the other size is for a cluster in clustering 2.
Furthermore, the map (4 above) of overlap sizes is larger. Each pair of cluster ids (i,j) 
maps to a tuple containing overlap size, and the sizes of clusters i and j.

Therefore, the adjusted mutual information calculation has the most space requirements. 
Enough compute ranks need to be provided for this. 

## Extra Calculations

In the `extras` directory, there are some additional programs to calculate quantities useful for cluster analysis. Use the `-h` help flag for more details.

`extra/cluster_incidence` outputs the an "incidence" file, which is like an weighted edge list of a graph where the vertices are cluster ids $i,j$ in clustering 2 and the edge weights are the number of cluster ids in clustering 1 the are present in both $i$ and $j$.

`graph_cluster_analysis` assumes that a clustering labels the nodes in a graph, for example the output of a graph community detection algorithm. The required inputs are a clustering and a graph edge list. The program outputs the counts of the number of inter- or intra- cluster edges for each combination of cluster ids. Each file line is white-space separated and contains cluster id $i$, cluster id $j$, number of edges in the graph that are between a vertex in cluster $i$ and a vertex in cluster $j$. There are also optional arguments to output other information including: edge probabilities for each cluster pair, neighborhood purities of vertices, cluster sizes, vertex degrees in the graph.

# License

This project is licensed under the BSD-Commercial license – see the [LICENSE](LICENSE) file for details.

# Release
LLNL-CODE-2015147
