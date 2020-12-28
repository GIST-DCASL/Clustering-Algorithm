# Clustering_Algorithm 
## Definition of the Cluster
> Ma, Jeong-Min, et al. "Clusters in multi-leader directed consensus networks." 2020 20th International Conference on Control, Automation and Systems (ICCAS). IEEE, 2020.

By the conensus protocol, the agents reach steady-state. If steady-states of every agent are same, the system reaches a consensus.

However, at the multi-leader directed consensus network, the system cannot reach a consensus.

At the above paper, *Cluster* is defined in this multi-leader directed consensus network.

The paper guarantees that the steady-state value of all nodes in same cluster converges to common value whatever intial values and weights are.

### Constraint of the clustering algorithm
There is an important constraint to use the clustering algorithm.

The input graph should not have the leader-like strongly connected component(LSCC) which is mentioned in the paper above.

The Fig 1 is an example of LSCC.

<p align="center"><img src="https://user-images.githubusercontent.com/39582428/103194986-a886f800-4924-11eb-9d6e-da170ef003c1.JPG" width="300" height="200">
<p align="center">Fig 1. Example of LSCC <p align="center">


## Use Clustering Algorithm
This algorithm finds all clusters by using only topological data. (The input is incidence matrix of the graph).

This algorithm works on any directed graph with multiple leaders.

You can test your multi-leader directed graph by making an incidence matrix in .xlsx file. (The node number can be setted arbitrary.)

The incidence matrix of your graph should be composed with 0, 1 and -1 like the Fig 2.

<p align="center"><img src="https://user-images.githubusercontent.com/39582428/103191180-5f7c7700-4917-11eb-8185-ec632fee019b.JPG" width="700" height="200">
<p align="center">Fig 2. The incidence matrix of the graph in Fig 3<p align="center">

<p align="center"><img src="https://user-images.githubusercontent.com/39582428/103191183-61463a80-4917-11eb-806a-884cbeae25dd.JPG" width="300" height="350">
<p align="center">Fig 3. An example of the multi-leader directed graph<p align="center">

At Fig 3, each circle of green dotted line means a cluster. The oragne nodes are CR node which is defined in the above paper. There is only one CR node in each cluster.

If you run the matlab code, you can take a result like the Fig 4.

<p align="center"><img src="https://user-images.githubusercontent.com/39582428/103192225-45dd2e80-491b-11eb-8887-440e647115ce.JPG" width="450" height="300">
<p align="center">Fig 4. Clustering result for the graph in Fig 3<p align="center">

*Clusternum* means the number of clusters.

And at the *Clusterset*, the nodes on the same row belong to same cluster.

Additionally, the first component of each row of clusterset is CR node.
