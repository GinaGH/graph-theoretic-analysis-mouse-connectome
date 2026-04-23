# Graph-Theoretic Analysis of the Mouse Connectome

Graph-theoretic analysis of the Allen mouse brain connectome, including sparsity, degree structure, centrality, modularity, and spectral properties.


## Motivation

I am interested in the intersection of biology and mathematics, particularly neural networks. Understanding the structure inherent in connectome data and what it reveals about the organization and dynamics of the brain is central to my work.


## Project Structure

- src/ — source code for analysis of normalized connection density and normalized connection strength
- data/ — dataset location (see Data section below)
- outputs/ — generated plots, `.csv` files, and other non-terminal outputs


## Data

This project uses aggregated mouse connectome data from the Allen Institute:

https://download.alleninstitute.org/publications/A_high_resolution_data-driven_model_of_the_mouse_connectome/

### How to obtain the data

Download the following data files and place them in the `data/` directory:

```
data/
  normalized_connection_density.csv
  normalized_connection_strength.csv
```


## Graph-Theoretic Analysis

The analysis focuses on structural properties of the connectome using graph-theoretic methods:

* **Heatmaps** — visualization of ipsilateral and contralateral connectivity patterns and strengths
* **Sparsity Analysis** — threshold-based filtering to isolate biologically relevant connections
* **Community Detection** — modular structure via `greedy_modularity_communities` (NetworkX) and Laplacian spectral analysis
* **Centrality Analysis** — in/out degree distributions, hub identification, betweenness, PageRank, and spectral properties of adjacency and Laplacian matrices


## Notes

This repository contains code corresponding to the structural connectome analysis presented in the associated work. Extensions to functional connectivity and dynamic modeling are ongoing.
