# Clustering Hi-C contact graphs using Graph Neural Networks
Bioinformatics Institute spring project, year 2022.
Student: Velikonivtsev Fyodor
Supervisors: Tolstoganov Ivan, Korobeynikov Anton

<img width="330" alt="logo-bi-18-5" src="https://user-images.githubusercontent.com/67659154/169664789-545303fd-91fd-4411-829b-88a17b1c0524.png">

## Goal & objectives
**Goal:** discover properties, opportunities and perspectives of clustering metagenome Hi-C data by graph deep learning methods. \
**Objectives**:
- [X] Develop core understanding of problem´s crucial concepts and current advances
- [X] Apply effective GNN clustering models - DMoN, GraphMB 
- [X] Apply VAE model - VAMB
- [X] Create interface & modify models API  for correctly processing Hi-C data formats
- [X] Explore models´ hyperparameters space and compare results
- [X] Compare tools efficiency with VAMB and Bin3C as baselines using AMBER & CheckM

## Key results:
Bin was considered as HQ (high-quality metagenome-assembled genome) if it had >95% completeness and <5% Contamination
3 datasets have been used:
1. Zymo dataset (**supervised**) [6625 contigs, 76799 Hi-C links]
2. IC9 dataset (**unsupervised**) size  [165712 contigs, 1150887 Hi-C links]
3. CAMI AIRWAYS synthetic dataset (**supervised**) [728682 contigs, 70405 Hi-C links]
### _DMoN_:
* Zymo dataset - 0 HQ genomes
* Wase taken out of experiment

### _GraphMB_:
* Restored 7-8/10 HQ MAGs vs. _VAMB’s_ 10/10 vs _bin3C`s_ (non-DL tool for clustering Hi-C data) 6/10 (a)
* Restored 98/600 HQ MAGs vs. _VAMB’s_ 93/600 in CAMI AIRWAYS (b)
* Restored 12 HQ bins vs _VAMB’s_ 12 in IC9 dataset (c)
![image](https://user-images.githubusercontent.com/67659154/169665489-a3b73409-a11c-4bb3-b7b6-d28afdde7cf0.png)

## Conclusions:

## References:

# Tool:
This repository contains tools for:
* Correctly preprocessing Hi-C contact map and related assembly graph with scaffolds to create compatible input formats for DMoN & GraphMB - `preprocess_files.py`
* Correctly transforming ground truth labels (if any) to AMBER specific input format & post-clustering output modification for VAMB output - `vamb2amber.py`
* Basic exploratory data analysis by providing length distribution for subsets of contigs with Hi-C links and without Hi-C links - might be helpful in choosing minimal length threshold while clustering - `preprocess_files.py`
* Conveniently comparing CheckM quality assessment summaries of several binning runs by providing metrics plot for each of the runs (_example can be seen at plot (c)_) - `compare_checkm_results.py`

## Requirements
### Software requirements
* Desired tested GNN - DMoN, GraphMB
* QC tools - AMBER, CHECKM
* Python 3.8+
* Python packages (see installation for details)
* UNIX command line

### Hardware requirements
* GPU is preferrable
* CPU - any, multicore is preferrable
* RAM - 16 Gb
* Disk space - depends on data, + 2 Gb for CheckM database

## Installation
1. Install [AMBER](https://github.com/CAMI-challenge/AMBER), [CheckM](https://github.com/Ecogenomics/CheckM) in separate environments
2. Install desired tool to separated environment
3. Install python packages:
```{bash}
pip install -U numpy scipy pandas sklearn tqdm plotly kaleido
4. Clone repository and add it to PATH:
```{bash}
git clone https://github.com/Abusagit/GNN_plus_HiC.git && 
```

## Workflow
1. Preprocess data for given GNN (_e.g. GraphMB: transform `contact_map.tsv, scaffolds.fasta`_):
```{bash}
py preprocess.py --graphmb -c contact_map.tsv --scaling log -f scaffolds.fasta --mimic-jgi -o graphmb_input/
```

2. Run GNN:

```{bash}

```
