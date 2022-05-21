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
* DMoN strongly depends from input node features and doesn’t support edge features, therefore it is incompatible for contact map clustering
* GraphMB showed better results on larger dataset where VAMB isn’t enough - it shows that some new information from contact graph was used to resolve bins
* GraphMB relies on the effective VAMB workflow - this explains close results on small Zymo dataset and its close to VAMB result on such small datasets
* GraphMB showed its competitive effectiveness in the problem of Hi-C contact map clustering  


## References:
1. Nissen, J.N., Johansen, J., Allesøe, R.L. et al. Improved metagenome binning and assembly using deep variational autoencoders. Nat Biotechnol 39, 555–560 (2021). https://doi.org/10.1038/s41587-020-00777-4
2. Metagenomic binning with assembly graph embeddings. Andre Lamurias, Mantas Sereika, Mads Albertsen, Katja Hose, Thomas Dyhre Nielsen. bioRxiv 2022.02.25.481923; doi: https://doi.org/10.1101/2022.02.25.481923
3. Tsitsulin, Anton & Palowitch, John & Perozzi, Bryan & Müller, Emmanuel. (2020). Graph Clustering with Graph Neural Networks. 
4. Fernando Meyer, Peter Hofmann, Peter Belmann, Ruben Garrido-Oter, Adrian Fritz, Alexander Sczyrba, Alice C McHardy, AMBER: Assessment of Metagenome BinnERs, GigaScience, Volume 7, Issue 6, June 2018, giy069, https://doi.org/10.1093/gigascience/giy069
5. Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Res. 2015 Jul;25(7):1043-55. doi: 10.1101/gr.186072.114. Epub 2015 May 14. PMID: 25977477; PMCID: PMC4484387.
6. DeMaere, M., Darling, A. bin3C: exploiting Hi-C sequencing data to accurately resolve metagenome-assembled genomes. Genome Biol 20, 46 (2019). https://doi.org/10.1186/s13059-019-1643-1


# __Toolkit__
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
1. Install [AMBER](https://github.com/CAMI-challenge/AMBER), [CheckM](https://github.com/Ecogenomics/CheckM) into separate environments
2. Install desired tool into the separated environment (crucial for GraphMB and AMBER) ([GraphMB - my modification](https://github.com/Abusagit/GraphMB), [DMoN - my modification](https://github.com/Abusagit/DMoN_for_HiC), [VAMB - my modification (minimal)](https://github.com/Abusagit/vamb))
3. Install python packages:
```{bash}
pip install -U numpy scipy pandas sklearn tqdm plotly kaleido
```
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
graphmb --assembly graphmb_input/ [--other-paraneters-for-graphmb]
```

3. Run AMBER or CheckM
4. You could also run VAMB - it produces clustering output incompatible for AMBER work but suitable for CheckM. You can use the followeing:

```{bash}
py vamb2amber.py -i amber_result.tsv -g golden_standard_with_amber_format.tsv -o outdir/vamb_for_amber.tsv
```


5. In the case of having labels you can directly compare binning results by AMBER - it provides comprehensive visualization plots. However, this is not the case for CheckM multiple study - here you can use `compare_checkm_results.py` and compare joint distribution of completeness and purity metrics among HQ genomes in m binning results (any number starting from 1):

```{bash}
py compare_checkm_results.py -i [checkm_input_1, ..., checkm_input_m] --labels [label_1, ..., label_m] --min-completeness 0.95 --min-purity 0.95
```
