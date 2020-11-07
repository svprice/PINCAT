# PINCAT: Protein Interaction Network Computational Analysis Tool
---
## Overview
PINCAT provides preprocessing and a number of ready-to-run computational tools for the analysis of protein-protein interaction networks. 

## Network Functions

  - [Guilia's network spatial entropy][l1]
  - [Guilia's network configuration entropy][l1]
  - [Marinka's network resilience][l2]

[//]: <# (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen.>

[l1]: <https://pubs.rsc.org/en/content/articlelanding/2015/MB/c5mb00143a#!divAbstract>
[l2]: <https://www.pnas.org/content/116/10/4426>

## Layout
- raw_data/
    - contains raw transcript expression data of the form GROUPID.txt
- preprocessing/
    - contains scipts for cleaning and manipulating the raw data into nodelists and edgelists
- network_data/
    - contains nodelists and edgelists 
- network_analysis/
    - contains a file network_features.py that has all the network functions
    - contains code to run experimeents on the network data and the corresponding outputs

## How to Perform Analysis on Partek Data
To run analysis, 
1. Clone the repository
2. Add your "GROUPID.txt" raw data file to the raw_data folder
3. Run partek_run.sh
    ```console
    bash preprocessing/partek_run.sh <GROUPID> partek_intermediate_files/ <CLASSIFICATION>
    ```
    where <CLASSIFICATION> is the label for the cells you want to select from the run.

    ```







