
# laveau 

<p align="center">
<img src="./Logo.png" width="300"/>
</p>

## Translation of RCLASS Atom-to-Atom information out of the KEGG database into DPO rules for the generation of Atom-to-Atom Maps

## Description

This repository is used to provide the laveau tool and as a database for the Atom-to-Atom maps and DPO rules generated from the KEGG database.

1. laveau tool

We developed laveau (can be found in the folder 'scripts'), a tool that computes explicit DPO rules from KEGG reactions and RCLASS data. The algorithm proceeds stepwise, starting from a translation of individual RDM codes into equivalent pattern graphs. Multiple RDM pattern graphs for the same RCLASS are then combined based on their embeddings into the educt and product molecules, observing certain consistency conditions. In the ﬁnal step, these combined pairwise patterns are merged into a pair of subgraphs of educts and products, respectively, which results in DPO graph rewrite rules. These rules are then used to generate Atom-to-Atom Maps for the reaction. Therefore, the Tool gives you an output with two lists of Atom-to-Atom map SMILES which are divided into a global (the complete reaction with all molecules) and partial (part of the molecules of the reaction) Atom-to-Atom Maps.

2. database

In the folder 'database' you can find all partial and global Atom-to-Atom Maps as well as DPO rules, which could be translated by laveau from the RCLASS RDM codes.

## Usage

Tool only works under python version > 3.6 

### Required packages

* networkx 2.8.8
* matplotlib 3.7.1
* rdkit 2023.9.2
* MØD version 0.16.0

### Step 0: Translate RDM Code in Graphs Fromat

RDM codes which are to be translated into DPO rules must be available as list in the following format as a .txt file like this example:

`RC00001	C1x-C8x:*-*:C2x+C2y-C8x+C8y N1y-N5y:*-*:C1y+C2x+C2x-C1y+C8x+C8x`

To process this list, the following call has to be run:

```console
python 00_rdm2graph.py <PATH:RCLASS_Data>
```

As output, a folder with the name 00_RCLASS_Graphs is created in which all graphs are saved as .xgef files

### Step 1: Translate RDM Code Graphs in RDM Graphs

This script requires the 'Labels' folder which can also be found in the 'scripts' folder. The path to the folder generated in step 0 must be specified as input.

To process this list, the following call has to be run:

```console
python 01_graph2pairwisegraph.py <PATH:00_RCLASS_Graphs>
```
As output a folder '01_RDM_Graphs' will be created which stores the RDM pattern graphs for each reaction side in named folders as .xgef format. The suffix number in the folder name indicates the variant resulting from the different possible type graphs of Step 0. The suffix number to the individual files, on the other hand, indicates the variants resulting from the merging of the RDM graphs. In addition, the individual files have the addition 'left or 'right' for the respective reaction side.


### Step 2: Translate RDM Graphs in Pairwise RDM Pattern

This script requires the 'RCLASS_RPAIR.txt' and 'List_UndefindAoms.txt' files which can also be found in the 'scripts/Additional Files' folder.

To process this list, the following call has to be run:

```console
python 02_pairwisegraph2reactionrules.py <PATH:01_RDM_Graphs>
```

As output, a folder '02_Reaction Rules' is generated with all reaction rule graphs in the respective separate folders. The suffix number to the individual files indicates the variants resulting from earlier steps. In addition, the individual files have the addition 'left or 'right' for the respective reaction side.

The files List_BigRCLASSES.txt and List_ErrorRCLASS.txt are also created. List_BigRCLASSES.tx lists all cases where the complexity of creating the rule was too great to execute it. List_ErrorRCLASS.txt, on the other hand, lists all RDM pattern graphs that could not be merged into a reaction rule.

### Step 3: Translate Pairwise RDM Pattern in Reaction Rules (DPO-Rules)

To process this list, the following call has to be run:

```console
python 03_reationrules2DPO.py <PATH:02_Pairwise_RDM_Pattern>
```
A folder '03_DPO_Rule' with all DPO rules is created as output.
The suffix number to the individual files, indicates the variants.

In addition, the following files are created in the folder '03_stats':
* log_unsucc: Lists all cases where no DPO rule could be created
* log_big: Lists all cases which were combinatorially too large
* log_missC: Lists all cases where structural information of single molecules was missing in the KEGG database.

### Step 4: Generate Atom-to-Atom Maps out of DPO-Rules

The MØD package is required for this step. Information on installation can be found at:
[https://jakobandersen.github.io/mod/](https://jakobandersen.github.io/mod/)

After activating the conda environment with the MØD package, step 4 can be called with the following command:

```console
mod -f 04_DPO2AAMs.py <PATH:03_Reaction_Rules>
```

The following files are created as output:

* output_PartialReactions.txt: Lists all discovered partial atom-to-atom maps in SMILES format
* output_GlobalReactions.txt: Lists all discovered global atom-to-atom maps in SMILES format
* output_FalseReactions.txt: Lists all cases where no atom-to-atom mpas could be found
