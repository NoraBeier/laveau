

# laveau 

<p align="center">
<img src="./Logo.png" width="300"/>
</p>

## Translation of RCLASS Atom-to-Atom Information out of the KEGG Database into DPO Rules for the Generation of Atom-to-Atom Maps

We developed laveau, a tool that computes explicit DPO rules from KEGG reactions and RCLASS data. The algorithm proceeds stepwise, starting from a translation of individual RDM codes into equivalent pattern graphs. Multiple RDM pattern graphs for the same RCLASS are then combined based on their embeddings into the educt and product molecules, observing certain consistency conditions. In the ﬁnal step, these combined pairwise patterns are merged into a pair of subgraphs of educts and products, respectively, which results in DPO graph rewrite rules. These rules are then used to generate Atom-to-Atom Maps for the reaction. Therefore, the Tool gives you an output with two lists of Atom-to-Atom map SMILES which are divided into a global (the complete reaction with all molecules) and partial (part of the molecules of the reaction) Atom-to-Atom Maps.

#### This repository serves to make the generated atom-to-atom maps and associated DPO rules available to the public. The scripts used to create them can be found under scripts.

Please note that all KEGG database reactions with status 06/04/2024 have already been processed by the tool. 

#### A demonstration of five sample reactions is available in demo/fun_demo.sh. For this purpose, simpler reactions were chosen so that the complete run finishes in about five minutes. Please note that runtime may increase considerably when using more complex examples.

Tool only works under python version > 3.6 

### Required packages

* networkx 2.8.8
* matplotlib 3.7.1
* rdkit 2023.9.2
* MØD version 0.16.0


