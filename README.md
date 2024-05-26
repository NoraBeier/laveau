<p style="color:#1e8a30">Thank you for your interest in laveau. The tool is currently still under development. If you would like to use this tool in the future, please put it on your watch list.</p>


# laveau 

<p align="center">
<img src="./Logo.png" width="300"/>
</p>

## Translation of RCLASS atom-to-atom information out of the KEGG database into DPO rules for the generation of atom-to-atom maps

## Usage

Tool only works under python version > 3.6 

### Required packages
* networkx 2.8.8
* matplotlib 3.7.1
* rdkit 2023.9.2

### step 1: RDM pattern in tree representation
RCLASSes which are to be translated into DPO rules must be available as list in the following format as a .txt file:
`RCLASS equation`
as an example:
`RC00001	C1x-C8x:*-*:C2x+C2y-C8x+C8y N1y-N5y:*-*:C1y+C2x+C2x-C1y+C8x+C8x`
To process this list, the following call has to be run:
```console
python 00_rdm2graph.py <PATH:RCLASS_Data>
```
As output, a folder with the name 00_RCLASS_Graphs is created in which all trees are saved as .xgef files

### step 2: Tree representation in molecular structure subgraphs 
```console
python 01_rdm2molecule.py <PATH:00_RCLASS_Graphs>
```
### step 3: Molecular structure subgraphs in RCLASS graphs
```console
python 02_moleculeSubgraphPuzzle.py <PATH:01_Molecule_Graphs>
```
### step 3: RCLASS graphs in DPO-Rules for the entire reaction
```console
python 03_reaction_DPOrules.py <PATH:01_Puzzled_Graphs>
```
### step 4:  generate Atom-to-Atom Maps

## Data 
## Citation

