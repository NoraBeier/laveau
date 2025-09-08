#!/bin/bash

DIRNAME=$(dirname "$0")

OUT=.
if [ "$#" -eq 1 ]; then
OUT=$1
mkdir -p $OUT
cd $OUT
fi

if [[ -z "${LAVEAU_PATH}" ]]; then
  PV="../scripts"
else
  PV="${LAVEAU_PATH}"
fi

cp $DIRNAME/rxn_demo.txt .

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~This is a demonstration for 4 Reactions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo
read -p "Press Enter to Start"
echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Step 00 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~ The RDM codes of the RCLASS for the corresponding reactions are converted into           ~~~" 
echo "~~~ the corresponding graph representation and can be found in the folder ‘00_RCLASS_Graphs’.~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python $PV/00_rdm2graph.py rxn_demo.txt
echo
echo "Step 01"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Step 01 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~ WARNING: this part maybe time expansive.                                                 ~~~"
echo "~~~ The RDM graphs are now assembled into the corresponding pairwise pattern graphs and can  ~~~" 
echo "~~~ be found in the folder ‘01_RDM_Graphs’.                                                  ~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python $PV/01_graph2pairwisegraph.py 00_RCLASS_Graphs
echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Step 02 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" 
echo "~~~ The RDM codes of the RCLASS for the corresponding reactions are converted into           ~~~" 
echo "~~~ the corresponding graph representation and can be found in the folder ‘02_ReactionRuless’.~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python $PV/02_pairwisegraph2reactionrules.py 01_RDM_Graphs
echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Step 03  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~ WARNING: this part maybe time expansive caused by high combinatorics.                    ~~~" 
echo "~~~ Now all RCLASS are combined into one reaction and the DPO rule formed                    ~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python $PV/03_reactionrules2DPO.py 02_ReactionRules rxn_demo.txt
echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Step 04 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~ The DPO rule is now applied to the substrates of the reaction, and the resulting hyper-edge~"         
echo "~~~ is used to derive the atom-to-atom map from it.’                                         ~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python $PV/preprocess_04.py 03_DPO_Rules
mod -f $PV/04_DPO2AAMs.py 
echo
echo "Fin."


