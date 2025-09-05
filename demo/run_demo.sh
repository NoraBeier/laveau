#!/bin/bash
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~This is a demonstration for 4 Reactions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo
read -p "Push Enter to Start"
echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Step 00 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~ The RDM codes of the RCLASS for the corresponding reactions are converted into           ~~~" 
echo "~~~ the corresponding graph representation and can be found in the folder ‘00_RCLASS_Graphs’.~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
cd ../scripts/
python 00_rdm2graph.py ../demo/rxn_demo.txt
echo
echo "Step 01"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Step 01 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~ WARNING: this part maybe time expansive.                                                 ~~~"
echo "~~~ The RDM graphs are now assembled into the corresponding pairwise pattern graphs and can  ~~~" 
echo "~~~ be found in the folder ‘01_RDM_Graphs’.                                                  ~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python 01_graph2pairwisegraph.py ../demo/00_RCLASS_Graphs
echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Step 02 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" 
echo "~~~ The RDM codes of the RCLASS for the corresponding reactions are converted into           ~~~" 
echo "~~~ the corresponding graph representation and can be found in the folder ‘02_ReactionRuless’.~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python 02_pairwisegraph2reactionrules.py ../demo/01_RDM_Graphs 
echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Step 03  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~ WARNING: this part maybe time expansive caused by high combinatorics.                    ~~~" 
echo "~~~ Now all RCLASS are combined into one reaction and the DPO rule formed                    ~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python 03_reactionrules2DPO.py ../demo/02_ReactionRules ../demo/rxn_demo.txt
echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Step 04 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~ The DPO rule is now applied to the substrates of the reaction, and the resulting hyper-edge~"         
echo "~~~ is used to derive the atom-to-atom map from it.’                                         ~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python ./Additional_Files/preprocess_04.py '../demo/03_DPO_Rules/'
mod -f 04_DPO2AAMs.py 
echo
echo "Fin."


