
### load SMILES of Compounds
with open('/homes/biertank/nora/Githubs/Copy/AAMRuleTool/laveau_demo_fst/scripts/Additional_Files/KEGG_MoleculeDB_updated.txt','r') as f:
    lines = f.readlines()

comp_smiles = {}
for line in lines:
    comp = line.split(',')
    comp_smiles[comp[0]] = comp[1]

### load Compounds of Reaction
with open('/homes/biertank/nora/Githubs/Copy/AAMRuleTool/laveau_demo_fst/scripts/Additional_Files/REACTION_RCLASS_DATA.txt','r') as f:
    lines = f.readlines()

rxns = {}
side = False 
left = []
right = []

for line in lines:   
    if line.startswith('Reaction'):
        rxn_id = line.split(':')[1].strip()
    elif line.startswith('Compound IDs'):
        line = line.split("'")
        for c in line:        
            if c.startswith('<'):           
                side = True
            elif c.startswith('C') and not c.startswith('Co'):
                if side == False:
                    left.append(c)
                else:
                    right.append(c)
    elif line.startswith('RCLASS'):
        rxns[rxn_id] = {'reactant':left,'product':right}
        side = False 
        left = []
        right = []

## Load Reactions
results = {}
with open('/homes/biertank/nora/Githubs/Copy/AAMRuleTool/laveau_demo/benchmark/globalAAMs.smiles','r') as f:
    lines = f.readlines()

for line in lines:
    print(line)
    rxn = line.split(' ')[0]    
    comp = rxns[rxn]
    smiles = ''
    for i in range(len(comp['reactant'])): 
        string_comp = comp_smiles[comp['reactant'][i]]
        string_comp = string_comp.strip()
        smiles = smiles + string_comp
        if i != len(comp['reactant'])-1:
            smiles = smiles +'.'

    smiles = smiles + '>>'    
    for i in range(len(comp['product'])):     
        string_comp = comp_smiles[comp['product'][i]]
        string_comp = string_comp.strip()
        smiles = smiles + string_comp
        if i != len(comp['product'])-1:
            smiles = smiles +'.'
        
    results[rxn] = smiles


with open('/homes/biertank/nora/Githubs/Copy/AAMRuleTool/laveau_demo/benchmark/globalAAMs_unmapped.txt','a') as f:
    for rxn in results.keys():
        line = rxn+','+results[rxn]+'\n'
        f.write(line)
    
        

