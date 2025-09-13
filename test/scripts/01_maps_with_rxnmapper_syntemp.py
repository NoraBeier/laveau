#conda activate syntemp-env

from syntemp.SynAAM.rxn_mapper_wrapper import map_with_rxn_mapper
import pandas as pd

df = pd.read_csv('globalAAMs_unmapped.txt', header=None)
df.rename(columns={0:'id', 1:'rsmi'}, inplace=True)
data = df.to_dict('records')

for i in data:
    i['rxn_mapper'] = map_with_rxn_mapper(i['rsmi'])
print(data)

with open('globalAAMs_rxnmapper.txt', 'w') as f:
    for i in data:
        f.write(i['id']+','+i['rxn_mapper']+'\n')
