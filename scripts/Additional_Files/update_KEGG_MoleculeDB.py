import requests
import pandas as pd
from rdkit import Chem
import os

# Optional: PubChemPy für Fallback
try:
    import pubchempy as pcp
except ImportError:
    pcp = None

# Optional: SPARQLWrapper für Nikkaji-Fallback
try:
    from SPARQLWrapper import SPARQLWrapper, JSON
except ImportError:
    SPARQLWrapper = None

def fetch_mol_from_kegg(compound_id):
    url = f"https://www.kegg.jp/dbget-bin/www_bget?-f+m+cpd:{compound_id}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200 and len(r.text.strip()) > 0:
            return r.text
        else:
            return None
    except Exception as e:
        print(f"Error while loading the file {compound_id}: {e}")
        return None

def fetch_all_kegg_ids():
    url = "http://rest.kegg.jp/list/compound"
    r = requests.get(url)
    r.raise_for_status()
    ids = []
    for line in r.text.strip().split("\n"):
        print(line)
        ids.append(line.split("\t")[0])
    return ids

def fetch_smiles_from_pubchem(compound_id):
    if pcp is None:
        return None
    try:
        results = pcp.get_compounds(compound_id, "name")
        if results:
            return results[0].canonical_smiles
    except Exception as e:
        print(f"PubChem-Suche für {compound_id} fehlgeschlagen: {e}")
    return None

def fetch_smiles_from_nikkaji(compound_id):
    if SPARQLWrapper is None:
        return None
    try:
        sparql = SPARQLWrapper("https://jglobal.jst.go.jp/sparql")
        sparql.setQuery(f"""
        PREFIX nc: <http://rdf.nikkaji.jp/cheminf#>
        SELECT ?smiles WHERE {{
          ?compound nc:hasKEGGCompoundID "{compound_id}" ;
                    nc:hasCanonicalSMILES ?smiles .
        }}""")
        sparql.setReturnFormat(JSON)
        res = sparql.query().convert()
        if res["results"]["bindings"]:
            return res["results"][0]["smiles"]["value"]
    except Exception as e:
        print(f"Nikkaji-Suche für {compound_id} fehlgeschlagen: {e}")
    return None

INPUT_TXT = "KEGG_MoleculeDB.txt"  
OUTPUT_TXT = "KEGG_MoleculeDB_updated.txt"
MOL_BACKUP_DIR = "mol_backup"  
os.makedirs(MOL_BACKUP_DIR, exist_ok=True)

try:
    df = pd.read_csv(INPUT_TXT, sep=",", header=None, names=["CompoundID", "SMILES"], dtype=str)
except Exception as e:
    raise RuntimeError(f"Error while reading the file {INPUT_TXT}: {e}")

existing_ids = set(df["CompoundID"])
all_kegg_ids = fetch_all_kegg_ids()

new_rows = []
updated = []
for compound_id in all_kegg_ids:
    if compound_id not in existing_ids or (df.loc[df["CompoundID"] == compound_id, "SMILES"] == "FALSE").any():
        smiles = "FALSE"

        # KEGG MOL → RDKit
        mol_text = fetch_mol_from_kegg(compound_id)
        if mol_text:
            mol = Chem.MolFromMolBlock(mol_text, sanitize=False)
            if mol:
                try:
                    Chem.SanitizeMol(mol)
                    smiles = Chem.MolToSmiles(mol)
                except Exception:
                    smiles = Chem.MolToSmiles(mol)
            else:
                with open(os.path.join(MOL_BACKUP_DIR, f"{compound_id}.mol"), "w") as f:
                    f.write(mol_text)

        # PubChem-Fallback
        if smiles == "FALSE":
            pubchem_smiles = fetch_smiles_from_pubchem(compound_id)
            if pubchem_smiles:
                smiles = pubchem_smiles

        # Nikkaji-Fallback
        if smiles == "FALSE":
            nikkaji_smiles = fetch_smiles_from_nikkaji(compound_id)
            if nikkaji_smiles:
                smiles = nikkaji_smiles

        if compound_id in existing_ids:
            df.loc[df["CompoundID"] == compound_id, "SMILES"] = smiles
            updated.append(compound_id)
        else:
            new_rows.append([compound_id, smiles])

if new_rows:
    df_new = pd.DataFrame(new_rows, columns=["CompoundID", "SMILES"])
    df = pd.concat([df, df_new], ignore_index=True)

df = df.sort_values(by="CompoundID").reset_index(drop=True)
df.to_csv(OUTPUT_TXT, sep=",", index=False, header=False)
print(f"Fertig! {len(new_rows)} neue IDs ergänzt, {len(updated)} ersetzt. Ergebnis in {OUTPUT_TXT} gespeichert.")

