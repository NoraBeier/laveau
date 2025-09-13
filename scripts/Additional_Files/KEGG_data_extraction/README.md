
# KEGG Data Extraction Protocols

This repository contains scripts to extract and preprocess data from the [KEGG REST API](https://www.kegg.jp/kegg/rest/keggapi.html). 
The workflow ensures reproducibility by documenting API endpoints, data formats, and transformation steps.

## Overview

The scripts in this repository support three main tasks:

1. **Extracting RCLASS definitions**

   * `get_list_RCLASS.py`
   * Downloads all KEGG `RCLASS` IDs and their textual definitions.
   * Output: `list_RCLASS.txt` (tab-separated file with ID and definition).

2. **Extracting REACTION information from the KEGG Database**

   * `get_REACTION_RCLASS_DATA.py`
   * Downloads all KEGG `REACTION` IDs and their associated `Compounds` IDs and `RCLASS` entries.
   * Output: `REACTION_RCLASS_DATA.txt` (text file with listed Reaction ID:R00001,Compound IDs,RCLASS).

2. **Extracting RCLASS–RPAIR relationships**

   * `get_RCLASS_RPAIR.py`
   * Downloads all KEGG `RCLASS` IDs and their associated `RPAIR` entries.
   * Output: `RCLASS_RPAIR.txt` (text file with mappings between RCLASS and RPAIRs).

3. **Updating molecule database with SMILES**

   * `update_KEGG_MoleculeDB.py`
   * Fetches KEGG compound structures and attempts to generate canonical SMILES strings.
   * Uses multiple sources in priority order:

     1. **KEGG MOL files** (converted via RDKit)
     2. **PubChem** (via PubChemPy, optional)
     3. **Nikkaji SPARQL endpoint** (optional)
   * Input: `KEGG_MoleculeDB.txt` (CSV, two columns: `CompoundID,SMILES`)
   * Output: `KEGG_MoleculeDB_updated.txt` (updated CSV with missing/fixed SMILES).
   * Fallback MOL files are saved in the `mol_backup/` directory.

---

## Data Sources & Versioning

* **KEGG REST API**

  * Endpoints used:

    * `/list/rclass`
    * `/get/{rclass_id}`
    * `/list/compound`
    * `/dbget-bin/www_bget?-f+m+cpd:{compound_id}`
  * Data retrieved reflects the KEGG version at the date of download.
  * It is recommended to record the retrieval date and KEGG release in publications for reproducibility.

* **PubChem (optional)**

  * Used as a fallback for missing SMILES representations.

* **Nikkaji (optional)**

  * Queried via SPARQL for KEGG compound IDs to SMILES mappings.

---

## File Formats

* **Text outputs** (`.txt`): tab- or line-based KEGG ID lists and mappings.
* **Molecule structures**: KEGG provides MOL blocks, which are converted into SMILES strings using RDKit.
* **CSV outputs**: KEGG MoleculeDB (`CompoundID,SMILES`) for structured storage and further processing.



