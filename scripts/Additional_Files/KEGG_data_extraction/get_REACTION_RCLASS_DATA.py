import requests
import time
import re

BASE_URL = "https://rest.kegg.jp"

def get_all_reaction_ids():
    url = f"{BASE_URL}/list/reaction"
    resp = requests.get(url)
    resp.raise_for_status()
    ids = [line.split("\t")[0] for line in resp.text.strip().split("\n")]
    return ids

def parse_reaction_entry(reaction_id):
    url = f"{BASE_URL}/get/{reaction_id}"
    resp = requests.get(url)
    resp.raise_for_status()
    lines = resp.text.splitlines()

    compound_ids = []
    rclass_list = []

    # parse compound equation
    for line in lines:
        if line.startswith("EQUATION"):
            eqn = line.replace("EQUATION", "").strip()
            # Tokenize: numbers, compound IDs (Cxxxxxx), operators (+, <=>), parentheses
            compound_ids = re.findall(r'\(n\+1\)|\(n\)|[A-Za-z0-9\+\<\=\>\(\)]+', eqn)
            break

    # parse RCLASS
    in_rclass = False
    for line in lines:
        if line.startswith("RCLASS"):
            in_rclass = True
            parts = line.replace("RCLASS", "").strip().split()
            if len(parts) >= 2:
                rclass_list.append([parts[0], parts[1]])
        elif in_rclass and line.startswith("            "):
            parts = line.strip().split()
            if len(parts) >= 2:
                rclass_list.append([parts[0], parts[1]])
        elif in_rclass and not line.startswith(" "):
            break

    return compound_ids, rclass_list

def main(output_file="REACTION_RCLASS_DATA.txt"):
    reaction_ids = get_all_reaction_ids()
    with open(output_file, "w", encoding="utf-8") as f:
        for i, rid in enumerate(reaction_ids, 1):
            try:
                compounds, rclasses = parse_reaction_entry(rid)
                f.write(f"Reaction ID:{rid}\n")
                f.write(f"Compound IDs:{compounds}\n")
                f.write(f"RCLASS:{rclasses}\n")
                print(f"[{i}/{len(reaction_ids)}] {rid} done")
                time.sleep(0.2)  # nicht zu schnell abfragen
            except Exception as e:
                print(f"Fehler bei {rid}: {e}")

if __name__ == "__main__":
    main()

