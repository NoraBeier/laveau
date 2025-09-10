import requests
import time

BASE_URL = "https://rest.kegg.jp"

def get_all_rclass_ids():
    url = f"{BASE_URL}/list/rclass"
    resp = requests.get(url)
    resp.raise_for_status()
    ids = [line.split("\t")[0] for line in resp.text.strip().split("\n")]
    return ids

def get_rclass_definition(rclass_id):
    url = f"{BASE_URL}/get/{rclass_id}"
    resp = requests.get(url)
    resp.raise_for_status()
    lines = resp.text.splitlines()
    defs = []
    for line in lines:
        if line.startswith("DEFINITION"):
            defs.append(line.replace("DEFINITION", "").strip())
        elif defs and line.startswith("            "):
            defs.append(line.strip())
        elif defs and not line.startswith(" "):
            break
    return " ".join(defs)

def main(output_file="list_RCLASS.txt"):
    rclass_ids = get_all_rclass_ids()
    with open(output_file, "w", encoding="utf-8") as f:
        for i, rid in enumerate(rclass_ids, 1):
            try:
                definition = get_rclass_definition(rid)
                f.write(f"{rid}\t{definition}\n")
                print(f"[{i}/{len(rclass_ids)}] {rid} done")
                time.sleep(0.2)  # KEGG mag es nicht zu schnell
            except Exception as e:
                print(f"Fehler bei {rid}: {e}")

if __name__ == "__main__":
    main()
