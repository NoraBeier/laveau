import requests
import time

BASE_URL = "https://rest.kegg.jp"

def get_all_rclass_ids():
    url = f"{BASE_URL}/list/rclass"
    resp = requests.get(url)
    resp.raise_for_status()
    ids = [line.split("\t")[0] for line in resp.text.strip().split("\n")]
    return ids

def get_rclass_rpairs(rclass_id):
    url = f"{BASE_URL}/get/{rclass_id}"
    resp = requests.get(url)
    resp.raise_for_status()
    lines = resp.text.splitlines()
    rpairs = []
    for line in lines:
        if line.startswith("RPAIR"):
            rpairs.append(line.replace("RPAIR", "").strip())
        elif rpairs and line.startswith("            "):
            rpairs.append(line.strip())
        elif rpairs and not line.startswith(" "):
            break
    return rpairs

def main(output_file="RCLASS_RPAIR.txt"):
    rclass_ids = get_all_rclass_ids()
    with open(output_file, "w", encoding="utf-8") as f:
        for i, rid in enumerate(rclass_ids, 1):
            try:
                rpairs = get_rclass_rpairs(rid)
                if rpairs:
                    rpairs_str = ", ".join(rpairs)
                    f.write(f"RCLASS: {rid}\n")
                    f.write(f"RPAIRs: {rpairs_str}\n")
                    f.write("--------------\n")
                print(f"[{i}/{len(rclass_ids)}] {rid} done")
                time.sleep(0.2)  # KEGG nicht überlasten
            except Exception as e:
                print(f"Fehler bei {rid}: {e}")

if __name__ == "__main__":
    main()

