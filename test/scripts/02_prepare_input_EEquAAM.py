def merge_files(file1, file2, output_file):
    def load_file(filename, separator=","):
        entries = {}
        with open(filename, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue  # leere Zeilen überspringen
                if separator not in line:
                    continue
                id_, content = line.split(separator, 1)
                entries[id_.strip()] = content.strip()
        return entries

    entries1 = load_file(file1, separator=",")
    entries2 = {}
    with open(file2, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(maxsplit=1)
            if len(parts) < 2:
                continue
            id_, content = parts
            entries2[id_.strip()] = content.strip()

    # Nur IDs, die in beiden Files vorkommen
    ids = sorted(set(entries1.keys()) & set(entries2.keys()))

    with open(output_file, "w", encoding="utf-8") as out:
        for id_ in ids:
            out.write(f"#,{id_}\n")
            out.write(f"RXNmap,{entries1[id_]}\n")
            out.write(f"LAVmap,{entries2[id_]}\n")

#merge_files(
    "../map_list/globalAAMs_rxnmapper.txt",
    "../map_list/globalAAMs.smiles",
    "../map_list/input_EEquAAM.smiles"
#    )


import re

def normalize_and_filter(input_file, output_file):
    """
    1. Removes hydrogen atoms from the SMILES (CH2->C, OH->O, NH2->N, H+ removed)
    2. Writes only those blocks where file1 and file2 lines
       have the same number of letters (A–Z, a–z).
    """

    hydrogen_pattern = re.compile(r"\[H[+]?(:\d+)?\]")
    carbon_h_pattern = re.compile(r"\[CH\d*:([0-9]+)\]")
    oxygen_h_pattern = re.compile(r"\[OH:([0-9]+)\]")
    nitrogen_h_pattern = re.compile(r"\[NH\d*:([0-9]+)\]")

    def clean_line(line: str) -> str:
        line = hydrogen_pattern.sub("", line)
        line = carbon_h_pattern.sub(r"[C:\1]", line)
        line = oxygen_h_pattern.sub(r"[O:\1]", line)
        line = nitrogen_h_pattern.sub(r"[N:\1]", line)

        line = re.sub(r"\.\.", ".", line)
        line = line.replace(">>.", ">>")
        line = line.replace(".>>", ">>")
        line = line.replace("()", "")
        line = line.replace("@", "")
        return line.strip()

    def count_letters(s: str) -> int:
        return sum(1 for ch in s if ch.isalpha())

    with open(input_file, "r", encoding="utf-8") as f_in, \
         open(output_file, "w", encoding="utf-8") as f_out:

        lines = f_in.readlines()
        i = 0
        while i < len(lines):
            if lines[i].startswith("#,"):
                block_id = lines[i]
                if i + 2 < len(lines):
                    file1_line = clean_line(lines[i+1])
                    file2_line = clean_line(lines[i+2])

                    file1_payload = file1_line.split(",", 1)[-1]
                    file2_payload = file2_line.split(",", 1)[-1]

                    if count_letters(file1_payload) == count_letters(file2_payload):
                        f_out.write(block_id)
                        f_out.write(file1_line + "\n")
                        f_out.write(file2_line + "\n")
                i += 3
            else:
                i += 1


normalize_and_filter("../map_list/input_EEquAAM.smiles", "../map_list/input_EEquAAM_noH.smiles")

