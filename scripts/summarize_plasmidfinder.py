#!/usr/bin/env python3
import os
import sys
import csv
from glob import glob

DEFAULT_EMPTY="none"

def find_results_file(sample_dir):
    """
    Try common PlasmidFinder output filenames.
    """
    candidates = [
        os.path.join(sample_dir, "results_tab.tsv"),
        os.path.join(sample_dir, "results_tab.txt"),
        os.path.join(sample_dir, "plasmidfinder_results.tsv"),
        os.path.join(sample_dir, "results.tsv"),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    # fallback: any tsv that contains 'Inc' in it
    for c in glob(os.path.join(sample_dir, "*.tsv")):
        return c
    return None

def main(root_dir, out_tsv):
    rows = []
    for sample in sorted(os.listdir(root_dir)):
        sample_dir = os.path.join(root_dir, sample)
        if not os.path.isdir(sample_dir):
            continue

        res = find_results_file(sample_dir)
        if res is None:
            rows.append((sample, DEFAULT_EMPTY))
            continue

        incs = []
        with open(res, newline="") as f:
            # Try tab-delimited, tolerate headers
            reader = csv.reader(f, delimiter="\t")
            header = None
            for r in reader:
                if not r:
                    continue
                if header is None:
                    header = [x.strip() for x in r]
                    continue

                # Heurística: replicón suele estar en una columna llamada "Plasmid" o "Replicon" o similar
                # Si no, toma la primera columna.
                if "Plasmid" in header:
                    inc = r[header.index("Plasmid")].strip()
                elif "Replicon" in header:
                    inc = r[header.index("Replicon")].strip()
                else:
                    inc = r[0].strip()

                if inc and inc.lower() not in ("plasmid", "replicon"):
                    incs.append(inc)

        incs = sorted(set(incs))

        if len(incs) == 0:
            inc_str = DEFAULT_EMPTY
        else:
            inc_str = ";".join(incs)

        rows.append((sample, inc_str))

    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "inc_types"])
        for sample, incs in rows:
            w.writerow([sample, incs])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: summarize_plasmidfinder.py <plasmidfinder_out_dir> <out.tsv>", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
