#!/usr/bin/env python3
"""
Annotate the CDR loops of a nanobody / antibody heavy chain in a PDB.

Extracts the chain's sequence, numbers it with ANARCI (IMGT scheme), maps the
IMGT CDR ranges back to the ACTUAL residue numbers in the PDB, and prints:
  1. A reference table: each CDR's residue range and sequence in this structure.
  2. Ready-to-paste PyMOL selection + colour commands.

Why map back to real residue numbers: RF2 / pipelines can renumber structures,
so fixed "CDR = resi 27-38" assumptions are unreliable. ANARCI numbers the
sequence; we then translate those positions to whatever residue IDs this PDB
actually uses.

CDR definitions (IMGT scheme, inclusive):
  CDR-H1: IMGT 27-38
  CDR-H2: IMGT 56-65
  CDR-H3: IMGT 105-117
Insertion-coded positions (e.g. 111A) inside a range are included.

Usage:
    python3 cdr_annotate.py <pdb_file> [--chain H]

Run inside the rfantibody-bw env (needs ANARCI + biopython).
"""

import sys
import argparse
import warnings

warnings.filterwarnings("ignore")

# IMGT CDR boundaries (inclusive)
CDR_RANGES = {
    "H1": (27, 38),
    "H2": (56, 65),
    "H3": (105, 117),
}


def extract_seq_and_resnums(pdb_path, chain):
    """Return (sequence, resnums) for a chain: parallel lists of one-letter
    amino acids and their actual PDB residue numbers, in chain order."""
    from Bio import PDB
    from Bio.Data.PDBData import protein_letters_3to1

    structure = PDB.PDBParser(QUIET=True).get_structure("s", pdb_path)
    model = structure[0]
    if chain not in model:
        sys.exit(f"Chain {chain} not found in {pdb_path}. "
                 f"Chains present: {[c.id for c in model]}")

    seq, resnums = "", []
    for res in model[chain]:
        if res.id[0] != " ":          # skip hetero/water
            continue
        aa = protein_letters_3to1.get(res.resname, "X")
        seq += aa
        resnums.append(res.id[1])     # the integer residue number
    return seq, resnums


def number_with_anarci(seq):
    """Return ANARCI IMGT numbering: list of ((imgt_pos, ins_code), aa),
    excluding gap positions (aa == '-')."""
    from anarci import anarci
    results, _, _ = anarci([("q", seq)], scheme="imgt", output=False)
    if not results or results[0] is None:
        sys.exit("ANARCI could not number this sequence (is it an antibody chain?).")
    numbering = results[0][0][0]
    # keep only occupied positions, in order
    return [((pos, ins), aa) for (pos, ins), aa in numbering if aa != "-"]


def map_cdrs(seq, resnums, numbering):
    """Walk the ANARCI numbering alongside the sequence to assign each
    occupied IMGT position to a real PDB residue number.

    ANARCI's occupied positions correspond 1:1, in order, to the input
    sequence residues -- so the Nth occupied ANARCI position maps to the Nth
    residue, i.e. resnums[N].
    """
    cdrs = {name: [] for name in CDR_RANGES}
    for i, ((imgt_pos, _ins), aa) in enumerate(numbering):
        if i >= len(resnums):
            break
        real_resnum = resnums[i]
        for name, (lo, hi) in CDR_RANGES.items():
            if lo <= imgt_pos <= hi:
                cdrs[name].append((real_resnum, aa))
    return cdrs


def main():
    ap = argparse.ArgumentParser(description="Annotate antibody CDR loops from a PDB.")
    ap.add_argument("pdb_file", help="Path to the PDB file")
    ap.add_argument("--chain", default="H", help="Chain ID to analyse (default: H)")
    args = ap.parse_args()

    seq, resnums = extract_seq_and_resnums(args.pdb_file, args.chain)
    numbering = number_with_anarci(seq)
    cdrs = map_cdrs(seq, resnums, numbering)

    print(f"\nCDR annotation for {args.pdb_file}  (chain {args.chain})")
    print("=" * 64)
    print(f"{'CDR':<5} {'residues':<14} {'length':<7} sequence")
    print("-" * 64)
    for name in ("H1", "H2", "H3"):
        residues = cdrs[name]
        if not residues:
            print(f"{name:<5} {'(none found)':<14}")
            continue
        nums = [rn for rn, _ in residues]
        loop_seq = "".join(aa for _, aa in residues)
        rng = f"{min(nums)}-{max(nums)}"
        print(f"{name:<5} {rng:<14} {len(residues):<7} {loop_seq}")

    # PyMOL commands (use the real residue numbers just computed)
    print("\n" + "=" * 64)
    print("PyMOL commands (paste into PyMOL after loading the structure):")
    print("-" * 64)
    ch = args.chain
    print(f"color grey70, chain {ch}")
    palette = {"H1": "yellow", "H2": "orange", "H3": "red"}
    for name in ("H1", "H2", "H3"):
        residues = cdrs[name]
        if not residues:
            continue
        nums = [rn for rn, _ in residues]
        sel = f"chain {ch} and resi {min(nums)}-{max(nums)}"
        print(f"select cdr{name.lower()}, {sel}")
        print(f"color {palette[name]}, cdr{name.lower()}")
    print()


if __name__ == "__main__":
    main()
