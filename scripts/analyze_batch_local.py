#!/usr/bin/env python3
"""
Analyze and rank RFantibody design batches (LOCAL / quiver-native version).

The local Blackwell build is quiver-native: metrics live inside the .qv files
as QV_SCORE lines, not in loose .trb / _best.pdb footers. This script reads
them via `qvscorefile`, which writes a tab-separated <name>.sc next to each
quiver.

Two scorefiles matter per batch:
  - <outdir>/1_rfdiffusion.sc : columns  mindist, averagemin, tag   (geometry)
  - <outdir>/4_rf2.sc         : columns  interaction_pae, pae, pred_lddt,
                                target_aligned_*_rmsd, framework_aligned_*_rmsd,
                                tag  (interface + fold quality)

We join the two on the design's backbone: an RF2 tag like
`samples_design_3_dldesign_1_best` traces back to RFdiffusion tag
`samples_design_3`, which carries the mindist geometry.

CDR-H3 length is computed from each design's final PDB by extracting the
heavy-chain sequence and numbering it with ANARCI (IMGT scheme). Under IMGT,
CDR-H3 is positions 105-117 inclusive; long (camelid) loops use insertion codes
(111A, 112A, ...) which ARE counted, so length is the number of occupied
positions in that range, NOT 117-105+1.

Usage:
    python3 analyze_batch_local.py <batch_dir> [<batch_dir> ...]
    # e.g. python3 analyze_batch_local.py outputs/batchB outputs/batchA

Requires the `qvscorefile` CLI and ANARCI (+ HMMER, biopython), i.e. run inside
the rfantibody-bw env.
"""

import sys
import os
import csv
import glob
import subprocess
import warnings

warnings.filterwarnings("ignore")

# ---- Thresholds for a "promising" design (unchanged from RunPod version) -----
PLDDT_MIN = 0.85          # pred_lddt above this
IPAE_MAX = 10.0           # interaction_pae below this (lower = better interface)
DIST_MAX = 8.0            # min hotspot-loop distance below this (Angstrom)

# ---- IMGT CDR-H3 boundary (inclusive). Strict IMGT definition. ---------------
H3_IMGT_START = 105
H3_IMGT_END = 117


def run_qvscorefile(qv_path):
    """Ensure the .sc scorefile exists for a given quiver, then return its path.

    qvscorefile writes <name>.sc next to <name>.qv and takes no output arg.
    Regenerated each run so a stale .sc never misleads us.
    """
    if not os.path.exists(qv_path):
        return None
    sc_path = qv_path[:-3] + ".sc" if qv_path.endswith(".qv") else qv_path + ".sc"
    subprocess.run(["qvscorefile", qv_path], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return sc_path


def parse_scorefile(sc_path):
    """Parse a tab-separated .sc file into a list of dict rows.

    Numeric values become floats; 'NaN'/blank become None. tag stays a string.
    """
    rows = []
    if not sc_path or not os.path.exists(sc_path):
        return rows
    with open(sc_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            parsed = {}
            for k, v in r.items():
                if k == "tag":
                    parsed[k] = v
                    continue
                try:
                    parsed[k] = float(v)
                except (ValueError, TypeError):
                    parsed[k] = None
            rows.append(parsed)
    return rows


def backbone_tag(rf2_tag):
    """Map an RF2 design tag back to its RFdiffusion backbone tag.

    e.g. 'samples_design_3_dldesign_1_best' -> 'samples_design_3'
    """
    stem = rf2_tag
    if stem.endswith("_best"):
        stem = stem[: -len("_best")]
    if "_dldesign_" in stem:
        stem = stem.split("_dldesign_")[0]
    return stem


def extract_heavy_seq(pdb_path, chain="H"):
    """Return the one-letter sequence of the given chain from a PDB, or None."""
    try:
        from Bio import PDB
        from Bio.Data.PDBData import protein_letters_3to1
    except ImportError:
        return None
    if not os.path.exists(pdb_path):
        return None
    structure = PDB.PDBParser(QUIET=True).get_structure("d", pdb_path)
    model = structure[0]
    if chain not in model:
        return None
    seq = "".join(
        protein_letters_3to1.get(res.resname, "X")
        for res in model[chain]
        if res.id[0] == " "
    )
    return seq or None


def h3_length_from_seq(seq):
    """CDR-H3 length via ANARCI IMGT numbering. Counts insertion codes.

    Returns an int, or None if ANARCI can't number the sequence.
    """
    if not seq:
        return None
    try:
        from anarci import anarci
    except ImportError:
        return None
    try:
        results, _, _ = anarci([("d", seq)], scheme="imgt", output=False)
    except Exception:
        return None
    if not results or results[0] is None:
        return None
    numbering = results[0][0][0]   # list of ((pos, ins), aa)
    count = 0
    for (pos, _ins), aa in numbering:
        if H3_IMGT_START <= pos <= H3_IMGT_END and aa != "-":
            count += 1
    return count if count > 0 else None


def get_h3_length(design_tag, batch_dir):
    """H3 length for a design: find its PDB, extract seq, count via ANARCI."""
    pdb_path = os.path.join(batch_dir, "final_pdbs", design_tag + ".pdb")
    if not os.path.exists(pdb_path):
        hits = glob.glob(os.path.join(batch_dir, "final_pdbs", design_tag + "*.pdb"))
        if not hits:
            return None
        pdb_path = hits[0]
    seq = extract_heavy_seq(pdb_path)
    return h3_length_from_seq(seq)


def analyze_batch(batch_dir):
    """Join RF2 interface metrics with RFdiffusion geometry + H3 length."""
    rf2_sc = run_qvscorefile(os.path.join(batch_dir, "4_rf2.qv"))
    diff_sc = run_qvscorefile(os.path.join(batch_dir, "1_rfdiffusion.qv"))

    rf2_rows = parse_scorefile(rf2_sc)
    diff_rows = parse_scorefile(diff_sc)

    mindist_by_backbone = {}
    for r in diff_rows:
        tag = r.get("tag")
        if tag is not None:
            mindist_by_backbone[tag] = r.get("mindist")

    rows = []
    for r in rf2_rows:
        tag = r.get("tag")
        if tag is None:
            continue
        bb = backbone_tag(tag)
        mindist = mindist_by_backbone.get(bb)
        rows.append({
            "batch": os.path.basename(os.path.normpath(batch_dir)),
            "design": tag,
            "pred_lddt": r.get("pred_lddt"),
            "interaction_pae": r.get("interaction_pae"),
            "pae": r.get("pae"),
            "target_aligned_cdr_rmsd": r.get("target_aligned_cdr_rmsd"),
            "min_hotspot_dist": round(mindist, 2) if mindist is not None else None,
            "h3_len": get_h3_length(tag, batch_dir),
        })
    return rows


def is_promising(r):
    return (
        r["pred_lddt"] is not None and r["pred_lddt"] >= PLDDT_MIN
        and r["interaction_pae"] is not None and r["interaction_pae"] <= IPAE_MAX
        and r["min_hotspot_dist"] is not None and r["min_hotspot_dist"] <= DIST_MAX
    )


def fmt(v, nd=2):
    return f"{v:.{nd}f}" if isinstance(v, float) else str(v)


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 analyze_batch_local.py <batch_dir> [<batch_dir> ...]")
        sys.exit(1)

    all_rows = []
    for batch_dir in sys.argv[1:]:
        rows = analyze_batch(batch_dir)
        all_rows.extend(rows)

        scored = [r for r in rows if r["pred_lddt"] is not None]
        ranked = sorted(
            scored,
            key=lambda r: (r["interaction_pae"] if r["interaction_pae"] is not None else 999,
                           -(r["pred_lddt"] or 0)),
        )
        promising = [r for r in rows if is_promising(r)]

        print("=" * 78)
        print(f"Batch: {os.path.basename(os.path.normpath(batch_dir))}")
        print(f"  Designs analysed: {len(rows)}")
        print(f"  Promising (pLDDT>={PLDDT_MIN}, iPAE<={IPAE_MAX}, dist<={DIST_MAX}A): "
              f"{len(promising)}")
        if rows:
            print(f"  Success rate:     {100*len(promising)/len(rows):.1f}%")
        print("-" * 78)
        print(f"  {'design':<42} {'pLDDT':>6} {'iPAE':>6} {'dist':>6} {'H3':>4}")
        for r in ranked[:15]:
            print(f"  {r['design']:<42} "
                  f"{fmt(r['pred_lddt'],3):>6} "
                  f"{fmt(r['interaction_pae']):>6} "
                  f"{fmt(r['min_hotspot_dist']):>6} "
                  f"{str(r['h3_len']):>4}")
        print()

    csv_path = "batch_analysis_local.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "batch", "design", "pred_lddt", "interaction_pae", "pae",
            "target_aligned_cdr_rmsd", "min_hotspot_dist", "h3_len"])
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"Wrote combined results to {csv_path}")

    print("\n" + "=" * 78)
    print("Batch comparison:")
    for batch_dir in sys.argv[1:]:
        bname = os.path.basename(os.path.normpath(batch_dir))
        brows = [r for r in all_rows if r["batch"] == bname]
        prom = [r for r in brows if is_promising(r)]
        if brows:
            ipaes = [r["interaction_pae"] for r in brows if r["interaction_pae"] is not None]
            mean_ipae = sum(ipaes) / len(ipaes) if ipaes else float("nan")
            median_ipae = sorted(ipaes)[len(ipaes)//2] if ipaes else float("nan")
            print(f"  {bname}: {len(prom)}/{len(brows)} promising "
                  f"({100*len(prom)/len(brows):.1f}%), "
                  f"mean iPAE {mean_ipae:.2f}, median iPAE {median_ipae:.2f}")


if __name__ == "__main__":
    main()
