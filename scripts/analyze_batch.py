#!/usr/bin/env python3
"""
Analyze and rank RFantibody design batches.

RF2 writes its metrics as `SCORE <name>: <value>` lines at the END of each
_best.pdb (after the coordinates). This script parses those, joins in the
RFdiffusion .trb geometry (min hotspot-loop distance, H3 length), ranks
designs, and reports per-batch success rates so you can compare H3-range
batches.

Usage:
    python3 analyze_batch.py outputs/batchA outputs/batchB

Key metrics used:
  - pred_lddt        : fold confidence (0-1, higher better)
  - interaction_pae  : interface predicted aligned error (lower better; the
                       real binding-confidence signal)
  - min_hotspot_dist : closest CDR-loop approach to a hotspot (from .trb)
  - H3 length        : from .trb (for the batch comparison)
"""

import sys
import os
import glob
import csv
import numpy as np

# ---- Thresholds for a "promising" design --------------------------------
PLDDT_MIN = 0.85          # pred_lddt above this
IPAE_MAX = 10.0           # interaction_pae below this (lower = better interface)
DIST_MAX = 8.0            # min hotspot-loop distance below this (Angstrom)


def parse_scores(pdb_path):
    """Parse all 'SCORE <name>: <value>' lines from an RF2 output PDB.

    IMPORTANT: SCORE lines appear at the END of the file (after coordinates),
    so we scan the whole file rather than stopping at the first ATOM line.
    """
    scores = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("SCORE "):
                parts = line.strip().split()
                if len(parts) >= 3:
                    name = parts[1].rstrip(":")
                    try:
                        scores[name] = float(parts[2])
                    except ValueError:
                        pass
    return scores


def find_trb_geometry(batch_dir, design_name):
    """Look up min hotspot-loop distance and H3 length from the backbone .trb."""
    stem = design_name
    if stem.endswith("_best"):
        stem = stem[:-len("_best")]
    if "_dldesign_" in stem:
        stem = stem.split("_dldesign_")[0]
    trb_path = os.path.join(batch_dir, "raw", stem + ".trb")
    if not os.path.exists(trb_path):
        return None, None
    d = np.load(trb_path, allow_pickle=True)
    return float(d["mindist"]), int(d["H3_len"])


def analyze_batch(batch_dir):
    rf2_dir = os.path.join(batch_dir, "rf2")
    rows = []
    for pdb_path in sorted(glob.glob(os.path.join(rf2_dir, "*_best.pdb"))):
        name = os.path.basename(pdb_path)[:-4]
        scores = parse_scores(pdb_path)
        mindist, h3len = find_trb_geometry(batch_dir, name)
        rows.append({
            "batch": os.path.basename(batch_dir),
            "design": name,
            "pred_lddt": scores.get("pred_lddt"),
            "interaction_pae": scores.get("interaction_pae"),
            "pae": scores.get("pae"),
            "target_aligned_cdr_rmsd": scores.get("target_aligned_cdr_rmsd"),
            "min_hotspot_dist": round(mindist, 2) if mindist is not None else None,
            "h3_len": h3len,
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
        print("Usage: python3 analyze_batch.py <batch_dir> [<batch_dir> ...]")
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
        print(f"Batch: {os.path.basename(batch_dir)}")
        print(f"  Designs analysed: {len(rows)}")
        print(f"  Promising (pLDDT>={PLDDT_MIN}, iPAE<={IPAE_MAX}, dist<={DIST_MAX}A): "
              f"{len(promising)}")
        if rows:
            print(f"  Success rate:     {100*len(promising)/len(rows):.1f}%")
        print("-" * 78)
        print(f"  {'design':<38} {'pLDDT':>6} {'iPAE':>6} {'dist':>6} {'H3':>4}")
        for r in ranked[:15]:
            print(f"  {r['design']:<38} "
                  f"{fmt(r['pred_lddt'],3):>6} "
                  f"{fmt(r['interaction_pae']):>6} "
                  f"{fmt(r['min_hotspot_dist']):>6} "
                  f"{str(r['h3_len']):>4}")
        print()

    csv_path = "batch_analysis.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "batch", "design", "pred_lddt", "interaction_pae", "pae",
            "target_aligned_cdr_rmsd", "min_hotspot_dist", "h3_len"])
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"Wrote combined results to {csv_path}")

    print("\n" + "=" * 78)
    print("Batch comparison (for the H3-length note):")
    for batch_dir in sys.argv[1:]:
        bname = os.path.basename(batch_dir)
        brows = [r for r in all_rows if r["batch"] == bname]
        prom = [r for r in brows if is_promising(r)]
        if brows:
            ipaes = [r["interaction_pae"] for r in brows if r["interaction_pae"] is not None]
            mean_ipae = np.mean(ipaes) if ipaes else float("nan")
            print(f"  {bname}: {len(prom)}/{len(brows)} promising "
                  f"({100*len(prom)/len(brows):.1f}%), mean iPAE {mean_ipae:.2f}")


if __name__ == "__main__":
    main()
