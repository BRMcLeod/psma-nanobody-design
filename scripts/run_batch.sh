#!/bin/bash
# =============================================================================
# RFantibody batch design funnel for PSMA nanobodies
#
# Runs: RFdiffusion (generate) -> filter backbones by hotspot distance ->
#       ProteinMPNN (sequences) -> RF2 (predict/validate)
#
# Usage:
#   bash run_batch.sh <batch_name> <h3_range> <num_designs>
#
# Examples:
#   bash run_batch.sh batchA "H3:5-15" 96
#   bash run_batch.sh batchB "H3:8-20" 96
#
# Everything except the H3 range is held constant so the two batches are
# directly comparable.
# =============================================================================

set -e  # stop if any command fails

# ---- Arguments -----------------------------------------------------------
BATCH_NAME=${1:?"Provide a batch name, e.g. batchA"}
H3_RANGE=${2:?"Provide an H3 range, e.g. H3:5-15"}
NUM_DESIGNS=${3:-96}

# ---- Fixed parameters (identical across batches) -------------------------
TARGET="inputs/psma_target_clean.pdb"
FRAMEWORK="scripts/examples/example_inputs/h-NbBCII10.pdb"
HOTSPOTS="T171,T173,T199,T223,T341"
LOOPS="H1:7,H2:6,${H3_RANGE}"
MPNN_SEQS=2                # sequences per surviving backbone
DIST_CUTOFF=8.0           # keep backbones with min hotspot-loop dist below this (Angstrom)

# ---- Output layout -------------------------------------------------------
OUTDIR="outputs/${BATCH_NAME}"
RAW="${OUTDIR}/raw"           # all RFdiffusion backbones
KEPT="${OUTDIR}/kept"         # backbones passing the distance filter
MPNN="${OUTDIR}/mpnn"         # ProteinMPNN sequence designs
RF2="${OUTDIR}/rf2"           # RF2 predictions
mkdir -p "$RAW" "$KEPT" "$MPNN" "$RF2"

echo "=========================================================="
echo "Batch:        $BATCH_NAME"
echo "H3 range:     $H3_RANGE"
echo "Num designs:  $NUM_DESIGNS"
echo "Loops:        $LOOPS"
echo "Dist cutoff:  ${DIST_CUTOFF} A"
echo "=========================================================="

# ---- Stage 1: RFdiffusion (generate backbones) ---------------------------
echo "[Stage 1] RFdiffusion: generating $NUM_DESIGNS backbones..."
rfdiffusion \
    -t "$TARGET" \
    -f "$FRAMEWORK" \
    -o "${RAW}/${BATCH_NAME}" \
    -n "$NUM_DESIGNS" \
    -l "$LOOPS" \
    -h "$HOTSPOTS"

echo "[Stage 1] Done. Backbones in $RAW"

# ---- Stage 2: filter backbones by hotspot-to-loop distance ---------------
# Reads the min-distance metric from each .trb and copies passing PDBs to KEPT.
echo "[Stage 2] Filtering backbones by min hotspot-loop distance < ${DIST_CUTOFF} A..."
python3 - "$RAW" "$KEPT" "$DIST_CUTOFF" <<'PYEOF'
import sys, glob, os, shutil
import numpy as np

raw_dir, kept_dir, cutoff = sys.argv[1], sys.argv[2], float(sys.argv[3])
kept, total = 0, 0
for trb_path in glob.glob(os.path.join(raw_dir, "*.trb")):
    total += 1
    d = np.load(trb_path, allow_pickle=True)
    mindist = float(d["mindist"])
    if mindist < cutoff:
        pdb_path = trb_path[:-4] + ".pdb"
        if os.path.exists(pdb_path):
            shutil.copy(pdb_path, kept_dir)
            kept += 1
print(f"  Kept {kept} of {total} backbones (min dist < {cutoff} A)")
PYEOF

NUM_KEPT=$(ls "$KEPT"/*.pdb 2>/dev/null | wc -l)
if [ "$NUM_KEPT" -eq 0 ]; then
    echo "[Stage 2] No backbones passed the filter. Stopping."
    exit 0
fi
echo "[Stage 2] Done. $NUM_KEPT backbones passed."

# ---- Stage 3: ProteinMPNN (sequence design on survivors) -----------------
echo "[Stage 3] ProteinMPNN: designing sequences for $NUM_KEPT backbones..."
proteinmpnn \
    -i "$KEPT" \
    -o "$MPNN" \
    -l "H1,H2,H3" \
    -n "$MPNN_SEQS"
echo "[Stage 3] Done. Sequences in $MPNN"

# ---- Stage 4: RF2 (structure prediction / validation) --------------------
echo "[Stage 4] RF2: predicting structures..."
rf2 \
    -i "$MPNN" \
    -o "$RF2"
echo "[Stage 4] Done. Predictions in $RF2"

echo "=========================================================="
echo "Batch $BATCH_NAME complete."
echo "  Raw backbones:   $RAW"
echo "  Passed filter:   $KEPT  ($NUM_KEPT)"
echo "  RF2 predictions: $RF2"
echo "Next: run analyze_batch.py to rank the RF2 designs."
echo "=========================================================="
