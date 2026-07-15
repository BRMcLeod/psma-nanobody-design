#!/bin/bash
# =============================================================================
# RFantibody batch design funnel for PSMA nanobodies  (LOCAL / Blackwell build)
#
# Runs: RFdiffusion (generate) -> filter backbones by hotspot distance ->
#       ProteinMPNN (sequences) -> RF2 (predict/validate) -> extract PDBs
#
# Quiver-native: each stage passes a single .qv file to the next, rather than
# directories of loose PDB/TRB files (that was the RunPod build's behaviour).
#
# Usage:
#   bash run_batch_local.sh <batch_name> <h3_range> <num_designs> <target.pdb> <hotspots>
#
# Examples:
#   bash run_batch_local.sh batchA "H3:5-15" 96 targets/psma_cropped.pdb "T171,T173,T199,T223,T341"
#   bash run_batch_local.sh batchB "H3:8-20" 96 targets/psma_cropped.pdb "T171,T173,T199,T223,T341"
#
# Everything except the H3 range is held constant so batches are comparable.
# =============================================================================

set -e  # stop if any command fails

# ---- Arguments -----------------------------------------------------------
BATCH_NAME=${1:?"Provide a batch name, e.g. batchA"}
H3_RANGE=${2:?"Provide an H3 range, e.g. H3:5-15"}
NUM_DESIGNS=${3:-96}
TARGET=${4:?"Provide a target PDB, e.g. targets/psma_cropped.pdb"}
HOTSPOTS=${5:?"Provide hotspots, e.g. T171,T173,T199,T223,T341"}

# ---- Fixed parameters (identical across batches) -------------------------
FRAMEWORK="scripts/examples/example_inputs/h-NbBCII10.pdb"
LOOPS="H1:7,H2:6,${H3_RANGE}"
MPNN_SEQS=2                # sequences per surviving backbone
DIST_CUTOFF=8.0           # keep backbones with min hotspot-loop dist below this (Angstrom)

# ---- Output layout (quivers are FILES, not directories) ------------------
OUTDIR="outputs/${BATCH_NAME}"
mkdir -p "$OUTDIR"

RAW_QV="${OUTDIR}/1_rfdiffusion.qv"       # all RFdiffusion backbones
SCOREFILE="${OUTDIR}/1_rfdiffusion.sc"    # qvscorefile writes here (same path, .sc)
KEPT_QV="${OUTDIR}/2_kept.qv"             # backbones passing the distance filter
MPNN_QV="${OUTDIR}/3_proteinmpnn.qv"      # ProteinMPNN sequence designs
RF2_QV="${OUTDIR}/4_rf2.qv"               # RF2 predictions
FINAL_PDBS="${OUTDIR}/final_pdbs"         # extracted PDBs for analysis/inspection
mkdir -p "$FINAL_PDBS"

echo "=========================================================="
echo "Batch:        $BATCH_NAME"
echo "H3 range:     $H3_RANGE"
echo "Num designs:  $NUM_DESIGNS"
echo "Target:       $TARGET"
echo "Hotspots:     $HOTSPOTS"
echo "Loops:        $LOOPS"
echo "Dist cutoff:  ${DIST_CUTOFF} A"
echo "=========================================================="

# ---- Stage 1: RFdiffusion (generate backbones) ---------------------------
echo "[Stage 1] RFdiffusion: generating $NUM_DESIGNS backbones..."
rfdiffusion \
    -t "$TARGET" \
    -f "$FRAMEWORK" \
    -q "$RAW_QV" \
    -n "$NUM_DESIGNS" \
    -l "$LOOPS" \
    -h "$HOTSPOTS" \
    --no-trajectory
echo "[Stage 1] Done. Backbones in $RAW_QV"

# ---- Stage 2: filter backbones by hotspot-to-loop distance ---------------
# qvscorefile dumps QV_SCORE lines to $SCOREFILE (cols: mindist, averagemin, tag).
# awk selects tags whose mindist is under the cutoff; qvslice builds a new quiver
# containing only those designs.
echo "[Stage 2] Filtering backbones by min hotspot-loop distance < ${DIST_CUTOFF} A..."
qvscorefile "$RAW_QV"
awk -F'\t' -v cut="$DIST_CUTOFF" 'NR>1 && $1 < cut {print $3}' "$SCOREFILE" \
    | qvslice "$RAW_QV" > "$KEPT_QV"

NUM_KEPT=$(qvls "$KEPT_QV" | wc -l)
if [ "$NUM_KEPT" -eq 0 ]; then
    echo "[Stage 2] No backbones passed the filter. Stopping."
    exit 0
fi
echo "[Stage 2] Done. $NUM_KEPT backbones passed."

# ---- Stage 3: ProteinMPNN (sequence design on survivors) -----------------
echo "[Stage 3] ProteinMPNN: designing sequences for $NUM_KEPT backbones..."
proteinmpnn \
    -q "$KEPT_QV" \
    --output-quiver "$MPNN_QV" \
    -l "H1,H2,H3" \
    -n "$MPNN_SEQS"
echo "[Stage 3] Done. Sequences in $MPNN_QV"

# ---- Stage 4: RF2 (structure prediction / validation) --------------------
echo "[Stage 4] RF2: predicting structures..."
rf2 \
    -q "$MPNN_QV" \
    --output-quiver "$RF2_QV" \
    --num-recycles 10
echo "[Stage 4] Done. Predictions in $RF2_QV"

# ---- Stage 4b: extract final PDBs ----------------------------------------
echo "[Stage 4b] Extracting final PDBs..."
qvextract "$RF2_QV" -o "$FINAL_PDBS"    # verify flag with: qvextract --help
echo "[Stage 4b] Done. PDBs in $FINAL_PDBS"

echo "=========================================================="
echo "Batch $BATCH_NAME complete."
echo "  RFdiffusion quiver: $RAW_QV"
echo "  Passed filter:      $KEPT_QV  ($NUM_KEPT)"
echo "  RF2 quiver:         $RF2_QV"
echo "  Final PDBs:         $FINAL_PDBS"
echo "Next: run analyze_batch.py on the RF2 output to rank the designs."
echo "=========================================================="