# PSMA Nanobody Design: Methods and Results

## Overview

A de novo nanobody (VHH) design campaign targeting the apical domain of
prostate-specific membrane antigen (PSMA), using the RFantibody pipeline
(RFdiffusion, ProteinMPNN, RoseTTAFold2). The work ran in two phases:

- **Phase 1** was a 4-design proof-of-concept to validate the full workflow end
  to end, shake out failure modes, and build the judgment needed to evaluate
  designs. It confirmed the pipeline runs, caught a target-renumbering bug, and
  produced a hint that CDR-H3 length matters.
- **Phase 2** scaled to a 192-design, two-batch experiment that turned that hint
  into a controlled test of CDR-H3 length, and produced a verified headline
  design.

The pipeline was subsequently rebuilt to run locally on a Blackwell GPU and the
Phase 2 campaign re-run to confirm the local build reproduces the cloud results;
that validation is documented at the end.

The target, epitope, and target-preparation steps below are shared across both
phases. The pipeline is shared too; only the scale and the H3 sampling range
differ between phases.

## Target and epitope

Target: PSMA (GCPII / FOLH1), extracellular domain, from PDB 4NGM (~1.8 A
resolution). PSMA is a validated prostate cancer target, overexpressed on
prostate cancer cells and the target of approved radioligand therapies.

Epitope: apical domain, the mapped epitope region of the clinical anti-PSMA
antibody J591. Chosen because it is a validated antibody-accessible surface,
projects away from the membrane on a live cell, and is distinct from the
active-site pocket (which is a small-molecule target and shared with the GCPIII
homolog, raising cross-reactivity risk).

Hotspot residues (original 4NGM numbering, chain A): Glu171, Asp173, Lys199,
Lys223, Lys341. Selected for solvent exposure, tight spatial clustering (~13 A
Ca-Ca across the set), and location on ordered secondary structure (Lys223 sits
in a beta turn between two helices). Four of the five are charged, favouring
specific surface contacts.

## Target preparation

1. Cropped the 4NGM extracellular domain to ~235 residues (residues 117-351)
   around the epitope, preserving local secondary structure and keeping all five
   hotspots well within the crop with margin.
2. Relabelled the target chain to T (the chain ID RFantibody expects for
   targets).
3. Removed alternate conformations (alt-locs), keeping only the primary
   conformer. This step was necessary: an alt-loc on a single residue produced a
   degenerate (all-zero) backbone frame that crashed RFdiffusion deterministically
   with a "Non-positive determinant in rotation matrix" error. Stripping alt-locs
   resolved it.

## Pipeline and parameters

Framework: h-NbBCII10 (the validated VHH/nanobody framework provided with
RFantibody).

The funnel, parameterized in `scripts/run_batch.sh`:

1. **RFdiffusion** (backbone design): heavy-chain loops only (H1, H2, H3), since
   a VHH has no light chain. Hotspots set to the five residues above on chain T.
2. **8 A hotspot-distance filter**: backbones whose designed loops come within
   8 A of the hotspots are kept; the rest are discarded before sequence design.
   This is the geometric gate that concentrates compute on backbones that
   actually reach the epitope.
3. **ProteinMPNN** (sequence design): H1, H2, H3 only; 2 sequences per surviving
   backbone.
4. **RoseTTAFold2** (structure prediction / validation): default 10 recycles.

`scripts/analyze_batch.py` parses the RF2 SCORE lines (which are footer lines at
the end of each `_best.pdb`, not the B-factor column), pulls mindist and H3
length from the `.trb` metadata, ranks designs by interaction_pae, and reports
per-batch pass rates.

Compute: run on a cloud RTX 4090 (Ada, sm_89) via RunPod. The local RTX 5060 Ti
(Blackwell, sm_120) is not supported by RFantibody's default CUDA 11.8 container,
so cloud compute on an Ada GPU was used for the original campaign. A later local
rebuild (below) removed this dependence.

---

# Phase 1: proof-of-concept (4 designs)

A minimal run to validate the pipeline end to end and develop evaluation
judgment before committing compute to scale.

Scale: 2 backbones x 2 sequences = 4 designs. H3 range initially set wide (6-20)
to reflect the longer CDR-H3 seen in camelid VHHs, then tightened to 6-16: the
longest H3 lengths (17-20) contributed to the degenerate-frame instability, and
tightening improved stability.

## Phase 1 results

RF2 confidence (pLDDT) for the four designs:

| Design               | pLDDT | Backbone reached epitope? |
|----------------------|-------|---------------------------|
| psma_nb_1_dldesign_0 | 0.921 | Yes (~8 A min dist)       |
| psma_nb_0_dldesign_0 | 0.913 | No (~18 A min dist)       |
| psma_nb_1_dldesign_1 | 0.912 | Yes (~8 A min dist)       |
| psma_nb_0_dldesign_1 | 0.905 | No (~18 A min dist)       |

All four designs fold confidently (pLDDT > 0.90), but pLDDT reflects fold
confidence, not binding. Hotspot-to-loop distance separates the two backbones:
psma_nb_1 (H3 length 11) placed its loops near the epitope (~8 A), while
psma_nb_0 (H3 length 20) missed (~18 A). This is a concrete illustration that
high pLDDT does not imply productive engagement.

Visual inspection of the two epitope-reaching designs (psma_nb_1):

- psma_nb_1_dldesign_0: dock is non-productive. The framework (FR2/FR3), not the
  CDRs, forms most of the target contact.
- psma_nb_1_dldesign_1: better. The CDRs approach the epitope, though the overall
  angle of approach is suboptimal. Plausibly weak but detectable binding at best.

Since dldesign_0 and dldesign_1 share the same backbone, the difference between a
framework-mediated dock and a CDR-mediated approach comes from the ProteinMPNN
sequence alone. This highlights that the sequence design step affects productive
engagement, not just the backbone.

## What Phase 1 established

- The pipeline runs end to end and produces confidently folded designs.
- None of the four is a convincing binder, which is expected at this scale.
- A target-renumbering bug exists in the outputs (see the QC section below),
  caught during Phase 1 analysis.
- A hint on H3 length: the one long-H3 backbone (length 20) missed while the
  shorter (length 11) reached. With only two backbones this is not a result, but
  it is a testable hypothesis, and it set up Phase 2.

---

# Phase 2: the two-batch H3-length experiment

Hypothesis: natural camelid VHHs carry long CDR-H3 loops, but does sampling
longer H3s actually improve de novo design success against this epitope, or does
the extra loop length just explore more space and miss more often? Phase 1
gestured at the second possibility with a single backbone; Phase 2 tests it
properly.

Design: two batches of 96 designs each, identical in every parameter except the
CDR-H3 length range sampled by RFdiffusion.

- batchA: H3 range 5-15
- batchB: H3 range 8-20

## Phase 2 results

| Batch  | H3 range | RF2 outputs (passed 8 A filter) | Best iPAE | Median iPAE |
|--------|----------|---------------------------------|-----------|-------------|
| batchA | 5-15     | 53 / 96                         | 13.89     | 18.65       |
| batchB | 8-20     | 64 / 96                         | 5.27      | 18.33       |

batchB won on both the geometric pass rate (64 vs 53 of 96 reaching the epitope)
and interface confidence at the top end. The medians are close (18.33 vs 18.65),
so most designs in both batches are poor, as expected for de novo design. The
difference is in the tail that matters: batchA produced zero designs under the
iPAE < 10 "promising" threshold, while batchB produced six, including a tight
cluster of five under iPAE 8.

### The batchB top cluster

| Design               | iPAE | pLDDT | target-aln CDR RMSD (A) | min hotspot dist (A) | H3 len |
|----------------------|------|-------|-------------------------|----------------------|--------|
| batchB_90_dldesign_1 | 5.27 | 0.91  | 4.06                    | 4.69                 | 16     |
| batchB_27_dldesign_0 | 6.97 | 0.90  | 9.41                    | 4.52                 | 12     |
| batchB_11_dldesign_0 | 7.10 | 0.91  | 9.92                    | 4.12                 | 13     |
| batchB_20_dldesign_0 | 7.41 | 0.91  | 17.53                   | 4.55                 | 10     |
| batchB_29_dldesign_1 | 7.66 | 0.91  | 1.28                    | 5.12                 | 11     |

For contrast, batchA's single best design across all 96:

| Design               | iPAE  | pLDDT | target-aln CDR RMSD (A) | min hotspot dist (A) | H3 len |
|----------------------|-------|-------|-------------------------|----------------------|--------|
| batchA_35_dldesign_1 | 13.89 | 0.91  | 4.89                    | 4.51                 | 12     |

All six curated designs (five batchB, one batchA) are in `designs/top/`.

## Interpretation

The distributional shift is the finding: widening the H3 range upward moved both
the pass rate and, more importantly, the best-case interface confidence in the
right direction. The structural rationale is that a longer H3 has more reach to
extend from the framework into the epitope.

Two honest caveats on over-reading this:

- Within the winning cluster, H3 lengths span 10-16, and the second-best design
  (batchB_27, iPAE 6.97) has an H3 of only 12. So H3 length does not cleanly
  predict quality at the single-design level. Combined with Phase 1, where the
  long-H3 backbone missed, the honest read is that the effect is real as a shift
  between sampling ranges, not as a per-design rule. This is exactly why a
  controlled batch comparison was needed rather than trusting a single backbone.
- target_aligned_cdr_rmsd varies widely across the cluster (1.28 to 17.53 A).
  This measures how far RF2's predicted CDR pose sits from the designed backbone.
  A low iPAE with a high CDR RMSD (e.g. batchB_20) means RF2 is confident about
  an interface it has repositioned, which is a weaker signal than a low iPAE with
  a low RMSD (e.g. batchB_29). The headline design sits in between at ~4 A.

## Headline design

**batchB_90_dldesign_1**: iPAE 5.27, pLDDT 0.91, min hotspot distance 4.69 A, H3
length 16, target-aligned CDR RMSD ~4 A.

Visual inspection in PyMOL confirms the entire CDR-H3 loop folds directly over
all five hotspot residues, reading Glu-Asp-Lys-Lys-Lys (EDKKK) along the
epitope. This is genuine CDR-mediated engagement of the intended epitope: the
loop is doing the binding, not the framework, and it is contacting the residues
that were actually targeted. Metrics and structure agree, which is the
combination that matters. The interface figure is at
`designs/top/tophit_batchB90_interface.png`.

The ~4 A target-aligned CDR RMSD is the caveat: RF2 predicts a pose slightly
shifted from the designed backbone while remaining confident in the interface.
This is a strong in silico hypothesis, not a validated binder.

One developability caveat sits alongside the interface result. Sequence-liability scanning (see "Sequence-liability triage" below) flags an NG motif in CDR2 of this design, at IMGT position 62. NG is the fastest-deamidating motif in antibodies and the one that dominates degradation in clinical-stage mAbs (Lu et al. 2018), and it sits in a CDR where modification is most likely to affect binding. It does not unseat batchB_90 as the best design by interface confidence, but in a real campaign this is exactly the liability to either engineer out (N->Q, or G->A at the +1 position) or carry forward with a stability assay attached. The strongest design by one axis is not automatically the cleanest by another, which is the tension a triage funnel exists to surface.       

---

## QC finding: target renumbering across outputs

Caught during Phase 1 analysis and handled throughout Phase 2. The RF2 output
structures renumber the target chain inconsistently between designs (e.g. chain
T starting at a different residue number in different outputs, versus 117 in the
input). The target sequence and fold are identical across all outputs (0
sequence differences; target Ca-Ca RMSD ~0.05 A between designs), but the
numbering offsets mean the original hotspot residue numbers point to the wrong
physical residues in the outputs.

This was caught by noticing that a single labelled residue number had a
different amino acid identity across two outputs, which is impossible for an
identical target unless the numbering had shifted. Hotspots were then re-mapped
onto each output by structural alignment to the correctly-numbered input and
proximity selection, and all hotspot-distance analysis was done on the correctly
mapped residues.

Lesson: residue numbering is not guaranteed to be preserved through this
pipeline. Any position-based analysis (hotspot contacts, interface residues)
must be mapped through structural alignment rather than trusting output residue
numbers. This is baked into the Phase 2 analysis approach: hotspots are located
by proximity and alignment, not by residue number.

---

# Local build reproduction (Blackwell / CUDA 12.8)

## Motivation

The original campaign ran on rented cloud GPUs because RFantibody's default
container targets CUDA 11.8, which does not support the local card's Blackwell
architecture (compute capability sm_120): designs failed with "CUDA error: no
kernel image is available for execution on the device." The stack was rebuilt
from source against CUDA 12.8 to run locally on an RTX 5060 Ti. The full build
account (dependency handling, the DGL from-source compile, an ANARCI install
requiring a numpy pin, the quiver-native pipeline adaptation) is kept in the
separate rfantibody-blackwell project notes.

The local build is quiver-native (structures and scores flow between stages in
`.qv` archive files rather than loose PDB/TRB files), so the pipeline and
analysis scripts were re-implemented accordingly:
`scripts/run_batch_local.sh`, `scripts/analyze_batch_local.py`, and
`scripts/cdr_annotate.py` (the last maps IMGT CDR ranges to real PDB residue
numbers via ANARCI).

## Validation approach

A rebuilt stack on a different CUDA numerical path and a newer RFantibody version
could in principle shift results, so the local build was not assumed equivalent.
The Phase 2 batchB campaign was re-run locally with identical parameters (H3
range 8-20, 96 designs, same target, same five hotspots, same 8 A filter, 2
sequences per backbone, RF2 at 10 recycles) and compared against the stored cloud
batchB. The comparison deliberately weights the stable high-N statistics (median
iPAE, pass rate) over the noisy best-of-N value.

## Validation results

| Metric               | RunPod (cu118) | Local (CUDA 12.8, sm_120) |
|----------------------|----------------|---------------------------|
| Median iPAE          | 18.33          | 18.50                     |
| Mean iPAE            | ~18            | 18.16                     |
| 8 A filter pass rate | 64 / 96        | 54 / 96                   |
| Best iPAE            | 5.27           | 4.78                      |

The median interface-pAE is essentially identical (18.50 vs 18.33), so the local
build reproduces the campaign's interface-quality distribution, which is the
statistic that actually matters. The pass-rate difference (54 vs 64 of 96) sits
within binomial sampling noise on the softer geometric filter and does not
reflect a quality shift. The target-renumbering behaviour persisted identically
on the local build and was handled the same way (structural mapping, not trusting
output residue numbers), which is itself a point of consistency between the two
builds.

## Structural triage of the local run

The same scrutiny applied to the original campaign was applied here, and it
mattered. The top design by iPAE (samples_design_17, iPAE 4.78) was a
framework-mediated dock, not a productive one: it carried a large target-aligned
CDR RMSD (27.8 A), the numerical fingerprint of an interface RF2 is confident
about but which is not CDR-mediated, and visual inspection confirmed the CDRs
sit perpendicular to the epitope with the framework making the contact. A pure
iPAE ranking would have mistaken it for the best design.

Ranking instead on iPAE together with pose agreement (low target-aligned CDR
RMSD) surfaced a genuine CDR-mediated candidate:

| Design                 | iPAE | pLDDT | target-aln CDR RMSD (A) | min hotspot dist (A) | H3 len |
|------------------------|------|-------|-------------------------|----------------------|--------|
| samples_design_0_dld_0 | 8.81 | 0.91  | 6.42                    | 4.47                 | 21     |

Visual inspection confirms this design's CDR-H3 (a long, 21-residue loop) engages
the epitope at a natural angle of approach, contacting the hotspots through the
CDRs rather than the framework. The interface figure is at
`designs/top/local_validation_design0_interface.png`.

This candidate is presented as evidence that the local build produces the same
kind of output as the cloud build and stands up to the same structural scrutiny,
not as a new scientific result. The Phase 2 finding (H3 length as a
distributional effect; batchB_90 as the campaign headline) is unchanged.

## Cost and practical outcome

The original 96-design campaign cost roughly a few USD of rented cloud time. The
equivalent local run costs a fraction of that in electricity, but the meaningful
gain is not the money: it is the removal of per-run cost and cloud-provisioning
friction, which allows unconstrained iteration and overnight runs. For very large
campaigns, rented cloud GPUs remain faster and are the sensible tool; the local
build is best for development, iteration, and modest batches. Knowing that split
is part of the outcome.

## Honest assessment and limitations

- This is an in silico campaign. Even the best design is a hypothesis, not a
  binder. de novo nanobody design against a defined epitope typically needs
  hundreds to thousands of designs and hard interface filtering to yield a real
  candidate; 192 designs producing one strong in silico hit is consistent with
  the method's known low per-design hit rate.
- iPAE is the primary interface-confidence signal here, but a single in silico
  metric is not validation, and, as the local triage shows, it can rank a
  framework dock above a productive one. Pose agreement and visual inspection are
  necessary alongside it.
- The headline design's ~4 A target-aligned CDR RMSD means RF2 repositions the
  CDR somewhat relative to the design; the interface is confident but not a
  perfect backbone match.
- The local reproduction validates the build against a matched PSMA benchmark at
  matched settings. It does not guarantee identical behaviour on an arbitrary new
  target, and it spans two changes at once (CUDA stack and RFantibody version),
  so it is a practical equivalence check, not a controlled single-variable one.
- No wet-lab validation. No developability assessment (aggregation, exposed
  hydrophobics, unpaired cysteines, glycosylation motifs). No screening against
  the GCPIII homolog for cross-reactivity.

## Sequence-liability triage

A per-design developability screen was added downstream of the interface funnel, implemented in `scripts/liability_scan.py`. The interface metrics (iPAE, pose agreement, hotspot contact) rank a design by whether it binds; they say nothing about whether the sequence is manufacturable and stable. Those are separate axes, and a design can score well on one while carrying real problems on the other.

The screen splits each design's heavy chain into IMGT regions with ANARCI, then scans for known chemical-liability motifs: N-linked glycosylation sequons (NxS/NxT), asparagine deamidation (NG high-risk, NS/NT/NN/NH moderate), aspartate isomerization (D[GSDNTH]), Asp-Pro fragmentation, and methionine/tryptophan oxidation. Motifs are scored ordinally and weighted 5x in CDRs versus framework, on the basis that a CDR liability is both more solvent-exposed and more likely to affect binding (motif set and weighting after Lu et al. 2018, mAbs 11:45-57, and the LAP/Adimab IMGT liability conventions). The tool reports both a flagged inventory (every motif, its position, and region) and per-region plus total scores.

Two design choices matter for interpretation. Oxidation-prone Met and Trp in the framework are flagged but scored zero: they are usually buried and often structurally essential (the conserved Trp41, for instance), so they are surfaced as reminders rather than penalties. And the canonical VHH disulfide cysteines (IMGT 23 and 104) are recognized as expected rather than counted as unpaired-cysteine liabilities.

The screen also tracks the camelid FR2 hallmark residues (IMGT 42, 49, 50, 52; Vincke et al. 2009), which govern VHH solubility and contribute to antigen binding and CDR3 pre-organization. Across all ten curated designs these read FGLA and are identical, confirming that the h-NbBCII10 scaffold retains its camelid hallmark character (F42 and A52 in particular) and that sequence design did not drift these binding-relevant positions.

Applied to the curated set, the screen makes the interface-versus-developability tension concrete: the two long-H3 Phase 1 designs (psma_nb_0, H3 length 22) carry the heaviest liability load, while the headline batchB_90 sits mid-pack, its main liability being the CDR2 NG motif noted above. This is the intended function of a triage layer: to rank on manufacturability independently of binding, so the two signals can be weighed against each other rather than conflated.

## What a real campaign would add

- Scale further and filter hard on interface pAE and CDR-mediated contact, then
  inspect survivors individually.
- Orthogonal structure prediction (e.g. ColabFold / AlphaFold) of the top
  candidates for a second, independent confidence signal.
- Rosetta interface scoring (dG_separated, shape complementarity) on the
  survivors.
- Cross-reactivity screening against GCPIII homoloh before any synthesis.
- Wet-lab expression and binding assays (e.g. SPR/BLI) as the only real
  validation.

## Files

- `targets/`: input structure (4NGM crop), cleaned target, hotspot PyMOL session
- `scripts/run_batch.sh` / `analyze_batch.py`: original RunPod funnel and analysis
- `scripts/run_batch_local.sh` / `analyze_batch_local.py`: local quiver-native
  funnel and analysis (ANARCI-based IMGT CDR-H3 counting)
- `scripts/cdr_annotate.py`: maps IMGT CDR ranges to real PDB residue numbers,
  emits a reference table and PyMOL selections
- `scripts/liability_scan.py`: sequence-liability triage (IMGT region split via ANARCI, motif scan, per-region and total scores, FR2 hallmark tracking)
- `designs/poc/`: the four Phase 1 proof-of-concept designs and figure
- `designs/top/`: curated Phase 2 designs, the headline interface figure, and the
  local-build validation figure
- `results/batch_analysis.csv`: full per-design metrics for both batches
- `docs/`: this document 


Raw per-design outputs (hundreds of PDBs per batch, ProteinMPNN intermediates,
trajectory folders) are kept locally and excluded from the repo via
`.gitignore`; the curated designs and the full metrics CSV are committed so the
results remain reproducible and inspectable without the bulk.
