# PSMA Nanobody Design: Methods and Results

## Overview

A proof-of-concept de novo nanobody (VHH) design campaign targeting the apical
domain of prostate-specific membrane antigen (PSMA), using the RFantibody
pipeline (RFdiffusion, ProteinMPNN, RoseTTAFold2). The goal was to build and run
the full computational workflow end to end against a chosen epitope, and to
develop the judgment needed to evaluate the resulting designs.

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

Step 1, RFdiffusion (backbone design):
- Loops designed: heavy-chain only (H1:7, H2:6, H3:6-20 then 6-16), since a VHH
  has no light chain.
- Hotspots: the five residues above, on chain T.
- H3 range set wider than a typical IgG (6-20) to reflect the longer CDR-H3 seen
  in camelid VHHs. Note: the longest H3 lengths (17-20) contributed to the
  degenerate-frame instability; tightening to 6-16 improved stability.

Step 2, ProteinMPNN (sequence design):
- Loops: H1, H2, H3 only.
- 2 sequences per backbone.

Step 3, RoseTTAFold2 (structure prediction / validation):
- Default 10 recycles.

Scale: 2 backbones x 2 sequences = 4 designs. This run was a pipeline
proof-of-concept, not a production campaign.

Compute: run on a cloud RTX 4090 (Ada, sm_89) via RunPod. The local RTX 5060 Ti
(Blackwell, sm_120) is not supported by RFantibody's CUDA 11.8 container, so
cloud compute on an Ampere/Ada GPU was used.

## Results

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
high pLDDT does not imply productive engagement, and that longer loops explore
more space at the cost of more frequent misses.

Visual inspection of the two epitope-reaching designs (psma_nb_1):
- psma_nb_1_dldesign_0: dock is non-productive. The framework (FR2/FR3), not the
  CDRs, forms most of the target contact.
- psma_nb_1_dldesign_1: better. The CDRs approach the epitope, though the overall
  angle of approach is suboptimal and remains too far from the hotspot residues. Very weak but potentially detectable binding at best.

Since dldesign_0 and dldesign_1 share the same backbone, the difference between a
framework-mediated dock and a CDR-mediated approach comes from the ProteinMPNN
sequence alone, highlighting that the sequence design step affects productive
engagement, not just the backbone.

## QC finding: target renumbering across outputs

During analysis, the RF2 output structures were found to renumber the target
chain inconsistently between designs (e.g. chain T started at residue 120 in one
output and 129 in another, versus 117 in the input). The target sequence and
fold were identical across all outputs (0 sequence differences; target Ca-Ca
RMSD ~0.05 A between designs), but the numbering offsets meant the original
hotspot residue numbers pointed to the wrong physical residues in the outputs.

This was caught by noticing that a single labelled residue (341) had a different
amino acid identity across two outputs, which is impossible for an identical
target unless the numbering had shifted. Hotspots were then re-mapped onto each
output by structural alignment to the correctly-numbered input and proximity
selection, and the design assessment was redone on the correct residues.

Lesson: residue numbering is not guaranteed to be preserved through this
pipeline. Any position-based analysis (hotspot contacts, interface residues)
must be mapped through structural alignment rather than trusting output residue
numbers.

## Honest assessment and limitations

- This was a 4-design proof-of-concept. Productive de novo nanobody designs
  against a defined epitope require generating hundreds to thousands of designs
  and filtering hard on interface quality. A handful rarely yields a good binder,
  and these results are consistent with the method's known low per-design hit
  rate.
- No formal interface pAE was available from this RF2 build; assessment used
  pLDDT plus hotspot-to-loop distance plus visual inspection of the interface.
- None of the four designs is a convincing binder. The best (psma_nb_1_dldesign_1)
  reaches the epitope with its CDRs but at a poor angle.
- No wet-lab validation. No developability assessment (aggregation, exposed
  hydrophobics, unpaired cysteines, glycosylation motifs). No screening against
  the GCPIII homolog for cross-reactivity.

## What a real campaign would add

- Generate hundreds+ of backbones, filter on interface pAE and CDR-mediated
  contact, and inspect survivors.
- Orthogonal structure prediction (e.g. ColabFold) of top candidates for a
  second confidence signal.
- Rosetta interface scoring (dG_separated, shape complementarity).
- Developability triage and cross-reactivity screening before any synthesis.

## Files

- targets/: input structure (4NGM crop), hotspot session
- designs/rf2/: the four RF2-predicted design structures
- designs/mpnn/: ProteinMPNN sequence-design intermediates
- designs/trb/: RFdiffusion metadata (pLDDT, hotspot distances, loop lengths)
- designs/*.png: interface and summary figures
