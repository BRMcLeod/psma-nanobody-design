# psma-nanobody-design
De novo nanobody design targeting PSMA using RFantibody
Remote dev environment confirmed working
## Target and epitope selection

### Target

Prostate-specific membrane antigen (PSMA, also GCPII / FOLH1) is a validated
prostate cancer target. It is overexpressed on prostate cancer cells and tumour
vasculature, and is the target of approved radioligand therapies
(177Lu-PSMA-617). PSMA is a homodimeric type II membrane glycoprotein with a
large extracellular domain, which is the region accessible to an antibody on a
live cell.

Structure: PDB 4NGM, the extracellular domain of human GCPII in complex with a
urea-based inhibitor, ~1.8 Å resolution. Chosen for high resolution, a complete
apical domain, and a clean single-ligand active site.

### Format

Single-domain antibody (VHH / nanobody) rather than a full IgG or Fab.
Nanobodies are ~110 residues, express well, are straightforward to engineer, and
are the format with the strongest published results in the RFantibody pipeline.
The smaller target footprint also suits a single 16 GB GPU compute budget.

### Epitope region

Apical domain, residues 153 to 347. This is the mapped epitope of the clinical
anti-PSMA antibody J591, giving direct precedent that a therapeutically useful
antibody binds here. It projects away from the membrane, so it is sterically
accessible on a live cell, and it is distinct from the active-site pocket.

Deliberately avoided:
- The active-site funnel (marked by two catalytic zincs): a deep pocket suited
  to small molecules, not antibodies, and nearly identical to the pocket in the
  GCPIII homolog, so targeting it risks cross-reactivity.
- Membrane-proximal and dimer-interface surfaces, which are occluded in the
  cellular context.

### Hotspot residues

Final hotspot set: Lys223, Glu171, Lys199, Asp173, Lys341 (chain A).

Reasoning:
- Clustered within ~13 Å (Cα-Cα across the widest pair), defining a single
  contiguous landing site sized for a nanobody paratope.
- Four charged residues plus Lys199, all solvent-exposed, giving specific
  surface contacts.
- Anchored on ordered secondary structure where possible. Lys223 sits in a
  beta turn between two alpha helices, presenting a defined, rigid shape rather
  than a mobile loop.
- Located on an open, outward-facing patch, confirmed by surface inspection in
  PyMOL.

Session: targets/psma_hotspots.pse (loadable in PyMOL to view target, epitope,
and hotspots).