# PyMOL Cheat Sheet

A practical reference for structural inspection, oriented toward antibody/epitope work. Commands are typed into the PyMOL command bar (the line at the top of the viewer, or the external GUI command line). Everything is case-sensitive.

---

## Loading and saving

| Command | What it does |
|---|---|
| `fetch 4ngm` | Download a structure from the PDB by ID and load it |
| `fetch 4ngm, type=pdb` | Force classic PDB format (default is mmCIF) |
| `load targets/4NGM.pdb` | Load a local file |
| `load targets/4NGM.pdb, psma` | Load and name the object `psma` |
| `save session.pse` | Save the full session (objects, view, colours) |
| `save out.pdb, chain A` | Export a selection's coordinates to a PDB file |
| `delete all` | Clear everything; `delete psma` removes one object |

PSE sessions are the right way to save work-in-progress; they preserve everything exactly as you left it.

---

## Mouse navigation (in the viewer)

| Action | Mouse |
|---|---|
| Rotate | Left-click drag |
| Pan / translate | Middle-click drag |
| Zoom | Right-click drag (vertical) or scroll wheel |
| Select an atom | Left-click on it |
| Clipping slab | Scroll wheel (adjusts near/far clip) |

If you ever lose the molecule off-screen, `orient` or `reset` brings it back.

---

## View commands

| Command | What it does |
|---|---|
| `orient` | Rotate/zoom to best-fit the whole structure (or a selection) |
| `orient chain A` | Best-fit just that selection |
| `zoom` | Fit everything in view |
| `zoom resi 153-347` | Zoom to a selection |
| `center selection` | Recenter rotation point on a selection |
| `reset` | Reset to default view |
| `turn y, 90` | Rotate the view 90° about the y-axis |
| `set_view (...)` | Restore an exact saved camera (copy from `get_view`) |

---

## Selections (PyMOL's real power)

Make a named selection: `select NAME, EXPRESSION`. Example: `select epitope, chain A and resi 153-347`.

**Building blocks:**

| Selector | Meaning | Example |
|---|---|---|
| `resi` | Residue number(s) | `resi 153-347`, `resi 700` |
| `resn` | Residue name | `resn HIS`, `resn ZN` |
| `chain` | Chain ID | `chain A` |
| `name` | Atom name | `name CA`, `name CA+CB` |
| `elem` | Element | `elem Zn`, `elem C` |
| `polymer` | Protein/nucleic atoms | `polymer` |
| `organic` | Small-molecule ligands | `organic` |
| `solvent` | Waters | `solvent` |
| `metals` | Metal ions | `metals` |
| `hetatm` | Non-standard (ligands, ions, water) | `hetatm` |

**Logic and proximity:**

| Operator | Meaning | Example |
|---|---|---|
| `and` / `or` / `not` | Boolean logic | `chain A and resi 200-250` |
| `within X of S` | Atoms within X Å of selection S (includes S) | `polymer within 5 of organic` |
| `around X` | Atoms within X Å of S, excluding S | `organic around 5` |
| `byres (...)` | Expand to whole residues | `byres (polymer within 5 of organic)` |
| `name CA and ...` | Restrict to specific atoms | `byres (... ) and name CA` |

**Handy ready-made selections:**
```
select pocket, byres (polymer within 5 of organic)   # residues lining the ligand site
select epitope, chain A and resi 153-347              # the J591 apical-domain region
select zinc, elem Zn                                  # the catalytic zincs
```

---

## Showing and hiding

| Command | What it does |
|---|---|
| `show cartoon` | Show cartoon (ribbon) representation |
| `hide everything` | Hide all representations |
| `as cartoon` | Replace all reps with cartoon only |
| `show sticks, resn HIS` | Add sticks for a selection |
| `show surface, chain A` | Molecular surface |
| `show spheres, elem Zn` | Spheres (good for ions) |
| `hide surface` | Hide one representation type |

Representation types: `lines`, `sticks`, `cartoon`, `ribbon`, `surface`, `spheres`, `dots`, `mesh`, `nonbonded`.

A clean starting point for a protein:
```
hide everything
show cartoon
show sticks, organic
show spheres, metals
```

---

## Colour

(The command is `color`, American spelling.)

| Command | What it does |
|---|---|
| `color cyan, chain A` | Colour a selection |
| `color red, epitope` | Colour a named selection |
| `util.cbc` | Colour each chain a different colour |
| `util.cbag` | Colour by element, green carbons (a, by atom; g, green) |
| `spectrum b` | Colour by B-factor across the structure |
| `spectrum b, blue_white_red, psma` | B-factor with a chosen palette |
| `bg_color white` | Set background colour |

`spectrum b` becomes very useful later: AlphaFold/ColabFold and RF2 store their per-residue confidence (pLDDT/pAE) in the B-factor column, so `spectrum b` instantly colours a predicted model by confidence.

---

## Measuring

| Command | What it does |
|---|---|
| `distance d1, sel1, sel2` | Create a distance object between two selections |
| `get_distance (resi 200 and name CA), (resi 250 and name CA)` | Print one distance to the log |
| `angle a1, s1, s2, s3` | Measure an angle |
| `dist contacts, chain A, chain B, mode=2` | Polar contacts (H-bonds) between chains |

The Wizard menu (Wizard → Measurement) also lets you click atoms to measure interactively.

---

## Cleaning up a structure

| Command | What it does |
|---|---|
| `remove solvent` | Strip waters |
| `remove not polymer` | Keep only protein/nucleic (drops ligands, ions, water) |
| `remove not chain A` | Keep only chain A |
| `create monomer, chain A` | Copy a selection into a new standalone object |
| `remove hydrogens` | Drop hydrogens |

For design prep you'll often want a single clean chain: `create target, chain A and polymer` then work on `target`.

---

## Making figures

```
bg_color white
set ray_opaque_background, 0     # transparent background (for slides/figures)
orient                            # frame the shot
ray 1600, 1200                    # ray-trace at this pixel size (slow but crisp)
png ~/Desktop/figure.png, dpi=300
```

Tips: `set cartoon_transparency, 0.4` to see ligands through cartoon; `set ray_shadows, 0` for flatter publication-style lighting; `bg_color white` before raytracing.

---

## Recipes for the PSMA project

**Look at the whole target, oriented and cleaned:**
```
fetch 4ngm
remove solvent
as cartoon
util.cbc
orient
```

**Highlight and zoom the candidate epitope (apical domain):**
```
select epitope, chain A and resi 153-347
color orange, epitope
show sticks, epitope and not name C+N+O
zoom epitope
```

**See what lines the active-site pocket (where small-molecule ligands bind):**
```
select pocket, byres (polymer within 5 of organic)
show sticks, pocket
color yellow, pocket
zoom pocket
```

**Isolate one monomer for design work** (PSMA crystallizes as a dimer):
```
create target, chain A and polymer
disable 4ngm
zoom target
```

---

## Getting information

| Command | What it does |
|---|---|
| `get_chains` | List chain IDs |
| `count_atoms chain A and polymer` | Count atoms in a selection |
| `iterate epitope and name CA, print(resi, resn)` | Loop over residues and print info |
| Display → Sequence | Show the sequence viewer (click residues to select) |

---

## Quick tips

- Press **ESC** to toggle between the graphics view and the text log.
- Object/selection names appear in the right-hand panel; click the toggles there to show/hide/colour without typing.
- `deselect` clears the pink selection dots.
- Commands can be abbreviated when unambiguous (`hide everything` works, so does most tab-completion; press Tab to autocomplete).
- If the view gets messy, `hide everything` then `show cartoon` is the fastest reset to a clean slate.
