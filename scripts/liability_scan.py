#!/usr/bin/env python3
"""
liability_scan.py - sequence liability triage for VHH / nanobody designs.

Splits each sequence into IMGT regions (FR1/CDR1/FR2/CDR2/FR3/CDR3/FR4) using
ANARCI, scans for developability liability motifs, and reports both a flagged
inventory (what is present, where) and a weighted score (per region and total).

Two output modes by design:
  - FLAGS: every motif occurrence with its position and region. Nothing is
    hidden, including occurrences that carry zero score.
  - SCORE: region-weighted sum. Only motifs whose region weight is non-zero
    contribute.

The distinction matters. Some liabilities are worth knowing about but not
worth penalising in a given context (a buried framework Met is a reminder,
not a risk). Those appear in the flag table with score 0. If flag counts and
the score appear to disagree, that is intended, not a bug.

Usage:
    python liability_scan.py designs.fasta
    python liability_scan.py designs/*.pdb --chain H
    python liability_scan.py designs.fasta --parent parent.fasta
    python liability_scan.py designs.fasta -o results.csv --flags flags.csv

Requires: ANARCI (conda install -c bioconda anarci "numpy<2.0").
Falls back to motif-based region assignment if ANARCI is unavailable; the
fallback is approximate and is reported as such in the output.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

# --------------------------------------------------------------------------
# Region definitions (IMGT)
# --------------------------------------------------------------------------
# IMGT boundaries. CDR3 is 105-117 inclusive; insertion codes (111A, 112A...)
# fall inside that range and are counted as occupied positions, so a long
# camelid H3 is measured by occupancy, not by 117-105+1.
IMGT_REGIONS = [
    ("FR1",  1,   26),
    ("CDR1", 27,  38),
    ("FR2",  39,  55),
    ("CDR2", 56,  65),
    ("FR3",  66,  104),
    ("CDR3", 105, 117),
    ("FR4",  118, 128),
]

REGION_ORDER = [r[0] for r in IMGT_REGIONS]
CDR_REGIONS = {"CDR1", "CDR2", "CDR3"}

# Canonical VHH disulfide, IMGT numbering. Cys23 and Cys104 are structural
# and expected in every domain; they are not liabilities.
CANONICAL_CYS_IMGT = {23, 104}

# Camelid VHH hallmark residues (IMGT FR2 positions 42, 49, 50, 52), per
# Vincke et al. (2009) JBC 284:3273-3284 and Fernandez-Quintero et al. (2024)
# Protein Sci 33:e5176. These positions substitute the conventional VH VGLW
# motif and both raise solubility (they replace the hydrophobic VL interface)
# and, at 42 and 52 especially, contribute directly to antigen binding and
# pre-organisation of the CDR3 loop. They are FLAGGED, never scored: on a
# humanized scaffold like h-NbBCII10 some are already human, so a change here
# is not inherently a defect - it is a position that couples to binding and
# stability and therefore warrants a look whenever it varies. 42 and 52 are
# the higher-risk pair.
HALLMARK_IMGT = {42: "high", 49: "moderate", 50: "moderate", 52: "high"}


# --------------------------------------------------------------------------
# Region weights
# --------------------------------------------------------------------------
# CDRs weighted 5x framework. Rationale: liabilities in CDRs sit at or near
# the paratope and are solvent-exposed, so a chemical modification there is
# far more likely to alter binding, and exposure raises the modification rate
# in the first place. Framework liabilities are more often buried and more
# often tolerated. The 5x figure is a triage convention, not a measured
# ratio - it is chosen to make any CDR liability outrank any single framework
# one, which is the ranking behaviour wanted from a funnel.
REGION_WEIGHTS = {
    "FR1": 1, "CDR1": 5,
    "FR2": 1, "CDR2": 5,
    "FR3": 1, "CDR3": 5,
    "FR4": 1,
}


# --------------------------------------------------------------------------
# Liability definitions
# --------------------------------------------------------------------------
# Each entry:
#   patterns      regex list, applied to the region sequence
#   severity      base score per occurrence, before region weighting
#   regions       which regions this liability is *scored* in.
#                 "all" = scored everywhere. A set = scored only there.
#                 Occurrences outside the scored set are still FLAGGED at
#                 score 0.
#   note          rationale, carried into the output so the scheme is
#                 self-documenting
#
# Severity is ordinal, not calibrated to any experimental rate. It encodes
# relative concern: 5 = would likely block a candidate, 3 = needs attention,
# 1-2 = note it.
#
# Primary literature anchors:
#   Lu et al. (2018) "Deamidation and isomerization liability analysis of 131
#     clinical-stage antibodies", mAbs 11(1):45-57. Empirical basis for
#     treating NG/DG as dominant and for the CDR-H2/H3 focus: across 131
#     clinical mAbs most deamidation/isomerization localised to three
#     positions (H54, H98, L30), with H54/H98 predominantly at NG and DG.
#   Jain et al. (2017) PNAS 114(5):944-949, the clinical-stage developability
#     dataset these conventions draw from.
#   Motif sets cross-checked against the Adimab/Geneious and LAP (Liability
#     Antibody Profiler, 2024) IMGT-region tables and against
#     developability-filtered VHH library patents (US 12,442,107 / 12,344,959).

LIABILITIES = {
    "N-glycosylation": {
        "patterns": [r"N[^P][ST]"],
        "severity": 5,
        "regions": "all",
        "note": "NxS/NxT sequon (x != P). Introduces heterogeneity and can "
                "block binding if in a CDR. High confidence motif.",
    },
    "Deamidation_high": {
        "patterns": [r"NG"],
        "severity": 5,
        "regions": "all",
        "note": "NG is the fastest-deamidating motif by a wide margin; "
                "the following Gly gives no steric hindrance to succinimide "
                "formation.",
    },
    "Deamidation_moderate": {
        "patterns": [r"N[STNH]"],
        "severity": 2,
        "regions": "all",
        "note": "NS/NT/NN/NH deamidate measurably but far slower than NG. "
                "Reverse motifs (SN, TN, KN) excluded: the effect is much "
                "weaker and including them floods the output.",
    },
    "Isomerization": {
        "patterns": [r"D[GSDNTH]"],
        "severity": 3,
        "regions": "all",
        "note": "Asp isomerisation via succinimide. DG highest risk, DS/DD/DN "
                "lower. Scored once here rather than double-counted under a "
                "separate deamidation heading.",
    },
    "Fragmentation": {
        "patterns": [r"DP"],
        "severity": 2,
        "regions": "all",
        "note": "Asp-Pro is the classic acid-labile peptide bond; a low-pH "
                "hold or formulation step can clip it.",
    },
    "Met_oxidation": {
        "patterns": [r"M"],
        "severity": 3,
        "regions": CDR_REGIONS,
        "note": "Met oxidation. Scored in CDRs only, where Met is solvent "
                "exposed and oxidation can abolish binding. Framework Met is "
                "flagged at score 0 as a reminder: it is usually buried, and "
                "without solvent-accessibility data region is the best "
                "available proxy for exposure.",
    },
    "Trp_oxidation": {
        "patterns": [r"W"],
        "severity": 2,
        "regions": CDR_REGIONS,
        "note": "Trp oxidation. Same logic as Met. Framework Trp is flagged "
                "at 0 - it includes the structurally essential conserved Trp "
                "positions, which must not be penalised.",
    },
    "Pyroglutamate": {
        "patterns": [r"^[QE]"],
        "severity": 1,
        "regions": {"FR1"},
        "note": "N-terminal Gln/Glu cyclisation. Anchored to the sequence "
                "start, so only an actual N-terminal residue matches. Very "
                "common and usually benign, hence severity 1.",
    },
    "Hydrophobic_patch": {
        "patterns": [r"[FILVWY]{3,}"],
        "severity": 3,
        "regions": CDR_REGIONS,
        "note": "Run of 3+ consecutive hydrophobics in a CDR. Aggregation and "
                "non-specific binding risk. Framework runs are normal "
                "secondary structure and are not flagged.",
    },
    "Aromatic_polyreactivity": {
        "patterns": [r"[WY]{2,}", r"W.W"],
        "severity": 2,
        "regions": CDR_REGIONS,
        "note": "Adjacent or near-adjacent aromatics in CDRs associate with "
                "polyreactivity. Note: the broad dipeptide lists sometimes "
                "used here (GG, VG, RR) are omitted - they occur too often in "
                "ordinary sequence to carry signal.",
    },
}


# --------------------------------------------------------------------------
# ANARCI numbering
# --------------------------------------------------------------------------

def run_anarci(sequences):
    """
    Number sequences with ANARCI (IMGT). Returns
    {name: {"numbered": [(num, ins, aa), ...], "status": "ok"|"failed"}}.
    Returns None only if ANARCI is not installed at all.

    Submits one sequence per call. Batch submission makes output-to-input
    matching fragile (a single mangled ID silently drops a design), and the
    cost of per-call invocation is trivial at panel scale. A design ANARCI
    cannot number is recorded as status "failed", never as an empty success.
    """
    if not _have_anarci():
        return None

    results = {}
    with tempfile.TemporaryDirectory() as td:
        for name, seq in sequences.items():
            fa = os.path.join(td, "in.fasta")
            with open(fa, "w") as fh:
                fh.write(">%s\n%s\n" % (name, seq))
            out = os.path.join(td, "out")
            for f in (out, out + "_H.csv"):
                if os.path.exists(f):
                    os.remove(f)
            try:
                subprocess.run(
                    ["ANARCI", "-i", fa, "-o", out, "--scheme", "imgt",
                     "--restrict", "H"],
                    check=True, capture_output=True, timeout=300,
                )
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired,
                    FileNotFoundError):
                results[name] = {"numbered": [], "status": "failed"}
                continue

            numbered = _parse_anarci_output(out)
            results[name] = {
                "numbered": numbered,
                "status": "ok" if numbered else "failed",
            }
    return results


def _have_anarci():
    try:
        subprocess.run(["ANARCI", "--help"], capture_output=True, timeout=30)
        return True
    except Exception:
        return False


def _parse_anarci_output(out_base):
    """
    Parse ANARCI output for a single sequence. This build writes the aligned
    text format (one residue per line: '<chain> <num> [ins] <aa>'), not a
    CSV, and does not append a '_H.csv' suffix. Handles both: prefers the
    plain file, falls back to a '_H.csv' if a CSV-writing build is used.

    Text format lines look like:
        H 42      F
        H 112   A N     (position 112, insertion code A, residue N)
        H 60      -      (gap; skipped)
    Comment lines start with '#', the record ends with '//'.
    """
    path = None
    if os.path.exists(out_base):
        path = out_base
    elif os.path.exists(out_base + "_H.csv"):
        return _parse_anarci_csv(out_base + "_H.csv")
    if path is None:
        return []

    numbered = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("//"):
                continue
            parts = line.split()
            # Expected: [chain, num, aa] or [chain, num, ins, aa]
            if len(parts) < 3 or parts[0] != "H":
                continue
            try:
                num = int(parts[1])
            except ValueError:
                continue
            if len(parts) == 3:
                ins, aa = "", parts[2]
            else:
                ins, aa = parts[2], parts[3]
            if aa == "-" or aa == ".":
                continue
            numbered.append((num, ins, aa))
    numbered.sort(key=lambda p: (p[0], p[1]))
    return numbered


def _parse_anarci_csv(path):
    """Fallback for ANARCI builds that write the CSV format instead."""
    numbered = []
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            for col, val in row.items():
                if not col or not val or val == "-":
                    continue
                m = re.match(r"^(\d+)([A-Z]?)$", col.strip())
                if not m:
                    continue
                numbered.append((int(m.group(1)), m.group(2), val))
            break  # single sequence per file
    numbered.sort(key=lambda p: (p[0], p[1]))
    return numbered


def split_regions_anarci(numbered_seq):
    """Split an ANARCI-numbered sequence into IMGT regions."""
    regions = {r: [] for r in REGION_ORDER}
    for num, ins, aa in numbered_seq:
        for rname, lo, hi in IMGT_REGIONS:
            if lo <= num <= hi:
                regions[rname].append((num, ins, aa))
                break
    return regions


# --------------------------------------------------------------------------
# Fallback region splitting
# --------------------------------------------------------------------------

def split_regions_fallback(seq):
    """
    Approximate region assignment using conserved VHH framework motifs.
    Used only when ANARCI is unavailable. Anchors:
      Cys22-ish (end FR1), W..Q (start FR2), second conserved Cys (end FR3),
      WGxG (start FR4).
    Returns regions with synthetic numbering; positions are sequence indices,
    not true IMGT numbers, and the output records that.
    """
    n = len(seq)

    c1 = seq.find("C", 15)
    if c1 == -1 or c1 > 30:
        c1 = 22
    w1 = seq.find("W", c1 + 5)
    if w1 == -1:
        w1 = c1 + 14

    c2 = seq.rfind("C")
    if c2 <= w1:
        c2 = int(n * 0.78)
    m = re.search(r"WG.G", seq[c2:])
    w2 = c2 + m.start() if m else min(c2 + 20, n - 10)

    cdr1_end = w1
    cdr2_start = w1 + 17
    cdr2_end = cdr2_start + 10

    bounds = {
        "FR1":  (0, c1 + 3),
        "CDR1": (c1 + 3, cdr1_end),
        "FR2":  (cdr1_end, cdr2_start),
        "CDR2": (cdr2_start, cdr2_end),
        "FR3":  (cdr2_end, c2 + 3),
        "CDR3": (c2 + 3, w2),
        "FR4":  (w2, n),
    }

    regions = {}
    for rname in REGION_ORDER:
        lo, hi = bounds[rname]
        lo, hi = max(0, lo), min(n, hi)
        regions[rname] = [(i + 1, "", seq[i]) for i in range(lo, hi)] if hi > lo else []
    return regions


# --------------------------------------------------------------------------
# Scanning
# --------------------------------------------------------------------------

def scan_regions(regions):
    """
    Scan each region for liability motifs.
    Returns (flags, scores_by_region) where flags is a list of dicts.
    """
    flags = []
    scores = {r: 0 for r in REGION_ORDER}

    for rname in REGION_ORDER:
        entries = regions.get(rname, [])
        if not entries:
            continue
        seq = "".join(e[2] for e in entries)
        weight = REGION_WEIGHTS[rname]

        for lname, ldef in LIABILITIES.items():
            scored_here = (ldef["regions"] == "all" or rname in ldef["regions"])

            for pattern in ldef["patterns"]:
                # lookahead keeps overlapping matches (YYY yields two YY)
                for m in re.finditer("(?=(%s))" % pattern, seq):
                    text = m.group(1)
                    start = m.start()
                    num, ins, _ = entries[start]
                    score = ldef["severity"] * weight if scored_here else 0
                    flags.append({
                        "region": rname,
                        "liability": lname,
                        "motif": text,
                        "imgt_position": "%d%s" % (num, ins),
                        "severity": ldef["severity"],
                        "region_weight": weight,
                        "score": score,
                        "scored": "yes" if scored_here else "flag_only",
                    })
                    scores[rname] += score

    return flags, scores


def scan_cysteines(regions, imgt=True):
    """
    Cysteine handling, separate from motif scanning.
    The canonical VHH disulfide (IMGT 23 and 104) is expected and is not a
    liability. Any additional Cys is flagged; an odd total guarantees at
    least one unpaired thiol, which is scored higher.

    Under real IMGT numbering the canonical pair is identified by position.
    Under the approximate fallback the numbers are sequence indices, so the
    pair is identified structurally instead: the first Cys in FR1 and the
    last Cys in FR3, which is where the canonical disulfide sits in every
    VHH. This is the same claim, made with the information available.
    """
    flags = []
    score = 0
    all_cys = []
    for rname in REGION_ORDER:
        for num, ins, aa in regions.get(rname, []):
            if aa == "C":
                all_cys.append((rname, num, ins))

    if imgt:
        canonical = [c for c in all_cys if c[1] in CANONICAL_CYS_IMGT]
    else:
        canonical = []
        fr1 = [c for c in all_cys if c[0] == "FR1"]
        fr3 = [c for c in all_cys if c[0] == "FR3"]
        if fr1:
            canonical.append(fr1[0])
        if fr3:
            canonical.append(fr3[-1])
    extra = [c for c in all_cys if c not in canonical]

    for rname, num, ins in canonical:
        flags.append({
            "region": rname, "liability": "Cys_canonical", "motif": "C",
            "imgt_position": "%d%s" % (num, ins), "severity": 0,
            "region_weight": REGION_WEIGHTS[rname], "score": 0,
            "scored": "flag_only",
        })

    odd = len(all_cys) % 2 == 1
    for rname, num, ins in extra:
        sev = 5 if odd else 3
        s = sev * REGION_WEIGHTS[rname]
        flags.append({
            "region": rname, "liability": "Cys_unpaired" if odd else "Cys_extra",
            "motif": "C", "imgt_position": "%d%s" % (num, ins),
            "severity": sev, "region_weight": REGION_WEIGHTS[rname],
            "score": s, "scored": "yes",
        })
        score += s

    return flags, score, len(all_cys), len(extra)


# --------------------------------------------------------------------------
# Physicochemical properties
# --------------------------------------------------------------------------

# Kyte-Doolittle hydropathy
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
    "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
    "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
    "Y": -1.3, "V": 4.2,
}


def net_charge(seq):
    """Approximate net charge at pH 7.4 from residue counts."""
    pos = seq.count("K") + seq.count("R")
    neg = seq.count("D") + seq.count("E")
    his = seq.count("H") * 0.1   # His mostly neutral at 7.4
    return pos + his - neg


def mean_hydropathy(seq):
    vals = [KD[a] for a in seq if a in KD]
    return sum(vals) / len(vals) if vals else 0.0


# --------------------------------------------------------------------------
# Parent comparison (optional, flag only)
# --------------------------------------------------------------------------

def compare_framework(regions, parent_regions):
    """
    Compare framework regions against a parent scaffold. Flag only, never
    scored: the consequence of a framework mutation depends entirely on
    whether the position is buried, part of the hydrophobic core, or a
    conserved functional residue, and none of that is knowable from sequence
    alone. A count is honest; a score would be false precision.
    """
    flags = []
    for rname in REGION_ORDER:
        if rname in CDR_REGIONS:
            continue
        design = "".join(e[2] for e in regions.get(rname, []))
        parent = "".join(e[2] for e in parent_regions.get(rname, []))
        if not design or not parent:
            continue

        n = min(len(design), len(parent))
        for i in range(n):
            if design[i] != parent[i]:
                num, ins, _ = regions[rname][i]
                flags.append({
                    "region": rname, "liability": "FW_deviation",
                    "motif": "%s>%s" % (parent[i], design[i]),
                    "imgt_position": "%d%s" % (num, ins),
                    "severity": 0, "region_weight": REGION_WEIGHTS[rname],
                    "score": 0, "scored": "flag_only",
                })
    return flags


def scan_hallmarks(regions):
    """
    Report the residue occupying each camelid FR2 hallmark position (IMGT 42,
    49, 50, 52). Flag only, never scored. Requires true IMGT numbering; under
    the fallback these positions can't be located reliably, so the caller
    should pass imgt=True only when ANARCI ran.

    Returns (flags, occupancy) where occupancy is {position: residue}.
    """
    flags = []
    occ = {}
    for rname in REGION_ORDER:
        for num, ins, aa in regions.get(rname, []):
            if num in HALLMARK_IMGT and ins == "":
                occ[num] = aa
                flags.append({
                    "region": rname, "liability": "Hallmark_FR2",
                    "motif": "%d:%s" % (num, aa),
                    "imgt_position": "%d" % num,
                    "severity": 0, "region_weight": REGION_WEIGHTS.get(rname, 1),
                    "score": 0, "scored": "flag_only",
                })
    return flags, occ


# --------------------------------------------------------------------------
# Input parsing
# --------------------------------------------------------------------------

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0].split("|")[0]
                buf = []
            elif line:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def read_pdb(path, chain="H"):
    """Extract a chain's sequence from CA records."""
    seq, seen = [], set()
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            if line[21] != chain:
                continue
            key = (line[22:27])
            if key in seen:
                continue
            seen.add(key)
            seq.append(AA3.get(line[17:20].strip(), "X"))
    return "".join(seq)


def strip_tag(seq):
    """
    Remove a trailing purification tag / linker from a construct sequence.
    Crystallised parent sequences (e.g. from the PDB) often carry these, and
    an uncleaved tag registers as a large FR4 deviation.
    """
    original = seq
    seq = re.sub(r"H{6,}$", "", seq)
    seq = re.sub(r"(GS|SG|RGR|GGGGS|AAA|EPEA|LE)+$", "", seq)
    return seq, (seq != original)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Sequence liability triage for VHH / nanobody designs.")
    ap.add_argument("inputs", nargs="+",
                    help="FASTA file(s) and/or PDB file(s)")
    ap.add_argument("--chain", default="H",
                    help="Chain to read from PDB inputs (default: H)")
    ap.add_argument("--parent", default=None,
                    help="Optional parent scaffold FASTA. Enables framework "
                         "deviation flagging (flag only, never scored).")
    ap.add_argument("-o", "--output", default="liability_summary.csv",
                    help="Per-design summary CSV")
    ap.add_argument("--flags", default="liability_flags.csv",
                    help="Per-occurrence flag CSV")
    ap.add_argument("--no-strip-tag", action="store_true",
                    help="Do not strip trailing His-tags / linkers")
    args = ap.parse_args()

    # ---- load sequences
    sequences = {}
    for path in args.inputs:
        if not os.path.exists(path):
            sys.stderr.write("skipping missing file: %s\n" % path)
            continue
        if path.lower().endswith((".fasta", ".fa", ".faa", ".fas")):
            sequences.update(read_fasta(path))
        elif path.lower().endswith(".pdb"):
            s = read_pdb(path, args.chain)
            if s:
                sequences[os.path.splitext(os.path.basename(path))[0]] = s
            else:
                sys.stderr.write("no chain %s in %s\n" % (args.chain, path))
        else:
            sys.stderr.write("unrecognised extension: %s\n" % path)

    if not sequences:
        sys.stderr.write("no sequences loaded\n")
        return 1

    # ---- parent
    parent_seq = None
    if args.parent:
        pseqs = read_fasta(args.parent)
        if not pseqs:
            sys.stderr.write("could not read parent: %s\n" % args.parent)
            return 1
        pname, parent_seq = next(iter(pseqs.items()))
        if not args.no_strip_tag:
            parent_seq, stripped = strip_tag(parent_seq)
            if stripped:
                sys.stderr.write(
                    "note: stripped trailing tag/linker from parent '%s'\n"
                    % pname)

    # ---- numbering
    to_number = dict(sequences)
    if parent_seq:
        to_number["__PARENT__"] = parent_seq

    numbered = run_anarci(to_number)
    anarci_available = numbered is not None
    if not anarci_available:
        sys.stderr.write(
            "WARNING: ANARCI not installed. Using approximate motif-based "
            "region assignment for all designs. Positions are sequence "
            "indices, not IMGT numbers, and hallmark / canonical-Cys checks "
            "are disabled.\n")

    parent_regions = None
    if parent_seq:
        pinfo = numbered.get("__PARENT__") if anarci_available else None
        if pinfo and pinfo["status"] == "ok":
            parent_regions = split_regions_anarci(pinfo["numbered"])
        else:
            parent_regions = split_regions_fallback(parent_seq)
            if anarci_available:
                sys.stderr.write(
                    "WARNING: ANARCI could not number the parent scaffold; "
                    "framework comparison uses approximate numbering.\n")

    # ---- scan
    summary_rows, flag_rows = [], []
    n_failed = 0

    for name, seq in sequences.items():
        info = numbered.get(name) if anarci_available else None
        if info and info["status"] == "ok":
            regions = split_regions_anarci(info["numbered"])
            method = "ANARCI/IMGT"
        elif anarci_available and info and info["status"] == "failed":
            # ANARCI is installed but could not number THIS design. Do not
            # silently fall back to approximate numbering and do not emit a
            # zero row that would rank as the safest design. Report it.
            method = "FAILED"
            n_failed += 1
        else:
            regions = split_regions_fallback(seq)
            method = "fallback/approximate"

        if method == "FAILED":
            row = {
                "design": name, "length": len(seq), "numbering": "FAILED",
                "total_score": "", "cdr_score": "", "framework_score": "",
                "cys_score": "", "n_flags_total": "", "n_flags_scored": "",
                "n_flags_only": "", "cdr3_length": "", "cdr3_long": "",
                "n_cys": "", "n_cys_noncanonical": "", "net_charge": "",
                "cdr_net_charge": "", "cdr_mean_hydropathy": "",
                "fw_deviations": "", "fr2_hallmark_42_49_50_52": "FAILED",
            }
            for r in REGION_ORDER:
                row["score_" + r] = ""
            for lname in list(LIABILITIES) + ["Cys_unpaired", "Cys_extra"]:
                row["n_" + lname] = ""
            summary_rows.append(row)
            continue

        flags, region_scores = scan_regions(regions)
        cys_flags, cys_score, n_cys, n_extra_cys = scan_cysteines(
            regions, imgt=(method == "ANARCI/IMGT"))
        flags.extend(cys_flags)

        if parent_regions:
            flags.extend(compare_framework(regions, parent_regions))

        hallmark_ran = (method == "ANARCI/IMGT")
        if hallmark_ran:
            hm_flags, hm_occ = scan_hallmarks(regions)
            flags.extend(hm_flags)
            hallmark_sig = "".join(
                hm_occ.get(p, "-") for p in (42, 49, 50, 52))
        else:
            hallmark_sig = "needs_anarci"

        cdr_seq = "".join(
            e[2] for r in CDR_REGIONS for e in regions.get(r, []))
        h3 = regions.get("CDR3", [])

        total = sum(region_scores.values()) + cys_score
        cdr_total = sum(region_scores[r] for r in CDR_REGIONS)

        fw_dev = sum(1 for f in flags if f["liability"] == "FW_deviation")
        flag_only = sum(1 for f in flags if f["scored"] == "flag_only")

        by_type = defaultdict(int)
        for f in flags:
            by_type[f["liability"]] += 1

        row = {
            "design": name,
            "length": len(seq),
            "numbering": method,
            "total_score": total,
            "cdr_score": cdr_total,
            "framework_score": total - cdr_total - cys_score,
            "cys_score": cys_score,
            "n_flags_total": len(flags),
            "n_flags_scored": len(flags) - flag_only,
            "n_flags_only": flag_only,
            "cdr3_length": len(h3),
            "cdr3_long": "yes" if len(h3) > 22 else "no",
            "n_cys": n_cys,
            "n_cys_noncanonical": n_extra_cys,
            "net_charge": round(net_charge(seq), 1),
            "cdr_net_charge": round(net_charge(cdr_seq), 1),
            "cdr_mean_hydropathy": round(mean_hydropathy(cdr_seq), 2),
            "fw_deviations": fw_dev if parent_regions else "not_run",
            "fr2_hallmark_42_49_50_52": hallmark_sig,
        }
        for r in REGION_ORDER:
            row["score_" + r] = region_scores[r]
        for lname in list(LIABILITIES) + ["Cys_unpaired", "Cys_extra"]:
            row["n_" + lname] = by_type.get(lname, 0)

        summary_rows.append(row)
        for f in flags:
            f2 = dict(f)
            f2["design"] = name
            flag_rows.append(f2)

    # ---- write
    # Failed rows (empty total_score) sort to the bottom, explicitly, rather
    # than crashing the numeric sort or masquerading as low-risk.
    def sort_key(r):
        t = r["total_score"]
        return (0, 0) if t == "" else (1, -t)
    summary_rows.sort(key=sort_key)

    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        w.writerows(summary_rows)

    if flag_rows:
        cols = ["design", "region", "liability", "motif", "imgt_position",
                "severity", "region_weight", "score", "scored"]
        with open(args.flags, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in flag_rows:
                w.writerow({c: r.get(c, "") for c in cols})

    # ---- console
    methods_seen = set(r["numbering"] for r in summary_rows)
    print("\nLiability scan: %d designs" % len(summary_rows))
    print("Region assignment: %s" % ", ".join(sorted(methods_seen)))
    print("Framework comparison: %s"
          % ("on" if parent_regions else "off (no --parent)"))
    if n_failed:
        print("ANARCI could not number %d design(s); reported as FAILED, "
              "not scored." % n_failed)
    print("\n%-32s %6s %6s %6s %5s %6s %6s"
          % ("design", "total", "CDR", "FW", "H3", "flags", "FWdev"))
    print("-" * 78)
    for r in summary_rows:
        if r["total_score"] == "":
            print("%-32s %6s %6s %6s %5s %6s %6s"
                  % (r["design"][:32], "FAIL", "-", "-", "-", "-", "-"))
            continue
        print("%-32s %6d %6d %6d %5d %6d %6s"
              % (r["design"][:32], r["total_score"], r["cdr_score"],
                 r["framework_score"], r["cdr3_length"], r["n_flags_total"],
                 str(r["fw_deviations"])))
    print("\nScores are ordinal triage values, not calibrated risk estimates.")
    print("'flag_only' entries appear in the flag table at score 0 by design.")
    if "fallback/approximate" in methods_seen:
        print("FR2 hallmark residues (IMGT 42/49/50/52) require ANARCI and were "
              "not resolved for fallback-numbered designs.")
    print("\nWrote %s and %s" % (args.output, args.flags))
    return 0


if __name__ == "__main__":
    sys.exit(main())
