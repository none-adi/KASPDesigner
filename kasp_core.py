#!/usr/bin/env python3
"""
kasp_core.py – KASP Primer Design Engine (all logic, no UI)

This module contains every function needed to design allele‑specific and common
primers for a list of loci. It is designed to be imported by both a command‑line
script and a Streamlit web application.

Changes from the original script:
- Added `reverse_complement()` function (was missing – caused NameError).
- Removed the artificial 'Q' padding; instead we validate that both flanks
  are ≥40 bp and reject loci that are too short.
- design_one_locus() now returns the original input row index (as `_index`)
  so the caller can preserve the original order.
- process_rows() accepts a list of rows (with index) instead of a file path.
- Full product size check implemented (≤100 bp) using binding positions
  returned by design_common_primer().
"""

import sys
import csv
import io
import math
import concurrent.futures
from typing import Dict, Optional, Tuple, List

# =============================================================================
# 1. Import Biopython if available – graceful fallback to empirical Tm
# =============================================================================
try:
    from Bio.SeqUtils import MeltingTemp as mt
    USE_BIOPYTHON = True
except ImportError:
    USE_BIOPYTHON = False
    # We print to stderr; in a Streamlit app this goes to the server log.
    print("Warning: Biopython not installed. Using empirical Tm formula (less accurate).",
          file=sys.stderr)

# =============================================================================
# 2. Constants (PCR tails, salt concentration)
# =============================================================================
FAM_TAIL = "GAAGGTGACCAAGTTCATGCT"      # FAM universal tail (5'→3')
HEX_TAIL = "GAAGGTCGGAGTCAACGGATT"      # HEX universal tail (5'→3')
SALT_CONC = 50.0                         # mM Na+ (default for Tm_NN)

# =============================================================================
# 3. DNA helper: reverse complement
# =============================================================================
def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence (5'→3').
    Handles only A, T, G, C; other characters are kept unchanged.
    """
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(comp.get(base, base) for base in reversed(seq.upper()))

# =============================================================================
# 4. Melting temperature calculations (two methods)
# =============================================================================
def tm_empirical(seq: str) -> float:
    """
    Original EasyKASP empirical formula: Tm = 72 + (37.4 * GC% - 747)/length
    Quick but less accurate for short oligos.
    """
    if not seq:
        return 0.0
    gc = seq.count('G') + seq.count('C')
    length = len(seq)
    return 72.0 + (37.4 * gc - 747.0) / length

def tm_nearest_neighbor(seq: str, na_conc: float = SALT_CONC) -> float:
    """
    Nearest‑neighbour Tm using Biopython's implementation of SantaLucia 1998.
    Falls back to empirical formula if an error occurs.
    """
    if not seq or len(seq) < 6:
        return 0.0
    try:
        return mt.Tm_NN(seq, Na=na_conc, DNA=True)
    except Exception:
        return tm_empirical(seq)

def tm_calc(seq: str) -> float:
    """
    Configurable Tm calculator – uses nearest‑neighbour if Biopython is available,
    otherwise the empirical formula.
    """
    if USE_BIOPYTHON:
        return tm_nearest_neighbor(seq)
    else:
        return tm_empirical(seq)

# =============================================================================
# 5. Sequence pattern checks (quality control)
# =============================================================================
def has_mono_run(seq: str, run_len: int = 6) -> bool:
    """
    Return True if the sequence contains a run of `run_len` identical bases
    (e.g. 'AAAAAA').
    """
    for i in range(len(seq) - run_len + 1):
        if all(seq[i+j] == seq[i] for j in range(1, run_len)):
            return True
    return False

def has_complementary_run(seq: str, run_len: int = 4) -> bool:
    """
    Return True if any window of length `run_len` consists of alternating
    complementary bases (e.g. 'ATAT', 'GCGC').
    """
    if len(seq) < run_len:
        return False
    for i in range(len(seq) - run_len + 1):
        window = seq[i:i+run_len]
        ok = True
        for j in range(1, run_len):
            pair = (window[j-1], window[j])
            if pair not in [('A','T'), ('T','A'), ('G','C'), ('C','G')]:
                ok = False
                break
        if ok:
            return True
    return False

def has_ssr(seq: str) -> bool:
    """
    Detect simple sequence repeats: di‑ or tri‑nucleotide units repeated
    at least 3 times (e.g. 'ATATAT', 'CAGCAGCAG').
    """
    for length in [2, 3]:
        for i in range(len(seq) - length*3 + 1):
            unit = seq[i:i+length]
            if seq[i:i+length*3] == unit * 3:
                return True
    return False

def all_bases_present(seq: str) -> bool:
    """True if the sequence contains at least one A, T, G, and C."""
    return all(base in seq for base in "ATGC")

def valid_3prime_end(seq: str) -> bool:
    """
    False if the last 4 bases are all G/C, or the last 5 are all A/T.
    (Primer 3' end should not be too stable or too weak.)
    """
    if len(seq) >= 4 and all(b in 'GC' for b in seq[-4:]):
        return False
    if len(seq) >= 5 and all(b in 'AT' for b in seq[-5:]):
        return False
    return True

# =============================================================================
# 6. Dimer and hairpin prediction
# =============================================================================
def check_primer_dimer(primer1: str, primer2: str) -> Tuple[bool, str]:
    """
    Check for 3' end complementarity between two primers.
    Returns (is_acceptable, reason).
    A primer pair is rejected if there are ≥4 consecutive complementary bases.
    Also checks self‑dimer (palindrome) of primer1.
    """
    max_match = 0
    # Compare 3' of primer1 with 5' of primer2 (both orientations)
    for i in range(1, min(len(primer1), len(primer2)) + 1):
        if primer1[-i:].upper() == primer2[:i].upper():
            max_match = i
        elif primer1[-i:].upper() == reverse_complement(primer2[:i]).upper():
            max_match = i
    if max_match >= 4:
        return False, f"Primer‑dimer: 3' end complementarity of {max_match} bp"

    # Self‑dimer check (primer1 vs its own reverse complement)
    rev_comp = reverse_complement(primer1)
    for i in range(1, min(len(primer1), len(rev_comp)) + 1):
        if primer1[-i:] == rev_comp[:i]:
            return False, f"Self‑dimer: {i} bp complementarity"
    return True, ""

def check_hairpin(primer: str, min_stem: int = 4) -> Tuple[bool, str]:
    """
    Simple hairpin detection: look for any substring of length ≥ min_stem
    that appears in the reverse complement (indicating an inverted repeat).
    Returns (is_acceptable, reason).
    """
    rev_comp = reverse_complement(primer)
    for length in range(min_stem, len(primer)//2 + 1):
        for i in range(len(primer) - length + 1):
            sub = primer[i:i+length]
            if sub in rev_comp:
                return False, f"Hairpin: {length} bp stem"
    return True, ""

def check_primer(seq: str, is_specific: bool, partner: Optional[str] = None) -> Tuple[bool, str]:
    """
    Run all quality checks on a primer sequence:
    - mono‑nucleotide runs
    - complementary runs
    - simple sequence repeats
    - all four bases present
    - valid 3' end
    - hairpin (self‑dimer)
    - cross‑dimer with a partner (if given)
    Returns (is_acceptable, error_message).
    """
    if not seq:
        return False, "Empty sequence"

    if has_mono_run(seq, 6):
        return False, "Mono‑nucleotide run (≥6 identical bases)"
    if has_complementary_run(seq, 4):
        return False, "Alternating complementary run (≥4 bp)"
    if has_ssr(seq):
        return False, "Simple sequence repeat (di‑/tri‑nucleotide)"
    if not all_bases_present(seq):
        return False, "Missing one or more base types (A,T,G,C)"
    if not valid_3prime_end(seq):
        return False, "Invalid 3' end (last 4 all G/C or last 5 all A/T)"

    ok, reason = check_hairpin(seq)
    if not ok:
        return False, reason

    if partner:
        ok, reason = check_primer_dimer(seq, partner)
        if not ok:
            return False, reason

    return True, ""

# =============================================================================
# 7. Design of allele‑specific primers (FAM and HEX)
# =============================================================================
def design_specific_primer(forward_direction: bool,
                           allele: str,
                           left_flank: str,
                           right_flank: str,
                           partner_seq: Optional[str] = None,
                           fam_length: Optional[int] = None,
                           m2: int = 3) -> Tuple[Optional[str], float, str]:
    """
    Design one allele‑specific primer.
    Parameters:
        forward_direction : True = forward, False = reverse
        allele            : the variant allele (e.g. 'A' or 'ATCG')
        left_flank        : left flank sequence (≥40 bp guaranteed)
        right_flank       : right flank sequence (≥40 bp guaranteed)
        partner_seq       : the other allele‑specific primer for dimer check
        fam_length        : length of the FAM primer (used to constrain HEX length)
        m2                : 0 if alleles are complementary, else 3 (for length matching)
    Returns:
        (primer_sequence, Tm, error_message) – error_message is empty on success.
    """
    if forward_direction:
        upstream = left_flank[-40:]   # last 40 bp of left flank
        candidate = upstream + allele
    else:
        downstream = right_flank[:40] # first 40 bp of right flank
        candidate = allele + downstream

    start_len, end_len = 27, 21

    # Determine allowed primer lengths
    if fam_length is not None:
        if m2 == 0:
            allowed_lengths = [fam_length]
        else:
            allowed_lengths = list(range(fam_length, min(fam_length+4, 28)))
        allowed_lengths = [l for l in allowed_lengths if end_len <= l <= start_len]
    else:
        allowed_lengths = list(range(start_len, end_len-1, -1))

    for k in allowed_lengths:
        if forward_direction:
            primer = candidate[-k:] if len(candidate) >= k else candidate
        else:
            primer = candidate[:k] if len(candidate) >= k else candidate
        if len(primer) < 21:
            continue

        tm = tm_calc(primer)
        if tm < 56.0 or tm > 59.0:
            continue

        ok, err = check_primer(primer, is_specific=True, partner=partner_seq)
        if not ok:
            continue

        if not forward_direction:
            primer = reverse_complement(primer)   # convert to 5'→3' orientation
        return primer, tm, ""

    return None, 0.0, "No valid specific primer (length 21‑27, Tm 56‑59°C, no repeats/dimer/hairpin)"

# =============================================================================
# 8. Design of the common primer
# =============================================================================
def design_common_primer(forward_direction: bool,
                         left_flank: str,
                         right_flank: str,
                         var_start: int,
                         var_end: int,
                         full_seq: str,
                         specific_fam: str,
                         specific_hex: str) -> Tuple[Optional[str], str, int, int]:
    """
    Design the common primer (binds opposite strand to the allele‑specific primers).
    Parameters:
        forward_direction : orientation of the allele‑specific primers
        left_flank, right_flank : unpadded flanks (≥40 bp)
        var_start, var_end : indices of '[' and ']' in full_seq
        full_seq          : complete padded? no – raw assembled sequence (no padding)
        specific_fam, specific_hex : the two allele‑specific primers (for dimer check)
    Returns:
        (primer_sequence (5'→3'), error_message, common_start, common_len)
        common_start and common_len refer to the extracted region in full_seq
        before any reverse complement (used later for product size calculation).
    """
    best_primer = None
    best_tm = 0.0
    best_start = 0
    best_len = 0

    # Paper: distance from variant = 4-30, primer length = 22-32
    for jj in range(4, 31):          # distance downstream (or upstream) of the variant
        for k in range(32, 21, -1):  # length (try longer first)
            if forward_direction:
                start_pos = var_end + jj
                if start_pos + k > len(full_seq):
                    continue
                candidate = full_seq[start_pos:start_pos + k]
            else:
                start_pos = var_start - k - jj
                if start_pos < 0:
                    continue
                candidate = full_seq[start_pos:start_pos + k]
            if len(candidate) < 22:
                continue

            tm = tm_calc(candidate)
            if tm < 59.0 or tm > 62.0:
                continue

            # Basic pattern checks on the candidate (before orientation)
            if has_mono_run(candidate, 6):
                continue
            if len(candidate) >= 4:
                if all(b in 'AT' for b in candidate[:4]) or all(b in 'GC' for b in candidate[:4]):
                    continue
                if all(b in 'AT' for b in candidate[-4:]) or all(b in 'GC' for b in candidate[-4:]):
                    continue

            ok, _ = check_hairpin(candidate)
            if not ok:
                continue

            # Dimer check against both allele‑specific primers.
            # The actual common primer in PCR is the reverse complement if forward_direction.
            common_for_dimer = reverse_complement(candidate) if forward_direction else candidate
            ok_fam, _ = check_primer_dimer(specific_fam, common_for_dimer)
            ok_hex, _ = check_primer_dimer(specific_hex, common_for_dimer)
            if not (ok_fam and ok_hex):
                continue

            # Prefer a Tm close to 60.4°C
            if abs(tm - 60.4) <= 0.57:
                best_primer = candidate
                best_tm = tm
                best_start = start_pos
                best_len = k
                break
            elif best_primer is None:
                best_primer = candidate
                best_tm = tm
                best_start = start_pos
                best_len = k
        if best_primer is not None:
            break

    if best_primer is None:
        return None, "No valid common primer (length 22‑32, Tm 59‑62°C, no repeats/dimer/hairpin)", 0, 0

    # Return primer in 5'→3' orientation as used in PCR
    if forward_direction:
        best_primer = reverse_complement(best_primer)
    return best_primer, "", best_start, best_len

# =============================================================================
# 9. Direction decision (forward or reverse design)
# =============================================================================
def decide_direction(allele1: str, allele2: str,
                     left_flank: str, right_flank: str) -> int:
    """
    Decide whether to design primers in forward (1) or reverse (-1) orientation.
    Based on GC content of the flanks and A/T runs, as in the original VBA.
    """
    left_win = left_flank[-26:] if len(left_flank) >= 26 else left_flank
    right_win = right_flank[:26] if len(right_flank) >= 26 else right_flank

    def gc_percent(seq):
        if not seq:
            return 50.0
        gc = seq.count('G') + seq.count('C')
        return 100.0 * gc / len(seq)

    gc_left = gc_percent(left_win)
    gc_right = gc_percent(right_win)

    if abs(gc_left - 45) <= abs(gc_right - 45):
        FR = 1
    else:
        FR = -1

    # Tie‑breaker: avoid A/T‑rich ends
    left_tail = left_flank[-4:] if len(left_flank) >= 4 else left_flank
    right_head = right_flank[:4] if len(right_flank) >= 4 else right_flank
    if left_tail.count('A') + left_tail.count('T') >= 3:
        FR = -1
    if right_head.count('A') + right_head.count('T') >= 3:
        FR = 1

    # For indels, prefer direction where the variable part is at the 3' end
    if len(allele1) != 1 or len(allele2) != 1:
        common_prefix = 0
        for i in range(min(len(allele1), len(allele2))):
            if allele1[i] == allele2[i]:
                common_prefix += 1
            else:
                break
        common_suffix = 0
        for i in range(1, min(len(allele1), len(allele2))+1):
            if allele1[-i] == allele2[-i]:
                common_suffix += 1
            else:
                break
        if common_prefix >= 3:
            FR = -1
        if common_suffix >= 3:
            FR = 1
    return FR

# =============================================================================
# 10. Main design function for a single locus
# =============================================================================
def design_one_locus(row_data: Tuple[int, List[str]]) -> Dict:
    """
    Process a single locus.
    Input:
        row_data = (original_index, [locus_name, allele1, allele2, sequence])
    Output:
        Dictionary with keys:
            locus      : locus name
            fam        : FAM primer (with universal tail)
            hex        : HEX primer (with universal tail)
            common     : common primer sequence
            status     : "Success" or "Failed"
            reason     : error message if failed
            _index     : the original row index (for sorting)
    """
    idx, row = row_data
    if len(row) < 4:
        return {"locus": f"Row_{idx}", "fam": "", "hex": "", "common": "",
                "status": "Failed", "reason": "Insufficient columns", "_index": idx}

    locus_name, allele1, allele2, raw_seq = row[0], row[1], row[2], row[3]

    # Clean sequence
    seq = raw_seq.upper().replace(" ", "").replace("-", "").replace("_", "")

    # Basic input validation
    if seq.count('[') != 1 or seq.count(']') != 1 or seq.count('/') != 1:
        return {"locus": locus_name, "fam": "", "hex": "", "common": "",
                "status": "Failed", "reason": "MALFORMED_INPUT: missing [ ] or /", "_index": idx}
    if allele1 == allele2:
        return {"locus": locus_name, "fam": "", "hex": "", "common": "",
                "status": "Failed", "reason": "IDENTICAL_ALLELES", "_index": idx}

    # Extract flanks and alleles
    pos_open = seq.find('[')
    pos_close = seq.find(']')
    pos_slash = seq.find('/')
    left_flank_raw = seq[:pos_open]
    right_flank_raw = seq[pos_close+1:]
    seq_allele1 = seq[pos_open+1:pos_slash]
    seq_allele2 = seq[pos_slash+1:pos_close]
    allele1 = seq_allele1
    allele2 = seq_allele2

    # Flank length check (must be ≥40 bp for the search windows)
    if len(left_flank_raw) < 40 or len(right_flank_raw) < 40:
        return {"locus": locus_name, "fam": "", "hex": "", "common": "",
                "status": "Failed",
                "reason": "FLANK_TOO_SHORT (need ≥40 bp on each side)", "_index": idx}

    # Check for ambiguous IUPAC codes in the 40‑bp search regions
    check_left = left_flank_raw[-40:]
    check_right = right_flank_raw[:40]
    ambiguous = set("RYSWKM")
    if any(c in ambiguous for c in check_left + check_right):
        return {"locus": locus_name, "fam": "", "hex": "", "common": "",
                "status": "Failed", "reason": "AMBIGUOUS_BASE_IN_FLANK (R,Y,S,W,K,M)", "_index": idx}

    # Assemble full sequence (without any padding) for common primer extraction
    full_seq = left_flank_raw + '[' + allele1 + '/' + allele2 + ']' + right_flank_raw
    var_start = len(left_flank_raw)                     # index of '['
    var_end = var_start + len(allele1) + len(allele2) + 3   # index of ']'

    # Decide direction and try both orientations
    FR = decide_direction(allele1, allele2, left_flank_raw, right_flank_raw)
    directions_to_try = [FR, -FR]
    last_error = ""

    for attempt_dir in directions_to_try:
        forward = (attempt_dir == 1)

        # Design FAM primer
        fam_seq, fam_tm, err_fam = design_specific_primer(forward, allele1,
                                                          left_flank_raw, right_flank_raw)
        if fam_seq is None:
            last_error = err_fam
            continue

        # Design HEX primer (with length constraints based on FAM)
        comp_pairs = [('A','T'), ('T','A'), ('G','C'), ('C','G')]
        m2_val = 0 if (allele1, allele2) in comp_pairs else 3
        hex_seq, hex_tm, err_hex = design_specific_primer(forward, allele2,
                                                          left_flank_raw, right_flank_raw,
                                                          partner_seq=fam_seq,
                                                          fam_length=len(fam_seq),
                                                          m2=m2_val)
        if hex_seq is None:
            last_error = err_hex
            continue

        # Design common primer
        common_seq, err_com, common_start, common_len = design_common_primer(
            forward, left_flank_raw, right_flank_raw,
            var_start, var_end, full_seq,
            specific_fam=fam_seq, specific_hex=hex_seq
        )
        if common_seq is None:
            last_error = err_com
            continue

        # Distinguishability: last 5 bases must differ
        if fam_seq[-5:] == hex_seq[-5:]:
            last_error = "FAM and HEX primers share last 5 bases – indistinguishable"
            continue

        # ---- Product size check (≤ 100 bp) ----
        if forward:
            # SNP position = end of allele1 inside the brackets
            snp_pos = pos_slash - 1
            product_end = common_start + common_len
            product_len = product_end - snp_pos
        else:
            # SNP position = start of allele1
            snp_pos = var_start + 1
            product_len = snp_pos - common_start + 1

        if product_len > 100:
            last_error = f"Product size too large ({product_len} bp > 100 bp maximum)"
            continue

        # Success!
        return {"locus": locus_name,
                "fam": FAM_TAIL + fam_seq,
                "hex": HEX_TAIL + hex_seq,
                "common": common_seq,
                "status": "Success",
                "reason": "",
                "_index": idx}

    # Both directions failed
    return {"locus": locus_name, "fam": "", "hex": "", "common": "",
            "status": "Failed",
            "reason": last_error or "No valid primer set after forward/reverse attempts",
            "_index": idx}

# =============================================================================
# 11. Batch processing function (used by Streamlit and CLI)
# =============================================================================
def process_rows(rows: List[Tuple[int, List[str]]], max_workers: int = 4) -> List[Dict]:
    """
    Process a list of rows (each row is (original_index, [col0, col1, col2, col3])).
    Returns a list of result dictionaries sorted by the original index.
    This function is thread‑safe and can be called from Streamlit.
    """
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_row = {executor.submit(design_one_locus, row): row for row in rows}
        for future in concurrent.futures.as_completed(future_to_row):
            try:
                res = future.result()
                results.append(res)
            except Exception as e:
                # Catch any unexpected exception from a locus
                row = future_to_row[future]
                idx = row[0]
                results.append({
                    "locus": f"Row_{idx}", "fam": "", "hex": "", "common": "",
                    "status": "Failed", "reason": f"Exception: {str(e)}",
                    "_index": idx
                })

    # Sort by original order (preserve the sequence as in the input CSV)
    results.sort(key=lambda r: r["_index"])
    # Remove the internal _index key before returning (optional – caller might want it)
    for r in results:
        del r["_index"]
    return results