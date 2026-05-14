"""
Tests for get_frequencies() — the mpileup base-call parser.

Regression coverage for the indel nucleotide bleed-through bug (issue #31):
  '-' and '+' are keys in fastadict, so without the explicit priority check
  the generic `if char in fastadict` branch would fire first, incrementing
  the '-'/'+' counter by 1 while leaving the subsequent nucleotide characters
  (e.g. the "ACG" in "-3ACG") to walk through the loop one-by-one and
  accumulate into A/T/G/C counts.
"""

import sys
from pathlib import Path

# Allow running from repo root without installing the package
sys.path.insert(0, str(Path(__file__).parent.parent))

from yleaf.Yleaf import get_frequencies


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _acgt(d):
    """Return (A, T, G, C) tuple from a get_frequencies result dict."""
    return d["A"], d["T"], d["G"], d["C"]


# ---------------------------------------------------------------------------
# Basic SNP-only sequences (no indels)
# ---------------------------------------------------------------------------

def test_simple_bases():
    result = get_frequencies("AATGC")
    assert result["A"] == 2
    assert result["T"] == 1
    assert result["G"] == 1
    assert result["C"] == 1
    assert result["-"] == 0
    assert result["+"] == 0


def test_case_insensitive():
    result = get_frequencies("aAtTgGcC")
    assert result["A"] == 2
    assert result["T"] == 2
    assert result["G"] == 2
    assert result["C"] == 2


def test_caret_skips_next_char():
    # '^' + mapping quality char should be consumed without counting either
    result = get_frequencies("A^]T")
    assert result["A"] == 1
    assert result["T"] == 1
    assert sum(_acgt(result)) == 2


# ---------------------------------------------------------------------------
# Deletion notation  (bug #31 regression)
# ---------------------------------------------------------------------------

def test_deletion_nucleotides_do_not_bleed():
    """
    '-3ACG' means a 3-bp deletion; 'ACG' must NOT be added to A/C/G counts.
    Before the fix: A+=1, C+=1, G+=1 (and -+=1 from the fastadict fast-path).
    After the fix:  no A/C/G increment; deletion event counted once in '-'.
    """
    result = get_frequencies("-3ACG")
    assert _acgt(result) == (0, 0, 0, 0), (
        "Deletion nucleotides bled into base counts — bug #31 not fixed"
    )
    assert result["-"] == 1  # one deletion event


def test_deletion_mixed_with_real_bases():
    """Real T reads plus a deletion; only the T count should increase."""
    # 4× T reads, 1× deletion of 2 bp (AT)
    result = get_frequencies("TTTT-2AT")
    assert result["T"] == 4, "Real T reads miscounted"
    assert result["A"] == 0, "Deletion 'A' bled into A count"
    assert result["-"] == 1, "Deletion event not counted"


def test_multi_digit_deletion():
    """Deletions with a 2-digit length: -12ACGTACGTACGT should skip 12 chars."""
    result = get_frequencies("A-12ACGTACGTACGTG")
    # Only the leading 'A' and trailing 'G' are real bases
    assert result["A"] == 1
    assert result["G"] == 1
    assert result["T"] == 0
    assert result["C"] == 0
    assert result["-"] == 1


# ---------------------------------------------------------------------------
# Insertion notation  (bug #31 regression)
# ---------------------------------------------------------------------------

def test_insertion_nucleotides_do_not_bleed():
    """+2AC means a 2-bp insertion; 'AC' must not bump A/C counts."""
    result = get_frequencies("+2AC")
    assert _acgt(result) == (0, 0, 0, 0)
    assert result["+"] == 1


def test_insertion_mixed_with_real_bases():
    """3× C reads, 1× insertion of 3 bp (TAG)."""
    result = get_frequencies("CCC+3TAG")
    assert result["C"] == 3
    assert result["T"] == 0, "Insertion 'T' bled into T count"
    assert result["A"] == 0, "Insertion 'A' bled into A count"
    assert result["G"] == 0, "Insertion 'G' bled into G count"
    assert result["+"] == 1


# ---------------------------------------------------------------------------
# Deletion-spanning '*' notation
# ---------------------------------------------------------------------------

def test_star_counted_as_deletion():
    """'*' at a position means the read has a deletion spanning this base."""
    result = get_frequencies("AA**")
    assert result["A"] == 2
    assert result["-"] == 2   # two spanning-deletion reads
    assert result["+"] == 0


def test_deletion_start_and_star_accumulate():
    """-2AT and two '*' at another position should both contribute to '-'."""
    result = get_frequencies("-2AT**")
    # deletion start: -1, two stars: +2  → total "-" == 3
    assert result["-"] == 3
    assert _acgt(result) == (0, 0, 0, 0)


# ---------------------------------------------------------------------------
# called_percentage denominator sanity (motivating case from issue #31)
# ---------------------------------------------------------------------------

def test_called_percentage_denominator():
    """
    Simulate the find_private_mutations scenario:
    10 reads call 'A', 5 reads have a 3-bp deletion '-3ACG'.
    Before fix: A=10+5=15, C=5, G=5, -=5 → sum=30, pct=50 % (wrong)
    After fix:  A=10,        -=5          → sum=15, pct=67 % (correct)
    """
    # 10× 'A' + 5× '-3ACG'
    sequence = "A" * 10 + "-3ACG" * 5
    result = get_frequencies(sequence)
    total = sum(result.values())
    pct_A = result["A"] / total * 100
    assert result["A"] == 10, f"Expected A=10, got {result['A']}"
    assert result["C"] == 0,  f"Deletion 'C' bled into C count"
    assert result["G"] == 0,  f"Deletion 'G' bled into G count"
    assert result["-"] == 5,  f"Expected -=5, got {result['-']}"
    assert abs(pct_A - 66.67) < 0.1, f"Unexpected called_percentage: {pct_A:.2f}"
