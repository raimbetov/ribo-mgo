#!/usr/bin/env python3
"""Measure the column geometry of a PDF table, so writing a new spec is a
measurement rather than a guess.

Writing a `pdf_anchor` source needs four things: which pages the table spans,
the x-range of the row-number margin, the x-ranges of each data column, and
whether a row's number sits on its first baseline or somewhere else. This prints
all four.

    # 1. find the pages a table spans (0-based indices, as specs use)
    python3 measure_columns.py doc.pdf --find "Table S3"

    # 2. dump one page's baselines to see how rows are laid out
    python3 measure_columns.py doc.pdf --page 12 --baselines

    # 3. histogram word x-positions to read off column boundaries
    python3 measure_columns.py doc.pdf --pages 12 18

Read the boundaries off the gaps in the histogram: a column spans a run of
occupied x-bins, and the empty bins between runs are where to cut.

CHOOSING `assign`. Look at the --baselines output. If every row's number shares a
baseline with the start of its protein name, use `assign: below`. If some rows
print the number on a LATER baseline than the name (Irshad's rows 1-4 do), use
`assign: nearest` -- `below` would attach those names to the previous row.
"""
import argparse
import re
import sys
from collections import Counter, defaultdict

try:
    import pdfplumber
except ImportError:
    sys.exit("pdfplumber is required: pip install -r ../requirements.txt")


def find_pages(pdf, needle):
    hits = []
    for i, page in enumerate(pdf.pages):
        text = page.extract_text() or ""
        if needle.lower() in text.lower():
            first = next((l.strip() for l in text.splitlines()
                          if needle.lower() in l.lower()), "")
            hits.append((i, first[:90]))
    return hits


def baselines(page, band, limit):
    bands = defaultdict(list)
    for w in page.extract_words():
        bands[round(w["top"] / band)].append(w)
    for k in sorted(bands)[:limit]:
        ws = sorted(bands[k], key=lambda w: w["x0"])
        line = " | ".join(f"{w['text']}@{w['x0']:.0f}" for w in ws)
        print(f"  y{k * band:7.0f}: {line[:190]}")


def histogram(pdf, lo, hi, binw):
    counts, samples = Counter(), defaultdict(list)
    for pno in range(lo, hi + 1):
        if pno >= len(pdf.pages):
            continue
        for w in pdf.pages[pno].extract_words():
            b = int(w["x0"] // binw) * binw
            counts[b] += 1
            if len(samples[b]) < 3:
                samples[b].append(w["text"][:18])
    if not counts:
        print("  no words found on those pages")
        return
    print(f"  x-bin  count  examples          (bin width {binw})")
    prev_empty = False
    for b in range(0, int(max(counts)) + binw, binw):
        n = counts.get(b, 0)
        if n == 0:
            if not prev_empty:
                print("  ----   ----   <gap: a column boundary belongs here>")
            prev_empty = True
            continue
        prev_empty = False
        bar = "#" * min(40, n // max(1, max(counts.values()) // 40 or 1))
        print(f"  {b:5d}  {n:5d}  {', '.join(samples[b]):34s} {bar}")


def anchors(pdf, lo, hi, xlo, xhi):
    """Report the row numbers found in a candidate anchor margin."""
    found = []
    for pno in range(lo, hi + 1):
        if pno >= len(pdf.pages):
            continue
        ns = sorted(int(w["text"]) for w in pdf.pages[pno].extract_words()
                    if xlo <= w["x0"] < xhi and re.fullmatch(r"\d{1,4}", w["text"]))
        found.append((pno, len(ns), ns[:1], ns[-1:]))
    print(f"  anchors in x=[{xlo},{xhi}):")
    tot = 0
    for pno, n, first, last in found:
        tot += n
        print(f"    page {pno}: {n:4d} numbers, range {first}..{last}")
    print(f"    total {tot}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdf")
    ap.add_argument("--find", help="locate pages containing this text")
    ap.add_argument("--page", type=int, help="single page (0-based) to dump")
    ap.add_argument("--pages", nargs=2, type=int, metavar=("LO", "HI"))
    ap.add_argument("--baselines", action="store_true",
                    help="print word baselines for --page")
    ap.add_argument("--band", type=float, default=2.0)
    ap.add_argument("--limit", type=int, default=25)
    ap.add_argument("--bin", type=int, default=10, help="histogram bin width")
    ap.add_argument("--anchor-x", nargs=2, type=float, metavar=("LO", "HI"),
                    help="test a candidate row-number margin")
    args = ap.parse_args()

    with pdfplumber.open(args.pdf) as pdf:
        print(f"{args.pdf}: {len(pdf.pages)} pages, "
              f"{pdf.pages[0].width:.0f}x{pdf.pages[0].height:.0f} pt\n")

        if args.find:
            for i, line in find_pages(pdf, args.find):
                print(f"  page {i} (0-based): {line}")
            return 0

        if args.page is not None and args.baselines:
            print(f"baselines on page {args.page} (band {args.band}):")
            baselines(pdf.pages[args.page], args.band, args.limit)
            return 0

        lo, hi = args.pages if args.pages else (args.page, args.page)
        if lo is None:
            ap.error("give --find, --page, or --pages")
        if args.anchor_x:
            anchors(pdf, lo, hi, *args.anchor_x)
            print()
        print(f"word x-position histogram, pages {lo}-{hi}:")
        histogram(pdf, lo, hi, args.bin)
    return 0


if __name__ == "__main__":
    sys.exit(main())
