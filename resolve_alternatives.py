"""
Resolve alternative nucleosome positions, then re-fill.

The problem
-----------
130 of the 3,420 linkers in the full-coverage map are shorter than 163 bp, which
the b2 zigzag cannot represent at all, and they are the entire cause of the
audit's 3% failure rate. 53 of them are measured-measured pairs sitting at or
near the 147 bp exclusion floor the peak caller enforced.

Those pairs are not two nucleosomes on one fibre. MNase occupancy is a
POPULATION average over millions of cells, so two dyads called 150 bp apart
almost always mean two mutually exclusive positions - the nucleosome sits here
in some cells and there in others. A single fibre can hold one or the other.
Placing both is a category error, not a geometry problem, and no amount of
fibre flexibility fixes it.

So each sub-163 pair is resolved to one position. Between two measured dyads the
stronger peak wins, because peak height is occupancy and the higher peak is the
position the larger fraction of cells actually used. A measured dyad always
beats a placeholder. Two placeholders cannot normally collide, since the fill
places them at 190 bp, but if it happens the one nearer a measurement is kept.

Dropping a nucleosome opens a gap, so the fill rule is re-applied afterwards and
the whole thing iterates until no linker is under the floor. Every dropped
measurement is recorded with the position that displaced it.
"""
import argparse
import numpy as np

MIN_NRL = 163            # smallest repeat the zigzag can represent
TARGET_NRL = 190
FILL_ABOVE = 285
CORE_BP, HALF = 147, 73
ROI_START, ROI_END = 12_771_273, 13_423_272
MAPPED_END = 12_954_295


def fill_between(p0, p1, nrl=TARGET_NRL):
    gap = p1 - p0
    if gap <= FILL_ABOVE:
        return []
    k = int(round(gap / nrl)) - 1
    while k > 0 and gap / (k + 1) < MIN_NRL:
        k -= 1
    if k <= 0:
        return []
    return [int(round(p0 + gap * (i + 1) / (k + 1))) for i in range(k)]


def load(path):
    out = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        sc = f[4]
        out.append(dict(pos=int(f[6]), src=f[7],
                        score=float(sc) if sc not in (".", "") else np.nan))
    return sorted(out, key=lambda r: r["pos"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inp", default="nucleosomes_GRCh38_full_coverage.bed")
    ap.add_argument("--out", default="nucleosomes_GRCh38_resolved.bed")
    ap.add_argument("--dropped", default="nucleosomes_alternative_positions.tsv")
    args = ap.parse_args()

    nuc = load(args.inp)
    n0 = len(nuc)
    r0 = sum(1 for r in nuc if r["src"] == "real")
    dropped = []

    for sweep in range(20):
        # ---- resolve every pair under the floor -------------------------
        changed = False
        i = 0
        while i < len(nuc) - 1:
            gap = nuc[i + 1]["pos"] - nuc[i]["pos"]
            if gap >= MIN_NRL:
                i += 1
                continue
            a, b = nuc[i], nuc[i + 1]
            if a["src"] != b["src"]:
                lose_i = i if a["src"] == "placeholder" else i + 1
            elif a["src"] == "real":
                # stronger peak wins; occupancy is what the score measures
                lose_i = i if (np.nan_to_num(a["score"]) <
                               np.nan_to_num(b["score"])) else i + 1
            else:
                lose_i = i + 1
            loser = nuc[lose_i]
            keeper = nuc[i + 1 if lose_i == i else i]
            if loser["src"] == "real":
                dropped.append((loser["pos"], loser["score"], keeper["pos"],
                                keeper["src"],
                                np.nan_to_num(keeper["score"]), gap))
            nuc.pop(lose_i)
            changed = True

        # ---- re-fill the gaps that opened --------------------------------
        new = []
        prev = ROI_START - TARGET_NRL
        for r in nuc:
            for p in fill_between(prev, r["pos"]):
                if ROI_START <= p <= ROI_END:
                    new.append(dict(pos=p, src="placeholder", score=np.nan))
            new.append(r)
            prev = r["pos"]
        for p in fill_between(prev, ROI_END + TARGET_NRL):
            if ROI_START <= p <= ROI_END:
                new.append(dict(pos=p, src="placeholder", score=np.nan))
        added = len(new) - len(nuc)
        nuc = sorted(new, key=lambda r: r["pos"])

        g = np.diff([r["pos"] for r in nuc])
        bad = int((g < MIN_NRL).sum())
        print(f"  sweep {sweep+1}: {len(nuc)} nucleosomes, "
              f"{added} re-filled, {bad} still under {MIN_NRL} bp")
        if not changed and bad == 0:
            break

    pos = np.array([r["pos"] for r in nuc])
    g = np.diff(pos)
    n_real = sum(1 for r in nuc if r["src"] == "real")

    with open(args.out, "w") as fh:
        fh.write(f"# GRCh38 chr1:{ROI_START}-{ROI_END}\n")
        fh.write(f"# {len(nuc)} nucleosomes: {n_real} measured, "
                 f"{len(nuc)-n_real} placeholder\n")
        fh.write(f"# alternative positions resolved; min NRL {int(g.min())} bp\n")
        fh.write("# chrom\tstart\tend\tname\tscore\tstrand\tcenter\tsource\tregion\n")
        for i, r in enumerate(nuc, start=1):
            nm = f"nuc_{i}" if r["src"] == "real" else f"gapfill_{i}"
            sc = f"{r['score']:.3f}" if r["src"] == "real" and not np.isnan(r["score"]) else "."
            reg = "mapped" if r["pos"] <= MAPPED_END else "unmapped_assembly"
            fh.write(f"chr1\t{r['pos']-1-HALF}\t{r['pos']-1+HALF+1}\t{nm}\t{sc}\t.\t"
                     f"{r['pos']}\t{r['src']}\t{reg}\n")

    with open(args.dropped, "w") as fh:
        fh.write("dropped_dyad\tdropped_score\tkept_dyad\tkept_source\t"
                 "kept_score\tseparation_bp\n")
        for d in sorted(dropped):
            fh.write("\t".join(str(x) for x in d) + "\n")

    print(f"\nwrote {args.out}")
    print(f"  nucleosomes {n0} -> {len(nuc)}   measured {r0} -> {n_real}")
    print(f"  measured dyads dropped as alternative positions: {len(dropped)}")
    print(f"  NRL: min {int(g.min())}, median {int(np.median(g))}, max {int(g.max())}")
    inw = int(((g >= MIN_NRL) & (g <= 205)).sum())
    print(f"  in 163-205 window : {inw} / {len(g)}  ({100*inw/len(g):.1f}%)")
    print(f"  below {MIN_NRL} bp      : {int((g < MIN_NRL).sum())}")
    print(f"  over 205 bp       : {int((g > 205).sum())}  (slack, still buildable)")
    print(f"\nwrote {args.dropped}")


if __name__ == "__main__":
    main()
