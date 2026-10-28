"""
Stage 2 driver: build all blocks in parallel. Resumable -- finished blocks are
skipped, so just run it again after an interruption.

    python run_stage2.py --workers 24                 # everything
    python run_stage2.py --workers 1 --only 0         # pilot: one block, timed
    python run_stage2.py --status                     # progress so far

Each block runs as its own process (so a crash in one cannot take down the
others) with its own log in blocks/logs/. BLAS threading is pinned to one
thread per worker; otherwise 24 workers x 32 threads would thrash the CPU.
"""
import argparse, csv, glob, hashlib, json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np


def block_sig(F, b, size, tethers, halo=6):
    """What a block's atoms depend on: frames of its nucleosomes +-halo (obstacles),
    their repeat lengths, and any tether waypoints on its linkers."""
    O, R, nrl = F["origins"], F["rotations"], F["nrl"]
    n = len(O); s0 = b * size; e0 = min(s0 + size, n) - 1
    lo, hi = max(0, s0 - halo), min(n, e0 + 2 + halo)
    h = hashlib.md5()
    h.update(np.round(O[lo:hi], 3).tobytes()); h.update(np.round(R[lo:hi], 6).tobytes())
    h.update(np.asarray(nrl[lo:min(hi, n - 1)]).tobytes())
    for t in tethers:
        for L in t["linkers"]:
            if s0 <= L["linker"] <= e0:
                h.update(json.dumps(L, sort_keys=True).encode())
                h.update(TETHER_RULES.encode())     # tether blocks rebuild when their build rules change
    return h.hexdigest()


TETHER_RULES = "v2: avoid histone CAs (tails) + partner's built P atoms"


def tether_waves(tethers, size):
    """{block of the higher-index tether linker: block of its partner}. The higher one is
    built second, so it can steer around the partner's finished atoms."""
    dep = {}
    for t in tethers:
        a_, b_ = sorted(int(L["linker"]) for L in t["linkers"])
        if a_ // size != b_ // size:
            dep[b_ // size] = a_ // size
    return dep


def status(outdir, nblocks):
    done = sorted(glob.glob(os.path.join(outdir, "block_??.cif")))
    rows = []
    for f in glob.glob(os.path.join(outdir, "block_??_qc.csv")):
        if os.path.exists(f.replace("_qc.csv", ".cif")):
            rows += list(csv.DictReader(open(f)))
    fail = [r for r in rows if r.get("ok") != "True"]
    bend = np.array([float(r["bend_max_deg"]) for r in rows]) if rows else np.array([0.0])
    secs = np.array([float(r["seconds"]) for r in rows]) if rows else np.array([0.0])
    print(f"blocks done {len(done)}/{nblocks}; linkers {len(rows)}; "
          f"failed closure {len(fail)}; linkers peaking >13 deg/step {(bend > 13).sum()}; "
          f"CPU time so far {secs.sum()/3600:.1f} h")
    if fail:
        print("  failed:", [r["linker"] for r in fail][:20])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", default="frames.npz")
    ap.add_argument("--size", type=int, default=50, help="nucleosomes per block")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 8))
    ap.add_argument("--outdir", default="blocks")
    ap.add_argument("--only", type=int, nargs="*", help="build just these block numbers")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--tethers", default="data/tethers.json")
    ap.add_argument("--adopt", help="frames file that existing unsigned blocks were built from")
    ap.add_argument("--adopt-tethers", help="tethers.json those blocks were built with (default: current)")
    a = ap.parse_args()

    n = len(np.load(a.frames)["origins"])
    nblocks = (n + a.size - 1) // a.size
    if a.status:
        status(a.outdir, nblocks); return

    os.makedirs(os.path.join(a.outdir, "logs"), exist_ok=True)
    # Per-block signatures: a block is rebuilt only if the frames it depends on
    # (its nucleosomes +-6, their repeat lengths, any tether waypoints) changed.
    F = np.load(a.frames)
    tjs = json.load(open(a.tethers)) if a.tethers and os.path.exists(a.tethers) else []
    sig = lambda Fr, b: block_sig(Fr, b, a.size, tjs)
    if a.adopt:
        Fa = np.load(a.adopt); tja = json.load(open(a.adopt_tethers)) if a.adopt_tethers else tjs
        k = 0
        for b in range(nblocks):
            cif = os.path.join(a.outdir, f"block_{b:02d}.cif")
            if os.path.exists(cif) and not os.path.exists(cif[:-4] + ".sig"):
                open(cif[:-4] + ".sig", "w").write(block_sig(Fa, b, a.size, tja)); k += 1
        print(f"adopted {k} existing blocks as built from {a.adopt}")
    stale, legacy = [], []
    for b in range(nblocks):
        cif = os.path.join(a.outdir, f"block_{b:02d}.cif")
        if os.path.exists(cif):
            sf = cif[:-4] + ".sig"
            if not os.path.exists(sf):
                legacy.append(b)
            elif open(sf).read().strip() != sig(F, b):
                stale.append(b)
    if legacy:
        sys.exit(f"{len(legacy)} blocks have no build signature (built by an older version). Rerun with "
                 f"--adopt <the frames file they were built from>, or delete '{a.outdir}'.")
    # a rebuilt first-wave tether block invalidates its partner's block (built around its atoms)
    dep = tether_waves(tjs, a.size)
    for b2, b1 in dep.items():
        c1 = os.path.join(a.outdir, f"block_{b1:02d}.cif")
        if (b1 in stale or not os.path.exists(c1)) and b2 not in stale:
            stale.append(b2)
    for b in stale:
        for ext in (".cif", "_qc.csv", "_seam.json", ".sig"):
            p = os.path.join(a.outdir, f"block_{b:02d}{ext}")
            if os.path.exists(p):
                os.remove(p)
    if stale:
        print(f"frames or tether rules changed for {len(stale)} blocks -> rebuilding: {sorted(stale)}")
    todo = a.only if a.only else list(range(nblocks))
    todo = [b for b in todo
            if not os.path.exists(os.path.join(a.outdir, f"block_{b:02d}.cif"))]
    print(f"{n} nucleosomes -> {nblocks} blocks of {a.size}; "
          f"{len(todo)} to build with {a.workers} workers", flush=True)
    if not todo:
        status(a.outdir, nblocks); return

    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", PYTHONUNBUFFERED="1")
    here = os.path.dirname(os.path.abspath(__file__))

    def run(b):
        t = time.time()
        log = os.path.join(a.outdir, "logs", f"block_{b:02d}.log")
        with open(log, "w") as lf:
            rc = subprocess.call([sys.executable, os.path.join(here, "stage2_block.py"),
                                  "--frames", a.frames, "--block", str(b),
                                  "--size", str(a.size), "--outdir", a.outdir,
                                  "--tethers", a.tethers or ""],
                                 stdout=lf, stderr=subprocess.STDOUT, env=env, cwd=here)
        if rc == 0:
            open(os.path.join(a.outdir, f"block_{b:02d}.sig"), "w").write(sig(F, b))
        return b, rc, time.time() - t

    t0 = time.time(); finished = 0
    wave2 = [b for b in todo if b in dep and dep[b] in todo]
    waves = [[b for b in todo if b not in wave2], wave2]
    if wave2:
        print(f"tether partners: blocks {wave2} wait for {[dep[b] for b in wave2]}", flush=True)
    for wave in waves:
        if not wave:
            continue
        with ThreadPoolExecutor(max_workers=a.workers) as ex:
            futs = [ex.submit(run, b) for b in wave]
            for fu in as_completed(futs):
                b, rc, dt = fu.result(); finished += 1
                el = time.time() - t0
                eta = el / finished * (len(todo) - finished)
                msg = "ok" if rc == 0 else f"FAILED (exit {rc}) -- see {a.outdir}/logs/block_{b:02d}.log"
                print(f"[{time.strftime('%H:%M')}] block {b:02d} {msg} in {dt/60:.1f} min | "
                      f"{finished}/{len(todo)} done, ~{eta/60:.0f} min left", flush=True)
    status(a.outdir, nblocks)


if __name__ == "__main__":
    main()
