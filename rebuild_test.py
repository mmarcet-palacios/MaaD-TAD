"""
Validation: can a base-pair-step builder reproduce a real nucleosome's DNA?

Rebuild A: real base-pair atoms replaced into the rebuilt frames.
           Tests the frame/step machinery only -- should be ~0.
Rebuild B: one canonical template per Watson-Crick pair type, placed into the
           rebuilt frames. This is what the actual builder does; it has no
           access to the real per-base-pair coordinates. RMSD here is the
           honest cost of ignoring intra-base-pair variation (propeller,
           buckle, opening) and of reusing one sugar conformation per type.
"""
import numpy as np, gemmi
import bpstep as bs

CIF = __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'data', '1KX5.cif')

# nucleoside atoms only -- the linking phosphate is placed analytically later
SUGAR = ["C1'", "C2'", "C3'", "C4'", "C5'", "O3'", "O4'", "O5'"]
BASE = {
    "A": ["N9","C8","N7","C5","C6","N6","N1","C2","N3","C4"],
    "G": ["N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"],
    "C": ["N1","C2","O2","N3","C4","N4","C5","C6"],
    "T": ["N1","C2","O2","N3","C4","O4","C5","C6","C7"],
}

def nucleoside_names(letter):
    return BASE[letter] + SUGAR


def load():
    st = gemmi.read_structure(CIF); st.setup_entities()
    return st[0]


def collect(model):
    idx = list(range(-73, 74))
    frames, seq, kept = bs.chain_frames(model, "I", "J", lambda i: -i, idx)
    ref = {r.seqid.num: r for r in model["I"] if bs.is_dna(r.name)}
    par = {r.seqid.num: r for r in model["J"] if bs.is_dna(r.name)}
    return frames, seq, kept, ref, par


def bp_local_atoms(frame, res_ref, res_par):
    """Atoms of a base pair expressed in its own frame, as an ordered dict."""
    out = {}
    for tag, res in (("R", res_ref), ("P", res_par)):
        atoms = bs.atom_dict(res)
        for n in nucleoside_names(bs.base_letter(res.name)):
            if n in atoms:
                out[(tag, n)] = bs.to_local(frame, atoms[n][None, :])[0]
    return out


def build_templates(frames, kept, ref, par):
    """One canonical nucleoside pair per WC type, chosen as the medoid of all
    real instances of that type expressed in their own bp frames."""
    buckets = {}
    for f, (i, j) in zip(frames, kept):
        letter = bs.base_letter(ref[i].name)
        buckets.setdefault(letter, []).append(bp_local_atoms(f, ref[i], par[j]))

    templates, report = {}, {}
    for letter, insts in buckets.items():
        common = set(insts[0])
        for d in insts[1:]:
            common &= set(d)
        common = sorted(common)
        arr = np.array([[d[k] for k in common] for d in insts])   # (n, atoms, 3)
        mean = arr.mean(axis=0)
        dists = np.sqrt(((arr - mean) ** 2).sum(-1).mean(-1))
        medoid = int(np.argmin(dists))
        templates[letter] = {k: arr[medoid, a] for a, k in enumerate(common)}
        report[letter] = (len(insts), float(dists.mean()), float(dists[medoid]))
    return templates, report


def place(frames, seq, templates):
    """Place template nucleosides into every frame -> world coordinates."""
    out = []
    for f, letter in zip(frames, seq):
        t = templates[letter]
        keys = sorted(t)
        local = np.array([t[k] for k in keys])
        world = bs.to_world(f, local)
        out.append({k: world[a] for a, k in enumerate(keys)})
    return out


def main():
    model = load()
    frames, seq, kept, ref, par = collect(model)
    params = np.array([bs.decompose_step(frames[i], frames[i+1])
                       for i in range(len(frames)-1)])

    rebuilt = bs.rebuild_frames(frames[0], params)
    fr_err = max(np.abs(a[0]-b[0]).max() for a, b in zip(rebuilt, frames))
    print(f"frame reconstruction max origin error: {fr_err:.2e} A")

    templates, report = build_templates(frames, kept, ref, par)
    print("\ncanonical templates (medoid of real instances, in own bp frame):")
    for letter in "ATGC":
        if letter in report:
            n, meanspread, medspread = report[letter]
            print(f"  {letter}: n={n:3d}  mean spread from centroid {meanspread:.2f} A"
                  f"   medoid spread {medspread:.2f} A")

    # --- Rebuild A: real atoms into rebuilt frames -------------------------
    realA, rebA = [], []
    for f_orig, f_new, (i, j) in zip(frames, rebuilt, kept):
        loc = bp_local_atoms(f_orig, ref[i], par[j])
        keys = sorted(loc)
        realA.append(bs.to_world(f_orig, np.array([loc[k] for k in keys])))
        rebA.append(bs.to_world(f_new,  np.array([loc[k] for k in keys])))
    realA, rebA = np.vstack(realA), np.vstack(rebA)
    print(f"\nRebuild A (real atoms, machinery test)   n={len(realA):5d} atoms")
    print(f"  RMSD (no superposition): {bs.kabsch_rmsd(rebA, realA, superpose=False):.4f} A")

    # --- Rebuild B: canonical templates into rebuilt frames ----------------
    placed = place(rebuilt, seq, templates)
    realB, rebB = [], []
    for pl, (i, j) in zip(placed, kept):
        a_ref, a_par = bs.atom_dict(ref[i]), bs.atom_dict(par[j])
        for (tag, n), xyz in pl.items():
            src = a_ref if tag == "R" else a_par
            if n in src:
                rebB.append(xyz); realB.append(src[n])
    realB, rebB = np.array(realB), np.array(rebB)
    print(f"\nRebuild B (canonical templates, real builder)  n={len(realB):5d} atoms")
    print(f"  RMSD (no superposition): {bs.kabsch_rmsd(rebB, realB, superpose=False):.3f} A")
    print(f"  RMSD (superposed)      : {bs.kabsch_rmsd(rebB, realB, superpose=True):.3f} A")

    # per-atom-class breakdown, to see where the error actually lives
    classes = {"base": [], "sugar": []}
    for pl, (i, j) in zip(placed, kept):
        a_ref, a_par = bs.atom_dict(ref[i]), bs.atom_dict(par[j])
        for (tag, n), xyz in pl.items():
            src = a_ref if tag == "R" else a_par
            if n in src:
                classes["sugar" if n in SUGAR else "base"].append(
                    np.linalg.norm(xyz - src[n]))
    for k, v in classes.items():
        v = np.array(v)
        print(f"    {k:5s}: mean dev {v.mean():.3f} A   n={len(v)}")

    np.save("rebuilt_frames.npy", np.array([np.vstack([o, R]) for o, R in rebuilt]))
    import pickle
    with open("templates.pkl", "wb") as fh:
        pickle.dump(templates, fh)
    with open("rebuild_context.pkl", "wb") as fh:
        pickle.dump({"seq": seq, "kept": kept, "params": params}, fh)
    print("\nsaved rebuilt_frames.npy, templates.pkl, rebuild_context.pkl")


if __name__ == "__main__":
    main()
