"""Regenerate the generated blocks of README.md from the installed files.

The flagged-element table, the per-set counts and the description of the
radial grid are read from the datasets' own records, never typed by hand:

    python tools/readme_tables.py            # print the blocks
    python tools/readme_tables.py --write    # rewrite them in README.md

Each block sits between ``<!-- generated: NAME -->`` and
``<!-- end generated: NAME -->`` in README.md.  Run it after any rebuild.
"""
import argparse
import os
import re
import warnings

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SETS = ("lda-sr", "lda-dirac", "pbe-sr")
LETTERS = "spdfgh"


def load(folder):
    from mandacaru.pseudopotentials.io import (available_elements,
                                               load_pseudopotential)
    from ase.data import atomic_numbers

    path = os.path.join(ROOT, folder)
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for symbol in sorted(available_elements(path),
                             key=lambda s: atomic_numbers[s]):
            out[symbol] = load_pseudopotential(
                os.path.join(path, f"{symbol}.parquet"))
    return out


def describe(defects):
    """One table cell: every recorded defect, by kind."""
    if not defects:
        return "—"
    parts = [f"{LETTERS[l]} ghost {e:+.3g} Ha"
             for l, e in sorted(defects.get("ghosts", {}).items())]
    if defects.get("s_miss") is not None:
        parts.append(f"1s miss {defects['s_miss']:.2f}")
    parts += [f"{LETTERS[l]} {near:.2f}/{far:.2f} rad"
              for l, (near, far) in sorted(defects.get("phases", {}).items())]
    parts += [f"{LETTERS[l]} (no projectors) {e:.2f} rad"
              for l, e in sorted(defects.get("unprojected", {}).items())]
    parts += [f"{LETTERS[l]} residual {e:.2g} Ha"
              for l, e in sorted(defects.get("residuals", {}).items())]
    return "; ".join(parts)


def flagged_block(sets):
    from ase.data import atomic_numbers

    lines = []
    for folder in SETS:
        data = sets[folder]
        flagged = [s for s, pp in data.items() if pp.defects]
        ghosts = [s for s in flagged if data[s].defects.get("ghosts")]
        lines.append(f"- `{folder}/`: {len(flagged)} of {len(data)} "
                     f"flagged, {len(ghosts) or 'none'} for a ghost state.")
    lines += ["", "| Element | " + " | ".join(f"`{f}/`" for f in SETS) + " |",
              "|---" * (len(SETS) + 1) + "|"]
    symbols = sorted({s for f in SETS for s, pp in sets[f].items()
                      if pp.defects}, key=lambda s: atomic_numbers[s])
    for symbol in symbols:
        cells = [describe(sets[f][symbol].defects) if symbol in sets[f]
                 else "not built" for f in SETS]
        lines.append(f"| {symbol} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def grid_block(sets):
    sizes = {}
    for folder in SETS:
        for symbol, pp in sets[folder].items():
            sizes[symbol] = (pp.r.size, float(pp.r[1] - pp.r[0]),
                             float(pp.r[0]), float(pp.r[-1]))
    first = min(sizes.values())
    last = max(sizes.values())
    return (f"- a uniform radial grid starting one step from the origin and "
            f"ending at {first[3]:.0f} Bohr, finer for heavier atoms: "
            f"{first[0]} points {first[1]:.4f} Bohr apart (hydrogen) to "
            f"{last[0]} points {last[1]:.5f} Bohr apart (uranium) -- every "
            f"fourth point of the grid the reference atom was solved on -- "
            f"with every quantity in atomic units.")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--write", action="store_true",
                        help="rewrite the blocks in README.md")
    args = parser.parse_args()
    sets = {folder: load(folder) for folder in SETS}
    blocks = {"flagged": flagged_block(sets), "grid": grid_block(sets)}
    if not args.write:
        for name, text in blocks.items():
            print(f"<!-- generated: {name} -->\n{text}\n"
                  f"<!-- end generated: {name} -->\n")
        return
    path = os.path.join(ROOT, "README.md")
    with open(path, encoding="utf-8") as handle:
        readme = handle.read()
    for name, text in blocks.items():
        pattern = re.compile(rf"(<!-- generated: {name} -->\n)(?:.*?\n)?"
                             rf"(<!-- end generated: {name} -->)", re.S)
        if not pattern.search(readme):
            raise SystemExit(f"README.md has no generated block {name!r}")
        readme = pattern.sub(lambda m: m.group(1) + text + "\n" + m.group(2),
                             readme)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(readme)


if __name__ == "__main__":
    main()
