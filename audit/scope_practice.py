#!/usr/bin/env python3
"""Scope-reporting practice among the audited DIA-NN studies (NEW, revision).

Answers the reviewer question "do DIA-NN users really build experiment-wide
lists from run-wise q-values?" from the scoring sheet: how many DIA-NN studies
state a scope at all, and how many of those that do filter only on run-wise
columns (Q.Value / Lib.Q.Value) while reporting experiment-wide counts.

    python scope_practice.py --sheet scoring.csv

Prints the numbers quoted in the revised Section 2 and writes
out/scope_practice.txt.
"""
import argparse, csv, os, re

RUNWISE = re.compile(r"\bQ\.Value\b|Lib\.Q\.Value|run-speci", re.I)
GLOBAL  = re.compile(r"Global|experiment-wide|matrix-qvalue|Global\.Q|global precursor|global protein", re.I)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sheet", default="scoring.csv")
    ap.add_argument("--out-dir", default="out")
    a = ap.parse_args()
    rows = [r for r in csv.DictReader(open(a.sheet, encoding="utf-8")) if r["included"] == "Yes"]
    diann = [r for r in rows if "DIA-NN" in r["software"]]
    scope_yes = [r for r in diann if r["scope"] == "Yes"]
    scope_any = [r for r in diann if r["scope"] in ("Yes", "Partial")]
    runwise_only = [r for r in scope_any
                    if RUNWISE.search(r["q_scope"]) and not GLOBAL.search(r["q_scope"])]
    lines = [
        f"audited studies: {len(rows)}",
        f"DIA-NN users: {len(diann)}",
        f"  scope stated (Yes): {len(scope_yes)}; Yes or Partial: {len(scope_any)}; "
        f"no scope stated: {len(diann) - len(scope_any)}",
        f"  of those stating a scope, filtering on run-wise columns only "
        f"(Q.Value / Lib.Q.Value, no global column): {len(runwise_only)}",
    ]
    for r in runwise_only:
        lines.append(f"    {r['id']} {r['first_author']} {r['year']} {r['journal']} | "
                     f"count_any={r['count_any']} | {r['q_scope'][:200]}")
    all_scope = sum(r["scope"] == "Yes" for r in rows)
    all_q = sum(r["qfield"] == "Yes" for r in rows)
    lines.append(f"all studies: scope stated {all_scope}/{len(rows)}, column named {all_q}/{len(rows)}")
    os.makedirs(a.out_dir, exist_ok=True)
    open(os.path.join(a.out_dir, "scope_practice.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))

if __name__ == "__main__":
    main()
