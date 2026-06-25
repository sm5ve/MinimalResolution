#!/usr/bin/env python3
"""Compute all Yoneda products in the given (MD, LEN) range.

Edit MD and LEN below, then:  python3 yoneda_table.py

Output: a dict  {((as,aa),(bs,bb)): result_string}  containing every
nontrivial product alpha*beta, printed to stdout and saved to
<MD>_yoneda_products.txt.

How generator counts are discovered:
  Run yoneda once with any valid beta to collect every {s-a} entry that
  appears in its table output; max(a)+1 for each s gives the number of
  generators at filtration s.  All (bs,bb) pairs with bs>=1 and
  bb < rank[bs] are then enumerated as beta.
"""

import subprocess, sys, os, re, json
from collections import defaultdict

# ── configure here ─────────────────────────────────────────────────────────
MD  = 40
LEN = 30
# ───────────────────────────────────────────────────────────────────────────

DIR = os.path.dirname(os.path.abspath(__file__))
pre = f"{MD}_"


# ── bootstrap ───────────────────────────────────────────────────────────────

def run_bootstrap(cmd, desc):
    print(f"  {desc} ...", flush=True)
    r = subprocess.run(cmd, cwd=DIR)
    if r.returncode != 0:
        sys.exit(f"FAILED (exit {r.returncode}): {' '.join(str(c) for c in cmd)}")

def need(path):
    return not os.path.isfile(os.path.join(DIR, path))

if need(f"{pre}ex2poly_index"):
    run_bootstrap([f"{DIR}/e2p", str(MD)], f"e2p {MD}")
if need(f"{pre}mot_deltas"):
    run_bootstrap([f"{DIR}/motTab", str(MD)], f"motTab {MD}")
if need(f"{pre}mot_res"):
    run_bootstrap([f"{DIR}/mr_ex", str(MD), str(LEN + 1)], f"mr_ex {MD} {LEN+1}")
if need(f"{pre}mot_gens{LEN}"):
    run_bootstrap([f"{DIR}/mr_mot", str(MD), str(LEN)], f"mr_mot {MD} {LEN}")


# ── helpers ─────────────────────────────────────────────────────────────────

TABLE_LINE = re.compile(r'^\{(\d+)-(\d+)\}\s*->\s*(.+)$')

def yoneda_table(bs, bb):
    """Run yoneda in table mode for beta={bs-bb}; return list of (s,a,result).
    Returns None only if yoneda printed nothing (truly invalid beta).
    A nonzero exit that still produced output is treated as a partial result
    (degree-boundary entries are marked '?' by the binary and skipped below).
    """
    r = subprocess.run(
        [f"{DIR}/yoneda", str(MD), str(LEN), str(bs), str(bb)],
        capture_output=True, text=True, cwd=DIR
    )
    entries = []
    for line in r.stdout.splitlines():
        m = TABLE_LINE.match(line.strip())
        if m:
            entries.append((int(m.group(1)), int(m.group(2)), m.group(3).strip()))
    if not entries:
        return None          # no output at all: invalid beta
    return entries


# ── discover generator ranks from one probe run ─────────────────────────────

print("Probing generator ranks ...", flush=True)
probe = yoneda_table(1, 0)
if probe is None:
    sys.exit("Probe run failed — check that bootstrap completed successfully.")

rank = defaultdict(int)          # rank[s] = number of generators at filtration s
for s, a, _ in probe:
    rank[s] = max(rank[s], a + 1)

# filtration 0 always has exactly 1 generator (the unit class)
rank[0] = max(rank[0], 1)

print(f"Generator counts per filtration:")
for s in sorted(rank):
    print(f"  s={s}: {rank[s]} generator(s)")


# ── enumerate all (bs,bb) pairs and collect nontrivial products ──────────────

products = {}          # {((as,aa),(bs,bb)): result_string}
total_betas = sum(rank[bs] for bs in range(1, LEN + 1) if rank[bs] > 0)
done = 0

print(f"\nComputing products for {total_betas} beta classes ...", flush=True)

for bs in range(1, LEN + 1):
    n = rank[bs]
    if n == 0:
        continue
    for bb in range(n):
        entries = yoneda_table(bs, bb)
        if entries is None:
            continue
        for (s, a, result) in entries:
            if result != "0" and result != "?":
                key = ((s, a), (bs, bb))
                products[key] = result
        done += 1
        print(f"  [{done}/{total_betas}] beta={{  {bs}-{bb}  }}: "
              f"{sum(1 for e in entries if e[2]!='0')} nontrivial product(s)",
              flush=True)


# ── output ──────────────────────────────────────────────────────────────────

# Enforce commutativity: if alpha*beta is known but beta*alpha is missing,
# fill in beta*alpha = alpha*beta (the motivic Ext ring is commutative over F_2).
sym_additions = {}
for (alpha, beta), result in list(products.items()):
    swapped = (beta, alpha)
    if swapped not in products:
        sym_additions[swapped] = result
products.update(sym_additions)
if sym_additions:
    print(f"Added {len(sym_additions)} product(s) by commutativity.")

print(f"\nTotal nontrivial products: {len(products)}")

outfile = os.path.join(DIR, f"{MD}_yoneda_products.txt")
with open(outfile, "w") as f:
    f.write(f"# Yoneda products  MD={MD}  LEN={LEN}\n")
    f.write(f"# format:  alpha * beta  =  result\n\n")
    for (alpha, beta), result in sorted(products.items()):
        f.write(f"{{{alpha[0]}-{alpha[1]}}} * {{{beta[0]}-{beta[1]}}}  =  {result}\n")

print(f"Results written to {outfile}")

# also available as a Python dict for import:
#   from yoneda_table import products
