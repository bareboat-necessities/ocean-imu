#!/usr/bin/env python3
"""Worst-case stack chains and reachability from GCC -fcallgraph-info=su.

Usage: callgraph.py <variant> <file.ci> <filter class>   (writes CSV rows)

For the sample-path entry (carma3_probe_step) and for the filter's
time_update / measurement_update_acc_only / measurement_update_mag_only,
reports the deepest static call chain (sum of the frames GCC reports; library
functions outside the translation unit -- libm, the soft-float helpers --
count as 0) and whether heap allocation, soft-float double arithmetic or a C++
throw is reachable through direct calls.  Static reachability: a reachable
call may sit on a branch that never executes.
"""
import csv
import re
import sys

sys.setrecursionlimit(100000)
ALLOC = r"^(malloc|calloc|realloc|_Znwj|_Znaj|_ZnwjSt11align_val_t)$"
DOUBLE = r"^__(add|sub|mul|div)df3$"
THROW = r"^__cxa_throw$"


def load(path):
    text = open(path).read()
    frames, labels, edges = {}, {}, {}
    for m in re.finditer(r'node: \{ title: "([^"]+)" label: "([^"]*)"', text):
        b = re.search(r"\\n(\d+) bytes \(", m.group(2))
        frames[m.group(1)] = int(b.group(1)) if b else 0
        labels[m.group(1)] = m.group(2).split("\\n")[0]
    for m in re.finditer(r'edge: \{ sourcename: "([^"]+)" targetname: "([^"]+)"', text):
        edges.setdefault(m.group(1), []).append(m.group(2))
    return frames, labels, edges


def deepest(t, frames, edges, memo, active):
    if t in memo:
        return memo[t]
    if t in active:   # recursion: count the cycle once
        return 0, [t]
    active.add(t)
    best = (0, [])
    for c in edges.get(t, []):
        d = deepest(c, frames, edges, memo, active)
        if d[0] > best[0]:
            best = d
    active.discard(t)
    memo[t] = (frames.get(t, 0) + best[0], [t] + best[1])
    return memo[t]


def reaches(t, pat, edges, seen):
    if t in seen:
        return False
    seen.add(t)
    if re.search(pat, t.split(":")[-1]):
        return True
    return any(reaches(c, pat, edges, seen) for c in edges.get(t, []))


def main():
    variant, path, cls = sys.argv[1:4]
    frames, labels, edges = load(path)
    roots = [("sample_path", r"(^|:)carma3_probe_step$")]
    for name, n in (("time_update", 11), ("measurement_update_acc_only", 27), ("measurement_update_mag_only", 27)):
        roots.append((name, rf"{cls}I.*{n}{name}E"))
    wr = csv.writer(sys.stdout, lineterminator="\n")
    memo = {}
    for name, pat in roots:
        for t in [t for t in frames if re.search(pat, t)]:
            d, chain = deepest(t, frames, edges, memo, set())
            top = max(chain, key=lambda x: frames.get(x, 0))
            wr.writerow([variant, name, d, len(chain), reaches(t, ALLOC, edges, set()),
                         reaches(t, DOUBLE, edges, set()), reaches(t, THROW, edges, set()),
                         labels.get(top, top)[:90]])


if __name__ == "__main__":
    main()
