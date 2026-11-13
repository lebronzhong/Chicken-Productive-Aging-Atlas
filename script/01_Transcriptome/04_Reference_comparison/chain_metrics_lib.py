#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Chain metrics lib.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import collections
from pipeline_7b import strict_sub


class Index:
    def __init__(self, keys):
        self.keys = set(keys)
        self.j = collections.defaultdict(set)
        self.start = collections.defaultdict(set)
        self.left = set()
        self.right = set()
        for k in self.keys:
            c, s, e = k
            self.start[c, s, e[0]].add(k)
            for a, b in e:
                self.j[c, s, (a, b)].add(k)
                self.left.add((c, s, a))
                self.right.add((c, s, b))

    def relations(self, k):
        if k in self.keys:
            return True, False, False
        c, s, e = k
        qsub = any(strict_sub(e, b[2]) for b in self.j.get((c, s, e[0]), ()))
        rsub = any(strict_sub(b[2], e) for x in e for b in self.start.get((c, s, x), ()))
        return False, qsub, rsub

    def classify(self, k):
        exact, qsub, rsub = self.relations(k)
        if exact:
            return "EXACT", rsub
        if qsub:
            return "QUERY_SUBCHAIN", rsub
        c, s, e = k
        known = all((c, s, a) in self.left and (c, s, b) in self.right for a, b in e)
        return ("NEW_COMBINATION_KNOWN_SITES" if known else "AT_LEAST_ONE_NEW_SPLICE_SITE"), rsub


def covered(refkeys, querykeys):
    idx = Index(querykeys)
    cats = collections.Counter()
    for k in refkeys:
        exact, refsub, querysub = idx.relations(k)
        cat = (
            "EXACT"
            if exact
            else (
                "BOTH_PARTIAL_DIRECTIONS"
                if refsub and querysub
                else (
                    "REFERENCE_SUBCHAIN_OF_QUERY"
                    if refsub
                    else ("QUERY_SUBCHAIN_OF_REFERENCE" if querysub else "NO_CHAIN_MATCH")
                )
            )
        )
        cats[cat] += 1
    n = len(refkeys)
    ex = cats["EXACT"]
    part = sum(
        cats[x]
        for x in (
            "BOTH_PARTIAL_DIRECTIONS",
            "REFERENCE_SUBCHAIN_OF_QUERY",
            "QUERY_SUBCHAIN_OF_REFERENCE",
        )
    )
    return {
        "reference_chains": n,
        "counts": dict(cats),
        "exact": ex,
        "partial_only": part,
        "exact_or_partial": ex + part,
        "exact_percent": 100 * ex / n if n else None,
        "exact_or_partial_percent": 100 * (ex + part) / n if n else None,
    }


def known_site(idx, k, edge, role):
    c, s, _ = k
    side = 0 if (role == "donor") == (s == "+") else 1
    return (c, s, edge[side]) in (idx.left if side == 0 else idx.right)
