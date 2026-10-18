"""Frozen exact controls for the explicit integer family; no statistical tests."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import resource
import sys
import time

from integer_family import make_matrix, parameter

ROOT = Path(__file__).resolve().parent
resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
START = time.process_time()


def dot(x, y):
    return sum((a * b for a, b in zip(x, y)), F(0))


def rational_reference(n, q):
    t = F(1, q)
    w = [t ** i for i in range(n)]
    directions = [w]
    h = t ** (4 * n)
    delta = t ** (4 * n * (n - 1))
    a = [[delta * int(i == j) for j in range(n)] for i in range(n)]
    eigenvalues = [F(1)]
    for m in range(1, n):
        v = [w[i] * w[m] if i < m else
             (-dot(w[:m], w[:m]) if i == m else F(0)) for i in range(n)]
        directions.append(v)
        eigenvalues.append(1 + h ** (m - 1) * dot(v, v) / delta)
        for i in range(n):
            for j in range(n):
                a[i][j] += h ** (m - 1) * v[i] * v[j]
    return [[x / delta for x in row] for row in a], directions, eigenvalues


def check_eigenbasis(a, directions, eigenvalues):
    assert eigenvalues[0] == 1 and all(x > 1 for x in eigenvalues[1:])
    for index, v in enumerate(directions):
        assert dot(v, v) > 0
        assert all(dot(v, u) == 0 for u in directions[:index])
        assert [dot(row, v) for row in a] == [eigenvalues[index] * x for x in v]


@lru_cache(None)
def schur_mean(state, steps, uniform=False):
    if not steps:
        return sum(state[i][i] for i in range(len(state)))
    trace = sum(state[i][i] for i in range(len(state)))
    positive = [i for i in range(len(state)) if state[i][i] > 0]
    if not positive:
        return F(0)
    probabilities = [F(1, len(positive)) if uniform else state[i][i] / trace for i in positive]
    assert sum(probabilities) == 1
    total = F(0)
    for i, probability in zip(positive, probabilities):
        keep = [j for j in range(len(state)) if j != i]
        residual = tuple(tuple(state[j][ell] - state[j][i] * state[i][ell] / state[i][i]
                               for ell in keep) for j in keep)
        total += probability * schur_mean(residual, steps - 1, uniform)
    return total


def state(a):
    return tuple(tuple(F(x) for x in row) for row in a)


def hex_fraction(x):
    return {"numerator_hex": hex(x.numerator), "denominator_hex": hex(x.denominator)}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


control = state([[1, F(1, 2), 0], [F(1, 2), 1, 0], [0, 0, 1]])
assert schur_mean(control, 2) == F(23, 28)
assert schur_mean(control, 2, True) != F(23, 28)
schur_mean.cache_clear()

identity_checks = []
for n in [2, 3, 4, 5, 6]:
    for q in [2, 3]:
        a = make_matrix(n, q)
        reference, directions, eigenvalues = rational_reference(n, q)
        assert a == reference
        check_eigenbasis(a, directions, eigenvalues)
        assert max(abs(x) for row in a for x in row) <= 3 * q ** (4 * n * (n - 1))
        identity_checks.append({"n": n, "q": q, "status": "PASS"})

rows = []
eta = F(1, 10)
for n in [2, 3, 4, 5]:
    k = n - 1
    q = parameter(n, eta)
    assert 4 * n * n <= eta * q * q
    assert q == 2 or eta * (q - 1) ** 2 < 4 * n * n
    a = make_matrix(n, q)
    reference, directions, eigenvalues = rational_reference(n, q)
    assert a == reference
    check_eigenbasis(a, directions, eigenvalues)
    mean = schur_mean(state(a), k)
    assert (1 - eta) * 2 ** k <= mean <= 2 ** k
    dpp = F(n) / sum((1 / x for x in eigenvalues), F(0))
    assert F(n) / (1 + F(n - 1, q ** (4 * n))) <= dpp <= n
    path = ROOT / f"integer-witness-k{k}.json"
    path.write_text(json.dumps({
        "n": n, "k": k, "q": q, "eta": str(eta),
        "matrix_decimal_strings": [[str(x) for x in row] for row in a],
        "tau_k": "1", "rp_expectation": hex_fraction(mean),
        "dpp_expectation": hex_fraction(dpp),
        "eigenvalues": [str(x) for x in eigenvalues],
        "method": "Exact full Schur-recursion expectation, no Monte Carlo"
    }, indent=2) + "\n")
    rows.append({"n": n, "k": k, "q": q, "status": "PASS",
                 "rp_expectation_display_rounded": float(mean),
                 "rp_over_2_to_k_display_rounded": float(mean / 2 ** k),
                 "maximum_entry_bits": max(abs(x).bit_length() for row in a for x in row),
                 "witness_sha256": sha(path)})
    schur_mean.cache_clear()

peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
peak_bytes = peak if sys.platform == "darwin" else peak * 1024
assert peak_bytes <= 512 * 2 ** 20
report = {"status": "PASS_FROZEN_EXACT_INTEGER_CONTROLS", "identity_cases": identity_checks,
          "exact_error_cases": rows, "known_control": "23/28", "uniform_mutation": "REJECTED",
          "cpu_seconds": time.process_time() - START, "peak_resident_bytes": peak_bytes,
          "memory_note": "Post-run observed peak check, not OS memory containment",
          "plan_sha256": sha(ROOT / "PLAN.json"),
          "generator_sha256": sha(ROOT / "integer_family.py"),
          "verifier_sha256": sha(Path(__file__)),
          "scope": "Finite exact checks support implementation; arbitrary-n claim rests on symbolic proof; novelty and human review unconfirmed."}
(ROOT / "EXACT-CONTROLS.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
