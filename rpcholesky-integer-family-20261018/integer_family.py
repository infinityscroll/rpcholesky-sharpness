"""Explicit integer SPD examples for RPCholesky sharpness. Python standard library.

This is a rescaling of SINGLE-PARAMETER-PROOF.md, not a floating-point
implementation. See README.md for attribution, claim limits, and exact checks.
"""
from fractions import Fraction
from math import isqrt


def parameter(n, eta=Fraction(1, 10)):
    """Smallest integer q >= 2 with 4*n*n/q**2 <= eta, exactly."""
    if not isinstance(n, int) or n < 2:
        raise ValueError("n must be an integer at least 2")
    eta = Fraction(eta)
    if eta <= 0:
        raise ValueError("eta must be positive")
    a, b = 4 * n * n * eta.denominator, eta.numerator
    q = isqrt(a // b)
    if q * q * b < a:
        q += 1
    return max(2, q)


def make_matrix(n, q):
    """Return the integer matrix with smallest eigenvalue exactly one."""
    if not isinstance(n, int) or n < 2 or not isinstance(q, int) or q < 2:
        raise ValueError("n and q must be integers at least 2")
    exponent = 4 * n * (n - 1)
    matrix = [[int(i == j) for j in range(n)] for i in range(n)]
    for m in range(1, n):
        z = [q ** (m - i) for i in range(1, m + 1)]
        z += [-sum(q ** (2 * m - 1 - 2 * ell) for ell in range(m))]
        z += [0] * (n - m - 1)
        scale = q ** (exponent - 4 * n * (m - 1) - 4 * m + 2)
        for i in range(n):
            for j in range(n):
                matrix[i][j] += scale * z[i] * z[j]
    return matrix


if __name__ == "__main__":
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int, help="matrix dimension; k=n-1 pivots")
    parser.add_argument("--eta", default="1/10", help="relative tolerance as a rational")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    eta = Fraction(args.eta)
    q = parameter(args.n, eta)
    matrix = make_matrix(args.n, q)
    args.output.write_text(json.dumps({
        "n": args.n, "k": args.n - 1, "q": q, "eta": str(eta),
        "encoding": "Matrix entries are hexadecimal integer strings; parse with int(value, 16).",
        "matrix": [[hex(x) for x in row] for row in matrix],
        "tau_k": "1",
        "guarantee": "(1-eta)*2**k <= exact expected RP residual <= 2**k",
        "status": "AI-assisted proof; priority and human review unconfirmed"
    }, indent=2) + "\n")
    print(f"Wrote {args.n}x{args.n} integer example with q={q} to {args.output}")
