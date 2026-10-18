# Sharpness of RPCholesky error doubling

**Public research draft — independent human review and novelty confirmation are pending.**

This repository proposes an explicit family of integer positive-definite matrices showing that the known worst-case expected trace-error factor for exactly **k ordinary RPCholesky pivots is sharp at 2^k**. The family has dimension k+1 and polynomial encoding length. The contribution is the proposed matching lower construction; the upper bound is due to Chen, Epperly, Tropp and Webber, [CPAM 2025, Lemma 5.5](https://doi.org/10.1002/cpa.22234).

**[Read the seven-page proof](sharpness-explicit-note.pdf)** · [LaTeX source](sharpness-explicit-note.tex) · [Exact verifier](rpcholesky-integer-family-20261018/verify_integer_family.py) · [Integer generator](rpcholesky-integer-family-20261018/integer_family.py)

## Claim and scope

Let R_k(A) be the residual after k pivots, each selected with probability proportional to the current residual diagonal. Let tau_k(A) be the sum of eigenvalues beyond the largest k. The proposed theorem is

```text
sup over d > k and A positive definite of E[trace(R_k(A))] / tau_k(A) = 2^k.
```

The supremum can already be approached by integer matrices of size n=k+1. For 0<eta<1, the generator chooses an integer q>=max(2,2n/sqrt(eta)), produces a matrix with tau_k=1, and gives the proposed finite guarantee

```text
(1-eta) * 2^k <= E[trace(R_k)] <= 2^k.
```

This is an expectation over the prescribed pivot law, not a claim about every path. Finite examples approach the supremum. **One additional pivot exhausts each matrix**, so this does not refute an oversampling conjecture. The matrices are deliberately ill-conditioned; no practical model or kernel-performance failure is claimed.

The mechanism tracks the omitted coordinate through a hierarchy of pivot directions. Its limiting survivor law produces a nearly doubling inverse-weight payoff. The actual residual trace still decreases when more pivots are taken. Uniform bounds control finite-scale and positive-definite perturbation errors.

The note also discusses [Epperly's thesis, Conjecture 11.2](https://doi.org/10.7907/pef3-mg80). Its displayed notation is ambiguous. The implication here concerns the expected-residual interpretation suggested by the adjacent citation; the author's intended interpretation has not been confirmed. The sharpness claim above is stated independently of that ambiguity.

## Reproduce the finite exact checks

Clone or download this repository, then run from its root:

```sh
python3 rpcholesky-integer-family-20261018/verify_integer_family.py
```

The verifier uses only Python's standard library. It requires a Unix-like system because it uses `resource`; run without `python -O` because checks use assertions. It performs ten fixed matrix-identity/eigenbasis checks and four complete exact Schur-recursion expectations, checks a known 23/28 control, and rejects a uniform-pivot mutation. It uses no random sampling, network access or GPU.

The command rewrites the four witness files and `EXACT-CONTROLS.json`; use a fresh clone if you want to preserve saved runtime metadata. A 60 CPU-second limit is enforced. Peak memory is checked against 512 MiB after execution; this is not an OS memory limit. Exact expectations enumerate branches and are intended for small dimensions. Finite checks support the implementation; the arbitrary-k claim rests on the symbolic proof.

At eta=1/10, the saved exact results are:

| k | Dimension | q | Expected residual (rounded) | Optimal error | Known upper bound | Largest entry, bits |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 13 | 1.9999999976 | 1 | 2 | 30 |
| 2 | 3 | 19 | 3.9944903578 | 1 | 4 | 102 |
| 3 | 4 | 26 | 7.9881918819 | 1 | 8 | 226 |
| 4 | 5 | 32 | 15.9765967765 | 1 | 16 | 401 |

The witness JSON files contain exact fractions, not just these rounded displays. Matrix entries are decimal strings; expectation numerators and denominators are hexadecimal strings.

Generate a separate example with:

```sh
python3 rpcholesky-integer-family-20261018/integer_family.py 5 --eta 1/10 --output example.json
```

The generator's CLI uses hexadecimal matrix strings, parsed with `int(value, 16)`. Generating a matrix does not compute its expected residual.

## Verification and attribution

The construction, proof and code underwent separately tasked AI reviews, and the frozen finite checks were replayed in a fresh directory with matching witness hashes. Those checks are not formal verification or independent human validation. Bounded primary-source searches found related simplex, scale-separation and inverse-diagonal ingredients but did not establish novelty or priority. The note credits those ingredients and the existing upper bound.

Jitendra Prajapati directed the project and set its research objective. OpenAI Codex produced the mathematical construction, derivation, implementation and source audit; separate AI agents reviewed them. This repository does not claim external endorsement, peer review, adoption, or hiring outcomes.

The full quantitative argument also appears in [the single-parameter proof](rpcholesky-single-parameter-family-20261018/SINGLE-PARAMETER-PROOF.md). `SCIENTIFIC-FILES.json` records file hashes for this release. Scientific code and witnesses match the reviewed internal package. Publication-status and attribution wording in the note were updated, and the supplementary proof's introductory link now points to the included note; the mathematical arguments are unchanged. No third-party papers or private correspondence are redistributed.

Corrections and precise prior-work pointers are welcome through repository issues.
