# An explicit one-parameter rational RPCholesky family

**Status:** symbolic strengthening of the local hierarchical-family proof. This document does not modify the frozen proof, note, or PDF. It supplies explicit integer exponents and uniform error bounds; it makes no novelty claim. No numerical experiment, external computation, paid run, outreach, or publication was used.

The family and tournament mechanism are included in the [self-contained research note](../sharpness-explicit-note.tex). The arguments below give the single-parameter bounds in detail. All expectations refer to ordinary diagonal-proportional RPCholesky or to fixed-cardinality determinant sampling, respectively.

## 1. Explicit statement

Fix an integer n >= 2 and put k = n - 1. For 0 < t <= 1/2, define

    w_i = t^(i-1),                 S_m = sum_{i=1}^m w_i^2,
    Z = S_n,                      p_i = w_i^2 / Z,
    h = t^(4n),                   delta = t^(4n(n-1)).

For m = 1,...,n-1 let v_m have entries

    (v_m)_i = w_i w_(m+1)  (i <= m),
    (v_m)_(m+1) = -S_m,
    (v_m)_i = 0           (i > m+1),

and set

    A0 = sum_{m=1}^{n-1} h^(m-1) v_m v_m^T,
    A  = A0 + delta I_n.

Thus explicit safe exponents are

    M(n) = 4n,             N(n) = 4n(n-1).

For every such t,

    | E_RP tr(R_k(A))/delta - 2^(n-1) |
      <= n 2^(n-1) t^2
         + [2^n n(n-1) + 2(1+(n-1)^2)] t^4
      <= 2^(n+1) n^2 t^2.                              (1)

Moreover,

    n / [1+(n-1)t^(4n)]
      <= E_DPP tr(R_k(A))/delta <= n,                  (2)

so in particular

    lim_{t -> 0} E_RP tr(R_k(A))/delta = 2^(n-1),
    lim_{t -> 0} E_DPP tr(R_k(A))/delta = n,
    lim_{t -> 0} E_RP tr(R_k(A)) / E_DPP tr(R_k(A))
      = 2^(n-1)/n.                                     (3)

The matrix is rational whenever t is rational. It is positive definite and has smallest eigenvalue exactly delta, hence tau_k(A) = delta.

For an explicit relative tolerance eta > 0, take any integer

    q >= max{2, 2n / sqrt(eta)},       t = 1/q.          (4)

Then (1) gives

    | E_RP tr(R_k(A))/delta - 2^(n-1) |
      <= eta 2^(n-1).

For an absolute tolerance epsilon > 0, it instead suffices that

    q >= max{2, sqrt(2^(n+1) n^2 / epsilon)}.

Every selected finite q specifies a finite rational SPD witness; no limiting matrix is used by the algorithm.

## 2. Spectral and minor identities

The v_m are mutually orthogonal, orthogonal to w, and satisfy

    ||v_m||^2 = S_m S_(m+1).

Consequently the nonzero eigenvalues of A0 are

    lambda_m = h^(m-1) S_m S_(m+1),   1 <= m <= n-1.

Since h <= 1 and S_m >= 1,

    lambda_* := min_m lambda_m >= h^(n-2).             (5)

Write u = w / sqrt(Z). Then

    A0 >= lambda_* (I - u u^T)

in the positive-semidefinite order. For every proper coordinate set T, compression gives

    A0[T,T] >= lambda_* (I_T - u_T u_T^T)
             >= lambda_* p_min I_T,

where p_min = p_n. Indeed the smallest eigenvalue of I_T - u_T u_T^T is 1-||u_T||^2 = sum_{j notin T} p_j >= p_min. We may therefore use the uniform proper-principal-block bound

    A0[T,T] >= mu I_T,       mu := h^(n-2) p_min.      (6)

Also t <= 1/2 implies

    1 <= Z <= 4/3,
    1/p_min <= 2 t^(-2(n-1)),
    1/p_min^2 <= 2 t^(-4(n-1)).                       (7)

The slack constants 2 in (7) are intentional.

For 1 <= r <= n-1, let H_r be the first r rows v_1,...,v_r restricted to coordinates 1,...,r+1. The rows are orthogonal and its kernel is spanned by (w_1,...,w_(r+1)). Thus its signed maximal-minor vector is proportional to that vector. Cauchy--Binet and the row norm identity determine the proportionality constant exactly:

    det(H_r with column j omitted)^2
       = w_j^2 (product_{ell=1}^r S_ell)^2.            (8)

To verify the constant, the sum of these squared minors is

    det(H_r H_r^T) = product_{ell=1}^r S_ell S_(ell+1).

Dividing this by sum_{j=1}^{r+1} w_j^2 = S_(r+1) yields (product_{ell=1}^r S_ell)^2 because S_1 = 1.

Finally put B = sum_{m=1}^{n-1} v_m v_m^T, so its nonzero eigenvalues are S_m S_(m+1). Since S_m <= 4/3,

    ||B|| <= 16/9 < 2.

For every size-r coordinate set T, another application of Cauchy--Binet gives

    sum_{I subset {1,...,n-1}, |I|=r} det(V[I,T])^2
      = det B[T,T] <= 2^r,                            (9)

where V has rows v_m. This bounds the sum of all minor coefficients, not merely an individual coefficient.

## 3. Uniform hierarchy error, including illegal paths

For any selected coordinate set S and unused i, the exact RPCholesky transition on A0 is

    P_0(i | S)
      = det A0[S union {i}, S union {i}]
        / sum_{j notin S} det A0[S union {j}, S union {j}].  (10)

The common denominator det A0[S,S] from the Schur-diagonal formula cancels. Every determinant in (10) is positive while |S| <= n-2.

Call a size-m selected set legal if

    S = {1,...,m+1} minus {J},      J <= m+1.

For m = 0 use S empty and J = 1. Let r = m+1 and e = r(r-1)/2. In the Cauchy--Binet expansion of each size-r determinant in (10), the unique row subset with minimum h-exponent e is {1,...,r}; every other row subset has exponent at least e+1.

For a legal S, only i = J and i = m+2 have a nonzero coefficient at exponent e. Formula (8) identifies those coefficients as

    c_J     = w_(m+2)^2 (product_{ell=1}^r S_ell)^2,
    c_(m+2) = w_J^2     (product_{ell=1}^r S_ell)^2.

For every unused i write its determinant weight as

    h^e c_i + E_i,       c_i = 0 for i > m+2.

All E_i are nonnegative. By (9), h <= 1, and the integer gap between exponents,

    0 <= E_i <= 2^r h^(e+1).

The total leading weight L and total remainder E therefore satisfy

    L = h^e (c_J+c_(m+2)) >= h^e t^(2(n-1)),
    E = sum_{i notin S} E_i <= n 2^(n-1) h^(e+1).

The first inequality uses J <= m+1 <= n-1 and product S_ell >= 1; it is deliberately weaker than the available t^(2(n-2)) lower bound.

Let Q be the ideal two-choice transition formed by normalizing the leading coefficients. The actual transition is exactly the mixture

    P_0(. | S) = [L/(L+E)] Q(. | S) + [E/(L+E)] U,

where U is the normalized remainder distribution if E > 0. Hence the total variation distance from Q, at every legal state, is at most

    alpha := n 2^(n-1) h t^(-2(n-1)).                 (11)

Here total variation means the maximal probability difference over events, equivalently the smallest disagreement probability in a coupling. Couple the full RP path to the ideal path step by step while they agree. The ideal path always remains legal. A union bound over k = n-1 transitions gives

    Pr{the two paths ever disagree} <= (n-1) alpha.    (12)

In particular (12) bounds the probability of every illegal path collectively. It does not require any control of transitions after the first disagreement.

Let P_0 also denote the A0 path law, let J be the final unselected coordinate, and define f(J) = 1/p_J. Since 0 < f <= 1/p_min, (7), (11), and (12) imply the required weighted bound

    | E_{P_0} f(J) - E_Q f(J) |
      <= 2^n n(n-1) h t^(-4(n-1)).                   (13)

This is uniform in t and explicitly includes the factor 1/p_min. In particular, convergence of transition probabilities alone is not being used to justify convergence of the diverging endpoint payoff.

Under Q, when the current representative is i and coordinate j is introduced, i survives with probability p_i/(p_i+p_j), and j survives otherwise. Thus Q is precisely the tournament in the starting proof. Denote its endpoint expectation by F_n(t).

## 4. Uniform small-noise error

Let P_delta be the RP path law on A = A0 + delta I. For every proper set T of size r, (6) yields

    1 <= det(A0[T,T] + delta I_T) / det A0[T,T]
      <= (1+x)^r,                x := delta/mu.

At any selected state S of size at most n-2, the determinant weights defining its transition are therefore tilted by factors between 1 and

    beta := (1+x)^(n-1).

For each unused i this implies

    P_delta(i | S) >= P_0(i | S)/beta.

The two transition laws consequently have total variation distance at most

    1 - beta^(-1) <= (n-1)x.

The last inequality follows from 1-(1+x)^(-k) <= kx for x >= 0, and needs no smallness assumption on x. Stepwise coupling over n-1 transitions gives

    TV(P_delta, P_0) <= (n-1)^2 delta / mu.            (14)

After n-1 pivots, the actual residual trace depends only on the omitted coordinate J and is 1/(A^(-1))_(J,J). By spectral decomposition,

    (A^(-1))_(j,j) = p_j/delta + b_j,
    0 <= b_j <= 1/lambda_* <= h^(-(n-2)).

Its normalized value is therefore

    g_delta(j) := tr(R_k(A) | J=j)/delta
                = 1/(p_j + delta b_j),

and uniformly in j,

    0 <= f(j)-g_delta(j)
      <= delta h^(-(n-2)) / p_min^2.                  (15)

Combining (14), (15), and mu = h^(n-2) p_min gives

    | E_{P_delta} g_delta(J) - E_{P_0} f(J) |
      <= [1+(n-1)^2] delta h^(-(n-2)) / p_min^2
      <= 2[1+(n-1)^2]
           delta h^(-(n-2)) t^(-4(n-1)).              (16)

The comparison random variable under P_0 is f(J), the limiting amplification weight. It is not the zero residual trace of the rank-(n-1) matrix A0 after n-1 pivots.

## 5. Substituting the exponents

With h = t^(4n) and delta = t^(4n(n-1)), both weighted error scales are exactly t^4:

    h t^(-4(n-1)) = t^4,

    delta h^(-(n-2)) t^(-4(n-1)) = t^4.               (17)

For the tournament, conditional on representative i before introducing j,

    E[1/p_(new representative) | i] = 2/(p_i+p_j).

Since p_j/p_i <= t^2 for i < j, its expectation satisfies

    Z [2/(1+t^2)]^(n-1) <= F_n(t) <= Z 2^(n-1).

Using Z >= 1, Z <= 1/(1-t^2), and t <= 1/2,

    2^(n-1)[1-(n-1)t^2] <= F_n(t)
      <= 2^(n-1)[1+(4/3)t^2].

It follows that

    |F_n(t)-2^(n-1)| <= n 2^(n-1) t^2.               (18)

The triangle inequality applied to (13), (16), (17), and (18) proves the first inequality in (1). For the simpler final constant, use t^4 <= t^2 and

    n 2^(n-1) + 2^n n(n-1) <= 2^n n^2,
    2[1+(n-1)^2] <= 2n^2 <= 2^n n^2.

These prove the second inequality in (1), and then the rational tolerance choices follow immediately.

## 6. Determinant-sampling comparator

For k = n-1 the exact fixed-size DPP residual expectation is

    E_DPP tr(R_k(A)) = n det(A) / e_(n-1)(A)
                    = n / tr(A^(-1)).

Indeed every omitted-coordinate contribution has determinant weight det(A[-j,-j]) and residual det(A)/det(A[-j,-j]), so all n numerator terms equal det(A). The eigenvalues above then give

    E_DPP tr(R_k(A))/delta
       = n / [1 + delta sum_{m=1}^{n-1} 1/(lambda_m+delta)].

By (5) and the chosen exponents,

    0 <= delta sum_m 1/(lambda_m+delta)
      <= (n-1) delta h^(-(n-2))
       = (n-1)t^(4n).

This proves (2), and also the explicit additive bound

    0 <= n - E_DPP tr(R_k(A))/delta
      <= n(n-1)t^(4n).

Together with (1), it proves (3). If 0 < eta < 1 and q obeys (4), the finite rational matrix also obeys the concrete ratio bound

    E_RP tr(R_k(A)) / E_DPP tr(R_k(A))
       >= (1-eta) 2^(n-1)/n.

This remains a same-pivot-count comparison in dimension n=k+1. One additional pivot exhausts every matrix in the family. No claim about an oversampled regime follows from these formulas.
