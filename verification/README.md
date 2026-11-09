# Supplementary material: Exact Survival Laws and Orbit Bounds for Line Units in Binary Group Algebras

These are finite numerical checks. They are not proofs; all theorems are proved in the manuscript.

## check_line_gates.py (standard library only)
Exhaustive enumeration with exact rational arithmetic, fields GF(4), GF(8), GF(16) given by x^2+x+1, x^3+x+1, x^4+x+1.

- Section 3 (q=4,8, rho=2): KD=K for every triple of distinct lines; no nonzero coefficient of support 1, 2 or 3 is fixed by D (all triples for q=4, one triple for q=8; supports through 0 by translation invariance).
- Section 4 (q=4,8): one-unit law (eq. 6) for every line; two-unit law (Theorem 4.2) for every pair of distinct lines; Theorem 4.3 for unit sequences with j=1,2,3 active lines, including repeated lines. Proposition 4.4: for every triple of lines (q=4), two general triples (q=8) and a four-line and a repeated-line family, the survival probability equals 1/|P| with |P| the number of reachable values of K*prod A_i, |P| is a power of 2, the probability is at least the product of the one-unit probabilities, and removing repeated lines does not change it.
- Section 5 (q=4 with rho=2,3; q=8 with rho=2): for every normalized nonparallel direction, every E with |E|<=q/2 and every nonzero relation: |F|>=R_M for the greedy exposed set, the exact formula (eq. 11), the bound (1-p)^{|F|}, the rank-deficiency bound (eq. 13), and the orbit bound (eq. 14) against the exact orbit-constancy probability.
- Section 6 (q=8,16, all gamma): H_Lambda = B_Lambda + diag(1_Lambda); rank q for odd |Lambda|; rank >= q/2 when |Lambda| even with nonzero sum; rank inequality (eq. 19). For q=8 and rho=3,4,5: the exact parity/sum formula (eq. 20) and both sides of (eq. 17) against the exact probability.
- Section 7 (q=4, rho=2 all lines; rho=3 subset of lines): eqs. (21)-(23) and the same-line equality of events (checked pointwise for every role-bit assignment, in addition to equality of probabilities).

Output: check_line_gates.log.

## f16-count/
C++17 program for the exact count of Section 8 and its complete output (q4.jsonl, q16.jsonl).
Build and run: `g++ -std=c++17 -O3 exposed_column_rank_counter.cpp -o rc && ./rc --q 4 > q4.jsonl && ./rc --q 16 > q16.jsonl`.
Rerun on 2026-10-26 (GCC, Linux x86_64, about 1 s): output identical to the files here apart from the timing field; summary line reports 672/672 cases validated.
