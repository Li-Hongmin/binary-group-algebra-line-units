# Exact Survival Laws and Orbit Bounds for Line Units in Binary Group Algebras

[![verify](https://github.com/Li-Hongmin/binary-group-algebra-line-units/actions/workflows/verify.yml/badge.svg)](https://github.com/Li-Hongmin/binary-group-algebra-line-units/actions/workflows/verify.yml) [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24054769-blue)](https://doi.org/10.5281/zenodo.24054769)

**Hongmin Li** (Institute of Science Tokyo; The University of Tokyo) · ORCID [0000-0003-0228-0600](https://orcid.org/0000-0003-0228-0600)

Main result of the manuscript: exact survival probabilities and orbit bounds for random products of line units in binary group algebras F₂[F_q^ρ], with an exact count over F₁₆.

- Manuscript on rStation: <https://rstation.ai/problems/math/papers/3/v2> (version 2)
- Archived on Zenodo: concept DOI [10.5281/zenodo.24054769](https://doi.org/10.5281/zenodo.24054769) (all versions); this repository mirrors version 2, DOI [10.5281/zenodo.24054776](https://doi.org/10.5281/zenodo.24054776)

## Contents

- `paper/paper.md`, `paper/binary-group-algebra-line-units.pdf` — the manuscript (version 2).
- `verification/` — the supplement of that version, unpacked unchanged; its own README describes each file.
- `.github/workflows/verify.yml` — reruns every check below on each push.

## What the programs check

| Check | Command (run inside `verification/`) |
|---|---|
| Sections 3-7: exhaustive exact checks over GF(4), GF(8), GF(16) | `python check_line_gates.py > out_line_gates.log` |
| Section 8: exact count (C++17) | `cd f16-count && g++ -std=c++17 -O3 exposed_column_rank_counter.cpp -o rc && ./rc --q 4 > out_q4.jsonl && ./rc --q 16 > out_q16.jsonl` |

The programs check the finite and exact computations used in the manuscript; the mathematical arguments are those written in the paper.

## Reproducing

```sh
cd verification
python check_line_gates.py > out_line_gates.log
cd f16-count && g++ -std=c++17 -O3 exposed_column_rank_counter.cpp -o rc && ./rc --q 4 > out_q4.jsonl && ./rc --q 16 > out_q16.jsonl
```

Requirements: Python 3.9+; the C++ step needs a C++17 compiler. Run the commands from a copy if you want to keep the recorded outputs untouched.

## Cite

```bibtex
@misc{li2026binarygroupalgebralineun,
  author    = {Li, Hongmin},
  title     = {Exact Survival Laws and Orbit Bounds for Line Units in Binary Group Algebras},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24054769},
  url       = {https://doi.org/10.5281/zenodo.24054769}
}
```

## License

Manuscript and data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: MIT (see `LICENSE`).
