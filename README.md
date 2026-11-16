# Maximal torus theory for crossed product quantum groups

TB 2.0 and GPT 6.2

For a discrete group acting by automorphisms on a compact connected Lie group, a single torus of the crossed product detects characters, amenability, and polynomial growth exactly when the action preserves a maximal torus of the compact factor. We give counterexamples to three earlier torus tests and formulate revised conjectures using the quantum center, quantum quotients, and cocycle twists.

## Files

- [008.pdf](008.pdf) is the paper. [008.tex](008.tex) is its self-contained LaTeX source.
- [tori_partitions.py](tori_partitions.py) computes the toral relations for rotations of the basis in rank two and checks ball volumes in free products of cyclic groups.
- [tori_twists.py](tori_twists.py) checks the Klein-four cocycle, its pullback to the free group, and the restrictions of character words for the twist of two free unitary quantum groups.

## Reproducing the computations

Run both scripts from the repository root with Python 3.8 or later. They use only the standard library.

```sh
python3 tori_partitions.py
python3 tori_twists.py
```

The first script finds 4 infinite dihedral, 4 Klein-four, and 632 order-two cases among the primitive signed pairs with both coordinates bounded by 16 in absolute value. The second checks 69,905 formal coefficient words and all 341 character words of length at most four, including the empty word in both counts.

## Building the paper

```sh
pdflatex 008.tex
pdflatex 008.tex
```
