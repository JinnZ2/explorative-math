# CLAUDE.md

Guidance for working in this repository. CC0. Co-authored.

## Purpose

Explorative geometry. Map where human curiosity and AI curiosity land on the
same objects, and record **all** outputs, including the ones that look
uninteresting. A run nobody found interesting is still a run, and it is
kept. The log is the product; no single result is.

## Rules

```text
license          CC0, everything
python           stdlib only: math, itertools, random, fractions for the maths;
                 json, pathlib, unittest for I/O and tests. No numpy, no pip.
                 phone-buildable: every script runs in seconds on a phone.
words            no word coinage. Use an existing name or leave the slot unnamed.
gaps             leave gaps. A slot with nothing in it stays void (null),
                 never dropped and never filled to look complete.
readings         "not interpretable atm" is a legal reading. Use it rather
                 than force one.
shape            keep shape: record the object as computed before any summary.
status tags      every claim is OBSERVED (computed here), DERIVED (follows from
                 computed values, reasoning shown) or PROPOSED (not yet run).
interest fields  human_interest is the human's to fill; the AI never writes it.
                 ai_interest is the AI's own reading, labelled as such.
co-authored      tag new files and commits co-authored.
```

## Method

The carrier's term for this method has no English name. It stays unnamed here.

```text
shape -> rotate (n-dim) -> project -> measure angles
      -> compose (phi-scale / second shape) -> compare nature catalog
```

Known related names, kept as pointers, not as a name for the method:
- grand tour
- projection pursuit
- Coxeter / Petrie projections
- cut-and-project
- principal angles

## Layout

```text
emlib.py               shared: run-file writer + schema check, LZ76 complexity,
                       small vector helpers
nature_catalog.py      entries a run can match against (all CARRIED, not verified)
experiments/
  near_orthogonality.py   a  random unit vectors: angle distribution vs dimension
  shell_concentration.py  b  ball volume fraction near the surface vs dimension
  fibonacci_chain.py      c  Z^2 cut at slope 1/phi, projected to a 1D chain
  golden_angle.py         d  phyllotaxis points: golden angle vs perturbed angles
  simplex_4d.py           e  5-cell: rotate in 4D planes, project to 2D
runs/*.json            one file per run, schema below
tests/                 unittest, one file per experiment plus emlib
DIVERGENCE_LOG.md      where human and AI interest differ (columns = schema)
```

## Run schema (`runs/*.json`)

```text
params            inputs, including the seed
outputs           everything computed, before any summary
angles            measured angles (degrees) or null when the run has none
human_interest    null until the human fills it
ai_interest       the AI's reading, or "not interpretable atm"
compressibility   LZ76 phrase count of a declared symbol string, with its definition
nature_match      a nature_catalog id, or null
divergence_note   null until the two interest fields differ
```

`emlib.write_run` refuses a run that is missing a key, or that carries a
non-null `human_interest`.

## Commands

```sh
python3 -m unittest discover -s tests          # all tests
python3 experiments/<name>.py                  # writes runs/<name>_seed<k>.json
```
