# Reproducibility contract

## Reference environment

- Python 3.12.3
- Exact Python packages: `requirements.lock`
- Source release: `v1.0.0`
- Random seeds: frozen in `config/*.yaml`
- Validation units: participant, population category, or cancer lineage; never random rows
- Public inputs: retrieved into ignored `data/raw/` paths and checksum verified

## Level 1: verify the archived evidence

This is the fast integrity path. It verifies every file named in each core analysis manifest and
checks the principal sample counts, decisions, null ensembles, sensitivity settings, and claim
boundaries.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.lock
.venv/bin/pip install --no-deps -e .
.venv/bin/python scripts/verify_release.py
.venv/bin/pytest -q
```

## Level 2: regenerate from public source data

`make reproduce` retrieves the source datasets, rebuilds the direct HumanGEM candidate catalogue,
runs the human and CCLE benchmarks, executes all hard-null and sensitivity analyses, regenerates
the mapping supplement and figures, and reruns the publication synthesis.

```bash
make env
make reproduce
```

The reference analyses use up to ten parallel workers. A 12-core, 64-GB machine is sufficient;
a GPU is not required. A complete clean-room run on 31 August 2026 took 37 minutes 34 seconds from
the first verified input through final synthesis on an Intel Xeon E-2176G host with 6 physical
cores, 12 threads and 62 GiB RAM. The dominant ST002081 reconstruction took about 20 minutes.
Runtime depends on network speed and BLAS implementation. The CCLE RNA file is approximately
144 MB, while the Workbench and Human-GEM inputs are smaller.

## Determinism and expected differences

Model splits, graph randomizations, bootstraps, and null draws use fixed seeds. All 83 files named
by core output manifests regenerate byte-identically under the locked environment. Run manifests
themselves record timestamps and absolute local paths and are therefore expected to differ; each
regenerated manifest recomputes checksums for its own outputs. The release verifier evaluates
declared output integrity and scientific invariants rather than requiring byte-identical run
metadata.

The clean-room validation record is in `docs/clean-room-validation.md`.


## Failure policy

Retrieval stops on a source-checksum mismatch. Analysis synthesis stops if any declared output is
missing or has changed. Such failures indicate upstream drift, an incomplete run, or a modified
artifact and must not be silently bypassed.
