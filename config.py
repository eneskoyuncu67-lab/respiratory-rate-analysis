"""
Central path configuration for the BART respiration pipeline.

Every script below reads its input/output locations from here instead of
hardcoding them, so the pipeline is portable across machines and safe to
publish without exposing any one person's folder layout.

Defaults assume a `data/edf/` folder (your raw EDF recordings) and an
`output/` folder, both placed next to this file - so the pipeline runs
out of the box on a fresh clone once you drop your own EDF files into
`data/edf/`.

To point at a different location instead (e.g. a data drive), either:
  - set the environment variables BART_EDF_DIR / BART_OUT_DIR, or
  - copy local_config.py.example to local_config.py (gitignored) and set
    EDF_DIR / OUT_DIR there - it is loaded last and wins over everything
    above, so your real paths never need to touch a tracked file.
"""
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

EDF_DIR = os.environ.get('BART_EDF_DIR', os.path.join(BASE_DIR, 'data', 'edf'))
OUT_DIR = os.environ.get('BART_OUT_DIR', os.path.join(BASE_DIR, 'output'))

try:
    from local_config import *  # noqa: F401,F403 - optional, gitignored machine-specific overrides
except ImportError:
    pass
