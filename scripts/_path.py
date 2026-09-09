"""Makes the esbm package importable when a script is run directly (python scripts/x.py)."""
import os, sys
_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if _root not in sys.path:
    sys.path.insert(0, _root)
