"""Make the repository root importable so tests can import the analysis scripts."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
