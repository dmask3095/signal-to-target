#!/usr/bin/env python
# /// script
# dependencies = ["pandas"]
# ///
import sys; sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__file__), "..", "..", ".."))
from s2t.faers import *  # noqa
if __name__ == "__main__":
    import runpy; runpy.run_module("s2t.faers", run_name="__main__")
