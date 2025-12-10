"""
SLAPMapper
==========

This package provides a Python implementation of the approximate graph matching algorithm based on sequential linear assignment problems (SLAPs) and its application for atom-to-atom mapping.

The public API of the package is intentionally kept small and consists of:

* **LabeledGraph**
    Lightweight representation of a labeled graph.

* **SlapMapper**
    Sequential LAP-based approximate mapper for a pair of labeled graphs, designed for general graph matching problems.

* **SlapAAM**
    SLAPMapper for atom-to-atom mapping.

Various helper routines are **not** part of the public API and may change without notice.

All classes listed above can be imported directly as follows:

.. code-block:: python

    from slapmapper import LabeledGraph, SlapMapper
    from slapmapper.aam import SlapAAM

See each API page for detailed documentation.
"""
from .core import LabeledGraph, SlapMapper

__all__ = ["LabeledGraph","SlapMapper"]
