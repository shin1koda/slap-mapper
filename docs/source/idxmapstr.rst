
.. _idxmapstr:

Index Mapping String Format
===========================

The *index mapping string* is a compact textual syntax for specifying
explicit correspondences between groups of nodes. It is used in the
``map_3d`` interface of :class:`SlapAAM` to impose mapping constraints
or to encode mapping results. In interactive mode, it is also used to
record the history of user choices.


An index mapping string consists of one or more *mapping blocks*,
separated by semicolons::

    <indices_left> >> <indices_right> ; <indices_left> >> <indices_right> ; ...

Each block defines a correspondence between a set of indices on the left
and a set of indices on the right.

The format is parsed by :func:`parse_index_mapping_string` in ``_idxmapstr.py``.
Below is the complete specification.

Mapping Block Structure
-----------------------

A mapping block has the form::

    1-3,5 >> 4-6

Here:

- ``1-3,5`` denotes a set of indices on the left
- ``4-6`` denotes a set of indices on the right

Both sides must contain the **same number of unique indices**. If they do
not match, a ``ValueError`` is raised.


Index Set Syntax
----------------

An index set is a comma-separated list of items. Each item can be:

- a single integer (e.g., ``3``), or
- a range written as ``start-end`` (e.g., ``3-6`` meaning ``3,4,5,6``)

Examples::

    5                → [5]
    1-3              → [1, 2, 3]
    1,4-7,10         → [1, 4, 5, 6, 7, 10]

Duplicate entries are removed automatically, and output indices are always
returned in sorted order.


Multiple Mappings
-----------------

Mapping blocks are separated by semicolons::

    1-3 >> 4-6 ; 7 >> 8

The parser returns a list of pairs::

    [
        ([1, 2, 3], [4, 5, 6]),
        ([7],       [8])
    ]


Index Base (0-based vs 1-based)
-------------------------------

The parser accepts a ``base`` argument to interpret the input indices:

- ``base = 0`` (default): indices are interpreted as zero-based.
- ``base = 1``: indices in the string are interpreted as one-based and
  converted internally to zero-based values.

Example::

    pair_str = "1-3 >> 4-6"
    base = 1

Output::

    [([0, 1, 2], [3, 4, 5])]


