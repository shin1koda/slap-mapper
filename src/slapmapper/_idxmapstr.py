def parse_index_string(s, base=0):
    """Parse a string like '1,3-5' into a list of unique integers."""
    indices = set()
    for part in s.split(','):
        if '-' in part:
            start, end = map(int, part.split('-'))
            indices.update(range(start, end + 1))
        else:
            indices.add(int(part))
    if base == 1:
        indices = {i - 1 for i in indices}
    return sorted(indices)

def parse_index_mapping_string(pair_str, base=0):
    """Parse a string of the form '1-3,5>>4-6;7>>8' with index base option."""
    pairs = pair_str.split(';')
    result = []

    for pair in pairs:
        if '>>' not in pair:
            raise ValueError(f"Missing '>>' in pair: {pair}")
        left_str, right_str = pair.split('>>')
        left_indices = parse_index_string(left_str.strip(), base=base)
        right_indices = parse_index_string(right_str.strip(), base=base)

        if len(left_indices) != len(right_indices):
            raise ValueError(
                f"Mismatch in number of unique elements: {left_indices} vs {right_indices}"
            )

        result.append((left_indices, right_indices))

    return result

def fmt_idxs(idxs, base=0):
    """
    Convert a sorted list of indices into a compact string such as '1-3,5,7'.

    Parameters
    ----------
    idxs : list[int]
        Sorted list of indices.
    base : int, optional (default=0)
        If 1, convert output to 1-based indices.

    Returns
    -------
    str
        Compact string representation.
    """
    if not idxs:
        return ""

    idxs = [i + base for i in idxs]

    parts = []
    start = prev = idxs[0]

    for x in idxs[1:]:
        if x == prev + 1:
            prev = x
            continue

        if start == prev:
            parts.append(str(start))
        else:
            parts.append(f"{start}-{prev}")

        start = prev = x

    if start == prev:
        parts.append(str(start))
    else:
        parts.append(f"{start}-{prev}")

    return ",".join(parts)

def lgp2idx_map_str(lgp,base=0):
    l2i0 = lgp[0].label2idxs
    l2i1 = lgp[1].label2idxs

    return ";".join(
        [">>".join( [fmt_idxs(l2i0[l],base),fmt_idxs(l2i1[l],base)] ) for l in l2i0]
    )

