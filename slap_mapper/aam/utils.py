
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

def parse_label_string(pair_str, base=0):
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
