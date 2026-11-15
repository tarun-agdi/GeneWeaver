
def find_off_targets(genome, guide_rna, max_mismatches=2):
    """
    Find approximate matches of a guide sequence
    inside a DNA sequence.
    """

    genome = genome.upper()
    guide_rna = guide_rna.upper().replace("U", "T")

    guide_length = len(guide_rna)
    results = []

    for i in range(len(genome) - guide_length + 1):

        segment = genome[i:i + guide_length]

        mismatches = sum(
            1 for a, b in zip(segment, guide_rna)
            if a != b
        )

        if mismatches <= max_mismatches:
            results.append({
                "position": i,
                "sequence": segment,
                "mismatches": mismatches
            })

    return results


if __name__ == "__main__":

    genome = "ATGCGTACGTAGCTAGCTAGGCTAACGT"
    guide = "GTACGTAGCTAG"

    matches = find_off_targets(genome, guide)

    for match in matches:
        print(match)
