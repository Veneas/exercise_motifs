# Lesson 03

# =-=-=-=-=-=-=-=-=-=-=-=-= Task 1 =-=-=-=-=-=-=-=-=-=-=-=-=
def count_matrix(motifs):
    """
    Builds a dict counting the occurences of each base at every
    position (column) across a list of motifs.
    :param motifs:
    :return:
    """

    # Finding the len() of single motif
    l = len(motifs[0])

    # Creating dict with zero-filled lists of the length[0] for each base
    counts = {
        "A": [0] * l,
        "C": [0] * l,
        "G": [0] * l,
        "T": [0] * l,
    }

    # Going through every motif seq
    for motif in motifs:
        # Going through every char and its index in the motif
        for i, base in enumerate(motif):
            # Icrement the count for that specific base at the specific position
            counts[base][i] += 1

    return counts


def score(motifs):
    """
    Calculates the score of a list of motifs. The score is the sum of the
    max counts for each columns in the count matrix.
    :param motifs:
    :return:
    """

    counts = count_matrix(motifs)
    l = len(motifs[0])
    total_score = 0

    # Looping through each column index
    for i in range(l):
        # Finding the highest count in this columns among all four bases
        max_in_column = max(counts["A"][i],
                            counts["C"][i],
                            counts["G"][i],
                            counts["T"][i])
        total_score += max_in_column

    return total_score


def consensus(motifs):
    """
    Returns the consensus string (the most freq base of
    every position in motif. Ties are resolved by picking
    the first base in the order A, C, G, T.
    :param motifs:
    :return:
    """

    counts = count_matrix(motifs)
    l = len(motifs[0])
    consensus_string = ""

    # Checking each column
    for i in range(l):
        best_base = ""
        highest_count = -1

        # Checking exactly in the A,C,G,T order -> tie breaker rule
        # Update best_base only if the count is GREATER.
        for base in "ACGT":
            if counts[base][i] > highest_count:
                highest_count = counts[base][i]
                best_base = base

        consensus_string += best_base

    return consensus_string


def hamming_distance(a, b):
    """
    Calculates the num of mismatched characters between two string of equal length.
    :param a:
    :param b:
    :return:
    """

    # Zipping (a,b) pairs up the chars of both strings side-by-side.
    # Adding 1 for every pair where the chars are not eq

    return sum(1 for char_a, char_b in zip(a,b) if char_a != char_b)


def total_distance(pattern, sequences):
    """
    For every seq, find the smallest Hamming dist between the pattern and any
    of its sliding windows, then summ minimums up.
    :param pattern:
    :param sequences:
    :return:
    """

    l = len(pattern)
    total_dist = 0

    # Processing one seq at a time
    for seq in sequences:
         # Start with infinitely high min dist
        min_dist_for_seq = float("inf")

        # Sliding the window of length 'l' across seq
        # range(len(seq) -l +1) to not go out of bounds at the end

        for i in range(len(seq) -l +1):
            window = seq[i: i+l]

            # Calculating dist between patter and the specific window
            current_dist = hamming_distance(pattern, window)

            # If the closes match so far, update record
            if current_dist < min_dist_for_seq:
                min_dist_for_seq = current_dist

        # Adding the best (minimum) dist from this eq to the total
        total_dist += min_dist_for_seq

    return total_dist


def main():
    lecture_dna = [
        "TGACGTATAAGTTGCGATGGACGAGATAGCAGAGAATAGGCAACGAGAGATAAGCAG",
        "GACGGTAGCAGATAGACAGATGAAGAGTATGAATTGCACAGATAGCAGATAGCAGAT",
        "GGAGTGTGACGTAGCAGAGACGAAAGACGTAGAGTAGCAGTAGCAGATAGAGGGAGT",
        "TAGACAGTATAGAGACAGCGAGTCGGATAGCACCCAGTATGACGATAGCAATGACAG",
        "GCAGTAGAGCAGATTAGCATTGACAGATAGACGATTGGAGAGATGTGTGGATGACGA",
        "GGCAGGTAGCACACTGGGTCGATAAAGAGTAGCATAGAGACATAGACATATTTTAGC",
    ]

    # =-=-=-=-=-=-=-=-=-=-=-=-= TASK 1 =-=-=-=-=-=-=-=-=-=-=-=-=
    red = ["TAAGTT", "TGAATT", "GGAGTG", "CGAGTC", "TGTGTG", "TGGGTC"]  # slide 19
    best = ["AGATAG", "AGATAG", "AGATAG", "AGACAG", "AGATAG", "AGGTAG"]

    print(score(red))  # 26
    print(consensus(best), score(best))  # AGATAG 34
    print(hamming_distance("TAAGTT", "TGAATT"))  # 2
    print(total_distance("TGCGTT", lecture_dna))  # 13

if __name__ == "__main__":
    main()

