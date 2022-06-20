#!/usr/bin/python3.10
import sys
from BoyerMoorePatternMatching.gusfields_z_alg import z_alg


def z_alg_hamming_distance_1(txt: str, pat: str) -> list[int | None]:
    """
    Finds all the occurrences of 'txt' in 'pat' with a Hamming distance ≤ 1
    :param txt: string to match 'pat' to
    :param pat: substring pattern
    :return: an array where the index is the Hamming distance (strictly ≤ 1)
        0: An exact match of 'txt' in 'pat' (i.e. Hamming distance of 0)
        1: A match of 'txt' in 'pat' but with one mismatch character (i.e. Hamming distance of 1)
        None: A Hamming distance > 1
    """
    if not pat or not txt: return [None] * max(len(txt), len(pat))  # empty string

    # special character must not be in the alphabet of the inputs
    str_std: str = pat + "\t" + txt
    str_rev: str = pat[::-1] + "\t" + txt[::-1]

    # removing everything before and including the special character
    z_pat_std = z_alg(str_std)[len(pat) + 1:len(str_std) - len(pat) + 1]
    z_pat_rev = z_alg(str_rev)[len(str_std) - len(pat):len(pat):-1]

    return [1 if sum(a) == len(pat) - 1 else
            0 if sum(a) == len(pat) * 2 else
            None
            for a in zip(z_pat_std, z_pat_rev)] + [None] * (len(pat) - 1)


def hd1_patmatch() -> None:
    """
    Command line version of z_alg_hamming_distance_1
    """
    txt_inp: str = open(sys.argv[1], "r").read()
    pat_inp: str = open(sys.argv[2], "r").read()
    hamm_dist_arr = z_alg_hamming_distance_1(txt_inp, pat_inp)
    with open("output_hd1_patmatch.txt", "w") as output_file:
        for i, el in enumerate(hamm_dist_arr):
            if el is not None:
                output_file.write(str(i + 1) + " " + str(el) +
                                  ("\n" if i != len(hamm_dist_arr) - 1 else ""))


if __name__ == "__main__":
    hd1_patmatch()

    # print(z_alg_hamming_distance_1("uyfmxcvhutualb", "h"))
