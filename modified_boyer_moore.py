#!/usr/bin/python3.10
import sys
from BoyerMoorePatternMatching.gusfields_z_alg import z_alg


def modified_boyer_moore_alg(txt: str, pat: str):
    """
    Sub-linear time pattern matching algorithm (without 'bad character rule', superseded by an improved 'good suffix
    rule')
    :param txt: base string to find instances of substring `pat` in
    :param pat: substring pattern of interest to find in `txt`
    :return: indices of exact matches of pat in txt (matching from left-to-right). Note: zero-indexed
    """
    if not pat or not txt: return []

    PAT_LEN = len(pat)
    PAT_LAST_IDX = PAT_LEN - 1

    ascii_pat: list[int] = list(map(ord, pat))
    max_alphabet_pat = max(ascii_pat)
    min_alphabet_pat = min(ascii_pat)
    alphabet_ord_delta = max_alphabet_pat - min_alphabet_pat + 1

    char_to_idx = lambda c: ord(c) - min_alphabet_pat
    is_idx_within_bounds = lambda e: 0 <= e < alphabet_ord_delta

    z_std = z_alg(pat)
    z_rev = z_alg(pat[::-1])

    z_suffix: list[int | None] = z_rev[::-1] + [None]
    matched_prefix: list[int] = [0] * (PAT_LEN + 1)
    spec_gs_table: list[list[int | None]] = [[None for _ in range(PAT_LEN + 1)]
                                             for _ in range(alphabet_ord_delta)]

    s_string_matches: list[int] = []  # zero-indexed txt characters

    # pre-processing modified good suffix table
    for p in range(PAT_LAST_IDX):
        if p - z_suffix[p] >= 0:
            spec_gs_table[char_to_idx(pat[p - z_suffix[p]])][PAT_LEN - z_suffix[p]] = p + 1
    spec_gs_table[char_to_idx(pat[PAT_LAST_IDX])][PAT_LEN] = PAT_LEN

    # pre-processing matched prefix values
    matched_prefix[0] = PAT_LEN
    for i in range(PAT_LEN - 1, 0, -1):
        matched_prefix[i] = z_std[i] if z_std[i] == PAT_LEN - i else matched_prefix[i + 1]

    n_gs_shift = 1
    main_idx = PAT_LAST_IDX
    stop_suf_repeat = start_suf_repeat = None
    while main_idx < len(txt):
        idx_pat = PAT_LAST_IDX

        while True:

            if stop_suf_repeat is not None and idx_pat == stop_suf_repeat:
                idx_pat = start_suf_repeat - 1
                stop_suf_repeat = start_suf_repeat = None
            idx_text = main_idx - PAT_LAST_IDX + idx_pat + (1 if idx_pat == -1 else 0)

            if idx_pat >= 0 and pat[idx_pat] != txt[idx_text]:

                # good suffix rule
                mismatch_char = char_to_idx(txt[idx_text])
                mismatch_char_in_bounds = is_idx_within_bounds(mismatch_char)
                good_suffix_v = spec_gs_table[mismatch_char][idx_pat + 1] if mismatch_char_in_bounds else None
                if good_suffix_v is None:
                    n_gs_shift = PAT_LEN - matched_prefix[idx_pat + 1]

                    start_suf_repeat = 0
                    stop_suf_repeat = matched_prefix[idx_pat + 1] - 1
                elif good_suffix_v > 0:
                    n_gs_shift = PAT_LEN - good_suffix_v

                    start_suf_repeat = good_suffix_v - PAT_LAST_IDX + idx_pat
                    stop_suf_repeat = good_suffix_v - 1
                break
            elif idx_pat <= 0:  # exact match
                n_gs_shift = PAT_LEN - matched_prefix[1]
                s_string_matches.append(idx_text)

                start_suf_repeat = 0
                stop_suf_repeat = matched_prefix[1] - 1
                break

            idx_pat -= 1

        main_idx += n_gs_shift

    return s_string_matches


def modified_BoyerMoore() -> None:
    """
    Command line version of boyer_moore_alg
    """
    txt_inp: str = open(sys.argv[1], "r").read()
    pat_inp: str = open(sys.argv[2], "r").read()
    b_more_arr = modified_boyer_moore_alg(txt_inp, pat_inp)
    with open("output_modified_BoyerMoore.txt", "w") as output_file:
        for i, el in enumerate(b_more_arr):
            output_file.write(str(el + 1) + ("\n" if i != len(b_more_arr) - 1 else ""))  # one-indexed txt characters


if __name__ == "__main__":
    modified_BoyerMoore()

    # print(modified_boyer_moore_alg("cccccbabaaacccccbaaaaa", "cccccbaaaaa"))
