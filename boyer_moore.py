#!/usr/bin/python3.10
from BoyerMoorePatternMatching.gusfields_z_alg import z_alg


def exact_match_bm(txt: str, pat: str):
    """
    Boyer-Moore sub-linear time (average-case) pattern matching algorithm
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

    cap_r_shift_table: list[list[int | None]] = [[None for _ in range(PAT_LEN)]
                                                 for _ in range(alphabet_ord_delta)]

    # set up r-shift bad character table
    for p_idx, p_alphabet_idx in enumerate(map(char_to_idx, pat)):
        if p_idx + 1 == PAT_LEN: continue
        if is_idx_within_bounds(p_alphabet_idx):
            cap_r_shift_table[p_alphabet_idx][p_idx + 1] = p_idx
    for x in range(alphabet_ord_delta):
        for r_k in range(1, PAT_LEN):
            if cap_r_shift_table[x][r_k] is None: cap_r_shift_table[x][r_k] = cap_r_shift_table[x][r_k - 1]

    z_std = z_alg(pat)
    z_rev = z_alg(pat[::-1])

    z_suffix: list[int | None] = z_rev[::-1] + [None]
    good_suffix: list[int] = [0] * (PAT_LEN + 1)
    matched_prefix: list[int] = [0] * (PAT_LEN + 1)

    s_string_matches: list[int] = []  # zero-indexed txt characters

    # pre-processing good suffix values
    for p in range(PAT_LEN - 1): good_suffix[PAT_LEN - z_suffix[p]] = p + 1

    # pre-processing matched prefix values
    matched_prefix[0] = PAT_LEN
    for i in range(PAT_LEN - 1, 0, -1):
        matched_prefix[i] = z_std[i] if z_std[i] == PAT_LEN - i else matched_prefix[i + 1]

    n_gs_shift = n_bc_shift = 1
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
                # bad character rule
                char_idx = char_to_idx(txt[idx_text])
                r_k_x = cap_r_shift_table[char_idx][idx_pat] if is_idx_within_bounds(char_idx) else None
                n_bc_shift = idx_pat - r_k_x if r_k_x is not None else 1

                # good suffix rule
                good_suffix_v = good_suffix[idx_pat + 1]
                if good_suffix_v > 0:
                    n_gs_shift = PAT_LEN - good_suffix_v

                    start_suf_repeat = good_suffix_v - PAT_LAST_IDX + idx_pat
                    stop_suf_repeat = good_suffix_v - 1
                elif good_suffix_v == 0:
                    n_gs_shift = PAT_LEN - matched_prefix[idx_pat + 1]

                    start_suf_repeat = 0
                    stop_suf_repeat = matched_prefix[idx_pat + 1] - 1
                break
            elif idx_pat <= 0:  # exact match
                n_gs_shift = PAT_LEN - matched_prefix[1]
                s_string_matches.append(idx_text)

                start_suf_repeat = 0
                stop_suf_repeat = matched_prefix[1] - 1
                break
            idx_pat -= 1

        main_idx += max(n_gs_shift, n_bc_shift)

    return s_string_matches


if __name__ == "__main__":
    print(exact_match_bm("cccccbaxaaacccccbaaaaa", "cccccbaaaaa"))
