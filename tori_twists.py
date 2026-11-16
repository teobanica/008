#!/usr/bin/env python3
"""Klein-four cocycle on the diagonal torus of U_2^+ * U_2^+.

Run: python3 tori_twists.py
All tests are exact, using integers, freely reduced words and group-ring
coefficients. The optional radius bounds the search over reduced words.
"""

import argparse
from collections import Counter
from itertools import product

ZERO = (0, 0)
LETTERS = (1, -1, 2, -2, 3, -3, 4, -4)
# We mark the first coordinate in each of the two fundamental blocks.
WEIGHTS = {1: (1, 0), 2: ZERO, 3: (0, 1), 4: ZERO}
BLOCKS = ((1, 2), (3, 4))


def add(p, q):
    return ((p[0] + q[0]) % 2, (p[1] + q[1]) % 2)


def sigma(p, q):
    return -1 if p[0] and q[1] else 1


def edge_sign(left, right):
    """A coefficient has (row weight, column weight)."""
    row, column = left
    other_row, other_column = right
    return sigma(row, other_row) * sigma(column, other_column)


def combined(left, right):
    return (add(left[0], right[0]), add(left[1], right[1]))


def weight(letter):
    return WEIGHTS[abs(letter)]


def append_letter(word, letter):
    if word and word[-1] == -letter:
        return word[:-1]
    return word + (letter,)


def degree_and_gauge(word):
    """Compute p(w) and b(w) from the ordered letters of a word."""
    degree = ZERO
    sign = 1
    for letter in word:
        step = weight(letter)
        sign *= sigma(degree, step)
        degree = add(degree, step)
    return degree, sign


def reduced_words(max_length, prefix=()):
    yield prefix
    if len(prefix) < max_length:
        for letter in LETTERS:
            if not prefix or prefix[-1] != -letter:
                yield from reduced_words(max_length, prefix + (letter,))


def test_cocycle():
    weights = list(product(range(2), repeat=2))
    for p in weights:
        assert sigma(ZERO, p) == sigma(p, ZERO) == 1
    for p, q, r in product(weights, repeat=3):
        assert sigma(p, q) * sigma(add(p, q), r) == sigma(q, r) * sigma(p, add(q, r))
    assert sigma((1, 0), (0, 1)) != sigma((0, 1), (1, 0))
    return len(weights) ** 3


def test_coefficients():
    def entries(block):
        return [(weight(i), weight(j)) for i, j in product(block, repeat=2)]

    first, second = map(entries, BLOCKS)
    assert edge_sign(first[1], second[0]) == -1  # u_12 followed by v_11
    assert edge_sign(second[0], first[1]) == 1
    forward = sum(edge_sign(a, b) < 0 for a, b in product(first, second))
    reverse = sum(edge_sign(b, a) < 0 for a, b in product(first, second))
    assert (forward, reverse) == (6, 0)
    all_entries = first + second
    for a, b, c in product(all_entries, repeat=3):
        assert (edge_sign(a, b) * edge_sign(combined(a, b), c)
                == edge_sign(b, c) * edge_sign(a, combined(b, c)))
    return forward, reverse, len(all_entries) ** 3


def test_word_gauge(radius):
    # Include inverse extensions, which cancel the last letter.
    checked = 0
    for word in reduced_words(radius):
        old_degree, old_sign = degree_and_gauge(word)
        for letter in LETTERS:
            new_degree, new_sign = degree_and_gauge(append_letter(word, letter))
            assert new_degree == add(old_degree, weight(letter))
            assert new_sign == old_sign * sigma(old_degree, weight(letter))
            checked += 1
    return checked


def test_coefficient_words(length):
    """Compare the twisted sign with the two free-word gauges, including adjoints."""
    generators = [(sign*i, sign*j) for block in BLOCKS
                  for i, j in product(block, repeat=2) for sign in (1, -1)]

    def check(row_word, col_word, degrees, sign, depth):
        _, row_gauge = degree_and_gauge(row_word)
        _, col_gauge = degree_and_gauge(col_word)
        assert row_gauge * col_gauge == sign
        checked = 1
        if depth < length:
            for i, j in generators:
                step = (weight(i), weight(j))
                checked += check(append_letter(row_word, i), append_letter(col_word, j),
                                 combined(degrees, step), sign * edge_sign(degrees, step), depth + 1)
        return checked

    return check((), (), (ZERO, ZERO), 1, 0)


def expand_character_product(labels, twisted):
    """Expand a toral character *-word in Z[F4], keeping multiplicities."""
    terms = Counter({(): 1})
    for label in labels:
        block = BLOCKS[abs(label) - 1]
        steps = block if label > 0 else tuple(-g for g in block)
        next_terms = Counter()
        for word, count in terms.items():
            for step in steps:
                sign = 1
                if twisted:
                    # Toral matrix coefficients have the same left and right weight.
                    p, _ = degree_and_gauge(word)
                    q = weight(step)
                    sign = edge_sign((p, p), (q, q))
                next_terms[append_letter(word, step)] += sign * count
        terms = next_terms
    return terms


def test_characters(length):
    checked = 0
    for n in range(length + 1):
        for labels in product((1, -1, 2, -2), repeat=n):
            assert expand_character_product(labels, True) == expand_character_product(labels, False)
            checked += 1
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius', type=int, default=4,
                        help='test all generator extensions of reduced words through this length')
    parser.add_argument('--character-length', type=int, default=4,
                        help='maximum length of character products (cost grows exponentially)')
    args = parser.parse_args()
    if args.radius < 0 or args.character_length < 0:
        parser.error('use --radius >= 0 and --character-length >= 0')
    print('Klein-four cocycle identities:', test_cocycle())
    forward, reverse, triples = test_coefficients()
    print(f'Negative mixed coefficient products: {forward}/16 forward, {reverse}/16 reverse')
    print('Coefficient associativity identities:', triples)
    print('Free-group coboundary identities on generator extensions:', test_word_gauge(args.radius))
    print('Signed coefficient words checked against the two gauges through length 4:',
          test_coefficient_words(4))
    checked = test_characters(args.character_length)
    print(f'Character *-words with identical toral expansion through length {args.character_length}: {checked}')


if __name__ == '__main__':
    main()
