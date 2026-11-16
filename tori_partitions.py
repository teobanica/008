#!/usr/bin/env python3
"""Exact partition and word calculations for Maximal torus theory for crossed product quantum groups.

Run with Python 3.8 or later. No external packages or input files are needed.
"""

import argparse
from collections import Counter
from itertools import product
from math import gcd, prod


def rotation_relators(a, b, arity):
    """Support of Q^tensor-arity applied to the one-block fixed tensor.

    Q is the rotation with rows (a,-b), (b,a), divided by sqrt(a*a+b*b).
    The common nonzero denominator is irrelevant for the support.
    """
    if a == b == 0:
        raise ValueError("the zero pair does not define a rotation")
    matrix = ((a, -b), (b, a))
    for word in product(range(2), repeat=arity):
        coefficient = sum(prod(matrix[i][j] for i in word) for j in range(2))
        if coefficient:
            yield word


def rotation_order(a, b):
    """Order of gh after the four-block relations, with 0 denoting infinity."""
    assert set(rotation_relators(a, b, 2)) == {(0, 0), (1, 1)}
    order = 0
    for word in rotation_relators(a, b, 4):
        reduced = []
        for letter in word:
            if reduced and reduced[-1] == letter:
                reduced.pop()
            else:
                reduced.append(letter)
        # An even reduced word in two involutions is (gh)^k or its inverse.
        assert len(reduced) % 2 == 0
        order = gcd(order, len(reduced) // 2)
    return order


def multiply_letter(word, factor, step, modulus):
    """Multiply a reduced syllable word in the free product by one generator."""
    if word and word[-1][0] == factor:
        exponent = (word[-1][1] + step) % modulus
        prefix = word[:-1]
    else:
        exponent = step % modulus
        prefix = word
    if exponent:
        return prefix + ((factor, exponent),)
    return prefix


def ball_sizes(factors, modulus, radius):
    """Breadth-first enumeration in the symmetric generating set."""
    visited = {()}
    frontier = {()}
    balls = [1]
    for _ in range(radius):
        next_frontier = {
            multiply_letter(word, i, step, modulus)
            for word in frontier for i in range(factors) for step in (-1, 1)
        } - visited
        visited.update(next_frontier)
        frontier = next_frontier
        balls.append(len(visited))
    return balls


def counted_balls(factors, modulus, radius):
    """Independent coefficient count by the length of the first syllable."""
    lengths = Counter(min(j, modulus-j) for j in range(1, modulus))
    # continuations[n] counts length-n continuations after a fixed factor.
    continuations = [1]
    balls = [1]
    for n in range(1, radius + 1):
        tails = sum(m * continuations[n-k] for k, m in lengths.items() if k <= n)
        continuations.append((factors-1) * tails)
        balls.append(balls[-1] + factors * tails)
    return balls


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bound', type=int, default=16,
                        help='maximum absolute numerator entries in the rotation scan')
    parser.add_argument('--radius', type=int, default=6,
                        help='word radius for the free products (memory grows exponentially)')
    args = parser.parse_args()
    if args.bound < 1 or args.radius < 0:
        parser.error('use --bound >= 1 and --radius >= 0')
    counts = Counter()
    for a, b in product(range(-args.bound, args.bound+1), repeat=2):
        if gcd(a, b) != 1:
            continue
        order = rotation_order(a, b)
        predicted = 0 if a*b == 0 else (2 if abs(a) == abs(b) else 1)
        assert order == predicted, (a, b, order)
        counts[order] += 1
    print('H_2^+ rotations: primitive signed pairs, bound', args.bound)
    for order, name in [(0, 'infinite dihedral'), (2, 'Klein four'), (1, 'cyclic of order two')]:
        print(f'  {name}: {counts[order]}')
    print('Balls in free products of N cyclic groups of order s, starting at radius 0')
    for factors, modulus in [(2, 2), (2, 3), (2, 4), (3, 3)]:
        balls = ball_sizes(factors, modulus, args.radius)
        assert balls == counted_balls(factors, modulus, args.radius)
        print(f'  N={factors}, s={modulus}:', *balls)


if __name__ == '__main__':
    main()
