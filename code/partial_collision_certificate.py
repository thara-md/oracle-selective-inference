"""Exact-rational illustration of partial observation-equality certificates.

This is independently authored demonstration code for the aggregate example
reported in Hara and Konishi's Oracle manuscript. It is NOT the archived
clinical analysis code, nor the simulation pipeline used in the manuscript.
No patient-level information or model outputs are included.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations


@dataclass(frozen=True)
class EqualityBlock:
    name: str
    n: int
    errors: int

    def __post_init__(self) -> None:
        if not (0 <= self.errors <= self.n):
            raise ValueError("Expected 0 <= errors <= block size")


def positive(x: Fraction) -> Fraction:
    return max(Fraction(0), x)


def candidate_multipliers(alpha: Fraction, blocks: tuple[EqualityBlock, ...]) -> tuple[Fraction, ...]:
    """The finite dual search set from manuscript Section 3.4."""
    if not (Fraction(0) <= alpha < Fraction(1)):
        raise ValueError("Risk alpha must be in [0, 1)")
    candidates = {Fraction(0), Fraction(1, 1) / (1 - alpha)}
    for block in blocks:
        denom = block.errors - alpha * block.n
        if denom > 0:
            candidates.add(Fraction(block.n) / denom)
    return tuple(sorted(candidates))


def dual_bound_at(
    n: int, total_errors: int, alpha: Fraction,
    blocks: tuple[EqualityBlock, ...], multiplier: Fraction,
) -> Fraction:
    """Upper coverage bound from enforcing the specified disjoint blocks."""
    if n <= 0 or not (0 <= total_errors <= n):
        raise ValueError("Invalid dataset totals")
    if multiplier < 0:
        raise ValueError("Multiplier must be nonnegative")
    if sum(b.n for b in blocks) > n or sum(b.errors for b in blocks) > total_errors:
        raise ValueError("Block counts exceed cohort totals")

    baseline = (n - total_errors) * (1 + alpha * multiplier)
    baseline += total_errors * positive(1 - (1 - alpha) * multiplier)
    reduction = Fraction(0)
    for block in blocks:
        casewise = (block.n - block.errors) * (1 + alpha * multiplier)
        casewise += block.errors * positive(1 - (1 - alpha) * multiplier)
        shared = positive(block.n - multiplier * (block.errors - alpha * block.n))
        reduction += casewise - shared
    return (baseline - reduction) / n


def optimal_upper_bound(
    n: int, total_errors: int, alpha: Fraction,
    blocks: tuple[EqualityBlock, ...], max_blocks: int | None = None,
) -> tuple[Fraction, tuple[str, ...], Fraction]:
    """Minimum optimistic bound using at most k prespecified disjoint blocks.

    For a small number of blocks, this enumerates subsets for clarity.
    It checks exact rational dual candidates and may be slow for many blocks.
    """
    if max_blocks is None:
        max_blocks = len(blocks)
    if not (0 <= max_blocks <= len(blocks)):
        raise ValueError("Invalid max_blocks")
    if len({b.name for b in blocks}) != len(blocks):
        raise ValueError("Block names must be unique")
    best = None
    for k in range(max_blocks + 1):
        for subset in combinations(blocks, k):
            for lam in candidate_multipliers(alpha, blocks):
                bound = dual_bound_at(n, total_errors, alpha, subset, lam)
                proposal = (bound, tuple(b.name for b in subset), lam)
                if best is None or (proposal[0], len(proposal[1]), proposal[1]) < (best[0], len(best[1]), best[1]):
                    best = proposal
    assert best is not None
    return best


def minimum_certifying_block_count(
    n: int, total_errors: int, alpha: Fraction,
    target_coverage: Fraction, blocks: tuple[EqualityBlock, ...],
) -> int | None:
    """First k with U_k(alpha) < target; a tie does NOT exclude the target."""
    if not (0 < target_coverage <= 1):
        raise ValueError("Target coverage must be in (0,1]")
    for k in range(len(blocks) + 1):
        bound, _, _ = optimal_upper_bound(n, total_errors, alpha, blocks, k)
        if bound < target_coverage:
            return k
    return None


def main() -> None:
    alpha = Fraction(1, 20)
    blocks = (EqualityBlock("PDACx20", 32, 6),
              EqualityBlock("NETx20", 2, 1),
              EqualityBlock("Otherx20", 2, 0))
    n, errors = 86, 31
    print("Published aggregate example: all-category GPT, n=86, errors=31, alpha=5%")
    for k in (0, 1, 2, 3):
        coverage, retained, lam = optimal_upper_bound(n, errors, alpha, blocks, k)
        print(f"k <= {k}: U={coverage} = {float(coverage)*100:.2f}% | blocks={retained} | lambda={lam}")
    for target in (Fraction(1, 2), Fraction(9, 20)):
        print(f"target={float(target)*100:g}% -> minimum blocks={minimum_certifying_block_count(n, errors, alpha, target, blocks)}")


if __name__ == "__main__":
    main()
