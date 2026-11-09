"""Tests for exact, published aggregate example only (not original studies)."""
import sys
import unittest
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from partial_collision_certificate import (EqualityBlock, dual_bound_at, optimal_upper_bound,
                                           minimum_certifying_block_count)


class PartialCertificateTests(unittest.TestCase):
    def setUp(self):
        self.alpha = Fraction(1, 20)
        self.blocks = (EqualityBlock('PDACx20', 32, 6),
                       EqualityBlock('NETx20', 2, 1),
                       EqualityBlock('Otherx20', 2, 0))

    def test_single_block_published_upper_bound(self):
        value = dual_bound_at(86, 31, self.alpha, self.blocks[:1], Fraction(80,11))
        self.assertEqual(value, Fraction(435, 946))

    def test_two_blocks_published_upper_bound(self):
        value = dual_bound_at(86, 31, self.alpha, self.blocks[:2], Fraction(80,11))
        self.assertEqual(value, Fraction(210, 473))

    def test_minimum_certificates(self):
        self.assertEqual(minimum_certifying_block_count(86,31,self.alpha,Fraction(1,2),self.blocks), 1)
        self.assertEqual(minimum_certifying_block_count(86,31,self.alpha,Fraction(9,20),self.blocks), 2)

    def test_casewise_upper_bound(self):
        value, _, _ = optimal_upper_bound(86,31,self.alpha,self.blocks,0)
        self.assertEqual(value, Fraction(550,817))


if __name__ == '__main__':
    unittest.main()
