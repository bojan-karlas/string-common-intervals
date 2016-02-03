#!/usr/bin/python -tt
# -*- coding: utf-8 -*-

# Test script for the didier.py file.

import unittest
import random

# Target of testing
import didier

class DidierAlgorithmsTests(unittest.TestCase):

  # Tests for the Problem 1 implementation
  def test_maximal_locations(self):
    string = 'aabbabcdbabdcaabab'
    charset = ['a', 'b']
    expected_locations = [(0, 6), (8, 11), (13, 18)]
    actual_locations = didier.maximal_locations(string, charset)
    self.assertEqual(expected_locations, actual_locations)

  def test_compute_occurlist(self):
    string = 'aabbabcdbabdcaabab'
    expected_occurlist = ['a', 'b', 'c', 'd']
    actual_occurlist = didier.compute_occurlist(string)
    self.assertEqual(expected_occurlist, actual_occurlist)

  def test_compute_ranks(self):
    alphabet = ['a', 'b', 'c', 'd', 'e']
    occurlist = ['b', 'c', 'd']
    inf = float('inf')
    expected_ranklist = { 'a' : inf, 'b' : 0, 'c' : 1, 'd' : 2, 'e' : inf }
    actual_ranklist = didier.compute_ranks(occurlist, alphabet)

    self.assertEqual(len(expected_ranklist), len(actual_ranklist))
    for c in expected_ranklist.keys():
      self.assertEqual(expected_ranklist[c], actual_ranklist[c])

  def test_compute_rank_intervals(self):
    # Example from the paper.
    string = 'abcdacadbab'
    alphabet = list(set(string))
    occurlist = ['a', 'b', 'c', 'd']
    ranks = didier.compute_ranks(occurlist, alphabet)

    expected_rank_intervals = {
      0 : (0, 1), 1 : (0, 2), 2 : (0, 3), 3 : (0, 11), 4 : (4, 5), 5 : (4, 7),
      6 : (6, 7), 7 : (0, 11), 8 : (8, 11), 9 : (9, 10), 10 : (8, 11)
    }
    actual_rank_intervals = didier.compute_rank_intervals(string, ranks)
    self.assertEqual(expected_rank_intervals, actual_rank_intervals)

  def test_compute_rank_successors(self):
    # Example from the paper.
    string = 'abcdacadbab'
    alphabet = list(set(string))
    occurlist = ['a', 'b', 'c', 'd']
    ranks = didier.compute_ranks(occurlist, alphabet)
    expected_succ = [1, 2, 3, None, 1, 3, 1, None, 5, 8, 5]
    actual_succ = didier.compute_rank_successors(string, ranks)
    self.assertEqual(expected_succ, actual_succ)

  def test_maximal_substrings(self):
    # Example from the paper.
    string = 'abcdacadbab'
    expected_intervals = [(0, 1), (4, 5), (6, 7), (9, 10), (1, 2), (8, 9),
      (10, 11), (2, 3), (5, 6), (3, 4), (7, 8), (0, 2), (8, 11), (4, 7), (3, 5),
      (6, 8), (1, 3), (7, 9), (2, 4), (0, 3), (6, 11), (2, 8), (1, 4), (0, 11)]
    actual_intervals = didier.maximal_substrings(string)
    self.assertItemsEqual(actual_intervals, expected_intervals)

def main():
  unittest.main()

if __name__ == '__main__':
  main()
