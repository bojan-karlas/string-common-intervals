#!/usr/bin/python -tt
# -*- coding: utf-8 -*-

# Test script for the rangemaxq.py file.

import unittest
import random

# Target of testing
import rangemaxq

class RangeMaxQueryTests(unittest.TestCase):

  # Tests that check if the RMQ functions properly handle a
  # single-element array.

  def common_test_singleton(self, rmq, rmq_pre):
    values = [1]
    rmqtable = rmq_pre(values)
    maxind_actual = rmq(0, 0, rmqtable)
    maxind_expected = 0
    self.assertEqual(maxind_actual, maxind_expected)

  def test_rmq1_singleton(self):
    self.common_test_singleton(rangemaxq.rmq1, rangemaxq.rmq1_pre)

  def test_rmq2_singleton(self):
    self.common_test_singleton(rangemaxq.rmq2, rangemaxq.rmq2_pre)

  def test_rmq3_singleton(self):
    self.common_test_singleton(rangemaxq.rmq3, rangemaxq.rmq3_pre)

  def test_rmq_singleton(self):
    self.common_test_singleton(rangemaxq.rmq, rangemaxq.rmq_pre)

  # Tests that check if the RMQ functions properly handla a 10-element
  # array made up of the same value.

  def common_test_same_values(self, rmq, rmq_pre, i, j):
    values = [1 for x in range(10)]
    rmqtable = rmq_pre(values)
    maxind_actual = rmq(i, j, rmqtable)
    maxind_expected = i if i <= j else -1
    self.assertEqual(maxind_actual, maxind_expected)

  def test_rmq1_same_values(self):
    self.common_test_same_values(rangemaxq.rmq1, rangemaxq.rmq1_pre, 0, 0)
    self.common_test_same_values(rangemaxq.rmq1, rangemaxq.rmq1_pre, 1, 0)
    self.common_test_same_values(rangemaxq.rmq1, rangemaxq.rmq1_pre, 0, 9)
    self.common_test_same_values(rangemaxq.rmq1, rangemaxq.rmq1_pre, 5, 9)
    self.common_test_same_values(rangemaxq.rmq1, rangemaxq.rmq1_pre, 9, 9)

  def test_rmq2_same_values(self):
    self.common_test_same_values(rangemaxq.rmq2, rangemaxq.rmq2_pre, 0, 0)
    self.common_test_same_values(rangemaxq.rmq2, rangemaxq.rmq2_pre, 1, 0)
    self.common_test_same_values(rangemaxq.rmq2, rangemaxq.rmq2_pre, 0, 9)
    self.common_test_same_values(rangemaxq.rmq2, rangemaxq.rmq2_pre, 5, 9)
    self.common_test_same_values(rangemaxq.rmq2, rangemaxq.rmq2_pre, 9, 9)

  def test_rmq3_same_values(self):
    self.common_test_same_values(rangemaxq.rmq3, rangemaxq.rmq3_pre, 0, 0)
    self.common_test_same_values(rangemaxq.rmq3, rangemaxq.rmq3_pre, 1, 0)
    self.common_test_same_values(rangemaxq.rmq3, rangemaxq.rmq3_pre, 0, 9)
    self.common_test_same_values(rangemaxq.rmq3, rangemaxq.rmq3_pre, 5, 9)
    self.common_test_same_values(rangemaxq.rmq3, rangemaxq.rmq3_pre, 9, 9)

  def test_rmq_same_values(self):
    self.common_test_same_values(rangemaxq.rmq, rangemaxq.rmq_pre, 0, 0)
    self.common_test_same_values(rangemaxq.rmq, rangemaxq.rmq_pre, 1, 0)
    self.common_test_same_values(rangemaxq.rmq, rangemaxq.rmq_pre, 0, 9)
    self.common_test_same_values(rangemaxq.rmq, rangemaxq.rmq_pre, 5, 9)
    self.common_test_same_values(rangemaxq.rmq, rangemaxq.rmq_pre, 9, 9)

  # Tests that check a random array of length 100.

  def common_test_random_values(self, rmq, rmq_pre, valuerange):
    N = 100
    random.seed(1)
    values = [random.randint(0,valuerange) for i in range(N)]
    rmqtable = rmq_pre(values)

    for i in range(N):
      for j in range(N):
        maxind_actual = rmq(i, j, rmqtable)
        maxind_expected = -1 if i > j else values[i:j+1].index(max(values[i:j+1])) + i
        self.assertEqual(maxind_actual, maxind_expected,
          'Expected: %d, Actual: %d - (i = %d; j = %d)'
          % (maxind_expected, maxind_actual, i, j))

  def test_rmq1_random_values(self):
    self.common_test_random_values(rangemaxq.rmq1, rangemaxq.rmq1_pre, 100)

  def test_rmq2_random_values(self):
    self.common_test_random_values(rangemaxq.rmq2, rangemaxq.rmq2_pre, 100)

  def test_rmq_random_values(self):
    self.common_test_random_values(rangemaxq.rmq, rangemaxq.rmq_pre, 100)

  # The test from rmq3 is a bit different because for this algorithm there is
  # a constraint saying that cosecutive values in the input array may differ
  # by either +1 or -1 (zero is not allowed).
  def test_rmq3_random_values(self):
    N = 100
    random.seed(1)
    valuediffs = [1 if random.randint(0,1) == 1 else -1 for i in range(N-1)]
    values = [N]
    for i in range(N-1): values.append(values[-1] + valuediffs[i])

    rmqtable = rangemaxq.rmq3_pre(values)

    for i in range(N):
      for j in range(N):
        maxind_actual = rangemaxq.rmq3(i, j, rmqtable)
        maxind_expected = -1
        if i <= j: maxind_expected = values[i:j+1].index(max(values[i:j+1])) + i

        self.assertEqual(maxind_actual, maxind_expected,
          'Expected: %d, Actual: %d - (i = %d; j = %d)'
          % (maxind_expected, maxind_actual, i, j))


  # Testing components of the final RMQ solution: Cartesian tree and Euler tour.

  # Check that the values of ancestor nodes are larger than their descendants.
  def test_cartesian_tree_values(self):
    N = 100
    random.seed(1)
    values = [random.randint(0, 100) for i in range(N)]
    (left, right, root) = rangemaxq.build_cartesian(values)

    maxind = values.index(max(values))
    self.assertEqual(maxind, root)

    for i in range(len(left)):
      if left[i] >= 0:
        self.assertGreaterEqual(values[i], values[left[i]])

      if right[i] >= 0:
        self.assertGreaterEqual(values[i], values[right[i]])

  # Check that for every tree node, the index of its left child is smaller and
  # that the index of its right child is greater.
  def test_cartesian_tree_indices(self):
    N = 100
    random.seed(1)
    values = [random.randint(0, 100) for i in range(N)]
    (left, right, root) = rangemaxq.build_cartesian(values)
    for i in range(len(left)):
      if left[i] >= 0:
        self.assertGreaterEqual(i, left[i])

      if right[i] >= 0:
        self.assertGreaterEqual(right[i], i)

  # Check that all consecutive nodes in the euler tour are in a parent-child
  # relationship, if they don't represent the same node.
  def test_euler_tour_parent_child(self):
    N = 100
    random.seed(1)
    values = [random.randint(0, 100) for i in range(N)]
    (left, right, root) = rangemaxq.build_cartesian(values)
    (tour, levels, represent) = rangemaxq.euler_tour(left, right, root)

    for i in range(len(tour) - 1):
      idx = tour[i]
      idxn = tour[i + 1]

      # Because of the constraint for RMQ3, the tour should not contain
      # the same element on two consecutive positions.
      self.assertNotEqual(idx, idxn)

      self.assertTrue(
        left[idx] == idxn or right[idx] == idxn or
        left[idxn] == idx or right[idxn] == idx,
        'Nodes %d and %d are not in a parent-child relationship.' % (idx, idxn))

  # Check that the represent array indeed contains, for each node, the position
  # of the first appearance of that node in the Euler tour.
  def test_euler_represent(self):
    N = 100
    random.seed(1)
    values = [random.randint(0, 100) for i in range(N)]
    (left, right, root) = rangemaxq.build_cartesian(values)
    (tour, levels, represent) = rangemaxq.euler_tour(left, right, root)

    self.assertEquals(len(values), len(represent))

    for i in range(len(values)):
      idx = tour.index(i)
      self.assertEquals(idx, represent[i])

  # Check that the for each appearance of a node in the tour, the level has
  # been correctly assigned.
  def test_euler_tour_levels(self):

    def check_level(node, level):
      indices = [i for i, x in enumerate(tour) if x == node]
      for i in indices:
        self.assertEquals(level, levels[i])

    N = 100
    random.seed(1)
    values = [random.randint(0, 100) for i in range(N)]
    (left, right, root) = rangemaxq.build_cartesian(values)
    (tour, levels, represent) = rangemaxq.euler_tour(left, right, root)

    stack = [root]
    check_level(root, max(levels))

    while len(stack) > 0:
      node = stack.pop()
      level = levels[represent[node]]

      if left[node] >= 0:
        check_level(left[node], level - 1)
        stack.append(left[node])

      if right[node] >= 0:
        check_level(right[node], level - 1)
        stack.append(right[node])


def main():
  unittest.main()

if __name__ == '__main__':
  main()
