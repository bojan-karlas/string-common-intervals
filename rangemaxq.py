#!/usr/bin/python -tt
# -*- coding: utf-8 -*-

import math

# Range Maximum Query (RMQ)
# Implementation of the various functions that solve the RMQ problem.
# Taken from the paper by Bender and Farach-Colton
# entitled 'The LCA Problem Revisited'

# Note: In the original paper, they solve the Range Minimum Query problem, which
#       equivalent to the Maximum version.

# TODO: In all rmq#() functions, the range should be expressed in the Python
#       fashion. I.e. [i,j] actually means from i, up to but not including j.

# Preprocessing step of the naïve solution for RMQ. Builds an O(n^2) table of
# maximal values for each pair of binding indices.
def rmq1_pre(values):
  N = len(values)
  rmqtable = [[-1 for x in range(N)] for x in range(N)]

  for i in range(N):
    # Init max with first elem of subarray starting at i.
    maxind = i

    for j in range(i, N):
      # Try to update the maxvalue.
      if values[j] > values[maxind]:
        maxind = j

      rmqtable[i][j] = maxind

  return rmqtable

# Naïve solution for RMQ. Uses an O(n^2) precomputed table.
def rmq1(i, j, rmqtable):
  return rmqtable[i][j]

# Preprocessing step for the Sparse Table (ST) solution for RMQ. Builds an
# O(n*log(n)) size lookup table in O(n*log(n)) time.
def rmq2_pre(values):
  N = len(values)
  M = int(math.ceil(math.log(N, 2)) + 1)
  # Init the table with every row initialized to the index of that row.
  # We use only the first column of this.
  lookup = [[i for j in range(M)] for i in range(N)]

  for j in range(1, M):
    for i in range(N):
      if i + 2**(j-1) >= N:
        lookup[i][j] = lookup[i][j-1]
      else:
        ind1 = lookup[i][j-1]
        ind2 = lookup[i + 2 ** (j-1)][j-1]

        if values[ind1] >= values[ind2]:
          lookup[i][j] = ind1
        else:
          lookup[i][j] = ind2

  rmqtable = {
    'lookup' : lookup,
    'values' : values
  }

  return rmqtable

# Sparse Table (ST) solution for RMQ. Uses an O(n*log(n)) size lookup table.
def rmq2(i, j, rmqtable):
  if i > j: return -1

  lookup = rmqtable['lookup']
  values = rmqtable['values']

  k = 0 if j==i else int(math.floor(math.log(j-i, 2)))
  ind1 = lookup[i][k]
  ind2 = lookup[j - 2**k + 1][k]

  if values[ind1] >= values[ind2]:
    return ind1
  else:
    return ind2

# Preprocessing step for the solution for plus/minus 1 RMQ. Builds an O(n) size
# data structure that is used for lookup. The only limitation is that the
# difference between any two consecutive elements in the input aray can be
# at most 1.
def rmq3_pre(values):
  # Split array into blocks of size log(N)/2 and find maximal values in
  # each block.
  N = len(values)
  M = int(math.ceil(math.log(N+1, 2) / 2))
  blocks = []
  blockmaxind = []
  blockmaxval = []

  for i in range(len(values)):
    j = i / M

    # Form blocks.
    if i % M == 0:
      blockmaxind.append(i)
      blockmaxval.append(values[i])
      if i + M < N:
        blocks.append(values[i:i+M])
      else:
        blocks.append(values[i:])

    # Keep track of maximal element in each block.
    if blockmaxval[j] < values[i]:
      blockmaxind[j] = i
      blockmaxval[j] = values[i]

  # Take the array of max values in each block and compute the ST-RMQ lookup.
  allblockslookup = rmq2_pre(blockmaxval)

  # For each block compute its stamp (binary diff over the elements).
  # For each distinct stamp compute compute a O(n^2) lookup table.
  inblocklookups = {}
  blockstamp = []
  for i in range(len(blocks)):
    # Compute block stamp.
    stamp = 0
    rblock = list(reversed(blocks[i]))
    for j in range(len(rblock) - 1):
      if rblock[j] < rblock[j+1]:
        stamp += 2**j
    blockstamp.append(stamp)

    # If stamp is a new one, then we haven't yet encountered a block of
    # this shape, so we need a new lookup table for it.
    if stamp not in inblocklookups:
      # We use the naïve algorithm to create a lookup table.
      inblocklookups[stamp] = rmq1_pre(blocks[i])

  # Finally we export all the data structures needed to do a lookup.
  rmqtable = {
    'values' : values,
    'blocks' : blocks,
    'blockmaxind' : blockmaxind,
    'blockmaxval' : blockmaxval,
    'blockstamp' : blockstamp,
    'inblocklookups' : inblocklookups,
    'allblockslookup' : allblockslookup
  }
  return rmqtable

# Solution for plus/minus 1 RMQ. The only limitation is that the difference
# between any two consecutive elements in the input aray can be at most 1.
def rmq3(i, j, rmqtable):
  if i > j: return -1

  values = rmqtable['values']
  N = len(values)
  M = int(math.ceil(math.log(N+1, 2) / 2))

  # Indices of blocks.
  block_i = i / M
  block_j = j / M

  # In case we are inside a single block, we just perform a single lookup.
  if block_i == block_j:
    stamp = rmqtable['blockstamp'][block_i]
    loockup = rmqtable['inblocklookups'][stamp]
    ind = rmq1(i % M, j % M, loockup)
    ind += block_i * M
    return ind

  # Otherwise, we need to do lookups from i to end of block_i and from beginning
  # of block_j to j.
  else:
    stamp_i = rmqtable['blockstamp'][block_i]
    loockup_i = rmqtable['inblocklookups'][stamp_i]
    ind1 = rmq1(i % M, M - 1, loockup_i)
    ind1 += block_i * M
    maxind = ind1

    stamp_j = rmqtable['blockstamp'][block_j]
    loockup_j = rmqtable['inblocklookups'][stamp_j]
    ind3 = rmq1(0, j % M, loockup_j)
    ind3 += block_j * M

    # In case block_i and block_j are not adjacent, we need to perform a third
    # lookup for the maximal element among all the blocks between them.
    if block_j - block_i > 1:
      maxblock = rmq2(block_i + 1, block_j - 1, rmqtable['allblockslookup'])
      ind2 = rmqtable['blockmaxind'][maxblock]
      if values[ind2] > values[ind1]:
        maxind = ind2

    if values[ind3] > values[maxind]:
      maxind = ind3

    return maxind

# Builds a Cartesian tree from a given array.
# Return a tuple (left_child, right_child, root)
#  - left_child, right_child - Array of same size as input array. Position i
#                              contains the index of the i'th node's
#                              left/right child, or -1 if there is no child.
#  - root                    - Index of the root node.
def build_cartesian(values):
  N = len(values)
  left = [-1 for x in range(N)]
  right = [-1 for x in range(N)]
  parent = [-1 for x in range(N)]

  for i in range(1, N):
    if values[i-1] >= values[i]:
      right[i-1] = i
      parent[i] = i-1
    else:
      v = i-1
      while values[v] < values[i] and parent[v] != -1:
        v = parent[v]

      if values[v] < values[i]:
        parent[v] = i
        left[i] = v

      else:
        if right[v] != -1:
          left[i] = right[v]
          parent[left[i]] = i

        right[v] = i
        parent[i] = v

  root = parent.index(-1)

  return (left, right, root)

# Builds an Euler tour of the given tree.
# Returns a tuple (tour, levels, represent):
#  - tour      - Indices of nodes given in the order of the tour. By definition
#                of the Euler tour, each node appears twice.
#  - levels    - Level within the tree of each node in the tour.
#                Note that here levels are given in a bottom-up order, so that
#                the parent node has a higher level than its children.
#  - represent - An array which in position i holds the first appearance of
#                the i'th node in the tour.
def euler_tour(left, right, root):
  stack = []
  curnode = root
  lastnode = -1
  tour = []
  levels = []
  represent = [-1 for i in range(len(left))]

  while len(stack) > 0 or curnode != -1:
    if curnode != -1:
      # First visit.
      represent[curnode] = len(tour)
      levels.append(len(stack))
      tour.append(curnode)

      stack.append(curnode)
      curnode = left[curnode]
    else:
      peeknode = stack[-1]
      if right[peeknode] != -1 and lastnode != right[peeknode]:
        # Second visit. Prevent repetition of nodes.
        if tour[-1] != peeknode:
          levels.append(len(stack) - 1)
          tour.append(peeknode)

        # If right child exists AND traversing node from left child, move right.
        curnode = right[peeknode]
      else:
        # Third visit. Prevent repetition of nodes.
        if tour[-1] != peeknode:
          levels.append(len(stack) - 1)
          tour.append(peeknode)

        lastnode = stack.pop()

  # Flip levels from top-down to bottom-up.
  maxlevel = max(levels)
  for i in range(len(levels)):
    levels[i] = maxlevel - levels[i]

  return (tour, levels, represent)

# Preprocessing step for the final solution for the general RMQ problem.
# Builds a O(n) size data structure which is then used to do lookups in
# constant time.
def rmq_pre(values):
  # Build a Cartesian tree of the N input values.
  (left, right, root) = build_cartesian(values)

  # Construct an Euler tour of the Cartesian tree. The tour of length 2N
  # contains pointers to nodes in the order of the tour. We also get a
  # 2N-sized array of levels of each node contained in the tour. Finally, we
  # get an N-sized array which at position i contains the index of the first
  # appearance of the i'th node within the tour.
  (tour, levels, represent) = euler_tour(left, right, root)

  # As every two consecutive nodes in the Euler tour have a parent-child
  # relationship, the respective levels differ by either +1 or -1. Therefore,
  # the levels array respects the plus/minus 1 constraint, so we can use the
  # plus/minus 1 RMQ algorithm.
  rmqtable = rmq3_pre(levels)

  # By reduction given in Lemma 1 of the paper, the result of the RMQ of range
  # [i,j] on the Euler tour is the Least Common Ancestor (LCA) between two
  # nodes i and j in a tree. Since our tree is a Cartesian tree constructed
  # on the input array, by Lemma 5 of the paper, the LCA of nodes i and j is
  # actually the maximal element of the range [i,j].

  rmqtable['tour'] = tour;
  rmqtable['levels'] = levels;
  rmqtable['represent'] = represent;

  return rmqtable

# Final solution for the general RMQ problem.
def rmq(i, j, rmqtable):
  if i > j: return -1

  # Find indices of i and j inside the Euler tour.
  represent = rmqtable['represent']
  itour = represent[i]
  jtour = represent[j]

  # It can be that i appears earlier in the tour than j. This is acceptable,
  # but we just need to swap the indices.
  if itour > jtour:
    (itour, jtour) = (jtour, itour)

  # Perform the RMQ on the Euler tour.
  maxind = rmq3(itour, jtour, rmqtable)

  # We got the maximal position within the Euler tour. Convert it back to
  # the position within the input array.
  tour = rmqtable['tour']
  return tour[maxind]
