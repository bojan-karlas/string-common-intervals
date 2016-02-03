#!/usr/bin/python -tt

# Implementation of various functions that are taken from the paper by
# Gilles Didier et al. entitled 'Character sets of strings'

# Algorithm for Problem 1:

# Takes a string and a character set, and constructs a list of maximal intervals
# such that the character set of each interval is equal to the given character
# set.
def maximal_locations(string, charset):
  # Check input parameters.
  if len(charset) < 1:
    raise ValueError('The charset argument cannot be empty.')

  locations = []
  register = {}
  count = 0
  start = -1

  # We use a hash table for faster lookup.
  charsethash = {}
  for i in range(len(charset)):
    charsethash[charset[i]] = None

  for i in range(len(string)):
    c = string[i]

    # Check if we have completed a 'location'
    if count >= len(charset) and c not in register:
      locations.append((start, i))
      register = {}
      count = 0
      start = -1

    # Check if the character is in charset. If not, reset.
    if c not in charsethash:
      register = {}
      count = 0
      start = -1
      continue
    elif start == -1:
      start = i

    # Check if current char has been encountered in current 'location'.
    if c not in register:
      register[c] = None
      count += 1

  # Check if the current location is tailing the string.
  if count >= len(charset):
      locations.append((start, len(string)))

  return locations

# Algorithm 1 for Problem 2

# Takes the input string and computes the list of chars in the string, ordered
# by their first occurance.
def compute_occurlist(string):
  occurset = {}
  occurlist = []
  for c in string:
    if c not in occurset:
      occurset[c] = None
      occurlist.append(c)

  return occurlist

# Takes the list of character occurences (L) and the alphabet (Sigma)
# and computes the rank of each character which is its index in the list L if
# it appears in the list, or +infinity otherwise.
def compute_ranks(occurlist, alphabet):
  ranks = {};

  for i in range(len(occurlist)): ranks[occurlist[i]] = i

  for c in alphabet:
    if c not in ranks:
      ranks[c] = float('inf')

  return ranks

# Computes the table of rank intervals described in Definition 7 of the paper.
# I.e. for every position k in the string S, a rank interval [a, b] satisfies:
#  (1) For every l in [a, b]: Rank[S[l]] <= Rank[S[k]]
#  (2) Rank[S[a-1]] > Rank[S[k]] AND Rank[S[b+1]] > Rank[S[k]]
def compute_rank_intervals(string, ranks):
  # Output is a hash table
  # with {key = position; value = (left_bound, right_bound)}.
  rankint = {};

  # Stack contains tuples: (rank, position, left_bound)
  # Init stack with string edge placeholder of infinite rank.
  stack = [(float('inf'), -1, -1)]

  for i in range(len(string)):
    c = string[i]

    # Pop from the stack all the positions of rank smaller than the current one.
    while stack[-1][0] < ranks[c]:
      rankint[stack[-1][1]] = (stack[-1][2], i)
      stack.pop()

    # Push the current position to the stack.
    if stack[-1][0] == ranks[c]:
      stack.append((ranks[c], i, stack[-1][2]))
    else:
      stack.append((ranks[c], i, stack[-1][1] + 1))

  # Flush the stack.
  while len(stack) > 1:
    rankint[stack[-1][1]] = (stack[-1][2], len(string))
    stack.pop()

  return rankint

# Build a rank lookup dictionary where the ranks are keys and values are
# lists of positions with the given rank. Also build a rank_table for each
# char in the input string to compute the rank distances.
def build_rank_table_and_lookup(string, ranks):
  rank_lookup = {}
  rank_table = [-1 for x in range(len(string))]

  for i in range(len(string)):
    rank = ranks[string[i]]
    rank_table[i] = rank

    if rank in rank_lookup:
      rank_lookup[rank].append(i)
    else:
      rank_lookup[rank] = [i]

  return (rank_table, rank_lookup)

# Utility method for computing the rank distance between two positions a and b.
#  - a, b        - Positions between which we want to calculate the distance.
#  - rmqtable    - The lookup table which allows us to perform the Range
#                  Maximum Query in O(1) time.
#  - rank_table  - Array of size the same as the input string S where each
#                  position i holds the rank of S[i].
#  - absolute    - If set to false and a > b then the function returns Inf.
#                  Otherwise, returns the absolute distance without paying
#                  attention to the inequality relation between a and b.
def compute_rank_distance(a, b, rmqtable, rank_table, absolute = True):
  if absolute == True and a > b: (a, b) = (b, a)
  maxind = rangemaxq.rmq(a, b, rmqtable)
  if maxind >= 0:
    return rank_table[maxind]
  else:
    return float('inf')

# Computes the table of rank successors described in Definition 10 of the paper.
# I.e. for every position k in the string S, is taken from the set P of
# positions in S such that their rank is equal to Rank[S[k]]+1. Find p1 and p2
# to be members of P such that p1 < k and p2 > k (if either p1 or p2 don't exist
# disregard it and set the other one to be the rank successor). Otherwise,
# the rank successor is p1 if rank_dist(p1, k) <= rank_dist(k, p2) is true
# and p2 if not.
# The rank distance (rank_dist(a, b)) is defined as the maximal rank of any
# character in the interval S[a, b].
import rangemaxq
def compute_rank_successors(string, ranks):

  # Build rank table and rank lookup.
  (rank_table, rank_lookup) = build_rank_table_and_lookup(string, ranks)

  # Build the successors table.
  succ = [None for x in range(len(string))]
  rmqtable = rangemaxq.rmq_pre(rank_table)

  for rank in rank_lookup:
    if (rank + 1) not in rank_lookup or rank == float('inf'): continue
    cur_rank = rank_lookup[rank]
    next_rank = rank_lookup[rank + 1]

    if len(next_rank) < 1:
      # The imput string must be complete.
      raise RuntimeError('Thre cannot be 0 chars with rank %d.' % rank)

    j1 = j2 = 0
    for i in range(len(cur_rank)):
      pos = cur_rank[i]
      while j2 < len(next_rank) and pos > next_rank[j2]:
        j1 = j2
        j2 += 1

      if j2 == len(next_rank): j2 = j1

      dist1 = compute_rank_distance(next_rank[j1], pos,
        rmqtable, rank_table, absolute = False)

      dist2 = compute_rank_distance(pos, next_rank[j2],
        rmqtable, rank_table, absolute = False)

      if dist1 <= dist2:
        succ[pos] = next_rank[j1]
      else:
        succ[pos] = next_rank[j2]

  return succ

# Takes a string and, for every possible character set that is a subset of the
# alphabet of the given string, computes all maximal locations of that
# character set.
def maximal_substrings(string):
  alphabet = list(set(string))

  # List of that collects the intervals of maximal substrings.
  intervals = []

  for i in range(len(string)):
    # Determine the maximal substring starting at i.
    j = len(string)
    if i > 0:
      j = i
      while j < len(string) and string[i-1] != string[j]: j += 1
    substring = string[i:j]

    # Initialize all required data structures.
    occurlist = compute_occurlist(substring)
    ranks = compute_ranks(occurlist, alphabet)
    rank_int = compute_rank_intervals(string, ranks)
    rank_succ = compute_rank_successors(string, ranks)
    (rank_table, rank_lookup) = build_rank_table_and_lookup(string, ranks)
    rmqtable = rangemaxq.rmq_pre(rank_table)

    # Initalize list with paths of rank 0 in increasing order.
    # Each element contains a tuple (position, (left, right))
    # Where position is the last position of the path, while left and right
    # are the bounds of the minimal interval that contains the path.
    path_list = []
    if 0 in rank_lookup:
      for pos in rank_lookup[0]:
        path_list.append((pos, (pos, pos + 1)))

    # Function that we will use to test if interval1 is contained
    # within interval2.
    def is_subset(interval1, interval2):
      return interval1[0] >= interval2[0] and interval1[1] <= interval2[1]

    while len(path_list) > 0:
      # ------------------------------------------------------------------------
      # Test if some of the paths in the list can result in locations to output.
      # ------------------------------------------------------------------------

      # Find first path with bounds contained within the rank interval of its
      # last position.
      k = 0
      while k < len(path_list) and not is_subset(path_list[k][1], rank_int[path_list[k][0]]):
        k += 1
      first = None if k == len(path_list) else path_list[k]

      if first != None and rank_int[first[0]][0] >= i:
        intervals.append(rank_int[first[0]])

        prev_rank_int = rank_int[first[0]]
        for cur in path_list[k:]:
          cur_rank_int = rank_int[cur[0]]
          if is_subset(cur[1], cur_rank_int) \
            and cur_rank_int != prev_rank_int:
            intervals.append(cur_rank_int)
            prev_rank_int = cur_rank_int

      # ------------------------------------------------------------------------
      # Compute the next level.
      # ------------------------------------------------------------------------

      # Build a list of paths that should be deleted, either because they don't
      # have a successor, or because their last position doesn't have the
      # smallest rank distance to the successor among all the last positions
      # that share the same successor.
      batch_succ = None
      batch_begin = 0
      nearest_in_batch = None
      min_dist = float('inf')
      mark_delete = []

      for k in range(len(path_list)):
        cur_path = path_list[k]
        cur_dist = float('inf') if rank_succ[cur_path[0]] == None \
          else compute_rank_distance(cur_path[0], rank_succ[cur_path[0]], rmqtable, rank_table)

        # If we've reached the end of a batch of paths with the same successor
        # or we reached the end of the list, we take the one that's nearest
        # to the successor and discard the rest.
        if batch_succ != rank_succ[cur_path[0]]:

          # Go through the batch and mark for deletion all except the nearest.
          for j in range(batch_begin, k):
            if rank_succ[path_list[j][0]] == None or j != nearest_in_batch:
              mark_delete.append(j)

          # Begin a new batch.
          batch_succ = rank_succ[cur_path[0]]
          batch_begin = nearest_in_batch = k
          min_dist = cur_dist

        elif cur_dist < min_dist:
          # If we are inside a batch, update the minimum.
          min_dist = cur_dist
          nearest_in_batch = k

      # Go through the batch and mark for deletion all except the nearest.
      for j in range(batch_begin, k + 1):
        if rank_succ[path_list[j][0]] == None or j != nearest_in_batch:
          mark_delete.append(j)

      # Delete all marked paths.
      for k in reversed(range(len(mark_delete))):
        path_list.pop(mark_delete[k])

      # Extend all the remaining paths with their successors and update their
      # bounds.
      for k in range(len(path_list)):
        cur_path = path_list[k]
        succ = rank_succ[cur_path[0]]
        left_bound = min(cur_path[1][0], succ)
        right_bound = max(cur_path[1][1], succ + 1)
        path_list[k] = (succ, (left_bound, right_bound))

  return intervals
