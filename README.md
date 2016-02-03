# String Common Intervals

Python modules that implement some algorithms for finding common intervals (i.e. maximal locations) in strings.

Currently, there are 2 modules:

* `rangemaxq.py` 

  Contains an implementation of the [Range Maximum Query (RMQ)](https://en.wikipedia.org/wiki/Range_minimum_query) problem. In short, the problem is to find the maximal (or minimal) value of a subarray, given the range of the subarray. The file follows the [paper](http://dl.acm.org/citation.cfm?id=690192) by Bender and Farach-Colton entitled 'The LCA Problem Revisited'. The module implementation follows the paper by giving 3 solutions with progressive improvement of complexity. The final fourth solution produces a O(n) data structure and is able to do a RMQ lookup in constant time. The preprocessing is done on the input array, and the queries are done by specifying the range of the subarray on which we perform the query. This module was implemented because it is used by `didier.py`.
  
  **Usage:**
  1. Call `rmq_pre#(values)` (where `#` is replaced by either 1, 2, 3 or nothing). The `values` argument represents the input array. The function returns a *helper data structure* `rmqtable` which is used in subsequent calls of the *query function*.
  2. Call `rmq#(i, j, rmqtable)` which is the query function. We pass it the helper data structure and it gives us the index of the maximal element in the [i,j] range of the input array.

* `didier.py`
  
  Implements two algorithms taken from a [paper](http://www.sciencedirect.com/science/article/pii/S1570866706000505) by Didier et al. entitled 'Character sets of strings'. Given a string, the module contains implementations of two possible queries that can be performed on that string:
  1. Given a string S and a character set C, find all maximal substrings of S that have C as their character set. Note that a substring is called maximal if it cannot be extended on either side without extending its character set. This query is performed by invoking `maximal_locations(string, charset)`.
  2. Given a string S, for each character set C that occurs at least once in S, all maximal substrings of S that have C as their character set. This query is performed by invoking `maximal_substrings(string)`.
