# Flatbush (Java Port)

A high-performance, static R-tree spatial index for 2D rectangles, ported from the popular JavaScript [Flatbush](https://github.com/mourner/flatbush) library.  
This Java version is a close line-by-line port for maximum speed and minimal footprint. No external dependencies.

## Features

- Bulk insertion of axis-aligned rectangles (or points).  
- Hilbert-curve sorting for optimal packing.  
- O(n) index construction time.  
- Very fast rectangle-intersection queries.  
- Efficient k-nearest-neighbors queries with optional distance limit and filtering.  
- Primitive arrays only; zero GC churn during queries.  
- All in a single small class (`Flatbush.java`).  
