package com.github.micycle1.flatbush;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.IntPredicate;

public class Flatbush {

	private static final int INITIAL_SEARCH_CAPACITY = 16;

	private final int numItems;
	private final int nodeSize;
	private final int[] levelBounds; // boundaries in the flat boxes array (in elements, not nodes)
	private final double[] boxes; // [minX,minY,maxX,maxY,...] for every node
	private final int[] indices; // maps node-id -> child node-element-index or leaf item-id
	private int pos; // next write position in boxes[]
	private double minX, minY, maxX, maxY;
	private final FlatQueue queue = new FlatQueue(); // for k-NN

	/**
	 * Create a Flatbush index that will hold {@code numItems} rectangles with
	 * default nodeSize=16.
	 */
	public Flatbush(final int numItems) {
		this(numItems, 16);
	}

	/**
	 * Create a Flatbush index that will hold {@code numItems} rectangles and pack
	 * {@code nodeSize} entries per node.
	 */
	public Flatbush(final int numItems, final int nodeSize) {
		if (numItems <= 0) {
			throw new IllegalArgumentException("Unexpected numItems value: " + numItems);
		}
		this.numItems = numItems;
		// clamp nodeSize to [2..65535]
		this.nodeSize = Math.min(Math.max(nodeSize, 2), 65535);

		// calculate total number of nodes and level boundaries
		int n = numItems;
		int numNodes = n;
		final List<Integer> lb = new ArrayList<>();
		lb.add(n * 4); // leaf-level boundary in box-elements
		do {
			n = (n + this.nodeSize - 1) / this.nodeSize;
			numNodes += n;
			lb.add(numNodes * 4); // cumulative element‐offset of next level
		} while (n != 1);

		// commit levelBounds[]
		this.levelBounds = new int[lb.size()];
		for (int i = 0; i < lb.size(); i++) {
			this.levelBounds[i] = lb.get(i);
		}

		// now allocate flat storage for numNodes * 4 doubles:
		this.boxes = new double[numNodes * 4];
		this.indices = new int[numNodes];
		this.pos = 0;
		this.minX = Double.POSITIVE_INFINITY;
		this.minY = Double.POSITIVE_INFINITY;
		this.maxX = Double.NEGATIVE_INFINITY;
		this.maxY = Double.NEGATIVE_INFINITY;
	}

	public double[] getBoxes() {
		return boxes;
	}

	public int[] getIndices() {
		return indices;
	}

	public int add(final double x, final double y) {
		return add(x, y, x, y);
	}

	/**
	 * Add a rectangle [minX,minY,maxX,maxY] to the index. Returns the leaf-id for
	 * this rectangle.
	 */
	public int add(final double minX, final double minY, final double maxX, final double maxY) {
		final int idx = pos >> 2;
		indices[idx] = idx;
		boxes[pos++] = minX;
		boxes[pos++] = minY;
		boxes[pos++] = maxX;
		boxes[pos++] = maxY;
		if (minX < this.minX) {
			this.minX = minX;
		}
		if (minY < this.minY) {
			this.minY = minY;
		}
		if (maxX > this.maxX) {
			this.maxX = maxX;
		}
		if (maxY > this.maxY) {
			this.maxY = maxY;
		}
		return idx;
	}

	/**
	 * Finalize and build the tree. Must be called after all add() calls.
	 */
	public void finish() {
		if ((pos >> 2) != numItems) {
			throw new IllegalStateException("Added " + (pos >> 2) + " items when expected " + numItems);
		}

		// single node?
		if (numItems <= nodeSize) {
			boxes[pos++] = minX;
			boxes[pos++] = minY;
			boxes[pos++] = maxX;
			boxes[pos++] = maxY;
			return;
		}

		// compute a Morton/Hilbert ordering on leaf centroids
		double w = (maxX - minX);
		double h = (maxY - minY);
		if (w == 0) {
			w = 1;
		}
		if (h == 0) {
			h = 1;
		}
		final long hilbertMax = (1 << 16) - 1;
		final long[] hilbertValues = new long[numItems];

		// map leaves into [0..hilbertMax]², compute Hilbert key
		int rp = 0;
		for (int i = 0; i < numItems; i++) {
			final double bx0 = boxes[rp++];
			final double by0 = boxes[rp++];
			final double bx1 = boxes[rp++];
			final double by1 = boxes[rp++];
			final double cx = (bx0 + bx1) * 0.5;
			final double cy = (by0 + by1) * 0.5;
			final int ix = (int) Math.floor(hilbertMax * (cx - minX) / w);
			final int iy = (int) Math.floor(hilbertMax * (cy - minY) / h);
			hilbertValues[i] = hilbert(ix, iy);
		}

		// sort leaves by Hilbert value (partial radix/quick trim by nodeSize)
		sort(hilbertValues, boxes, indices, 0, numItems - 1, nodeSize);

		// build parent nodes, bottom-up
		int readPos = 0;
		for (int lvl = 0; lvl < levelBounds.length - 1; lvl++) {
			final int end = levelBounds[lvl];
			while (readPos < end) {
				// start a new parent under construction
				final int childStart = readPos; // remember element‐offset
				double nminX = boxes[readPos];
				double nminY = boxes[readPos + 1];
				double nmaxX = boxes[readPos + 2];
				double nmaxY = boxes[readPos + 3];
				readPos += 4;
				// absorb up to nodeSize children
				for (int j = 1; j < nodeSize && readPos < end; j++) {
					final double bx0 = boxes[readPos++];
					final double by0 = boxes[readPos++];
					final double bx1 = boxes[readPos++];
					final double by1 = boxes[readPos++];
					if (bx0 < nminX) {
						nminX = bx0;
					}
					if (by0 < nminY) {
						nminY = by0;
					}
					if (bx1 > nmaxX) {
						nmaxX = bx1;
					}
					if (by1 > nmaxY) {
						nmaxY = by1;
					}
				}
				// write the new parent node
				final int parentId = pos >> 2;
				indices[parentId] = childStart; // element‐offset of first child
				boxes[pos++] = nminX;
				boxes[pos++] = nminY;
				boxes[pos++] = nmaxX;
				boxes[pos++] = nmaxY;
			}
		}
	}

	/**
	 * Search all leaf rectangles intersecting [minX,minY,maxX,maxY]. If
	 * filter==null, all hits are returned. Returns a packed int[] of leaf-ids.
	 */
	public int[] search(final double minX, final double minY, final double maxX, final double maxY, final IntPredicate filter) {
		if (pos != boxes.length) {
			throw new IllegalStateException("Data not indexed. Call finish().");
		}
		// a small stack of nodes to visit
		int[] stack = new int[16];
		int sp = 0;
		// results
		int[] result = new int[INITIAL_SEARCH_CAPACITY];
		int rc = 0;

		// start at root
		int nodeIndex = boxes.length - 4;
		while (nodeIndex >= 0) {
			// find child range
			final int end = Math.min(nodeIndex + nodeSize * 4, upperBound(nodeIndex, levelBounds));
			// scan children
			for (int p = nodeIndex; p < end; p += 4) {
				if ((maxX < boxes[p]) || (maxY < boxes[p + 1]) || (minX > boxes[p + 2]) || (minY > boxes[p + 3])) {
					continue;
				}
				final int id = indices[p >> 2];
				if (nodeIndex >= numItems * 4) {
					// internal node -> push its element-index
					if (sp == stack.length) {
						stack = Arrays.copyOf(stack, stack.length * 2);
					}
					stack[sp++] = id;
				} else {
					// leaf
					if (filter == null || filter.test(id)) {
						// grow array if needed
						if (rc == result.length) {
							result = Arrays.copyOf(result, result.length * 2);
						}
						result[rc++] = id;
					}
				}
			}
			// pop next
			nodeIndex = (sp > 0 ? stack[--sp] : -1);
		}
		return Arrays.copyOf(result, rc);
	}

	public int[] neighbors(final double x, final double y) {
		return neighbors(x, y, numItems, Double.POSITIVE_INFINITY, null);
	}

	public int[] neighbors(final double x, final double y, final int maxResults) {
		return neighbors(x, y, maxResults, Double.POSITIVE_INFINITY, null);
	}

	/**
	 * Find nearest neighbors to (x,y), up to maxResults and within maxDistance. If
	 * filter==null, all leaves qualify. Returns leaf-ids in ascending distance
	 * order.
	 */
	public int[] neighbors(final double x, final double y, final int maxResults, final double maxDistance, final IntPredicate filter) {
		if (pos != boxes.length) {
			throw new IllegalStateException("Data not indexed. Call finish().");
		}
		final double maxDist2 = maxDistance * maxDistance;
		final FlatQueue q = queue;
		q.clear();
		final int[] result = new int[Math.min(maxResults, numItems)];
		int rc = 0;

		int nodeIndex = boxes.length - 4;
		searchLoop: while (nodeIndex >= 0) {
			final int end = Math.min(nodeIndex + nodeSize * 4, upperBound(nodeIndex, levelBounds));
			// push all children that intersect the maxDistance circle
			for (int p = nodeIndex; p < end; p += 4) {
				final int id = indices[p >> 2];
				final double bx0 = boxes[p], by0 = boxes[p + 1];
				final double bx1 = boxes[p + 2], by1 = boxes[p + 3];
				final double dx = x < bx0 ? bx0 - x : (x > bx1 ? x - bx1 : 0);
				final double dy = y < by0 ? by0 - y : (y > by1 ? y - by1 : 0);
				final double dist2 = dx * dx + dy * dy;
				if (dist2 > maxDist2) {
					continue;
				}
				if (nodeIndex >= numItems * 4) {
					// internal: wrap id<<1 even
					q.push(id << 1, dist2);
				} else {
					// leaf: wrap id<<1|1 odd
					if (filter == null || filter.test(id)) {
						q.push((id << 1) | 1, dist2);
					}
				}
			}
			// drain all leaf items at the top of the queue
			while (q.length() > 0 && ((q.peek() & 1) == 1)) {
				final double d = q.peekValue();
				if (d > maxDist2) {
					break searchLoop;
				}
				result[rc++] = q.pop() >> 1;
				if (rc == maxResults) {
					break searchLoop;
				}
			}
			// next node
			if (q.length() > 0) {
				nodeIndex = q.pop() >> 1;
			} else {
				nodeIndex = -1;
			}
		}
		q.clear();
		return Arrays.copyOf(result, rc);
	}

	// ------------------------ static helpers ------------------------

	/** binary search in levelBounds for the first element > value */
	private static int upperBound(final int value, final int[] arr) {
		int lo = 0, hi = arr.length - 1;
		while (lo < hi) {
			final int mid = (lo + hi) >>> 1;
			if (arr[mid] > value) {
				hi = mid;
			} else {
				lo = mid + 1;
			}
		}
		return arr[lo];
	}

	/** partial quicksort by Hilbert value in blocks of nodeSize */
	private static void sort(final long[] values, final double[] boxes, final int[] indices, final int left, final int right, final int nodeSize) {
		if ((left / nodeSize) >= (right / nodeSize)) {
			return;
		}
		final long start = values[left];
		final long mid = values[(left + right) >>> 1];
		final long end = values[right];
		// median-of-three pivot
		long pivot = end;
		final long x = Math.max(start, mid);
		if (end > x) {
			pivot = x;
		} else if (x == start) {
			pivot = Math.max(mid, end);
		} else if (x == mid) {
			pivot = Math.max(start, end);
		}
		int i = left - 1, j = right + 1;
		while (true) {
			do {
				i++;
			} while (values[i] < pivot);
			do {
				j--;
			} while (values[j] > pivot);
			if (i >= j) {
				break;
			}
			swap(values, boxes, indices, i, j);
		}
		sort(values, boxes, indices, left, j, nodeSize);
		sort(values, boxes, indices, j + 1, right, nodeSize);
	}

	/** swap two leaf records in values, boxes and indices */
	private static void swap(final long[] values, final double[] boxes, final int[] indices, final int i, final int j) {
		final long tmpV = values[i];
		values[i] = values[j];
		values[j] = tmpV;
		final int bi = i << 2, bj = j << 2;
		final double a0 = boxes[bi], a1 = boxes[bi + 1], a2 = boxes[bi + 2], a3 = boxes[bi + 3];
		boxes[bi] = boxes[bj];
		boxes[bi + 1] = boxes[bj + 1];
		boxes[bi + 2] = boxes[bj + 2];
		boxes[bi + 3] = boxes[bj + 3];
		boxes[bj] = a0;
		boxes[bj + 1] = a1;
		boxes[bj + 2] = a2;
		boxes[bj + 3] = a3;
		final int t = indices[i];
		indices[i] = indices[j];
		indices[j] = t;
	}

	/**
	 * Fast 16-bit Hilbert curve from http://threadlocalmutex.com/ (public domain
	 * C++)
	 */
	private static final long hilbert(final int x, final int y) {
		int a = x ^ y;
		int b = 0xFFFF ^ a;
		int c = 0xFFFF ^ (x | y);
		int d = x & (y ^ 0xFFFF);

		int A = a | (b >>> 1);
		int B = (a >>> 1) ^ a;
		int C = ((c >>> 1) ^ (b & (d >>> 1))) ^ c;
		int D = ((a & (c >>> 1)) ^ (d >>> 1)) ^ d;

		a = A;
		b = B;
		c = C;
		d = D;
		A = ((a & (a >>> 2)) ^ (b & (b >>> 2)));
		B = ((a & (b >>> 2)) ^ (b & ((a ^ b) >>> 2)));
		C ^= ((a & (c >>> 2)) ^ (b & (d >>> 2)));
		D ^= ((b & (c >>> 2)) ^ ((a ^ b) & (d >>> 2)));

		a = A;
		b = B;
		c = C;
		d = D;
		A = ((a & (a >>> 4)) ^ (b & (b >>> 4)));
		B = ((a & (b >>> 4)) ^ (b & ((a ^ b) >>> 4)));
		C ^= ((a & (c >>> 4)) ^ (b & (d >>> 4)));
		D ^= ((b & (c >>> 4)) ^ ((a ^ b) & (d >>> 4)));

		a = A;
		b = B;
		c = C;
		d = D;
		C ^= ((a & (c >>> 8)) ^ (b & (d >>> 8)));
		D ^= ((b & (c >>> 8)) ^ ((a ^ b) & (d >>> 8)));

		a = C ^ (C >>> 1);
		b = D ^ (D >>> 1);

		int i0 = x ^ y;
		int i1 = b | (0xFFFF ^ (i0 | a));

		// interleave
		i0 = (i0 | (i0 << 8)) & 0x00FF00FF;
		i0 = (i0 | (i0 << 4)) & 0x0F0F0F0F;
		i0 = (i0 | (i0 << 2)) & 0x33333333;
		i0 = (i0 | (i0 << 1)) & 0x55555555;

		i1 = (i1 | (i1 << 8)) & 0x00FF00FF;
		i1 = (i1 | (i1 << 4)) & 0x0F0F0F0F;
		i1 = (i1 | (i1 << 2)) & 0x33333333;
		i1 = (i1 | (i1 << 1)) & 0x55555555;

		// combine to a 32-bit unsigned
		return (((long) i1 << 1) | i0) & 0xFFFFFFFFL;
	}

	/**
	 * Minimal binary heap of (int id, double priority) used by neighbors(...)
	 */
	private static class FlatQueue {
		private int[] ids = new int[16];
		private double[] vals = new double[16];
		private int size = 0;

		void clear() {
			size = 0;
		}

		int length() {
			return size;
		}

		void push(final int id, final double v) {
			if (size + 1 > ids.length) {
				final int newCap = ids.length * 2;
				ids = Arrays.copyOf(ids, newCap);
				vals = Arrays.copyOf(vals, newCap);
			}
			int pos = size++;
			while (pos > 0) {
				final int parent = (pos - 1) >>> 1;
				final double pv = vals[parent];
				if (v >= pv) {
					break;
				}
				ids[pos] = ids[parent];
				vals[pos] = pv;
				pos = parent;
			}
			ids[pos] = id;
			vals[pos] = v;
		}

		int peek() {
			return size > 0 ? ids[0] : -1;
		}

		double peekValue() {
			return size > 0 ? vals[0] : Double.NaN;
		}

		final int pop() {
			if (size == 0) {
				return -1;
			}
			final int ret = ids[0];
			final int lastIndex = --size;
			if (lastIndex > 0) {
				// take the last element up to root
				final int nodeId = ids[lastIndex];
				final double nodeVal = vals[lastIndex];
				ids[0] = nodeId;
				vals[0] = nodeVal;
				// bubble down
				int pos = 0;
				final int half = lastIndex >>> 1;
				while (pos < half) {
					final int left = (pos << 1) + 1;
					final int right = left + 1;
					int best = left;
					double bestVal = vals[left];
					if (right < lastIndex && vals[right] < bestVal) {
						best = right;
						bestVal = vals[right];
					}
					if (bestVal >= nodeVal) {
						break;
					}
					ids[pos] = ids[best];
					vals[pos] = bestVal;
					pos = best;
				}
				ids[pos] = nodeId;
				vals[pos] = nodeVal;
			}
			return ret;
		}
	}
}