package com.github.micycle1.flatbush;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.Arrays;
import java.util.function.IntPredicate;

import org.junit.jupiter.api.Test;

class FlatbushTest {

	static final double[] data = { 8, 62, 11, 66, 57, 17, 57, 19, 76, 26, 79, 29, 36, 56, 38, 56, 92, 77, 96, 80, 87, 70, 90, 74, 43, 41, 47, 43, 0, 58, 2, 62,
			76, 86, 80, 89, 27, 13, 27, 15, 71, 63, 75, 67, 25, 2, 27, 2, 87, 6, 88, 6, 22, 90, 23, 93, 22, 89, 22, 93, 57, 11, 61, 13, 61, 55, 63, 56, 17, 85,
			21, 87, 33, 43, 37, 43, 6, 1, 7, 3, 80, 87, 80, 87, 23, 50, 26, 52, 58, 89, 58, 89, 12, 30, 15, 34, 32, 58, 36, 61, 41, 84, 44, 87, 44, 18, 44, 19,
			13, 63, 15, 67, 52, 70, 54, 74, 57, 59, 58, 59, 17, 90, 20, 92, 48, 53, 52, 56, 92, 68, 92, 72, 26, 52, 30, 52, 56, 23, 57, 26, 88, 48, 88, 48, 66,
			13, 67, 15, 7, 82, 8, 86, 46, 68, 50, 68, 37, 33, 38, 36, 6, 15, 8, 18, 85, 36, 89, 38, 82, 45, 84, 48, 12, 2, 16, 3, 26, 15, 26, 16, 55, 23, 59,
			26, 76, 37, 79, 39, 86, 74, 90, 77, 16, 75, 18, 78, 44, 18, 45, 21, 52, 67, 54, 71, 59, 78, 62, 78, 24, 5, 24, 8, 64, 80, 64, 83, 66, 55, 70, 55, 0,
			17, 2, 19, 15, 71, 18, 74, 87, 57, 87, 59, 6, 34, 7, 37, 34, 30, 37, 32, 51, 19, 53, 19, 72, 51, 73, 55, 29, 45, 30, 45, 94, 94, 96, 95, 7, 22, 11,
			24, 86, 45, 87, 48, 33, 62, 34, 65, 18, 10, 21, 14, 64, 66, 67, 67, 64, 25, 65, 28, 27, 4, 31, 6, 84, 4, 85, 5, 48, 80, 50, 81, 1, 61, 3, 61, 71,
			89, 74, 92, 40, 42, 43, 43, 27, 64, 28, 66, 46, 26, 50, 26, 53, 83, 57, 87, 14, 75, 15, 79, 31, 45, 34, 45, 89, 84, 92, 88, 84, 51, 85, 53, 67, 87,
			67, 89, 39, 26, 43, 27, 47, 61, 47, 63, 23, 49, 25, 53, 12, 3, 14, 5, 16, 50, 19, 53, 63, 80, 64, 84, 22, 63, 22, 64, 26, 66, 29, 66, 2, 15, 3, 15,
			74, 77, 77, 79, 64, 11, 68, 11, 38, 4, 39, 8, 83, 73, 87, 77, 85, 52, 89, 56, 74, 60, 76, 63, 62, 66, 65, 67 };

	private Flatbush createIndex() {
		int n = data.length / 4;
		Flatbush idx = new Flatbush(n);
		for (int i = 0; i < data.length; i += 4) {
			idx.add(data[i], data[i + 1], data[i + 2], data[i + 3]);
		}
		idx.finish();
		return idx;
	}

	private Flatbush createSmallIndex(int numItems, int nodeSize) {
		Flatbush idx = new Flatbush(numItems, nodeSize);
		for (int i = 0; i < 4 * numItems; i += 4) {
			idx.add(data[i], data[i + 1], data[i + 2], data[i + 3]);
		}
		idx.finish();
		return idx;
	}

	@Test
	void testIndexesABunchOfRectangles() {
		Flatbush idx = createIndex();
		double[] boxes = idx.getBoxes();
		int[] inds = idx.getIndices();

		assertEquals(540, boxes.length + inds.length);

		int len = boxes.length;
		double[] rootBox = Arrays.copyOfRange(boxes, len - 4, len);
		assertArrayEquals(new double[] { 0, 1, 96, 95 }, rootBox, 0.0);

		// the very last node’s child‐pointer must be 400
		assertEquals(400, inds[len / 4 - 1]);
	}

	@Test
	void testSkipsSortingWhenNumItemsLessThanNodeSize() {
		int numItems = 14, nodeSize = 16;
		Flatbush idx = createSmallIndex(numItems, nodeSize);

		int[] inds = idx.getIndices();
		// expected = [0,1,2,...,13,0]
		int[] exp = new int[numItems + 1];
		for (int i = 0; i < numItems; i++) {
			exp[i] = i;
		}
		exp[numItems] = 0;
		assertArrayEquals(exp, inds);

		// check root box extent
		double rx0 = Double.POSITIVE_INFINITY, ry0 = Double.POSITIVE_INFINITY;
		double rx1 = Double.NEGATIVE_INFINITY, ry1 = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < 4 * numItems; i += 4) {
			rx0 = Math.min(rx0, data[i]);
			ry0 = Math.min(ry0, data[i + 1]);
			rx1 = Math.max(rx1, data[i + 2]);
			ry1 = Math.max(ry1, data[i + 3]);
		}
		double[] boxes = idx.getBoxes();
		double[] rootBox = Arrays.copyOfRange(boxes, boxes.length - 4, boxes.length);
		assertArrayEquals(new double[] { rx0, ry0, rx1, ry1 }, rootBox, 0.0);
		assertEquals((numItems + 1) * 4, boxes.length);
	}

	@Test
	void testPerformsBBoxSearch() {
		Flatbush idx = createIndex();
		int[] ids = idx.search(40, 40, 60, 60, null);

		// reconstruct rectangles from data
		int[] found = new int[ids.length * 4];
		int p = 0;
		for (int id : ids) {
			for (int k = 0; k < 4; k++) {
				found[p++] = (int) data[4 * id + k];
			}
		}
		int[] exp = { 57, 59, 58, 59, 48, 53, 52, 56, 40, 42, 43, 43, 43, 41, 47, 43 };
		Arrays.sort(found);
		Arrays.sort(exp);
		assertArrayEquals(exp, found);
	}

	@Test
	void testDefaultsToAddingAPoint() {
		Flatbush idx = new Flatbush(1);
		idx.add(10, 10); // uses the 2-arg overload
		idx.finish();
		assertArrayEquals(new int[] { 0 }, idx.search(0, 0, 20, 20, null));
	}

	@Test
	void testThrowsIfLessItemsThanSize() {
		Flatbush idx = new Flatbush(data.length / 4);
		assertThrows(IllegalStateException.class, idx::finish);
	}

	@Test
	void testThrowsIfSearchBeforeIndexing() {
		Flatbush idx = new Flatbush(data.length / 4);
		assertThrows(IllegalStateException.class, () -> idx.search(0, 0, 20, 20, null));
	}

	@Test
	void testConstructorZeroItemsThrows() {
		assertThrows(IllegalArgumentException.class, () -> new Flatbush(0));
	}

	@Test
	void testKNearestNeighborsK3() {
		Flatbush idx = createIndex();
		int[] ids = idx.neighbors(50, 50, 3, Double.POSITIVE_INFINITY, null);
		Arrays.sort(ids);
		int[] exp = { 6, 31, 75 };
		Arrays.sort(exp);
		assertArrayEquals(exp, ids);
	}

	@Test
	void testKNearestWithMaxDistance() {
		Flatbush idx = createIndex();
		int[] ids = idx.neighbors(50, 50, Integer.MAX_VALUE, 12, null);
		Arrays.sort(ids);
		int[] exp = { 6, 29, 31, 75, 85 };
		Arrays.sort(exp);
		assertArrayEquals(exp, ids);
	}

	@Test
	void testKNearestWithFilter() {
		Flatbush idx = createIndex();
		IntPredicate even = i -> (i & 1) == 0;
		int[] ids = idx.neighbors(50, 50, 6, Double.POSITIVE_INFINITY, even);
		Arrays.sort(ids);
		int[] exp = { 6, 16, 18, 24, 54, 80 };
		Arrays.sort(exp);
		assertArrayEquals(exp, ids);
	}

	@Test
	void testKNearestAllItems() {
		Flatbush idx = createIndex();
		int[] ids = idx.neighbors(50, 50);
		assertEquals(data.length / 4, ids.length);
	}

	@Test
	void testAddReturnsIncrementalIndices() {
		int count = 5;
		Flatbush idx = new Flatbush(count);
		int[] ids = new int[count];
		for (int i = 0; i < count; i++) {
			ids[i] = idx.add(data[4 * i], data[4 * i + 1], data[4 * i + 2], data[4 * i + 3]);
		}
		assertArrayEquals(new int[] { 0, 1, 2, 3, 4 }, ids);
	}

	@Test
	void testQuicksortOnInbalancedDatasetDoesNotBlowUp() {
		int n = 15_000;
		Flatbush idx = new Flatbush(2 * n);

		double[] a1 = linspace(0, 1000, n, true);
		double[] a2 = linspace(0, 1000, n, true);
		for (double x : a1) {
			idx.add(x, 0, x, 0);
		}
		for (double x : a2) {
			idx.add(x, 0, x, 0);
		}
		idx.finish();

		// should not throw
		assertDoesNotThrow(() -> idx.search(-100, -1, 15_000, 1, null));
	}

	// helper for the inbalanced data test
	private static double[] linspace(double start, double stop, int num, boolean endpoint) {
		int div = endpoint ? (num - 1) : num;
		double step = (stop - start) / div;
		double[] out = new double[num];
		for (int i = 0; i < num; i++) {
			out[i] = start + step * i;
		}
		return out;
	}
}