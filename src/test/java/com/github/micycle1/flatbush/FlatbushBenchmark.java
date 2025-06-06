package com.github.micycle1.flatbush;

import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

/**
 * Benchmark for Flatbush.
 */
@State(Scope.Benchmark)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(time = 1)
@Measurement(time = 1)
public class FlatbushBenchmark {

	@Param({ "1000000" })
	private int N; // number of rectangles

	// match JS's 0..100 range
	private static final double REGION = 100.0;

	// 4*N array of [x0,y0,x1,y1, ... ]
	private double[] rects;
	private Flatbush bush;

	// three sets of 1000 query windows (10%, 1%, 0.01% of the area)
	private double[][] search10, search1, search001;

	// use the same rect‐minima as neighbor‐query centers
	private double[] center100X, center100Y;
	private double[] center1X, center1Y;

	@Setup(Level.Trial)
	public void setup() {
		final Random rnd = new Random(123);

		// 1) generate N random rectangles in [0..REGION]² with
		// width and height up to 1 (i.e. 1% of 100)
		rects = new double[N * 4];
		for (int i = 0; i < N; i++) {
			final double x0 = rnd.nextDouble() * REGION;
			final double y0 = rnd.nextDouble() * REGION;
			final double w = rnd.nextDouble() * 1.0; // ≤1
			final double h = rnd.nextDouble() * 1.0; // ≤1
			rects[4 * i] = x0;
			rects[4 * i + 1] = y0;
			rects[4 * i + 2] = x0 + w;
			rects[4 * i + 3] = y0 + h;
		}

		// 2) build the Flatbush index
		bush = new Flatbush(N);
		for (int i = 0; i < N; i++) {
			final int j = 4 * i;
			bush.add(rects[j], rects[j + 1], rects[j + 2], rects[j + 3]);
		}
		bush.finish();

		// 3) prepare 1000 random search windows
		search10 = new double[1000][4];
		search1 = new double[1000][4];
		search001 = new double[1000][4];
		final double side10 = REGION * Math.sqrt(0.10); // ≈ 31.62
		final double side1 = REGION * Math.sqrt(0.01); // = 10
		final double side001 = REGION * Math.sqrt(0.0001);// = 1

		search10 = new double[1000][4];
		search1 = new double[1000][4];
		search001 = new double[1000][4];

		for (int i = 0; i < 1000; i++) {
			// exactly like JS: origin in [0 .. REGION - side]
			double x0 = rnd.nextDouble() * (REGION - side10);
			double y0 = rnd.nextDouble() * (REGION - side10);
			double w = rnd.nextDouble() * side10; // random width ≤ side10
			double h = rnd.nextDouble() * side10; // random height ≤ side10
			search10[i][0] = x0;
			search10[i][1] = y0;
			search10[i][2] = x0 + w;
			search10[i][3] = y0 + h;

			// same for 1% windows
			x0 = rnd.nextDouble() * (REGION - side1);
			y0 = rnd.nextDouble() * (REGION - side1);
			w = rnd.nextDouble() * side1;
			h = rnd.nextDouble() * side1;
			search1[i][0] = x0;
			search1[i][1] = y0;
			search1[i][2] = x0 + w;
			search1[i][3] = y0 + h;

			// and for 0.01%
			x0 = rnd.nextDouble() * (REGION - side001);
			y0 = rnd.nextDouble() * (REGION - side001);
			w = rnd.nextDouble() * side001;
			h = rnd.nextDouble() * side001;
			search001[i][0] = x0;
			search001[i][1] = y0;
			search001[i][2] = x0 + w;
			search001[i][3] = y0 + h;
		}

		// 4) pick the *same* N rectangles as neighbor-query centers
		// first 1000 for k=100, next 100_000 for k=1
		center100X = new double[1000];
		center100Y = new double[1000];
		for (int i = 0; i < 1000; i++) {
			center100X[i] = rects[4 * i];
			center100Y[i] = rects[4 * i + 1];
		}

		center1X = new double[100000];
		center1Y = new double[100000];
		for (int i = 0; i < 100000; i++) {
			center1X[i] = rects[4 * i];
			center1Y[i] = rects[4 * i + 1];
		}
	}

	@Benchmark
	public void buildIndex(final Blackhole bh) {
		final Flatbush idx = new Flatbush(N);
		for (int i = 0; i < N; i++) {
			final int j = 4 * i;
			idx.add(rects[j], rects[j + 1], rects[j + 2], rects[j + 3]);
		}
		idx.finish();
		bh.consume(idx);
	}

	@Benchmark
	public void search10pct(final Blackhole bh) {
		for (int i = 0; i < 1000; i++) {
			final double[] r = search10[i];
			bh.consume(bush.search(r[0], r[1], r[2], r[3], null));
		}
	}

	@Benchmark
	public void search1pct(final Blackhole bh) {
		for (int i = 0; i < 1000; i++) {
			final double[] r = search1[i];
			bh.consume(bush.search(r[0], r[1], r[2], r[3], null));
		}
	}

	@Benchmark
	public void search001pct(final Blackhole bh) {
		for (int i = 0; i < 1000; i++) {
			final double[] r = search001[i];
			bh.consume(bush.search(r[0], r[1], r[2], r[3], null));
		}
	}

	@Benchmark
	public void neighbors100(final Blackhole bh) {
		// exactly 1000 queries of k=100
		for (int i = 0; i < 1000; i++) {
			bh.consume(bush.neighbors(center100X[i], center100Y[i], 100, Double.POSITIVE_INFINITY, null));
		}
	}

	@Benchmark
	public void neighborsAll(final Blackhole bh) {
		// single query requesting ALL N neighbors
		bh.consume(bush.neighbors(REGION / 2, REGION / 2, N, Double.POSITIVE_INFINITY, null));
	}

	@Benchmark
	public void neighbors1_100k(final Blackhole bh) {
		// 100k queries of k=1
		for (int i = 0; i < 100000; i++) {
			bh.consume(bush.neighbors(center1X[i], center1Y[i], 1, Double.POSITIVE_INFINITY, null));
		}
	}
}