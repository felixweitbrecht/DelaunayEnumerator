package delaunayKD.misc;

import static delaunayKD.AllSimplicesFinder.DIM;

import java.util.ArrayList;
import java.util.Random;
import java.util.Stack;

import delaunayKD.AllSimplicesFinder;
import delaunayKD.geometry.AbstractSimplex;
import delaunayKD.geometry.Point;
import delaunayKD.triangulator.HoleTriangulator;
import delaunayKD.triangulator.IncrementalTriangulator;
import delaunayKD.triangulator.Star;

public class Benchmarks {

	public static void main(String[] args) {
		// args: n, d, repetitions, [unit, moment]
		if (args.length < 4) {
			System.out.println("need more args: n, d, repetitions, [unit, moment]");
			System.exit(0);
		}

		int pointCount = Integer.parseInt(args[0]);
		int dim = Integer.parseInt(args[1]);
		int repetitions = Integer.parseInt(args[2]);
		boolean unitCube = args[3].equals("unit");
		AllSimplicesFinder.DIM = dim;

//		System.out.println("Benchmarking with " + args[3]);
		for (int iter = 0; iter < repetitions; iter++) {
			ArrayList<Point> points = unitCube ? UtilityMethods.generatePoints(pointCount)
					: generatePointsMomentCurve(pointCount);
			doBenchmarks(points, pointCount, dim);
		}
	}

	private static void doBenchmarks(ArrayList<Point> points, int n, int d) {
		long start = System.currentTimeMillis();

		IncrementalTriangulator incTriangulator = new IncrementalTriangulator();
		ArrayList<AbstractSimplex> allSimplices = new ArrayList<AbstractSimplex>();
		// stack of simplices that need to trigger updates
		Stack<AbstractSimplex> simplexStack = new Stack<AbstractSimplex>();

		for (int pIdx = 0; pIdx < points.size(); pIdx++) {
			if (pIdx % 100 == 0) {
//				System.out.println("inserting point " + pIdx);
			}
			Point pNew = points.get(pIdx);
			new HoleTriangulator(pNew); // hole triangulator for new point
			new Star(pNew); // star for new point

			// insert point into incremental construction (row 0)
			AbstractSimplex loc = pIdx >= DIM ? AllSimplicesFinder.locate(pNew, points.get(pIdx - DIM).star) : null;
			ArrayList<AbstractSimplex> incrementalNewSimplices = incTriangulator.addPoint(pNew, loc);
			simplexStack.addAll(incrementalNewSimplices);
			// work off stack, trigger updates for hole triangulations (rows >0)
			while (!simplexStack.isEmpty()) {
				AbstractSimplex simplex = simplexStack.pop();
				allSimplices.add(simplex);
				// trigger star/hole triangulation update
				ArrayList<AbstractSimplex> moreSimplices = simplex.minPoint().star.registerSimplex(simplex, pNew);
				simplexStack.addAll(moreSimplices);
			}

			if (pIdx % 100 == 99) {
				// remember stats
				long currTime = System.currentTimeMillis();
				long duration = currTime - start;
				// print stats: dimension, points inserted so far, simplices and
				// facets found so far, time in ms spent so far
				String statLine = d + "\t" + (pIdx + 1) + "\t" + allSimplices.size() + "\t" + duration;
				System.out.println("\t" + statLine);

				// adjust timing so time spent on printing stats is omitted
				long timeLost = System.currentTimeMillis() - currTime;
				start += timeLost;
			}
		}

//		System.out.println("creating simplices finished, got " + allSimplices.size() + " simplices/facets.");
	}

	// noisy moment curve to avoid precision errors
	private static ArrayList<Point> generatePointsMomentCurve(int count) {
		Random random = new Random(42L);
		ArrayList<Point> points = new ArrayList<Point>(count);
		for (int i = 0; i < count; i++) {
			double[] vals = new double[DIM];
			double val = i;
			for (int dimIdx = 0; dimIdx < DIM; dimIdx++) {
				vals[dimIdx] = val + 10 * random.nextDouble();
				val *= i;
			}
			points.add(new Point(vals, i));
		}
		return points;
	}

}
