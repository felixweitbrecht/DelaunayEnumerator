package delaunayKD.misc;

import java.util.ArrayList;

import delaunayKD.AllSimplicesFinder;
import delaunayKD.geometry.AbstractSimplex;
import delaunayKD.geometry.Point;
import delaunayKD.geometry.Simplex;

public class DoSomething {

	public static void main(String[] args) {
		// some code examples

		System.out.println("Code with 3D examples...");
		// generate some points (point at points.get(i) must have "time stamp"
		// set to i, i.e. points.get(i).i==i for all i)
		ArrayList<Point> points = UtilityMethods.generatePointsSphere(1 << 10);
		// compute all Delaunay simplices and facets over all subsequences of
		// the point set
		ArrayList<AbstractSimplex> result = AllSimplicesFinder.findAllSimplices(points);

		// do some things with the result
		countSimplicesAndFacets(result);
		identifySimplicesAndFacetsOfSubsequence(result, 123, 456);

		System.out.println("\n" + "Code with 2D examples...");
		AllSimplicesFinder.DIM = 2; // change dimension

		// run another computation with 2D point set
		ArrayList<Point> points2D = UtilityMethods.generatePointsSphere(1 << 12);
		AllSimplicesFinder.findAllSimplices(points2D);

		System.out.println("\n" + "Code with benchmarks...");
		System.out.println("Data will be printed with the following columns:" + "\n\t"
				+ "dimension, points inserted so far, simplices and facets found so far, time in ms spent so far");
		// run some benchmarks
		System.out.println("4D moment curve benchmark:");
		Benchmarks.main(new String[] { "1000", "4", "1", "moment" });

		System.out.println("5D unit hypercube benchmark:");
		Benchmarks.main(new String[] { "1000", "5", "1", "unit" });
	}

	private static void countSimplicesAndFacets(ArrayList<AbstractSimplex> simplices) {
		int simplexCount = 0;
		int facetCount = 0;
		for (AbstractSimplex s : simplices) {
			if (s instanceof Simplex) { // s is a proper simplex
				simplexCount++;
			} else { // s is a facet
				facetCount++;
			}
		}
		System.out.println("\t" + "The result contains " + simplexCount + " simplices and " + facetCount + " facets.");
	}

	// identify which simplices and facets are part of the Delaunay
	// triangulation of the subsequence [startIndex, endIndex].
	// startIndex and endIndex are inclusive
	private static void identifySimplicesAndFacetsOfSubsequence(ArrayList<AbstractSimplex> simplices, int startIndex,
			int endIndex) {
		// s is part of that triangulation iff the following conditions are met:
		// 1. all vertices of s are part of the subsequence
		// 2. no points of the subsequence are in the circumsphere of s, i.e.
		// the last previous killer point and the first subsequent killer point
		// are outside the subsequence.

		int count = 0;
		for (AbstractSimplex s : simplices) {
			if (s.minPoint().i >= startIndex && s.maxPoint().i <= endIndex && s.lastPreviousKillerIndex < startIndex
					&& s.firstSubsequentKillerIndex > endIndex) {
				count++;
			}
		}
		System.out.println("\t" + "The Delaunay triangulation of the subsequence [" + startIndex + ", " + endIndex
				+ "] contains " + count + " simplices and facets.");
	}

}
