package delaunayKD;

import java.util.ArrayList;
import java.util.Stack;

import delaunayKD.geometry.AbstractSimplex;
import delaunayKD.geometry.Facet;
import delaunayKD.geometry.Point;
import delaunayKD.geometry.Face;
import delaunayKD.triangulator.HoleTriangulator;
import delaunayKD.triangulator.IncrementalTriangulator;
import delaunayKD.triangulator.Star;

public class AllSimplicesFinder {
	public static int DIM = 3;

	// whether to do bookkeeping for computation of temporal alpha-shape
	public static boolean doAlphaBookkeeping = false;

	public static ArrayList<AbstractSimplex> findAllSimplices(ArrayList<Point> points) {
		IncrementalTriangulator incTriangulator = new IncrementalTriangulator();
		return findAllSimplices(points, incTriangulator);
	}

	// allow supplying incTriangulator so its data structures remain accessible
	public static ArrayList<AbstractSimplex> findAllSimplices(ArrayList<Point> points,
			IncrementalTriangulator incTriangulator) {
		ArrayList<AbstractSimplex> allSimplices = new ArrayList<AbstractSimplex>();
		// stack of simplices that need to be registered with stars
		Stack<AbstractSimplex> simplexStack = new Stack<AbstractSimplex>();

		for (int pIdx = 0; pIdx < points.size(); pIdx++) {
			if (pIdx % 1000 == 0) {
				System.out.println("inserting point " + pIdx);
			}

			Point pNew = points.get(pIdx);
			new HoleTriangulator(pNew); // hole triangulator for new point
			new Star(pNew); // star for new point

			// insert point into incremental construction (row 0)
			AbstractSimplex location = pIdx >= DIM ? locate(pNew, points.get(pIdx - DIM).star) : null;
			ArrayList<AbstractSimplex> incrementalNewSimplices = incTriangulator.addPoint(pNew, location);
			simplexStack.addAll(incrementalNewSimplices);
			// work off stack, register simplices with stars and trigger
			// updates for hole triangulations (rows >0)
			while (!simplexStack.isEmpty()) {
				AbstractSimplex simplex = simplexStack.pop();
				allSimplices.add(simplex);
				if (AllSimplicesFinder.doAlphaBookkeeping) {
					// store simplex with its faces' lists of known simplices
					for (Face f : simplex.faces) {
						f.knownSimplices.add(simplex);
					}
				}
				// trigger star/hole triangulation update
				ArrayList<AbstractSimplex> moreSimplices = simplex.minPoint().star.registerSimplex(simplex, pNew);
				simplexStack.addAll(moreSimplices);
			}
		}

		System.out.println("creating simplices finished, got " + allSimplices.size() + " simplices/facets.");
		return allSimplices;
	}

	// given the star of p_(new-DIM), locate pNew in the incremental
	// construction. returns a simplex of the incremental construction which
	// contains pNew in its circumsphere.
	public static AbstractSimplex locate(Point pNew, Star star) {
		AbstractSimplex destroyedSimplex = (star.faceLatest.facesPoint(pNew) ? star.faceLatest
				: star.faceLatest.r).simplex.original;
		// loop: based on a destroyed hole triangulation simplex, find a
		// destroyed simplex in a row above. once we reach the first row, i.e.
		// the incremental triangulation, we are finished.
		while (destroyedSimplex.lastPreviousKillerIndex != Integer.MIN_VALUE) {
			// find a destroyed simplex in row above

			// we are in a hole triangulation
			// BFS to find more destroyed simplices in hole triangulation
			Stack<AbstractSimplex> simplicesToExplore = new Stack<AbstractSimplex>();
			ArrayList<AbstractSimplex> destroyedSimplices = new ArrayList<AbstractSimplex>();
			simplicesToExplore.add(destroyedSimplex);
			while (!simplicesToExplore.isEmpty()) {
				AbstractSimplex simplex = simplicesToExplore.pop();
				if (simplex != null && !simplex.isMarked() && simplex.containsPointInCircumsphere(pNew)) {
					simplex.mark();
					destroyedSimplices.add(simplex);
					for (Face face : simplex.faces) {
						simplicesToExplore.add(face.r.simplex);
					}
					if (simplex instanceof Facet) {
						for (Face face : simplex.faces[0].hNeighbors) {
							simplicesToExplore.add(face.simplex);
						}
					}
				}
			}
			AbstractSimplex.unmarkAll();

			// find a destroyed simplex in the corresponding star
			AbstractSimplex destroyedSimplexStar = null;
			// try finding an old boundary face (that would give us a destroyed
			// star simplex)
			for (int i = 0; i < destroyedSimplices.size() && destroyedSimplexStar == null; i++) {
				AbstractSimplex simplex = destroyedSimplices.get(i);
				for (Face face : simplex.faces) {
					if (face.r.simplex == null && face.r.faceBoundary.r.simplex.containsPointInCircumsphere(pNew)) {
						// this face is an old boundary face
						destroyedSimplexStar = face.r.faceBoundary.r.simplex;
						break;
					}
				}
			}
			// no success? look at facets and see if they have a corresponding
			// destroyed facet in the star
			for (int i = 0; i < destroyedSimplices.size() && destroyedSimplexStar == null; i++) {
				AbstractSimplex simplex = destroyedSimplices.get(i);
				if (simplex instanceof Facet) {
					Face face = simplex.faces[0];
					Point[] facePoints = face.points();
					for (int faceIdx = 0; faceIdx < DIM; faceIdx++) {
						Point pOpposite = facePoints[faceIdx];
						Face faceNeighbor = face.hNeighbors[faceIdx];
						if (faceNeighbor.simplex == null) {
							// faceNeighbor is a boundary face
							// look at adjacent facet in the star
							AbstractSimplex starFacet = faceNeighbor.faceBoundary.hNeighbors[faceNeighbor
									.pointOppositeIndex(face, pOpposite)].simplex;
							if (starFacet.containsPointInCircumsphere(pNew)) {
								// this facet faces pNew
								destroyedSimplexStar = starFacet;
								break;
							}
						}
					}
				}
			}
			destroyedSimplex = destroyedSimplexStar.original;
		}
		return destroyedSimplex;
	}

}
