package delaunayKD.geometry;

import java.util.ArrayList;
import java.util.Arrays;

import delaunayKD.misc.UtilityMethods;

import static delaunayKD.AllSimplicesFinder.DIM;

public class Face {
	// Every instance of a face (i.e. instance of a (DIM-1)-simplex) is created
	// twice: once for each direction induced by the normals of its supporting
	// hyperplane. The two objects are the reverse of each other, and they each
	// store a reference to the simplex attached to the respective side of that
	// face, allowing navigation between neighboring simplices.

	// In hole triangulations and star shapes, we call boundary faces those
	// sides of hull faces which do not face the point in the middle, i.e. the
	// outside sides of hull faces.

	// whether this is a reverse face (matters for orientation/insphere tests)
	public boolean isReverse = false;

	// backside
	public Face r;

	// the simplex bounded by this face, always null on boundary faces
	public AbstractSimplex simplex;

	// If face is part of the hull (i.e. not part of a proper simplex), these
	// are its neighboring faces on the hull. Old hull pointers aren't removed.
	// Convention: neighbor across ridge opposite point i is stored at index i.
	public Face hNeighbors[] = new Face[DIM];

	// Pointers on boundary faces between corresponding instances in star shapes
	// and hole triangulations, used for point location
	public Face faceBoundary;

	// Unique face: for any other instance of this face in other data
	// structures, the unique face will be the same for both instances
	public UniqueFace uniqueFace;

	// for search algorithms
	private boolean marked = false;
	private static ArrayList<Face> markedFaces = new ArrayList<Face>();

	// constructors ensuring uniqueFace is always passed on properly
	public Face(Point[] points) {
		this.uniqueFace = new UniqueFace(points);
		r = new Face(uniqueFace, true);
		r.r = this;
	}

	// constructor which doesn't also create the reverse object
	private Face(UniqueFace uniqueFace, boolean isReverse) {
		this.uniqueFace = uniqueFace;
		this.isReverse = isReverse;
	}

	public Face clone() {
		Face faceCloned = new Face(uniqueFace, isReverse);
		Face faceClonedR = new Face(uniqueFace, !isReverse);
		faceCloned.r = faceClonedR;
		faceClonedR.r = faceCloned;
		return faceCloned;
	}

	// requires that pNew faces this face.
	// creates a new face with pNew, sharing the ridge opposite pFacing with
	// this face. the returned face is oriented such that it faces pFacing.
	public Face createFaceFacing(Point pNew, int pFacingIdx) {
		// copy all points of faceNeighbor, but replace pFacing with pNew
		Point[] points = Arrays.copyOf(points(), DIM);
		points[pFacingIdx] = pNew;
		Face faceNew = new Face(points);
		// return reversed so new face faces this face
		return isReverse ? faceNew : faceNew.r;
	}

	// returns the proper simplex attached to this face, null if there's none
	public Simplex simplex() {
		return simplex instanceof Simplex ? (Simplex) simplex : null;
	}

	public Point[] points() {
		return uniqueFace.points;
	}

	// returns the hull face linked across the ridge opposite the given point.
	// no guarantees made for non-hull faces and non-vertex points.
	public Face hNeighborOpposite(Point pOpposite) {
		Point[] points = points();
		for (int i = 0; i < DIM; i++) {
			if (points[i] == pOpposite) {
				return hNeighbors[i];
			}
		}
		return null;
	}

	// sets back-and-forth hull neighbor pointers between this face and
	// hNeighbor over the shared ridge opposite the given point
	public void hLinkTo(Face hNeighbor, Point pOpposite) {
		Point[] points = points();
		for (int i = 0; i < DIM; i++) {
			if (points[i] == pOpposite) {
				hLinkTo(hNeighbor, i);
				return;
			}
		}
	}

	// sets back-and-forth hull neighbor pointers between this face and
	// hNeighbor over the shared ridge opposite point pIdx
	public void hLinkTo(Face hNeighbor, int pIdx) {
		hNeighbors[pIdx] = hNeighbor;
		int pNeighborOppositeIdx = hNeighbor.pointOppositeIndex(this, points()[pIdx]);
		hNeighbor.hNeighbors[pNeighborOppositeIdx] = this;
	}

	// orientation test
	public boolean facesPoint(Point q) {
		// www.cs.cmu.edu/~quake/robust.html
		Point[] points = points();
		double[] qVals = q.v;
		double[] matrix = new double[DIM * DIM];
		for (int row = 0; row < DIM; row++) {
			double[] vals = points[row].v;
			for (int col = 0; col < DIM; col++) {
				matrix[row * DIM + col] = vals[col] - qVals[col];
			}
		}
		return isReverse ^ (UtilityMethods.det(matrix, DIM) > 0.0);
	}

	// whether the simplex created from pTop and this face (must be facing pTop)
	// contains q in its circumsphere
	public boolean simplexContainsPointInCircumsphere(Point pTop, Point q) {
		// www.cs.cmu.edu/~quake/robust.html
		Point[] points = points();
		double[] qVals = q.v;
		double[] matrix = new double[(DIM + 1) * (DIM + 1)];
		for (int row = 0; row < DIM + 1; row++) {
			double[] vals = row < DIM ? points[row].v : pTop.v;
			double squareSum = 0.0;
			for (int col = 0; col < DIM; col++) {
				double diff = vals[col] - qVals[col];
				matrix[row * (DIM + 1) + col] = diff;
				squareSum += diff * diff;
			}
			matrix[row * (DIM + 1) + DIM] = squareSum;
		}
		return isReverse ^ (UtilityMethods.det(matrix, DIM + 1) > 0.0);
	}

	public boolean hasVertex(Point q) {
		Point[] points = points();
		for (Point point : points) {
			if (point == q) {
				return true;
			}
		}
		return false;
	}

	// finds this face's point opposite a given ridge (identified by the point
	// opposite of that ridge in another face sharing that ridge with this face)
	public Point pointOpposite(Face faceNeighbor, Point pOppositeNeighbor) {
		return points()[pointOppositeIndex(faceNeighbor, pOppositeNeighbor)];
	}

	// finds the index of this face's point opposite a given ridge (identified
	// by the point
	// opposite of that ridge in another face sharing that ridge with this face)
	public int pointOppositeIndex(Face faceNeighbor, Point pOppositeNeighbor) {
		Point[] points = points();
		for (int pIdx = 0; pIdx < DIM; pIdx++) {
			Point point = points[pIdx];
			if (point == pOppositeNeighbor || !faceNeighbor.hasVertex(point)) {
				// first condition is fallback for when this.r === faceNeighbor
				return pIdx;
			}
		}
		return -1;
	}

	// returns this face's corresponding instance in the appropriate star
	// shape, null if it doesn't exist
	public Face faceStar() {
		Face faceStar = uniqueFace.faceStar;
		if (faceStar != null && isReverse) {
			faceStar = faceStar.r;
		}
		return faceStar;
	}

	public Point maxPoint() {
		Point[] points = points();
		Point maxPoint = points[0];
		for (Point pCompare : points) {
			if (pCompare.i > maxPoint.i) {
				maxPoint = pCompare;
			}
		}
		return maxPoint;
	}

	public Point minPoint() {
		Point[] points = points();
		Point minPoint = points[0];
		for (Point pCompare : points) {
			if (pCompare.i < minPoint.i) {
				minPoint = pCompare;
			}
		}
		return minPoint;
	}

	@Override
	public String toString() {
		Point[] points = points();
		String s = points[0].i + "";
		if (isReverse) {
			s = "r: " + s;
		}
		for (int i = 1; i < DIM; i++) {
			s += " - " + points[i].i;
		}
		return s;
	}

	public void mark() {
		marked = true;
		markedFaces.add(this);
	}

	public boolean isMarked() {
		return marked;
	}

	static public void unmarkAll() {
		for (Face face : markedFaces) {
			face.marked = false;
		}
		markedFaces.clear();
	}

}
