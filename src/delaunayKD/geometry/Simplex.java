package delaunayKD.geometry;

import static delaunayKD.AllSimplicesFinder.DIM;

public class Simplex extends AbstractSimplex {

	// convention: simplex has face 0 as base, with newest point as top point
	public Simplex(Face[] faces, int lastPreviousKillerIndex) {
		super(lastPreviousKillerIndex);
		this.faces = faces;
		introduceSelfToFaces();
//		verify();
	}

	private void verify() {
		// determine involved points
		Point[] points = new Point[DIM + 1];
		for (int idx = 0; idx < DIM; idx++) {
			points[idx] = faces[0].points()[idx];
		}
		for (Point p : faces[1].points()) {
			if (!faces[0].hasVertex(p)) {
				points[DIM] = p;
				break;
			}
		}
		if (points[DIM] == null) {
			throw new RuntimeException("duplicate face supplied to simplex");
		}

		// sanity checks
		for (Face face : faces) {
			// ensure there are no more than DIM+1 points in total
			for (Point p : face.points()) {
				boolean found = false;
				for (Point pCompare : points) {
					if (pCompare == p) {
						found = true;
						break;
					}
				}
				if (!found) {
					throw new RuntimeException("simplex doesn't have exactly " + (DIM + 1) + " points!");
				}
			}
			// ensure each face faces the point that's not a vertex of that face
			for (Point p : points) {
				if (!face.hasVertex(p)) {
					if (!face.facesPoint(p)) {
						throw new RuntimeException("face of simplex faces the wrong way!");
					}
					break;
				}
			}
		}
	}

	// p must be vertex of this simplex
	public Face faceOpposite(Point p) {
		for (Face face : faces) {
			if (!face.hasVertex(p)) {
				return face;
			}
		}
		return null;
	}

	public Point pointOpposite(Face face) {
		for (Face face2 : faces) {
			if (face2 != face) {
				for (Point p : face2.points()) {
					if (!face.hasVertex(p)) {
						return p;
					}
				}
			}
		}
		return null;
	}

	@Override
	public boolean containsPointInCircumsphere(Point q) {
		return faces[0].simplexContainsPointInCircumsphere(maxPoint(), q);
	}

	@Override
	public String toString() {
		return "(" + faces[0].toString() + ") + " + pointOpposite(faces[0]).i;
	}

}
