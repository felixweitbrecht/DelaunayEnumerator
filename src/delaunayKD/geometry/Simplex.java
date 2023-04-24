package delaunayKD.geometry;

import static delaunayKD.AllSimplicesFinder.DIM;

import javax.management.RuntimeErrorException;

import delaunayKD.AllSimplicesFinder;
import delaunayKD.misc.UtilityMethods;

public class Simplex extends AbstractSimplex {

	public double[] circumcenter = new double[DIM];
	public double circumradius;

	// convention: simplex has face 0 as base, with newest point as top point
	public Simplex(Face[] faces, int lastPreviousKillerIndex) {
		super(lastPreviousKillerIndex);
		this.faces = faces;
		introduceSelfToFaces();
		if (AllSimplicesFinder.doAlphaBookkeeping) {
			computeCircum();
		}
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

	private void computeCircum() {
		// TODO general implementation
		if (DIM == 2) {
			// https://en.wikipedia.org/wiki/Circumscribed_circle#Circumcenter_coordinates
			Point a = faces[0].uniqueFace.points[0];
			Point b = faces[0].uniqueFace.points[1];
			Point c = pointOpposite(faces[0]);
			double d = 2 * (a.v[0] * (b.v[1] - c.v[1]) + b.v[0] * (c.v[1] - a.v[1]) + c.v[0] * (a.v[1] - b.v[1]));
			double ccX = ((a.v[0] * a.v[0] + a.v[1] * a.v[1]) * (b.v[1] - c.v[1])
					+ (b.v[0] * b.v[0] + b.v[1] * b.v[1]) * (c.v[1] - a.v[1])
					+ (c.v[0] * c.v[0] + c.v[1] * c.v[1]) * (a.v[1] - b.v[1])) / d;
			double ccY = ((a.v[0] * a.v[0] + a.v[1] * a.v[1]) * (c.v[0] - b.v[0])
					+ (b.v[0] * b.v[0] + b.v[1] * b.v[1]) * (a.v[0] - c.v[0])
					+ (c.v[0] * c.v[0] + c.v[1] * c.v[1]) * (b.v[0] - a.v[0])) / d;
			circumcenter[0] = ccX;
			circumcenter[1] = ccY;
			circumradius = Math.sqrt((a.v[0] - ccX) * (a.v[0] - ccX) + (a.v[1] - ccY) * (a.v[1] - ccY));
		} else if (DIM == 3) {
			// https://mathworld.wolfram.com/Circumsphere.html
			Face f = faces[0];
			Point[] points = new Point[DIM + 1];
			for (int i = 0; i < DIM; i++) {
				points[i] = f.uniqueFace.points[i];
			}
			points[DIM] = pointOpposite(f);
			double a = UtilityMethods.det(new double[] { points[0].v[0], points[0].v[1], points[0].v[2], 1,
					points[1].v[0], points[1].v[1], points[1].v[2], 1, points[2].v[0], points[2].v[1], points[2].v[2],
					1, points[3].v[0], points[3].v[1], points[3].v[2], 1, }, 4);
			double[] matrix = new double[] {
					points[0].v[0] * points[0].v[0] + points[0].v[1] * points[0].v[1] + points[0].v[2] * points[0].v[2],
					points[0].v[1], points[0].v[2], 1,
					points[1].v[0] * points[1].v[0] + points[1].v[1] * points[1].v[1] + points[1].v[2] * points[1].v[2],
					points[1].v[1], points[1].v[2], 1,
					points[2].v[0] * points[2].v[0] + points[2].v[1] * points[2].v[1] + points[2].v[2] * points[2].v[2],
					points[2].v[1], points[2].v[2], 1,
					points[3].v[0] * points[3].v[0] + points[3].v[1] * points[3].v[1] + points[3].v[2] * points[3].v[2],
					points[3].v[1], points[3].v[2], 1 };
			double dx = UtilityMethods.det(matrix, 4);
			matrix[1] = points[0].v[0];
			matrix[5] = points[1].v[0];
			matrix[9] = points[2].v[0];
			matrix[13] = points[3].v[0];
			double dy = -1 * UtilityMethods.det(matrix, 4);
			matrix[2] = points[0].v[1];
			matrix[6] = points[1].v[1];
			matrix[10] = points[2].v[1];
			matrix[14] = points[3].v[1];
			double dz = UtilityMethods.det(matrix, 4);
			matrix[3] = points[0].v[2];
			matrix[7] = points[1].v[2];
			matrix[11] = points[2].v[2];
			matrix[15] = points[3].v[2];
			double c = UtilityMethods.det(matrix, 4);

			circumcenter[0] = dx / (2 * a);
			circumcenter[1] = dy / (2 * a);
			circumcenter[2] = dz / (2 * a);

			circumradius = Math.sqrt(dx * dx + dy * dy + dz * dz - 4 * a * c) / (2 * Math.abs(a));
		} else {
			throw new RuntimeErrorException(new Error("Circumcenter computation not implemented for d = " + DIM));
		}
	}

}
