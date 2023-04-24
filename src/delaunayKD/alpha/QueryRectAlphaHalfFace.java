package delaunayKD.alpha;

import static delaunayKD.AllSimplicesFinder.DIM;

import javax.management.RuntimeErrorException;

import delaunayKD.geometry.Face;
import delaunayKD.geometry.Facet;
import delaunayKD.geometry.Point;
import delaunayKD.geometry.Simplex;
import delaunayKD.misc.UtilityMethods;

public class QueryRectAlphaHalfFace extends QueryRect {

	// boundaries describe where the edge ceases to be Delaunay with a non-empty
	// radius range, and where its radius range changes

	public Face f;

	// radius range in which this half edge is alpha, both inclusive
	public double radiusMin;
	public double radiusMax;

	// create alpha half edge for intersection of two adjacent faces' lifetimes
	// if back face is a triangle, its circumcenter must be in front of the edge
	public QueryRectAlphaHalfFace(Face f, QueryRectSimplex front, QueryRectSimplex back) {
		super(Math.max(front.lowerMin, back.lowerMin), Math.min(front.lowerMax, back.lowerMax),
				Math.max(front.upperMin, back.upperMin), Math.min(front.upperMax, back.upperMax));

		this.f = f;

		if (front.simplex instanceof Facet) {
			radiusMax = Double.POSITIVE_INFINITY;
		} else {
			radiusMax = ((Simplex) front.simplex).circumradius;
		}
		if (back.simplex instanceof Facet) {
			radiusMin = radiusOfSmallestSphereThroughVertices();

		} else {
			radiusMin = ((Simplex) back.simplex).circumradius;
		}

		// TODO better handling of precision errors
		if (radiusMax < radiusMin && radiusMin - radiusMax > 0.00001) {
			throw new RuntimeException("radius range is empty");
		}
	}

	private double radiusOfSmallestSphereThroughVertices() {
		// TODO general implementation
		if (DIM == 2) { // half the edge length
			return 0.5
					* Math.sqrt(((f.points()[0].v[0] - f.points()[1].v[0]) * (f.points()[0].v[0] - f.points()[1].v[0]))
							+ ((f.points()[0].v[1] - f.points()[1].v[1]) * (f.points()[0].v[1] - f.points()[1].v[1])));
		} else if (DIM == 3) { // hacky
			// pick 2 other points, compute circumcenters with each point as
			// fourth point together with the face, compute distance of one
			// of f's points to the line through the circumcenters
			Point p41 = new Point(new double[] { 1, 2, 3 }, -1);
			Point p42 = new Point(new double[] { 0, -3, 14 }, -1);

			// https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
			double[] x0 = f.uniqueFace.points[0].v;
			double[] x1 = circumcenterWithFourthPoint3D(p41);
			double[] x2 = circumcenterWithFourthPoint3D(p42);

			double[] diff1 = new double[] { x0[0] - x1[0], x0[1] - x1[1], x0[2] - x1[2] };
			double[] diff2 = new double[] { x0[0] - x2[0], x0[1] - x2[1], x0[2] - x2[2] };
			double[] diff3 = new double[] { x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2] };
			double[] cross = new double[] { diff1[1] * diff2[2] - diff1[2] * diff2[1],
					diff1[2] * diff2[0] - diff1[0] * diff2[2], diff1[0] * diff2[1] - diff1[1] * diff2[0] };
			double len1 = Math.sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
			double len2 = Math.sqrt(diff3[0] * diff3[0] + diff3[1] * diff3[1] + diff3[2] * diff3[2]);
			return len1 / len2;
		} else {
			throw new RuntimeErrorException(new Error(
					"Computation of smallest sphere through vertices of face is not implemented for d = " + DIM));
		}
	}

	private double[] circumcenterWithFourthPoint3D(Point fourth) {
		double[] circumcenter = new double[3];

		// https://mathworld.wolfram.com/Circumsphere.html
		Point[] points = new Point[DIM + 1];
		for (int i = 0; i < DIM; i++) {
			points[i] = f.uniqueFace.points[i];
		}
		points[DIM] = fourth;

		double a = UtilityMethods.det(new double[] { points[0].v[0], points[0].v[1], points[0].v[2], 1, points[1].v[0],
				points[1].v[1], points[1].v[2], 1, points[2].v[0], points[2].v[1], points[2].v[2], 1, points[3].v[0],
				points[3].v[1], points[3].v[2], 1, }, 4);
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

		circumcenter[0] = dx / (2 * a);
		circumcenter[1] = dy / (2 * a);
		circumcenter[2] = dz / (2 * a);

		return circumcenter;
	}

	@Override
	public String toString() {
		return super.toString() + "x[" + radiusMin + "," + radiusMax + "[";
	}

}
