package delaunayKD.geometry;

import static delaunayKD.AllSimplicesFinder.DIM;

import delaunayKD.triangulator.HoleTriangulator;
import delaunayKD.triangulator.Star;

public class Point {

	// for search algorithms
	public boolean marked = false;

	// coordinates
	public double v[];

	// index/time stamp
	public int i;

	// hole triangulator for the cavity resulting from the removal of this point
	public HoleTriangulator ht;
	public Star star;

	public Point(double[] v, int i) {
		this.v = v;
		this.i = i;
	}

	@Override
	public String toString() {
		String s = i + " (";
		for (int idx = 0; idx < DIM; idx++) {
			String val = v[idx] + "";
			if (val.length() > 5) {
				val = val.substring(0, 5);
			}
			s += val;
			if (idx < DIM - 1) {
				s += ", ";
			}
		}
		s += ")";
		return s;
	}

}
