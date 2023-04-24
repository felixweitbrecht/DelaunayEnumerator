package delaunayKD.alpha;

import java.util.ArrayList;

import delaunayKD.geometry.AbstractSimplex;

//holds a sorted list of query rectangles with a common left boundary,
//touching their neighbors at their top/bottom boundaries
public class AlignedQueryRects extends ArrayList<QueryRectSimplex> {

	private static final long serialVersionUID = -3560134750319703332L;

	// shared lowerMin boundary
	public int lowerMin;

	// create rects from non-empty list of simplices
	public AlignedQueryRects(ArrayList<AbstractSimplex> simplices) {
		super(simplices.size());
		if (simplices.isEmpty()) {
			throw new RuntimeException("Simplex list must not be empty!");
		}

		lowerMin = simplices.get(0).lastPreviousKillerIndex;

		// create QueryRectSimplex for each simplex and store it in list
		for (int i = 0; i < simplices.size(); i++) {
			add(new QueryRectSimplex(simplices.get(i)));
		}
	}
}
