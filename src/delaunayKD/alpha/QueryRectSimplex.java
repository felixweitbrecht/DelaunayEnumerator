package delaunayKD.alpha;

import delaunayKD.geometry.AbstractSimplex;

public class QueryRectSimplex extends QueryRect {
	// left boundary is last previous killer index
	// right boundary is lowest index point of simplex
	// bottom boundary is highest index point of simplex
	// top boundary is first subsequent killer index

	public AbstractSimplex simplex;

	public QueryRectSimplex(AbstractSimplex simplex) {
		super(simplex.lastPreviousKillerIndex, simplex.minPoint().i, simplex.maxPoint().i,
				simplex.firstSubsequentKillerIndex);
		this.simplex = simplex;
	}

}
