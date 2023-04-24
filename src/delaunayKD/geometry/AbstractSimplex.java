package delaunayKD.geometry;

import java.util.ArrayList;

public abstract class AbstractSimplex {

	// for search algorithms
	private boolean marked = false;
	private static ArrayList<AbstractSimplex> markedSimplices = new ArrayList<AbstractSimplex>();

	// the bounding faces of this simplex
	// convention: the pre-existing face of new simplices must be at index 0
	public Face[] faces;

	// the indices of the last killer with index < smallest point index of
	// simplex and first killer with index > largest point index of simplex
	public int lastPreviousKillerIndex; // known at time of creation
	public int firstSubsequentKillerIndex = Integer.MAX_VALUE;

	// in stars: reference to original simplex copied into the star
	public AbstractSimplex original;

	public AbstractSimplex(int lastPreviousKillerIndex) {
		this.lastPreviousKillerIndex = lastPreviousKillerIndex;
	}

	abstract public boolean containsPointInCircumsphere(Point p);

	public boolean isDead() {
		return firstSubsequentKillerIndex < Integer.MAX_VALUE;
	}

	public boolean isAlive() {
		return firstSubsequentKillerIndex == Integer.MAX_VALUE;
	}

	public void mark() {
		marked = true;
		markedSimplices.add(this);
	}

	public boolean isMarked() {
		return marked;
	}

	static public void unmarkAll() {
		for (AbstractSimplex simplex : markedSimplices) {
			simplex.marked = false;
		}
		markedSimplices.clear();
	}

	public Point maxPoint() {
		return faces[faces.length - 1].maxPoint();
	}

	public Point minPoint() {
		return faces[0].minPoint();
	}

	protected void introduceSelfToFaces() {
		for (Face face : faces) {
			face.simplex = this;
		}
	}

}
