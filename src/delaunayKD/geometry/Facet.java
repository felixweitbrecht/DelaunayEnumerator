package delaunayKD.geometry;

public class Facet extends AbstractSimplex {

	public Facet(Face face, int lastPreviousKillerIndex) {
		super(lastPreviousKillerIndex);
		this.faces = new Face[] { face };
		introduceSelfToFaces();
	}

	// insphere test degenerates to an orientation test
	@Override
	public boolean containsPointInCircumsphere(Point q) {
		return faces[0].facesPoint(q);
	}

	@Override
	public String toString() {
		return "(" + faces[0].toString() + ")";
	}

}
