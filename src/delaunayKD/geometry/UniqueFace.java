package delaunayKD.geometry;

public class UniqueFace {

	public Point[] points;

	// this face's instance in minPoint's starShape (front face!)
	public Face faceStar;

	public UniqueFace(Point[] points) {
		this.points = points;
	}

	public void setFaceStar(Face faceStar) {
		this.faceStar = faceStar.isReverse ? faceStar.r : faceStar;
	}

}
