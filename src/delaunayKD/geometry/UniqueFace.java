package delaunayKD.geometry;

import java.util.ArrayList;

import delaunayKD.alpha.AlignedQueryRects;

public class UniqueFace {

	public Point[] points;

	// this face's instance in minPoint's star (front face!)
	public Face faceStar;

	// for alpha face extraction:
	// whether this face is already in the list of unique faces
	public boolean picked = false;
	// for each orientation of this face: list of simplex query rects, nicely
	// sorted
	public ArrayList<AlignedQueryRects> rectsFront = new ArrayList<AlignedQueryRects>();
	public ArrayList<AlignedQueryRects> rectsBack = new ArrayList<AlignedQueryRects>();

	public UniqueFace(Point[] points) {
		this.points = points;
	}

	public void setFaceStar(Face faceStar) {
		this.faceStar = faceStar.isReverse ? faceStar.r : faceStar;
	}

}
