package delaunayKD.triangulator;

import java.util.ArrayList;

import delaunayKD.geometry.AbstractSimplex;
import delaunayKD.geometry.Facet;
import delaunayKD.geometry.Point;
import delaunayKD.geometry.Face;
import delaunayKD.geometry.Simplex;

import static delaunayKD.AllSimplicesFinder.DIM;

public class IncrementalTriangulator extends Triangulator {

	// collect first DIM points (as they appear one by one)
	private Point[] firstFacePoints = new Point[DIM];

	// insert a point given its location (in the form of a simplex or facet
	// destroyed by it)
	// no location is necessary for the first DIM points
	public ArrayList<AbstractSimplex> addPoint(Point pNew, AbstractSimplex location) {
		if (location != null) {
			// destroy simplices that have pNew in their circumsphere
			excavate(pNew, location.faces[0]);
			// determine with which faces pNew creates new simplices
			findAttachingFaces(pNew);
			return createNewSimplices(pNew);
		} else {
			// remember first DIM points to create first face with
			for (int i = 0; i < DIM; i++) {
				if (firstFacePoints[i] == null) {
					firstFacePoints[i] = pNew;
					if (i == DIM - 1) {
						// initialize with first face and 2 facets
						Face faceNew = new Face(firstFacePoints);
						ArrayList<AbstractSimplex> newSimplices = new ArrayList<AbstractSimplex>(2);
						newSimplices.add(new Facet(faceNew, getLastPreviousKillerIndex()));
						newSimplices.add(new Facet(faceNew.r, getLastPreviousKillerIndex()));
						for (int faceIdx = 0; faceIdx < DIM; faceIdx++) {
							faceNew.hNeighbors[faceIdx] = faceNew.r;
							faceNew.r.hNeighbors[faceIdx] = faceNew;
						}
						return newSimplices;
					}
					break;
				}
			}
			return new ArrayList<AbstractSimplex>(0);
		}
	}

	// create new simplices with all attachingFaces
	private ArrayList<AbstractSimplex> createNewSimplices(Point pNew) {
		ArrayList<AbstractSimplex> newSimplices = new ArrayList<AbstractSimplex>();
		for (Face faceBase : attachingFaces) {
			// find or create the DIM new faces of the new simplex
			Face[] faces = new Face[DIM + 1];
			faces[0] = faceBase;
			Point[] faceBasePoints = faceBase.points();
			for (int pIdx = 0; pIdx < DIM; pIdx++) {
				Point pOpposite = faceBasePoints[pIdx];
				Face faceNeighboringBase = findNeighboringAttachingFace(faceBase, pOpposite);
				if (faceNeighboringBase != null) {
					// face is shared with neighboring new simplex
					if (faceNeighboringBase.simplex.maxPoint() == pNew) {
						// neighboring new simplex was already created, so the
						// shared face already exists
						faces[pIdx + 1] = faceNeighboringBase.simplex().faces[1
								+ faceNeighboringBase.pointOppositeIndex(faceBase, pOpposite)].r;
					} else {
						// neighboring new simplex is yet to be created, so we
						// create the shared face now
						faces[pIdx + 1] = faceBase.createFaceFacing(pNew, pIdx);
					}
				} else {
					// face is shared with new facet, which we create now
					Face faceNew = faceBase.createFaceFacing(pNew, pIdx);
					faces[pIdx + 1] = faceNew;
					newSimplices.add(new Facet(faceNew.r, getLastPreviousKillerIndex()));
					// find and set hull link
					Face faceReplacedHull = rotateThroughDeadSimplices(faceBase.r, pOpposite).r;
					faceNew.r.hLinkTo(
							faceReplacedHull.hNeighbors[faceReplacedHull.pointOppositeIndex(faceBase, pOpposite)],
							pIdx);
				}
			}
			newSimplices.add(new Simplex(faces, getLastPreviousKillerIndex()));
		}

		// set hull links between new outside faces by rotating through new
		// simplices around new outside ridges to find hull neighbor
		for (AbstractSimplex simplex : newSimplices) {
			if (simplex instanceof Facet) {
				Face face = simplex.faces[0];
				Point[] facePoints = face.points();
				for (int pIdx = 0; pIdx < DIM; pIdx++) {
					// link opposite pNew already set above
					if (face.hNeighbors[pIdx] == null) {
						Point pOpposite = facePoints[pIdx];
						Face faceHullNeighbor = rotateThroughSimplices(face, pOpposite).r;
						face.hLinkTo(faceHullNeighbor, pIdx);
					}
				}
			}
		}

		attachingFaces.clear();
		return newSimplices;
	}

	// given an attaching face, find the other attaching face using the ridge
	// opposite pOpposite. if none exists (happens with faces next to convex
	// hull boundary), return null.
	private Face findNeighboringAttachingFace(Face face, Point pOpposite) {
		Face faceRotate = rotateThroughDeadSimplices(face.r, pOpposite);
		if (faceRotate.r.simplex.isAlive()) {
			// we already found the excavation boundary
			return faceRotate;
		} else {
			// new point is outside convex hull. rotation hit the hull.
			// on the hull, step across the ridge being rotated around.
			Face faceHullNeighbor = faceRotate.r.hNeighbors[faceRotate.pointOppositeIndex(face, pOpposite)];
			if (faceHullNeighbor.simplex instanceof Simplex) {
				// attaching to the neighboring hull face. the new simplex there
				// was already created.
				return faceHullNeighbor;
			} else if (faceHullNeighbor.simplex.isDead()) {
				// not attaching to the neighboring hull face. we keep rotating.
				return rotateThroughDeadSimplices(faceHullNeighbor, faceHullNeighbor.pointOpposite(face, pOpposite));
			} else {
				// fallback for invalid cases
				return null;
			}
		}
	}

	@Override
	protected int getLastPreviousKillerIndex() {
		return Integer.MIN_VALUE;
	}

	@Override
	protected void excavate(Point pNew, Face faceInit) {
		simplicesToExplore.add(faceInit.simplex);
		exploreAndDestroySimplices(pNew);
	}

	@Override
	protected void findAttachingFaces(Point pNew) {
		// attach to cavity faces if simplex behind is alive
		for (AbstractSimplex simplex : destroyedSimplices) {
			for (Face face : simplex.faces) {
				if (face.r.simplex.isAlive()) {
					attachingFaces.add(face);
				}
			}
		}
		destroyedSimplices.clear();
	}

}
