package delaunayKD.triangulator;

import java.util.ArrayList;
import java.util.Stack;

import delaunayKD.geometry.AbstractSimplex;
import delaunayKD.geometry.Facet;
import delaunayKD.geometry.Point;
import delaunayKD.geometry.Face;
import delaunayKD.geometry.Simplex;

import static delaunayKD.AllSimplicesFinder.DIM;

// terminology used in this class:
// boundary faces are that side of boundary faces which does not face the middle point
// front hull are those outside faces which are not boundary faces
// back hull are those outside faces which are also boundary faces
public class HoleTriangulator extends Triangulator {

	// this triangulation maintains the hole left after removing pMid
	public Point pMid;

	public HoleTriangulator(Point pMid) {
		this.pMid = pMid;
		pMid.ht = this;
	}

	public ArrayList<AbstractSimplex> update(Point pNew, ArrayList<Face> oldBoundaryFaces,
			ArrayList<Face> newBoundaryFaces) {
		// determine destroyed simplices, and faces to attach to
		destroySimplicesAndFindAttachingFaces(pNew, oldBoundaryFaces, newBoundaryFaces);
		// create all boundary faces with preliminary hull links
		createAndLinkNewBoundary(pNew, oldBoundaryFaces, newBoundaryFaces);
		// create new simplices and set remaining hull links
		return createNewSimplices(pNew, newBoundaryFaces);
	}

	// destroys all simplices that have pNew in their circumsphere and
	// determines with which faces new simplices are to be created
	private void destroySimplicesAndFindAttachingFaces(Point pNew, ArrayList<Face> oldBoundaryFaces,
			ArrayList<Face> newBoundaryFaces) {
		// gather faces to explore current hole triangulation from...
		ArrayList<Face> facesToExplore = new ArrayList<Face>();
		// old boundary faces from the face facing pNew
		for (Face faceStarOld : oldBoundaryFaces) {
			facesToExplore.add(faceStarOld.facesPoint(pNew) ? faceStarOld.faceBoundary : faceStarOld.faceBoundary.r);
		}
		// visible covering faces on the front hull
		// (those on the back hull are already considered as old star faces)
		for (Face faceStarNew : newBoundaryFaces) {
			Face faceBoundaryNeighbor = faceStarNew.hNeighborOpposite(pNew).faceBoundary;
			if (faceBoundaryNeighbor != null) {
				Face faceCovering = faceBoundaryNeighbor.hNeighbors[faceBoundaryNeighbor.pointOppositeIndex(faceStarNew,
						pNew)];
				if (faceCovering.simplex != null && faceCovering.facesPoint(pNew)) {
					facesToExplore.add(faceCovering);
				}
			}
		}
		// excavate from multiple faces because cavity may not be fully
		// connected
		for (Face face : facesToExplore) {
			excavate(pNew, face);
		}
		Face.unmarkAll();
		findAttachingFaces(pNew);
	}

	// creates all boundary faces with boundary pointers, sets some hull links:
	// boundary faces (old and new): all links inbetween
	// internally exposed boundary faces: link opposite pNew
	// boundary faces: link opposite pNew if linked face already exists
	private void createAndLinkNewBoundary(Point pNew, ArrayList<Face> oldBoundaryFaces,
			ArrayList<Face> newBoundaryFaces) {
		for (Face faceStar : newBoundaryFaces) {
			// copy new boundary faces into star shape
			Face faceHole = faceStar.clone();
			faceStar.faceBoundary = faceHole;
			faceHole.faceBoundary = faceStar;

			// for the ridges incident to pNew, if the adjacent boundary face
			// already exists (if it will at all), link to it (back hull)
			Point[] faceHolePoints = faceHole.points();
			for (int faceIdx = 0; faceIdx < DIM; faceIdx++) {
				Point pOpposite = faceHolePoints[faceIdx];
				if (pOpposite != pNew && faceHole.hNeighbors[faceIdx] == null) {
					Face faceBoundaryNeighbor = faceStar.hNeighbors[faceIdx].faceBoundary;
					if (faceBoundaryNeighbor != null) { // boundary continues?
						faceHole.hLinkTo(faceBoundaryNeighbor, faceIdx);
					}
				}
			}

			// for the ridge opposite pNew, if the boundary extends past it,
			// compute the hull link of faceHole.r if it stays exposed
			Face faceBoundaryNeighborOpposite = faceStar.hNeighborOpposite(pNew).faceBoundary;
			if (faceBoundaryNeighborOpposite != null) {
				Face faceBoundaryNeighborNeighbor = faceBoundaryNeighborOpposite.hNeighbors[faceBoundaryNeighborOpposite
						.pointOppositeIndex(faceStar, pNew)];
				if (faceBoundaryNeighborNeighbor.simplex instanceof Facet) {
					if (faceBoundaryNeighborNeighbor.simplex.isAlive()) {
						// new boundary faces continue where old boundary ended
						// and pNew doesn't destroy the facet next to the ridge
						faceHole.r.hLinkTo(faceBoundaryNeighborNeighbor, pNew);
					}
				} else {
					// faceHole replaces some boundary face. rotate through dead
					// simplices from there. if we hit the front hull and it's
					// dead, but its neighbor is alive, link to it.
					Face facePotentialFrontHull = rotateThroughDeadSimplices(faceBoundaryNeighborNeighbor,
							faceBoundaryNeighborNeighbor.pointOpposite(faceHole, pNew)).r;
					if (facePotentialFrontHull.simplex instanceof Facet && facePotentialFrontHull.simplex.isDead()) {
						Face facePotentialFrontHullNeighbor = facePotentialFrontHull.hNeighbors[facePotentialFrontHull
								.pointOppositeIndex(faceHole, pNew)];
						if (facePotentialFrontHullNeighbor.simplex.isAlive()) {
							faceHole.r.hLinkTo(facePotentialFrontHullNeighbor, pNew);
						}
					}
				}
				// finally link back hull if necessary and possible
				faceHole.hLinkTo(faceBoundaryNeighborOpposite, pNew);
			}
		}

		// set hull links on ridges opposite pNew where star shape ends
		for (Face faceStarOld : oldBoundaryFaces) {
			for (int faceIdx = 0; faceIdx < DIM; faceIdx++) {
				Face faceNeighborStar = faceStarOld.hNeighbors[faceIdx];
				// did the star shape end here, and does it end here still?
				if (faceNeighborStar.hasVertex(pMid) && faceNeighborStar.r.simplex != faceStarOld.r.simplex) {
					// (if it continued now, faceNeighborStar.r.simplex wouldn't
					// have been updated to a new simplex involving pNew)
					Face faceHoleReplacing = faceNeighborStar.hNeighborOpposite(pMid).faceBoundary;
					// faceHoleReplacing is the new boundary face replacing
					// faceStarOld
					Face faceFrontHullOld = faceStarOld.faceBoundary.hNeighbors[faceStarOld
							.pointOppositeIndex(faceHoleReplacing, pNew)];
					if (faceFrontHullOld.simplex.isDead()) {
						// current front hull survives, link to it
						faceHoleReplacing.hLinkTo(faceHoleReplacing.r, pNew);
					} else {
						// reverse side of new boundary face is new front hull
						faceHoleReplacing.hLinkTo(faceFrontHullOld, pNew);
					}
				}
			}
		}

		// special case: link over ridge opposite pNew for first face
		if (pNew.i == pMid.i + DIM) {
			Face faceBoundary = newBoundaryFaces.get(0).faceBoundary;
			faceBoundary.hLinkTo(faceBoundary.r, pNew);
		}
	}

	// creates all new simplices, sets the remaining hull links
	private ArrayList<AbstractSimplex> createNewSimplices(Point pNew, ArrayList<Face> newBoundaryFaces) {
		// temporarily mark all attaching faces so they're easy to recognize
		for (Face face : attachingFaces) {
			face.mark();
		}

		// create all simplices, and all facets on non-boundary faces
		ArrayList<AbstractSimplex> newSimplices = new ArrayList<AbstractSimplex>();
		for (Face faceAttaching : attachingFaces) {
			Face[] faces = new Face[DIM + 1];
			faces[0] = faceAttaching;
			Point[] facePoints = faceAttaching.points();

			// find or create the DIM new faces of the new simplex
			for (int pIdx = 0; pIdx < DIM; pIdx++) {
				Point pOpposite = facePoints[pIdx];
				Face faceCurr = faceAttaching.r;
				// rotate around the ridge of faceAttaching, towards pNew, to
				// find the neighboring attaching face, or determine that no
				// such face exists.
				// loop can take up to three iterations to terminate because
				// each iteration can only do one hull traversal step. there are
				// at most two hulls traversals necessary (front and back hull),
				// and one more rotation through dead simplices afterwards, so
				// three iterations is the maximum number of iterations.
				while (true) {
					faceCurr = rotateThroughDeadSimplices(faceCurr, faceCurr.pointOpposite(faceAttaching, pOpposite));

					// did we find the neighboring attaching face?
					if (faceCurr.isMarked()) {
						Simplex simplexNeighbor = faceCurr.simplex();
						if (simplexNeighbor != null && simplexNeighbor.maxPoint() == pNew) {
							// neighboring simplex was created already, so the
							// shared face already exists
							faces[pIdx + 1] = simplexNeighbor
									.faceOpposite(faceCurr.pointOpposite(faceAttaching, pOpposite)).r;
						} else {
							// neighboring simplex was not created yet
							// create the shared face now
							faces[pIdx + 1] = faceAttaching.createFaceFacing(pNew, pIdx);
						}
						break;
					}

					// we've hit a hull, step over ridge on hull
					Face faceHullNeighbor = faceCurr.r.hNeighbors[faceCurr.pointOppositeIndex(faceAttaching,
							pOpposite)];
					Face faceHullNeighborHullNeighbor = faceHullNeighbor.hNeighbors[faceHullNeighbor
							.pointOppositeIndex(faceAttaching, pOpposite)];
					if (faceHullNeighborHullNeighbor != faceCurr.r) {
						// new boundary begins here, use the new boundary face
						// (i.e. faceHullNeighborHullNeighbor.r)
						faces[pIdx + 1] = faceHullNeighborHullNeighbor.r;
						break;
					} else if ((faceCurr.r.simplex != null || faceHullNeighbor.simplex != null)
							&& !faceHullNeighbor.facesPoint(pNew)) {
						// reached front hull and pNew can not see the face
						// behind the ridge on the front hull
						// create face shared with new facet on the new front
						// hull
						Face faceNew = faceAttaching.createFaceFacing(pNew, pIdx);
						faces[pIdx + 1] = faceNew;
						newSimplices.add(new Facet(faceNew.r, getLastPreviousKillerIndex()));
						faceNew.r.hLinkTo(faceHullNeighbor, pIdx);
						break;
					}
					faceCurr = faceHullNeighbor;
				}
			}
			newSimplices.add(new Simplex(faces, getLastPreviousKillerIndex()));
		}
		attachingFaces.clear();
		Face.unmarkAll();

		// create all facets on the backside of boundary faces
		for (Face face : newBoundaryFaces) {
			if (face.faceBoundary.r.simplex == null) {
				// no simplex created here, so this must be on the front hull
				newSimplices.add(new Facet(face.faceBoundary.r, getLastPreviousKillerIndex()));
			}
		}

		// add missing front hull links by rotating through new simplices
		for (AbstractSimplex simplex : newSimplices) {
			if (simplex instanceof Facet) {
				Face face = simplex.faces[0];
				Point[] facePoints = face.points();
				for (int faceIdx = 0; faceIdx < DIM; faceIdx++) {
					Point pOpposite = facePoints[faceIdx];
					if (pOpposite != pNew) {
						// all links opposite pNew are already set at this point
						if (face.hNeighbors[faceIdx] == null) {
							Face faceCurr = face.r;
							while (!(faceCurr.simplex instanceof Facet)) { // rotate
								Simplex simplexNext = faceCurr.simplex();
								if (simplexNext != null) {
									faceCurr = simplexNext.faceOpposite(faceCurr.pointOpposite(face, pOpposite)).r;
								} else {
									// hit the back hull. does it continue?
									Face faceBoundaryNeighbor = faceCurr.hNeighbors[faceCurr.pointOppositeIndex(face,
											pOpposite)];
									if (faceBoundaryNeighbor != null) {
										faceCurr = faceBoundaryNeighbor.r;
									} else {
										// hull ends here, so we link to here
										break;
									}
								}
							}
							face.hLinkTo(faceCurr, faceIdx);
						}
					}
				}
			}
		}

		return newSimplices;
	}

	@Override
	protected void findAttachingFaces(Point pNew) {
		// attach to cavity faces if resulting circumsphere contains pMid
		// note that some attaching faces were already found in excavate()
		for (AbstractSimplex simplex : destroyedSimplices) {
			for (Face face : simplex.faces) {
				if ((face.r.simplex == null || face.r.simplex.isAlive())
						&& face.simplexContainsPointInCircumsphere(pNew, pMid)) {
					attachingFaces.add(face);
				}
			}
		}
		destroyedSimplices.clear();
	}

	@Override
	protected void excavate(Point pNew, Face faceInit) {
		if (faceInit.simplex instanceof Simplex) {
			simplicesToExplore.add(faceInit.simplex);
		} else { // hull face
			// manually explore hull because back hull doesn't have simplices
			Stack<Face> facesToExplore = new Stack<Face>();
			facesToExplore.add(faceInit);
			while (!facesToExplore.isEmpty()) {
				Face face = facesToExplore.pop();
				if (!face.isMarked()
						// ^^^ not yet explored
						&& ((face.simplex != null && face.simplex.isAlive())
								// ^^^ front hull visible to pNew
								|| (face.simplex == null && face.faceBoundary.r.simplex.isDead()))
						// ^^^ old boundary face
						&& face.facesPoint(pNew)) {
					// ^^^ visible from pNew
					face.mark(); // "whatever simplex would be here in the full
									// triangulation is dead"
					if (face.simplex != null) { // destroy face
						destroy(face.simplex, pNew);
					}
					if (face.r.simplex != null) { // keep excavating behind face
						if (face.r.simplex.containsPointInCircumsphere(pNew)) {
							simplicesToExplore.add(face.r.simplex);
						} else if (face.simplex == null) {
							// attach to old boundary face from outside the old
							// boundary because the cavity already ends here
							attachingFaces.add(face);
						}
					}
					// keep exploring hull
					for (Face faceNeighbor : face.hNeighbors) {
						facesToExplore.add(faceNeighbor);
					}
				}
			}
		}
		exploreAndDestroySimplices(pNew);
	}

	@Override
	protected int getLastPreviousKillerIndex() {
		return pMid.i;
	}

	@Override
	public String toString() {
		return "ht-" + pMid.i;
	}

}