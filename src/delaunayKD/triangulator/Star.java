package delaunayKD.triangulator;

import java.util.ArrayList;
import java.util.Stack;

import delaunayKD.geometry.AbstractSimplex;
import delaunayKD.geometry.Facet;
import delaunayKD.geometry.Point;
import delaunayKD.geometry.Face;
import delaunayKD.geometry.Simplex;

import static delaunayKD.AllSimplicesFinder.DIM;

public class Star {
	// point in the middle
	public Point pMid;

	// state information of ongoing star updates
	private ArrayList<AbstractSimplex> registeredSimplices = new ArrayList<AbstractSimplex>();
	private int unmatchedFaces = 0; // currently known unmatched faces
	public Face faceLatest = null; // the newest star face
	private Simplex simplexLatest = null; // some simplex from the last update
											// (original instance)

	public Star(Point pMid) {
		this.pMid = pMid;
		pMid.star = this;
	}

	// registers a new simplex of the star. if now all new simplices are
	// known, the star and hole triangulation are updated and new
	// simplices found in the hole triangulation are returned
	public ArrayList<AbstractSimplex> registerSimplex(AbstractSimplex newSimplex, Point pNew) {
		registeredSimplices.add(newSimplex);
		simplexLatest = newSimplex instanceof Simplex ? (Simplex) newSimplex : simplexLatest;
		ArrayList<AbstractSimplex> holeTriangulationSimplices;
		// ensure new faces also exist as star faces, and track how many faces
		// don't have 2 simplices yet (update once all have 2 simplices)
		unmatchedFaces += ensureStarFaceExistence(pNew, newSimplex);
		if (unmatchedFaces == 0) {
			if (pNew.i == pMid.i + DIM - 1) {
				initWithFirstFace();
				holeTriangulationSimplices = new ArrayList<AbstractSimplex>(0);
			} else {
				ArrayList<AbstractSimplex> destroyedSimplices = findDestroyedSimplices(simplexLatest, pNew);
				ArrayList<Face> oldBoundaryFaces = findOldBoundaryFaces(destroyedSimplices);
				ArrayList<Face> newBoundaryFaces = updateStar(pNew, destroyedSimplices);
				holeTriangulationSimplices = pMid.ht.update(pNew, oldBoundaryFaces, newBoundaryFaces);
			}
			registeredSimplices.clear();
		} else {
			// no update yet
			holeTriangulationSimplices = new ArrayList<AbstractSimplex>(0);
		}
		return holeTriangulationSimplices;
	}

	// returns (matched star faces count) - (unmatched star faces count) of
	// newSimplex
	private int ensureStarFaceExistence(Point pNew, AbstractSimplex newSimplex) {
		int existingFaces = 0;
		int newFaces = 0;
		// ensure new faces also exist as star faces, and track how many faces
		// don't have 2 simplices yet (update once all have 2 simplices)
		for (Face face : newSimplex.faces) {
			if (face.hasVertex(pMid) && face.hasVertex(pNew)) {
				// face incident to pMid and pNew, so it must have a new simplex
				// on both sides eventually
				if (face.faceStar() == null) {
					Face faceStar = face.clone();
					face.uniqueFace.setFaceStar(faceStar);
					faceLatest = faceStar;
					newFaces++;
				} else {
					existingFaces++;
				}
			}
		}
		return existingFaces - newFaces;
	}

	// inserts all new simplices into the star
	private ArrayList<Face> updateStar(Point pNew, ArrayList<AbstractSimplex> destroyedSimplices) {
		ArrayList<Facet> newFacets = new ArrayList<Facet>();
		ArrayList<Face> newBoundaryFaces = new ArrayList<Face>();
		// clone all new simplices into star
		for (AbstractSimplex simplex : registeredSimplices) {
			if (simplex instanceof Simplex) {
				Face[] faces = new Face[DIM + 1];
				// find or create star face instances of faces
				for (int i = 0; i < DIM + 1; i++) {
					Face faceOrig = simplex.faces[i];
					if (faceOrig.hasVertex(pMid)) {
						// star face (incident to pMid)
						faces[i] = faceOrig.faceStar();
					} else {
						// boundary face
						Face faceBoundaryInside = faceOrig.clone();
						faces[i] = faceBoundaryInside;
						newBoundaryFaces.add(faceBoundaryInside.r);
					}
				}
				Simplex simplexClone = new Simplex(faces, -1);
				simplexClone.original = simplex;
			} else {
				Face faceStar = simplex.faces[0].faceStar();
				Facet facetClone = new Facet(faceStar, -1);
				facetClone.original = simplex;
				newFacets.add(facetClone);
			}
		}

		// set hull links of boundary faces
		for (Face face : newBoundaryFaces) {
			Point[] facePoints = face.points();
			for (int faceIdx = 0; faceIdx < DIM; faceIdx++) {
				Point pOpposite = facePoints[faceIdx];
				if (face.hNeighbors[faceIdx] == null) {
					Face faceOutsideSimplex = face.r.simplex().faceOpposite(pOpposite).r;
					Simplex simplexNeighbor = faceOutsideSimplex.simplex();
					if (simplexNeighbor != null) {
						// boundary continues, link to neighboring boundary face
						face.hLinkTo(simplexNeighbor.faceOpposite(pMid).r, faceIdx);
					} else {
						// boundary ends, link to facet incident to pMid
						face.hLinkTo(faceOutsideSimplex, faceIdx);
					}
				}
			}
		}

		// set remaining hull pointers
		// 1. links over ridges with pNew: rotate through new simplices
		for (Facet facet : newFacets) {
			Face face = facet.faces[0];
			Point[] facePoints = face.points();
			for (int faceIdx = 0; faceIdx < DIM; faceIdx++) {
				Point pOpposite = facePoints[faceIdx];
				if (pOpposite != pNew && face.hNeighbors[faceIdx] == null) {
					face.hLinkTo(Triangulator.rotateThroughSimplices(face, pOpposite).r, faceIdx);
				}
			}
		}
		// 2. links over ridges opposite pNew: from every ridge of destroyed
		// facets which stays on the hull, rotate through dead simplices, then
		// rotate back through new simplices, to find the new face that attached
		// at this ridge
		for (AbstractSimplex destroyedSimplex : destroyedSimplices) {
			if (destroyedSimplex instanceof Facet) {
				Face face = destroyedSimplex.faces[0];
				Point[] facePoints = face.points();
				for (int faceIdx = 0; faceIdx < DIM; faceIdx++) {
					Point pOpposite = facePoints[faceIdx];
					if (pOpposite == pMid) {
						// link to boundary face - already set
						continue;
					}
					Face faceNeighbor = face.hNeighbors[faceIdx];
					if (faceNeighbor.simplex.isDead() || faceNeighbor.simplex.maxPoint() == pNew) {
						// ridge opposite pOpposite isn't on the hull
						continue;
					}
					Face faceRotate = Triangulator.rotateThroughDeadSimplices(face, pOpposite);
					Face faceReplacing = Triangulator.rotateThroughSimplices(faceRotate.r,
							faceRotate.r.pointOpposite(face, pOpposite)).r;
					// faceReplacing is the face of a new facet
					faceReplacing.hLinkTo(faceNeighbor, pNew);
				}
			}
		}

		return newBoundaryFaces;
	}

	private ArrayList<Face> findOldBoundaryFaces(ArrayList<AbstractSimplex> destroyedSimplices) {
		ArrayList<Face> oldBoundaryFaces = new ArrayList<Face>();
		for (AbstractSimplex simplex : destroyedSimplices) {
			if (simplex instanceof Simplex) {
				oldBoundaryFaces.add(((Simplex) simplex).faceOpposite(pMid).r);
			}
		}
		return oldBoundaryFaces;
	}

	private ArrayList<AbstractSimplex> findDestroyedSimplices(Simplex simplexNew, Point pNew) {
		ArrayList<AbstractSimplex> destroyedSimplices = new ArrayList<AbstractSimplex>();
		// BFS to find all simplices destroyed upon insertion of pNew
		// (uses firstSubsequentKillerIndex as visited marker)
		Stack<AbstractSimplex> simplicesToDestroy = new Stack<AbstractSimplex>();
		simplicesToDestroy.add(simplexNew.faceOpposite(pNew).faceStar().simplex);
		while (!simplicesToDestroy.isEmpty()) {
			AbstractSimplex simplex = simplicesToDestroy.pop();
			if (simplex.isAlive() && simplex.containsPointInCircumsphere(pNew)) {
				simplex.firstSubsequentKillerIndex = pNew.i;
				destroyedSimplices.add(simplex);
				for (Face face : simplex.faces) {
					if (face.r.simplex != null) { // nothing behind boundary
						simplicesToDestroy.add(face.r.simplex);
					}
				}
				if (simplex instanceof Facet) {
					Face face = simplex.faces[0];
					Point[] facePoints = face.points();
					for (int faceIdx = 0; faceIdx < DIM; faceIdx++) {
						if (facePoints[faceIdx] != pMid) {
							// there are no simplices behind boundary!
							simplicesToDestroy.add(face.hNeighbors[faceIdx].simplex);
						}
					}
				}
			}
		}
		return destroyedSimplices;
	}

	private void initWithFirstFace() {
		AbstractSimplex facet1 = registeredSimplices.get(0);
		AbstractSimplex facet2 = registeredSimplices.get(1);

		Face faceStar = facet1.faces[0].faceStar();
		for (int faceIdx = 0; faceIdx < DIM; faceIdx++) {
			faceStar.hNeighbors[faceIdx] = faceStar.r;
			faceStar.r.hNeighbors[faceIdx] = faceStar;
		}

		Facet facet1clone = new Facet(faceStar, -1);
		facet1clone.original = facet1;
		Facet facet2clone = new Facet(faceStar.r, -1);
		facet2clone.original = facet2;
	}

	@Override
	public String toString() {
		return "star-" + pMid.i;
	}

}
