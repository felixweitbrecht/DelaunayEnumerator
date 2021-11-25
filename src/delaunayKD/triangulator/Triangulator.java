package delaunayKD.triangulator;

import java.util.ArrayList;
import java.util.Stack;

import delaunayKD.geometry.AbstractSimplex;
import delaunayKD.geometry.Facet;
import delaunayKD.geometry.Point;
import delaunayKD.geometry.Face;
import delaunayKD.geometry.Simplex;

public abstract class Triangulator {
	// faces to attach to during point insertion
	protected ArrayList<Face> attachingFaces = new ArrayList<Face>();

	// simplices destroyed during point insertion
	protected ArrayList<AbstractSimplex> destroyedSimplices = new ArrayList<AbstractSimplex>();

	// for exploration of destroyed simplices during point insertion
	protected Stack<AbstractSimplex> simplicesToExplore = new Stack<AbstractSimplex>();

	// destroys a simplex which has a new point in its circumsphere
	protected void destroy(AbstractSimplex simplex, Point pNew) {
		simplex.firstSubsequentKillerIndex = pNew.i;
		destroyedSimplices.add(simplex);
	}

	// the index set as time of birth on newly created faces
	abstract protected int getLastPreviousKillerIndex();

	// explore and destroy simplices with pNew in their circumsphere, starting
	// at faceInit
	abstract protected void excavate(Point pNew, Face faceInit);

	// determine faces of the cavity to attach to
	abstract protected void findAttachingFaces(Point pNew);

	// explore and destroy simplices with pNew in their circumsphere
	// needs initialization via simplicesToExplore
	protected void exploreAndDestroySimplices(Point pNew) {
		while (!simplicesToExplore.isEmpty()) {
			AbstractSimplex simplex = simplicesToExplore.pop();
			if (simplex != null && simplex.isAlive() && simplex.containsPointInCircumsphere(pNew)) {
				destroy(simplex, pNew);
				for (Face face : simplex.faces) {
					simplicesToExplore.add(face.r.simplex);
				}
				if (simplex instanceof Facet) {
					for (Face face : simplex.faces[0].hNeighbors) {
						simplicesToExplore.add(face.simplex);
					}
				}
			}
		}
	}

	// rotate around the ridge of faceFirst which is opposite pOpposite, against
	// the direction faceFirst is facing, while the encountered simplices are
	// proper simplices and dead
	public static Face rotateThroughDeadSimplices(Face faceFirst, Point pOpposite) {
		return rotateThroughSimplicesGeneral(faceFirst, pOpposite, true);
	}

	// rotate around the ridge of faceFirst which is opposite pOpposite, against
	// the direction faceFirst is facing, while we encounter proper simplices
	public static Face rotateThroughSimplices(Face faceFirst, Point pOpposite) {
		return rotateThroughSimplicesGeneral(faceFirst, pOpposite, false);
	}

	// rotates around the ridge opposite the specified point behind faceFirst
	// while proper simplices are found. if onlyDeadSimplices, rotation stops
	// when encountering an alive simplex.
	private static Face rotateThroughSimplicesGeneral(Face face, Point pOpposite, boolean onlyDeadSimplices) {
		while (true) {
			// rotate one step around the ridge, if possible
			Simplex simplexNext = face.r.simplex();
			if (simplexNext == null || (onlyDeadSimplices && simplexNext.isAlive())) {
				return face;
			}
			Face faceNext = simplexNext.faceOpposite(pOpposite);
			pOpposite = faceNext.pointOpposite(face, pOpposite);
			face = faceNext;
		}
	}

}
