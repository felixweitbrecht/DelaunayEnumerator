package delaunayKD;

import static delaunayKD.AllSimplicesFinder.DIM;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import delaunayKD.alpha.AlignedQueryRects;
import delaunayKD.alpha.QueryRectAlphaHalfFace;
import delaunayKD.alpha.QueryRectSimplex;
import delaunayKD.geometry.Face;
import delaunayKD.geometry.Facet;
import delaunayKD.geometry.Point;
import delaunayKD.geometry.Simplex;
import delaunayKD.geometry.UniqueFace;
import delaunayKD.triangulator.IncrementalTriangulator;

//extractor for the temporal alpha-shape
public class AlphaFaceExtractor {

	// requires the incremental triangulator that was used to find all simplices
	public static ArrayList<QueryRectAlphaHalfFace> extractAlphaFaces(IncrementalTriangulator incTriangulator,
			ArrayList<Point> points) {
		// 1. have all faces register their simplices with the corresponding
		// unique face and fetch all UniqueFace instances
		ArrayList<UniqueFace> uniqueFaces = registerAllSimplicesAndGatherUniqueFaces(incTriangulator, points);

		// 2. combine pairs of simplices with overlapping lifetimes to find the
		// radius range for that face in the overlapping lifetime interval
		ArrayList<QueryRectAlphaHalfFace> allAlphaHalfFaces = new ArrayList<QueryRectAlphaHalfFace>();
		for (int i = 0; i < uniqueFaces.size(); i++) {
			UniqueFace u = uniqueFaces.get(i);
			Face testFace = new Face(u.points);
			new Facet(testFace, -1);
			new Facet(testFace.r, -1);
			int bottom = testFace.maxPoint().i;

			allAlphaHalfFaces.addAll(intersectQueryRects(testFace, u.rectsFront, u.rectsBack, bottom, u));
			allAlphaHalfFaces.addAll(intersectQueryRects(testFace.r, u.rectsBack, u.rectsFront, bottom, u));
		}

		return allAlphaHalfFaces;
	}

	// front are the rects of this half face, back are the rects of the reverse.
	// testFace is a face with facets on both sides, representing the face whose
	// intervals are to be extracted.
	public static ArrayList<QueryRectAlphaHalfFace> intersectQueryRects(Face testFace,
			ArrayList<AlignedQueryRects> front, ArrayList<AlignedQueryRects> back, int bottom, UniqueFace uniqueFace) {
		// right/bottom-aligned lists of query rects for both sides of the edge,
		// simplified to avoid overhead when intersecting
		ArrayList<List<QueryRectSimplex>> lists = constructSimplifiedRectLists(testFace, front, back, bottom);
		List<QueryRectSimplex> bottomListFront = lists.get(0);
		List<QueryRectSimplex> rightListFront = lists.get(1);
		List<QueryRectSimplex> bottomListBack = lists.get(2);
		List<QueryRectSimplex> rightListBack = lists.get(3);

		// neighbor indices:
		// for every rectangle in right lists, the index of its left neighbor
		// for every rectangle in bottom lists, the index of its top neighbor
		int[][] neighborsFront = constructNeighborIndices(bottomListFront, rightListFront);
		int[][] neighborsBack = constructNeighborIndices(bottomListBack, rightListBack);
		int[] topNeighborsFront = neighborsFront[0];
		int[] leftNeighborsFront = neighborsFront[1];
		int[] leftNeighborsBack = neighborsBack[1];

		// Deconstruct the back staircase by taking away one rectangle at a time
		// from the left or the top, only choosing rectangles which divide the
		// (remaining) staircase into two, guaranteeing that the remaining
		// rectangles on the (remaining) front staircase are entirely contained
		// in the chosen rectangle, or cut into two rectangles (but never a
		// rectangle and an L-shape!). For the front staircase, we keep
		// traversal indices for the right and bottom to indicate which rects
		// have not been covered entirely by back rects yet. These indices form
		// the starting point when intersecting with a back rectangle, and
		// neighbor indices as well as neighboring rectangles within lists allow
		// finding all other front rects intersecting the back rect.
		ArrayList<QueryRectAlphaHalfFace> alphaFaces = new ArrayList<QueryRectAlphaHalfFace>();

		int iRightFront = rightListFront.size() - 1;
		int iBottomFront = 0;

		int iRightBack = rightListBack.size() - 1;
		for (int iBottomBack = 0; iBottomBack < bottomListBack.size(); iBottomBack++) {
			QueryRectSimplex qrsBottomBack = bottomListBack.get(iBottomBack);
			// take away rects from rightListBack until qrsBottomBack is exposed
			while (iRightBack >= 0) {
				QueryRectSimplex qrsRightBack = rightListBack.get(iRightBack);
				int leftNeighborIndexBack = leftNeighborsBack[iRightBack];
				if (leftNeighborIndexBack == -1 || leftNeighborIndexBack < iBottomBack) {
					// qrsRightBack exposed, explore and intersect front rects

					// explore frontRightList while qrsRightBack intersects
					while (iRightFront >= 0) {
						QueryRectSimplex qrsRightFront = rightListFront.get(iRightFront);
						if (qrsRightFront.upperMax > qrsRightBack.upperMin) {
							// intersect (first figure out which of
							// qrsRightFront.simplex' faces we're working on
							// currently)
							Face face = null;
							for (int i = 0; i < DIM + 1; i++) {
								// for facets the first face must be a match so
								// we don't need to limit the dim index
								if (qrsRightFront.simplex.faces[i].uniqueFace == uniqueFace) {
									face = qrsRightFront.simplex.faces[i];
									break;
								}
							}
							alphaFaces.add(new QueryRectAlphaHalfFace(face, qrsRightFront, qrsRightBack));
						}
						if (qrsRightFront.upperMin > qrsRightBack.upperMin) {
							// done with this rect
							iRightFront--;
						} else {
							break;
						}
						// qrsRightFront.upperMin == qrsRightBack.upperMin never
						// happens because inserting a point once the face
						// exists can not create new simplices on both sides of
						// the face simultaneously, so after the first
						// intersection we are guaranteed that iRightFront is an
						// index of a rect that intersects qrsRightBack, or 0
					}

					// Explore remaining part of frontBottomList left-wards
					// starting at the left neighbor of the last intersected
					// rect (or the very end if none were intersected). This
					// works even if frontRightList is empty, and once all its
					// rects are worked off. Don't advance iBottomFront because
					// qrsRightBack doesn't reach the bottom boundary.
					int leftNeighborIndexFront = iRightFront >= 0 ? leftNeighborsFront[iRightFront]
							: bottomListFront.size() - 1;
					// leftNeighborIndexFront may be -1 for "no neighbor", but
					// the while condition catches that already
					while (leftNeighborIndexFront >= iBottomFront) {
						QueryRectSimplex qrsBottomFront = bottomListFront.get(leftNeighborIndexFront);
						if (qrsBottomFront.upperMax > qrsRightBack.upperMin) {
							// still overlaps, intersect and keep exploring
							// (but first figure out which of
							// qrsBottomFront.simplex'
							// faces we're working on currently)
							Face face = null;
							for (int i = 0; i < DIM + 1; i++) {
								// for facets the first face must be a match so
								// we don't need to limit the dim index
								if (qrsBottomFront.simplex.faces[i].uniqueFace == uniqueFace) {
									face = qrsBottomFront.simplex.faces[i];
									break;
								}
							}
							alphaFaces.add(new QueryRectAlphaHalfFace(face, qrsBottomFront, qrsRightBack));
							leftNeighborIndexFront--;
						} else {
							// no more overlap
							break;
						}
					}

					iRightBack--;
				} else {
					// explored all exposed rects of rightListBack
					break;
				}
			}

			// qrsBottomBack now exposed, explore and intersect front rects

			// explore frontBottomList while qrsBottomBack intersects
			while (iBottomFront < bottomListFront.size()) {
				QueryRectSimplex qrsBottomFront = bottomListFront.get(iBottomFront);
				if (qrsBottomFront.lowerMin < qrsBottomBack.lowerMax) {
					// intersect (first figure out which of
					// qrsBottomFront.simplex' faces we're working on
					// currently)
					Face face = null;
					for (int i = 0; i < DIM + 1; i++) {
						// for facets the first face must be a match so
						// we don't need to limit the dim index
						if (qrsBottomFront.simplex.faces[i].uniqueFace == uniqueFace) {
							face = qrsBottomFront.simplex.faces[i];
							break;
						}
					}

					alphaFaces.add(new QueryRectAlphaHalfFace(face, qrsBottomFront, qrsBottomBack));
				}
				if (qrsBottomFront.lowerMax < qrsBottomBack.lowerMax || iBottomBack == bottomListBack.size() - 1) {
					// done with this rect
					iBottomFront++;
				} else {
					break;
				}
				// qrsBottomFront.lowerMax == qrsBottomBack.lowerMax only
				// happens at the right boundary because inserting a point once
				// the face already exists can not create new simplices on both
				// sides of the face simultaneously, so after the first
				// intersection we are guaranteed that iBottomFront is an index
				// of a rect that intersects qrsRightBack
			}

			// Explore remaining part of frontRightList towards the top starting
			// at the top neighbor of the last intersected rect. This also works
			// once all rects of frontBottomList are worked off. No need to
			// advance iRightFront because once qrsBottomBack reaches the right
			// boundary, no more back rects to intersect with remain.
			int topNeighborIndexFront = iBottomFront < bottomListFront.size() ? topNeighborsFront[iBottomFront] : 0;
			if (topNeighborIndexFront != -1) {
				while (topNeighborIndexFront <= iRightFront) {
					QueryRectSimplex qrsRightFront = rightListFront.get(topNeighborIndexFront);
					if (qrsRightFront.lowerMin < qrsBottomBack.lowerMax) {
						// still overlaps, intersect and keep exploring
						// (but first figure out which of qrsRightFront.simplex'
						// faces we're working on currently)
						Face face = null;
						for (int i = 0; i < DIM + 1; i++) {
							// for facets the first face must be a match so
							// we don't need to limit the dim index
							if (qrsRightFront.simplex.faces[i].uniqueFace==uniqueFace) {
								face = qrsRightFront.simplex.faces[i];
								break;
							}
						}
						alphaFaces.add(new QueryRectAlphaHalfFace(face, qrsRightFront, qrsBottomBack));
						topNeighborIndexFront++;
					} else {
						// no more overlap
						break;
					}
				}
			}
		}

		return alphaFaces;
	}

	// for every rectangle in right list: a pointer to its left neighbor
	// for every rectangle in bottom list: a pointer to its top neighbor
	// returns [topNeighbors, leftNeighbors]
	private static int[][] constructNeighborIndices(List<QueryRectSimplex> bottomList,
			List<QueryRectSimplex> rightList) {
		int[] topNeighbors = new int[bottomList.size()];
		int[] leftNeighbors = new int[rightList.size()];
		// neighbor indices are -1 where no neighbor exists
		Arrays.fill(topNeighbors, -1);
		Arrays.fill(leftNeighbors, -1);

		// starting at bottom-right walk along the dividing boundaries between
		// bottom and right rects, matching neighbors to each other
		int rightIndex = 0;
		for (int bottomIndex = bottomList.size() - 1; bottomIndex >= 0; bottomIndex--) {
			QueryRectSimplex qrsBottom = bottomList.get(bottomIndex);
			// find top neighbor of current bottom rect by advancing in the
			// right rect list until the top neighbor is found. all right rects
			// advanced over in this step have that bottom rect as left
			// neighbor.
			while (rightIndex < rightList.size()) {
				QueryRectSimplex qrsRight = rightList.get(rightIndex);
				if (qrsRight.upperMin == qrsBottom.upperMax) {
					// found top neighbor
					topNeighbors[bottomIndex] = rightIndex;
					break;
				} else {
					// found left neighbor of right rect, keep advancing
					leftNeighbors[rightIndex] = bottomIndex;
					rightIndex++;
				}
			}
		}
		return new int[][] { topNeighbors, leftNeighbors };
	}

	// Constructs sorted lists of right-aligned/bottom-aligned front/back
	// simplices. Front lists exclude rects of proper simplices with
	// circumcenter on back side. Back lists merge all rects of proper simplices
	// with circumcenter on back side into a fake rect belonging to testFaces's
	// reverse's facet. Returns in this order: bottom front, right front, bottom
	// back, right back
	private static ArrayList<List<QueryRectSimplex>> constructSimplifiedRectLists(Face testFace,
			ArrayList<AlignedQueryRects> front, ArrayList<AlignedQueryRects> back, int bottom) {
		// build sorted lists of right-aligned and bottom-aligned rects (front):
		List<QueryRectSimplex> rightListFront = constructRightList(front, bottom);
		List<QueryRectSimplex> bottomListFront = constructBottomList(front, bottom);
		// discard rects of proper simplices with circumcenter on back side
		for (int i = 0; i < rightListFront.size(); i++) {
			QueryRectSimplex qrs = rightListFront.get(i);
			if (qrs.simplex instanceof Simplex) {
				if (!testFace.facesPoint(((Simplex) qrs.simplex).circumcenter)) {
					// discard from here
					rightListFront = rightListFront.subList(0, i);
					break;
				}
			}
		}
		for (int i = bottomListFront.size() - 1; i >= 0; i--) {
			QueryRectSimplex qrs = bottomListFront.get(i);
			if (qrs.simplex instanceof Simplex) {
				if (!testFace.facesPoint(((Simplex) qrs.simplex).circumcenter)) {
					// discard up until here
					bottomListFront = bottomListFront.subList(i + 1, bottomListFront.size());
					break;
				}
			}
		}

		// build sorted lists of right-aligned and bottom-aligned rects (back):
		List<QueryRectSimplex> rightListBack = constructRightList(back, bottom);
		List<QueryRectSimplex> bottomListBack = constructBottomList(back, bottom);
		// merge rects of proper simplices with circumcenter on back side, if
		// any
		// find left/right boundaries and potential top boundary in bottom
		// list
		int left = -1;
		int top = -1;
		for (int i = bottomListBack.size() - 1; i >= 0; i--) {
			boolean merge = true;
			QueryRectSimplex qrs = bottomListBack.get(i);
			if (qrs.simplex instanceof Simplex) {
				if (testFace.facesPoint(((Simplex) qrs.simplex).circumcenter)) {
					merge = false;
				}
			}
			if (merge) {
				left = qrs.lowerMin;
				top = qrs.upperMax;
				if (i == 0) {
					// entire list will be merged
					bottomListBack.clear();
				}
			} else {
				// remove merged rects from list
				bottomListBack = bottomListBack.subList(0, i + 1);
				break;
			}
		}
		// correct/confirm potential top boundary in right list
		for (int i = 0; i < rightListBack.size(); i++) {
			QueryRectSimplex qrs = rightListBack.get(i);
			if (!testFace.facesPoint(((Simplex) qrs.simplex).circumcenter)) {
				// right rects' top can't be below bottom rects' top, this maxes
				top = qrs.upperMax;
				if (i == rightListBack.size() - 1) {
					// entire list will be merged
					rightListBack.clear();
				}
			} else {
				// remove merged rects from list
				rightListBack = rightListBack.subList(i, rightListBack.size());
				break;
			}
		}
		if (left != -1) {
			// there are rects to merge
			// create fake query rect simplex using the facet on testFace.r
			Facet facet = (Facet) testFace.r.simplex;
			facet.lastPreviousKillerIndex = left;
			facet.firstSubsequentKillerIndex = top;
			QueryRectSimplex qrs = new QueryRectSimplex(facet);
			bottomListBack.add(qrs);
		}

		ArrayList<List<QueryRectSimplex>> returnList = new ArrayList<List<QueryRectSimplex>>(4);
		returnList.add(bottomListFront);
		returnList.add(rightListFront);
		returnList.add(bottomListBack);
		returnList.add(rightListBack);
		return returnList;
	}

	// constructs a sorted list of query rects touching the right boundary,
	// excluding the bottom-right rect
	private static ArrayList<QueryRectSimplex> constructRightList(ArrayList<AlignedQueryRects> columns, int bottom) {
		ArrayList<QueryRectSimplex> rightList = new ArrayList<QueryRectSimplex>();
		// get rects below (lowest rect of columns that don't align to bottom)
		for (int colIndex = columns.size() - 1; colIndex >= 0; colIndex--) {
			AlignedQueryRects column = columns.get(colIndex);
			if (column.get(0).upperMin == bottom) {
				rightList.addAll(column.subList(1, column.size()));
			}
		}
		// get rects of columns that don't align to bottom
		for (int colIndex = 0; colIndex < columns.size(); colIndex++) {
			AlignedQueryRects column = columns.get(colIndex);
			if (column.get(0).upperMin > bottom) {
				rightList.addAll(column);
			}
		}
		return rightList;
	}

	// constructs a sorted list of query rects touching the bottom boundary
	private static ArrayList<QueryRectSimplex> constructBottomList(ArrayList<AlignedQueryRects> columns, int bottom) {
		ArrayList<QueryRectSimplex> bottomList = new ArrayList<QueryRectSimplex>();
		// get bottom rects of columns that align to bottom
		for (int colIndex = 0; colIndex < columns.size(); colIndex++) {
			QueryRectSimplex qrs = columns.get(colIndex).get(0);
			if (qrs.upperMin == bottom) {
				bottomList.add(qrs);
			}
		}
		return bottomList;
	}

	// registers all simplices of all faces (and their reverse counterparts) of
	// all triangulators with the corresponding unique face. returns all
	// UniqueFace instances
	public static ArrayList<UniqueFace> registerAllSimplicesAndGatherUniqueFaces(
			IncrementalTriangulator incTriangulator, ArrayList<Point> points) {
		ArrayList<UniqueFace> uniqueFaces = new ArrayList<UniqueFace>();
		for (Face f : incTriangulator.knownFaces) {
			registerSimplicesAndGatherUniqueFaces(f, uniqueFaces);
		}
		for (int i = 0; i < points.size(); i++) {
			for (Face f : points.get(i).ht.knownFaces) {
				registerSimplicesAndGatherUniqueFaces(f, uniqueFaces);
			}
		}
		return uniqueFaces;
	}

	// converts all face's known simplices to query rectangles and inserts them
	// into a sorted list of (sorted) AlignedQueryRects
	private static void registerSimplicesAndGatherUniqueFaces(Face f, ArrayList<UniqueFace> uniqueFaces) {
		if (f.knownSimplices.isEmpty() && f.r.knownSimplices.isEmpty()) {
			System.out.println("this shouldn't happen");
		}

		UniqueFace u = f.uniqueFace;

		ArrayList<AlignedQueryRects> ownList = f.isReverse ? u.rectsBack : u.rectsFront;
		ArrayList<AlignedQueryRects> otherList = f.isReverse ? u.rectsFront : u.rectsBack;

		if (!f.knownSimplices.isEmpty()) {
			ownList.add(new AlignedQueryRects(f.knownSimplices));
		}

		if (!f.r.knownSimplices.isEmpty()) {
			otherList.add(new AlignedQueryRects(f.r.knownSimplices));
		}

		// collect unique faces
		if (!u.picked) {
			u.picked = true;
			uniqueFaces.add(u);
		}
	}

}
