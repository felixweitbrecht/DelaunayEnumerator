package delaunayKD.misc;

import static delaunayKD.AllSimplicesFinder.DIM;

import java.util.ArrayList;
import java.util.Random;

import delaunayKD.geometry.Point;

public class UtilityMethods {

	public static ArrayList<Point> generatePoints(int count) {
		return generatePoints(new Random(42L), count);
	}

	public static ArrayList<Point> generatePointsSphere(int count) {
		return generatePointsSphere(new Random(42L), count);
	}

	// generates points u.a.r. in the unit hypercube
	public static ArrayList<Point> generatePoints(Random random, int count) {
		ArrayList<Point> points = new ArrayList<Point>(count);
		for (int i = 0; i < count; i++) {
			double[] vals = new double[DIM];
			for (int dimIdx = 0; dimIdx < DIM; dimIdx++) {
				vals[dimIdx] = random.nextDouble();
			}
			points.add(new Point(vals, i));
		}
		return points;
	}

	public static ArrayList<Point> generatePointsSphere(Random random, int count) {
		int counter = 0;
		ArrayList<Point> points = new ArrayList<Point>(count);
		while (counter < count) {
			double[] vals = new double[DIM];
			double squareSum = 0.0;
			for (int dimIdx = 0; dimIdx < DIM; dimIdx++) {
				vals[dimIdx] = 2 * random.nextDouble() - 1;
				squareSum += vals[dimIdx] * vals[dimIdx];
			}
			if (squareSum < 1.0) {
				for (int dimIdx = 0; dimIdx < DIM; dimIdx++) {
					vals[dimIdx] = 5 + 5 * vals[dimIdx];
				}
				points.add(new Point(vals, counter));
				counter++;
			}
		}
		return points;
	}

	// computes the determinant of a DIM x DIM matrix
	// vals contains the rows of the matrix concatenated
	public static double det(double[] vals, int size) {
		if (size == 4) {
			return vals[0]
					* (vals[5] * (vals[10] * vals[15] - vals[11] * vals[14])
							- vals[6] * (vals[9] * vals[15] - vals[11] * vals[13])
							+ vals[7] * (vals[9] * vals[14] - vals[10] * vals[13]))
					- vals[1] * (vals[4] * (vals[10] * vals[15] - vals[11] * vals[14])
							- vals[6] * (vals[8] * vals[15] - vals[11] * vals[12])
							+ vals[7] * (vals[8] * vals[14] - vals[10] * vals[12]))
					+ vals[2] * (vals[4] * (vals[9] * vals[15] - vals[11] * vals[13])
							- vals[5] * (vals[8] * vals[15] - vals[11] * vals[12])
							+ vals[7] * (vals[8] * vals[13] - vals[9] * vals[12]))
					- vals[3] * (vals[4] * (vals[9] * vals[14] - vals[10] * vals[13])
							- vals[5] * (vals[8] * vals[14] - vals[10] * vals[12])
							+ vals[6] * (vals[8] * vals[13] - vals[9] * vals[12]));
		} else if (size == 3) {
			return vals[0] * (vals[4] * vals[8] - vals[5] * vals[7]) - vals[1] * (vals[3] * vals[8] - vals[5] * vals[6])
					+ vals[2] * (vals[3] * vals[7] - vals[4] * vals[6]);
		} else if (size == 2) {
			return vals[0] * vals[3] - vals[1] * vals[2];
		} else if (size == 1) {
			return vals[0];
		} else if (size < 2) {
			return -1.0;
		} else {
			double sum = 0.0;
			int factor = 1;
			for (int col = 0; col < size; col++) {
				double[] subMatrix = new double[(size - 1) * (size - 1)];
				for (int row = 0; row < size - 1; row++) {
					// copy row, leaving out the current column
					int colIndex = 0;
					for (int col2 = 0; col2 < size; col2++) {
						if (col2 != col) {
							subMatrix[row * (size - 1) + colIndex] = vals[(row + 1) * size + col2];
							colIndex++;
						}
					}
				}
				sum += factor * vals[col] * det(subMatrix, size - 1);
				factor *= -1;
			}
			return sum;
		}
	}

}
