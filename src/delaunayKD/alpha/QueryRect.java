package delaunayKD.alpha;

public class QueryRect {

	// Class to hold lifetime parameters of a query rectangle: ranges for upper
	// and lower index of interval of interest. lifetime here can be visualized
	// as an axis-aligned rectangle: x-axis: lower index, y-axis: upper index.

	// lowerMin and upperMax are exclusive, lowerMax and upperMin are inclusive.
	public int lowerMin;
	public int lowerMax;
	public int upperMin;
	public int upperMax;

	public QueryRect(int lowerMin, int lowerMax, int upperMin, int upperMax) {
		this.lowerMin = lowerMin;
		this.lowerMax = lowerMax;
		this.upperMin = upperMin;
		this.upperMax = upperMax;
		if (lowerMin >= lowerMax || upperMin >= upperMax) {
			throw new RuntimeException("illegal boundaries");
		}
	}

	@Override
	public String toString() {
		return "]" + lowerMin + ", " + lowerMax + "] x [" + upperMin + ", " + upperMax + "[";
	}
}
