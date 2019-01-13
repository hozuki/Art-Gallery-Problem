package artgallery;

public final class Util {

	private Util() {
	}

	public static boolean equals(double d1, double d2) {
		return equals(d1, d2, DOUBLE_DELTA);
	}

	public static boolean notEquals(double d1, double d2) {
		return !equals(d1, d2);
	}

	public static boolean equals(double d1, double d2, double delta) {
		double d = Math.abs(d2 - d1);
		return d <= delta;
	}

	public static boolean notEquals(double d1, double d2, double delta) {
		return !equals(d1, d2, delta);
	}

	public static final double DOUBLE_DELTA = 0.000000001;

}
