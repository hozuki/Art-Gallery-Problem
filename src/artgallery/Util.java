package artgallery;

import java.util.List;
import java.util.Objects;

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

	public static <E> int indexOfReference(final List<E> list, final E o) {
		final int listSize = list.size();

		for (int i = 0; i < listSize; ++i) {
			if (list.get(i) == o) {
				return i;
			}
		}

		return -1;
	}

	public static final double DOUBLE_DELTA = 0.000000001;

}
