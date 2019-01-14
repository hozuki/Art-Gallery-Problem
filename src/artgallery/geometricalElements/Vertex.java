package artgallery.geometricalElements;

import artgallery.Util;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Objects;

public class Vertex implements Comparable<Vertex> {

	private int id;
	private double X;
	private double Y;
	private boolean art;
	private boolean exit;
	private Polygon polygon;
	private ArrayList<Vertex> neighbors = new ArrayList<Vertex>();
	private ArrayList<Polygon> visibilityPolygon = new ArrayList<Polygon>();
	private boolean onBoundary;

	private Edge inEdge;
	private Edge outEdge;

	public Vertex(Polygon polygon) {
		this.polygon = polygon;
	}

	public Vertex(Polygon polygon, double x, double y) {
		this(polygon);
		this.setId(-1);
		this.X = x;
		this.Y = y;
		this.setArt(false);
		this.setExit(false);
	}

	public Vertex(Polygon polygon, double x, double y, int id, int artFlag, int exitFlag) {
		this(polygon, x, y);
		this.setId(id);
		this.setArt(artFlag == 1);
		this.setExit(exitFlag == 1);
	}

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	public Polygon getPolygon() {
		return polygon;
	}

	public double getX() {
		return X;
	}

	public void setX(double x) {
		X = x;
	}

	public double getY() {
		return Y;
	}

	public void setY(double y) {
		Y = y;
	}

	public boolean isArt() {
		return art;
	}

	public void setArt(boolean art) {
		this.art = art;
	}

	public boolean isExit() {
		return exit;
	}

	public void setExit(boolean exit) {
		this.exit = exit;
	}

	public void addNeighbor(Vertex n) {
		if (!neighbors.contains(n)) {
			neighbors.add(n);
		}
	}

	public void removeNeighbor(Vertex v) {
		if (neighbors.contains(v)) {
			v.neighbors.remove(this);
			neighbors.remove(v);
		}
	}

	public void swapNeighbors() {
		assert neighbors.size() == 2;
		Vertex n = neighbors.get(0);
		neighbors.set(0, neighbors.get(1));
		neighbors.set(1, n);
	}

	public ArrayList<Vertex> getNeighbors() {
		return neighbors;
	}

	public boolean isNeighbor(Vertex v) {
		return neighbors.contains(v);
	}

	public boolean isOnBoundary() {
		return onBoundary;
	}

	public void setOnBoundary(boolean onBoundary) {
		this.onBoundary = onBoundary;
	}

	public ArrayList<Polygon> getVisibilityPolygon() {
		return visibilityPolygon;
	}

	public void setVisibilityPolygon(ArrayList<Polygon> visibilityPolygon) {
		this.visibilityPolygon = visibilityPolygon;
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) {
			return true;
		}

		if (o instanceof Vertex) {
			final Vertex v = (Vertex) o;
			return Util.equals(getX(), v.getX()) && Util.equals(getY(), v.getY());
		} else {
			return false;
		}
	}

	@Override
	public int hashCode() {
		return Objects.hash(X, Y);
	}

	@Override
	public int compareTo(Vertex o) {
		if (this.equals(o)) {
			return 0;
		} else {
			return Double.compare(this.X, o.getX());
		}
	}

	@Override
	public Vertex clone() {
		return clone(polygon);
	}

	public Vertex clone(Polygon newPolygon) {
		return new Vertex(newPolygon, this.X, this.Y, this.id, this.art ? 1 : 0, this.exit ? 1 : 0);
	}

	@Override
	public String toString() {
		return "(" + this.getX() + ":" + this.getY() + ")";
	}

	public Edge getOutEdge() {
		return outEdge;
	}

	public Edge getInEdge() {
		return inEdge;
	}

	void setInEdge(Edge e) {
		inEdge = e;
	}

	void setOutEdge(Edge e) {
		outEdge = e;
	}

	public static final Comparator<Vertex> TOP_TO_BOTTOM_LEFT_TO_RIGHT = (final Vertex v1, final Vertex v2) -> {
		int result = -Double.compare(v1.getY(), v2.getY());

		if (result != 0) {
			return result;
		} else {
			return Double.compare(v1.getX(), v2.getX());
		}
	};

}
