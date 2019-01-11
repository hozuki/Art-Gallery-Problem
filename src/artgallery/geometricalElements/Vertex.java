package artgallery.geometricalElements;

import java.util.ArrayList;

public class Vertex implements Comparable<Vertex> {

	private int id;
	private double X;
	private double Y;
	private boolean art;
	private boolean exit;
	private ArrayList<Vertex> neighbors = new ArrayList<Vertex>();
	private ArrayList<Polygon> visibilityPolygon = new ArrayList<Polygon>();

	public Vertex() {
	}

	public Vertex(double x, double y) {
		this.setId(-1);
		this.X = x;
		this.Y = y;
		this.setArt(false);
		this.setExit(false);
	}

	public Vertex(double x, double y, int id, int artFlag, int exitFlag) {
		this.X = x;
		this.Y = y;
		this.setId(id);
		this.setArt(artFlag == 1 ? true : false);
		this.setExit(exitFlag == 1 ? true : false);
	}

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
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

	public ArrayList<Vertex> getNeighbors() {
		return neighbors;
	}

	public boolean isNeighbor(Vertex v) {
		if (neighbors.contains(v)) {
			return true;
		}
		return false;
	}

	public ArrayList<Polygon> getVisibilityPolygon() {
		return visibilityPolygon;
	}

	public void setVisibilityPolygon(ArrayList<Polygon> visibilityPolygon) {
		this.visibilityPolygon = visibilityPolygon;
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof Vertex) {
			if (Math.abs(this.X - ((Vertex) o).getX()) < 0.01 && Math.abs(this.Y - ((Vertex) o).getY()) < 0.01) {
				return true;
			}
		}
		return false;
	}

	@Override
	public int compareTo(Vertex o) {
		if (this.equals(o)) {
			return 0;
		} else {
			return Double.compare(this.X, ((Vertex) o).getX());
		}
	}

	@Override
	public Vertex clone() {
		Vertex clone = new Vertex(this.X, this.Y, this.id, this.art == true ? 1 : 0, this.exit == true ? 1 : 0);
		return clone;
	}

	@Override
	public String toString() {
		return "(" + this.getX() + ":" + this.getY() + ")";
	}

	/**
	 * Get edge e_i from vertex v_i
	 *
	 * @return
	 */
	public Edge getAssociatedEdge() {
		// TODO
		throw new UnsupportedOperationException("Not implemented");
	}

	/**
	 * Get e_{i-1} from vertex v_i
	 *
	 * @return
	 */
	public Edge getPreviousAssociatedEdge() {
		throw new UnsupportedOperationException("Not implemented");
	}

}
