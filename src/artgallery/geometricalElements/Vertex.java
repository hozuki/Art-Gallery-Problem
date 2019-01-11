package artgallery.geometricalElements;

import java.util.ArrayList;

public class Vertex {

	private int X;
	private int Y;
	private boolean art;
	private boolean exit;
	private ArrayList<Edge> inEdges = new ArrayList<Edge>();
	private ArrayList<Polygon> visibilityPolygon = new ArrayList<Polygon>();


	public Vertex() {
	}

	public Vertex(int x, int y) {
		this.X = x;
		this.Y = y;
	}

	public Vertex(int x, int y, int artFlag, int exitFlag) {
		this.X = x;
		this.Y = y;
		this.setArt(artFlag == 1 ? true : false);
		this.setExit(exitFlag == 1 ? true : false);
	}

	public int getX() {
		return X;
	}

	public void setX(int x) {
		X = x;
	}

	public int getY() {
		return Y;
	}

	public void setY(int y) {
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

	public void addInEdge(Edge e) {
		inEdges.add(e);
	}

	public ArrayList<Edge> getInEdges() {
		return inEdges;
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
		assert inEdges.size() == 2;

		Edge next = getAssociatedEdge();

		Edge e1 = inEdges.get(0);
		Edge e2 = inEdges.get(1);

		if (e1 == next) {
			return e2;
		} else {
			return e1;
		}
	}

	public boolean isNeighbor(Vertex v) {
		boolean neighbor = false;
		for (Edge edge : inEdges) {
			if (edge.containsVertex(v)) {
				neighbor = true;
			}
		}
		return neighbor;
	}

	public float computeClockwiseAngle(Vertex referenceVertex) {
		float angle = (float) Math.toDegrees(Math.atan2(this.getY() - referenceVertex.getY(), this.getX() - referenceVertex.getX()));
		if (angle < 0) angle += 360;
		return -angle;
	}

	public void removeNeighbor(Vertex v) {
	}

	public ArrayList<Polygon> getVisibilityPolygon() {
		return visibilityPolygon;
	}

	public void setVisibilityPolygon(ArrayList<Polygon> visibilityPolygon) {
		this.visibilityPolygon = visibilityPolygon;
	}
}
