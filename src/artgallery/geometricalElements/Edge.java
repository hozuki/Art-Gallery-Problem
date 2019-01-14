package artgallery.geometricalElements;

import artgallery.Util;

import java.util.ArrayList;
import java.util.Objects;

public class Edge implements Comparable<Edge> {

	private String id;
	private Vertex firstVertex;
	private Vertex secondVertex;
	private ArrayList<Edge> neighbors = new ArrayList<Edge>();
	private Polygon polygon;

	public Edge(Polygon polygon) {
		this.polygon = polygon;
		setFirstVertex(new Vertex(polygon));
		setSecondVertex(new Vertex(polygon));
	}

	public Edge(Polygon polygon, Vertex v1, Vertex v2) {
		this.polygon = polygon;
		setFirstVertex(v1);
		setSecondVertex(v2);
		setId(v1.getId() + "-" + v2.getId());
		v1.addNeighbor(v2);
		v2.addNeighbor(v1);
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public double length() {
		return Math.sqrt(
			Math.pow(firstVertex.getY() - secondVertex.getY(), 2) +
				Math.pow((firstVertex.getX() - secondVertex.getX()), 2));
	}

	public Vertex getFirstVertex() {
		return firstVertex;
	}

	public void setFirstVertex(Vertex firstVertex) {
		assert firstVertex.getPolygon() == getPolygon();
		this.firstVertex = firstVertex;
	}

	public Vertex getSecondVertex() {
		return secondVertex;
	}

	public void setSecondVertex(Vertex secondVertex) {
		assert secondVertex.getPolygon() == getPolygon();
		this.secondVertex = secondVertex;
	}

	public Vertex getOtherVertex(Vertex v) {
		if (this.containsVertex(v)) {
			return firstVertex.equals(v) ? secondVertex : firstVertex;
		}
		return null;
	}

	public boolean containsVertex(Vertex v) {
		if (firstVertex.equals(v) || secondVertex.equals(v)) {
			return true;
		}
		return false;
	}

	public Polygon getPolygon() {
		return polygon;
	}

	public boolean isHorizontal() {
		return Util.equals(firstVertex.getY(), secondVertex.getY());
	}

	public boolean isVertical() {
		return Util.equals(firstVertex.getX(), secondVertex.getX());
	}

	public double getSlope() {
		if (isVertical()) {
			return firstVertex.getY() < secondVertex.getY() ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
		} else {
			return (secondVertex.getY() - firstVertex.getY()) / (secondVertex.getX() - firstVertex.getX());
		}
	}

	public void addNeighbor(Edge e) {
		if (!neighbors.contains(e)) {
			neighbors.add(e);
		}
	}

	public void removeNeighbor(Edge e) {
		if (neighbors.contains(e)) {
			e.neighbors.remove(this);
			neighbors.remove(e);
		}
	}

	public ArrayList<Edge> getNeighbors() {
		return neighbors;
	}

	public boolean isNeighbor(Edge e) {
		if (neighbors.contains(e)) {
			return true;
		}
		return false;
	}

	public Vertex getMidpoint() {
		double midX = firstVertex.getX() + (secondVertex.getX() - firstVertex.getX()) / 2;
		double midY = firstVertex.getY() + (secondVertex.getY() - firstVertex.getY()) / 2;
		return new Vertex(polygon, midX, midY);
	}

	public Vertex getUpperVertex() {
		return firstVertex.getY() > secondVertex.getY() ? firstVertex : secondVertex;
	}

	public Vertex getLowerVertex() {
		return firstVertex.getY() < secondVertex.getY() ? firstVertex : secondVertex;
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof Edge)) {
			return false;
		}
		if (o instanceof Edge) {
			if (this.containsVertex(((Edge) o).getFirstVertex()) && this.containsVertex(((Edge) o).getSecondVertex())) {
				return true;
			}
		}
		return false;
	}

	@Override
	public int hashCode() {
		return Objects.hash(firstVertex, secondVertex);
	}

	@Override
	public int compareTo(Edge edge) {
		if (this.equals(edge)) {
			return 0;
		} else {
			int result = Integer.compare(firstVertex.getId(), edge.firstVertex.getId());

			if (result == 0) {
				return result;
			} else {
				return Integer.compare(secondVertex.getId(), edge.secondVertex.getId());
			}
		}
	}

	@Override
	public String toString() {
		return String.format("%s-%s", firstVertex.toString(), secondVertex.toString());
	}

}
