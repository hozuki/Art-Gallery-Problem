package artgallery.geometricalElements;

import java.util.ArrayList;

public class Edge implements Comparable<Edge> {
	private String id;
	private Vertex firstVertex;
	private Vertex secondVertex;

	public Edge() {
	}

	public Edge(Vertex v1, Vertex v2) {
		setFirstVertex(v1);
		setSecondVertex(v2);
		setId(v1 + "-" + v2);
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public Vertex getFirstVertex() {
		return firstVertex;
	}

	public void setFirstVertex(Vertex firstVertex) {
		this.firstVertex = firstVertex;
		this.firstVertex.addInEdge(this);
	}

	public Vertex getSecondVertex() {
		return secondVertex;
	}

	public void setSecondVertex(Vertex secondVertex) {
		this.secondVertex = secondVertex;
		this.secondVertex.addInEdge(this);
	}

	public boolean containsVertex(Vertex v) {
		if (firstVertex.equals(v) || secondVertex.equals(v)) {
			return true;
		}
		return false;
	}

	public Vertex getTheOtherVertexOf(Vertex v) {
		if (!containsVertex(v)) {
			throw new IllegalArgumentException("The vertex is not in the edge");
		}

		if (firstVertex.equals(v)) {
			return secondVertex;
		} else if (secondVertex.equals(v)) {
			return firstVertex;
		}

		throw new RuntimeException("Oops");
	}

	public Edge getPreviousEdgeLinkedBy(Vertex v) {
		if (!containsVertex(v)) {
			throw new IllegalArgumentException("The vertex is not in the edge");
		}

		ArrayList<Edge> edges = v.getInEdges();

		assert edges.contains(this);
		assert edges.size() == 2;

		Edge e1 = edges.get(0);
		Edge e2 = edges.get(1);

		if (this == e1) {
			return e2;
		} else if (this == e2) {
			return e1;
		}

		throw new RuntimeException("Oops");
	}

	// TODO: compare according to index in polygon
	@Override
	public int compareTo(Edge edge) {
		return this.getId().compareTo(edge.getId());
	}

}
