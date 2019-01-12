package artgallery.geometricalElements;

import java.util.ArrayList;
import java.util.Collections;

public class Hole {
	private ArrayList<Vertex> vertices;
	private ArrayList<Edge> edges;
	private int size;
	private Polygon polygon;

	public Hole(Polygon polygon) {
		this.polygon = polygon;
		this.vertices = new ArrayList<Vertex>();
		edges = new ArrayList<Edge>();
	}

	public Hole(Polygon polygon, ArrayList<Vertex> vertices) {
		this.polygon = polygon;
		this.setVertices(vertices);
		edges = new ArrayList<Edge>();
		for (int i = 0; i < vertices.size(); ++i) {
			Vertex v1 = vertices.get(i);
			Vertex v2 = vertices.get((i + 1) % vertices.size());
			edges.add(new Edge(polygon, v1, v2));
		}
	}

	public Polygon getPolygon() {
		return polygon;
	}

	public ArrayList<Vertex> getVertices() {
		return vertices;
	}

	public void setVertices(ArrayList<Vertex> vertices) {
		this.vertices = new ArrayList<>(vertices);

		for (Vertex v : vertices) {
			v.setOnBoundary(false);
		}
	}

	public void addVertex(Vertex v) {
		if (this.vertices.size() < this.size) {
			this.vertices.add(v);
			v.setOnBoundary(false);
		}
		if (this.vertices.size() == this.size) {
			for (int i = 0; i < vertices.size(); ++i) {
				Vertex v1 = vertices.get(i);
				Vertex v2 = vertices.get((i + 1) % vertices.size());
				Edge edge = new Edge(polygon, v1, v2);
				edges.add(edge);
			}
		}
	}

	public int getSize() {
		return size;
	}

	public void setSize(int size) {
		this.size = size;
	}

	public ArrayList<Edge> getEdges() {
		return edges;
	}

	public void setEdges(ArrayList<Edge> edges) {
		this.edges = edges;
	}

	public void reverseVerticesAndEdges() {
		Collections.reverse(vertices);
		Collections.reverse(edges);
	}

	public void generateEdges() {
		for (int i = 0; i < edges.size(); ++i) {
			vertices.get(i).setOutEdge(edges.get(i));
			vertices.get(i).setInEdge(edges.get((i + 1 + edges.size()) % edges.size()));
		}
	}

}
