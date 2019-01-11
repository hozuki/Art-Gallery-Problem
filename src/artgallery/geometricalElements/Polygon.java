package artgallery.geometricalElements;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import artgallery.GeometricAlgorithms;

public class Polygon {

	private ArrayList<Vertex> vertices = new ArrayList<Vertex>();
	private ArrayList<Edge> edges = new ArrayList<Edge>();
	private ArrayList<Hole> holes = new ArrayList<Hole>();
	private ArrayList<Polygon> triangulation = null;
	private Map<Vertex, Integer> chainMap = new HashMap<Vertex, Integer>();
	private boolean visibilityComputed;

	public Polygon() {
	}

	@SuppressWarnings("unchecked")
	public Polygon(ArrayList<?> elements) {
		if (elements.get(0) instanceof Vertex) {
			this.vertices = (ArrayList<Vertex>) elements;
			generateEdges();
		} else if (elements.get(0) instanceof Edge) {
			this.edges = (ArrayList<Edge>) elements;
			for (Edge e : this.edges) {
				if (!this.vertices.contains(e.getFirstVertex())) {
					this.vertices.add(e.getFirstVertex());
				}
				if (!this.vertices.contains(e.getSecondVertex())) {
					this.vertices.add(e.getSecondVertex());
				}
			}
		}
	}

	public Polygon(ArrayList<Vertex> vertices, ArrayList<Hole> holes) {
		this.vertices = vertices;
		this.holes = holes;
		generateEdges();
	}

	public Polygon(ArrayList<Vertex> vertices, ArrayList<Edge> edges, ArrayList<Hole> holes) {
		this.vertices = vertices;
		this.edges = edges;
		this.holes = holes;
	}

	public ArrayList<Vertex> getVertices() {
		return vertices;
	}

	public ArrayList<Vertex> getVerticesAndHoles() {
		ArrayList<Vertex> allVertices = new ArrayList<Vertex>(vertices);
		for (Hole h : holes) {
			allVertices.addAll(new ArrayList<Vertex>(h.getVertices()));
		}
		return allVertices;
	}

	public void setVertices(ArrayList<Vertex> vertices) {
		this.vertices = vertices;
	}

	public ArrayList<Edge> getEdges() {
		return edges;
	}

	public void generateEdges() {
		for (int i = 0; i < vertices.size(); ++i) {
			Vertex v1 = vertices.get(i);
			Vertex v2 = vertices.get((i + 1) % vertices.size());
			Edge edge = new Edge(v1, v2);
			this.edges.add(edge);
			v1.addInEdge(edge);
			v2.addInEdge(edge);
		}
	}

	public void computeTriangulation() {
		if (this.triangulation == null) {
			GeometricAlgorithms GA = new GeometricAlgorithms();
			triangulation = GA.triangulateMonotonePolygon(this);
		}
	}

	public void computeVisibility() {
		if (!visibilityComputed) {
			GeometricAlgorithms GA = new GeometricAlgorithms();

			for (Vertex v : vertices) {
				ArrayList<Polygon> visibilityPolygon = GA.computeVisibilityPolygon(v, this);
				v.setVisibilityPolygon(visibilityPolygon);
			}

			visibilityComputed = true;
		}
	}

	public ArrayList<Polygon> getTriangulation() {
		if (triangulation == null) {
			computeTriangulation();
		}

		return triangulation;
	}

	public double getMaxDistance() {
		double maxDistance = 0;
		for (Vertex v1 : vertices) {
			for (Vertex v2 : vertices) {
				double tempDistance = Math.sqrt(Math.pow(v1.getY() - v2.getY(), 2) + Math.pow((v1.getX() - v2.getX()), 2));
				if (tempDistance > maxDistance) {
					maxDistance = tempDistance;
				}
			}
		}
		return maxDistance;
	}

	// Bad implementation of monotone-chains.
	public void constructChains() {
		tiltHorizontals();

		int top = 0;
		for (int i = 0; i < vertices.size() - 1; ++i) {
			if (vertices.get(i).getY() > vertices.get(top).getY()) {
				top = i;
			}
		}
	}

	// Modified the vertices positions to avoid horizontal lines in the polygon.
	public void tiltHorizontals() {
		int lastY = Integer.MAX_VALUE;
		for (Vertex v : vertices) {
			if (lastY == Integer.MAX_VALUE) {
				lastY = v.getY();
			} else if (lastY == v.getY()) {
				v.setY(v.getY() + 1);
				lastY = v.getY();
			}
		}
	}

	public int getChain(Vertex v) {
		return chainMap.get(v);
	}

	public ArrayList<Hole> getHoles() {
		return holes;
	}

	public void setHoles(ArrayList<Hole> holes) {
		this.holes = holes;
	}

	public int[] getCentroid() {
		int[] centroid = {0, 0};
		for (Vertex v : vertices) {
			centroid[0] += v.getX();
			centroid[1] += v.getY();
		}
		centroid[0] /= vertices.size();
		centroid[1] /= vertices.size();
		return centroid;
	}

	public boolean edgeMatch(Polygon p) {
		if (this.getEdges().containsAll(p.getEdges()) && p.getEdges().containsAll(this.getEdges())) {
			return true;
		} else {
			return false;
		}
	}


}