package artgallery.geometricalElements;

import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import artgallery.GeometricAlgorithms;

public class Polygon {

	private ArrayList<Vertex> vertices = new ArrayList<Vertex>();
	private ArrayList<Edge> edges = new ArrayList<Edge>();
	private ArrayList<Hole> holes = new ArrayList<Hole>();
	private ArrayList<Polygon> triangulation = new ArrayList<>();
	private Map<Vertex, Integer> chainMap = new HashMap<Vertex, Integer>();

	public Polygon() {
	}

	@SuppressWarnings("unchecked")
	public Polygon(ArrayList<?> elements) {
		if (elements.get(0) instanceof Vertex) {
			this.vertices = (ArrayList<Vertex>) elements;
			for (Vertex v : (ArrayList<Vertex>) elements) {
				v.setOnBoundary(true);
			}
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
		for (Vertex v : vertices) {
			v.setOnBoundary(true);
		}
		this.holes = holes;
		generateEdges();
	}

	public Polygon(ArrayList<Vertex> vertices, ArrayList<Edge> edges, ArrayList<Hole> holes) {
		this.vertices = vertices;
		for (Vertex v : vertices) {
			v.setOnBoundary(true);
		}
		this.edges = edges;
		this.holes = holes;
	}

	public ArrayList<Vertex> getVertices() {
		return vertices;
	}

	public ArrayList<Vertex> getVerticesAndHoles() {
		ArrayList<Vertex> allVertices = new ArrayList<Vertex>(vertices);

		for (Hole h : holes) {
			allVertices.addAll(h.getVertices());
		}

		return allVertices;
	}

	public void setVertices(ArrayList<Vertex> vertices) {
		this.vertices = vertices;

		for (Vertex v : vertices) {
			v.setOnBoundary(true);
		}
	}

	public ArrayList<Edge> getEdges() {
		return edges;
	}

	public ArrayList<Edge> getEdgesAndHoles() {
		ArrayList<Edge> allEdges = new ArrayList<Edge>(edges);
		for (Hole h : holes) {
			allEdges.addAll(new ArrayList<Edge>(h.getEdges()));
		}
		return allEdges;
	}

	public Vertex getMatchingVertex(Vertex vertex) {
		for (Vertex localVertex : this.getVerticesAndHoles()) {
			if (localVertex.equals(vertex)) {
				return localVertex;
			}
		}
		return null;

	}

	public void generateEdges() {
		//Construct the edges for the gallery bounds by connecting points counter-clockwise.
		for (int i = 0; i < vertices.size(); ++i) {
			Vertex v1 = vertices.get(i);
			Vertex v2 = vertices.get((i + 1) % vertices.size());
			this.edges.add(new Edge(this, v1, v2));
			v1.addNeighbor(v2);
			v2.addNeighbor(v1);
		}

		//Connect edges to their neighbors
		for (int i = 0; i < edges.size(); ++i) {
			Edge e1 = edges.get(i);
			Edge e2 = edges.get((i + 1) % edges.size());
			e1.addNeighbor(e2);
			e2.addNeighbor(e1);
		}

		for (int i = 0; i < edges.size(); ++i) {
			vertices.get(i).setOutEdge(edges.get(i));
			vertices.get(i).setInEdge(edges.get((i - 1 + edges.size()) % edges.size()));
		}

		for (Hole h : holes) {
			h.generateEdges();
		}
	}

	public void fixVertexNeighbors() {
		vertices.get(0).swapNeighbors();

		for (Hole h : holes) {
			ArrayList<Vertex> vs = h.getVertices();

			for (int i = 0; i < vs.size() - 1; i++) {
				vs.get(i).swapNeighbors();
			}
		}
	}

	public void computeTriangulation() {
		if (this.triangulation.isEmpty()) {
			GeometricAlgorithms GA = new GeometricAlgorithms();
			triangulation = GA.triangulateSimplePolygon(this);
		}
	}

	public void computeVisibility(Vertex v) {
		GeometricAlgorithms GA = new GeometricAlgorithms();
		v.setVisibilityPolygon(GA.computeVisibilityPolygon(v, this));
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
				double tempDistance = Math
						.sqrt(Math.pow(v1.getY() - v2.getY(), 2) + Math.pow((v1.getX() - v2.getX()), 2));
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
		double lastY = Integer.MAX_VALUE;
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

	public Rectangle2D getBoundingBox() {
		double minX = Double.MAX_VALUE, minY = Double.MAX_VALUE;
		double maxX = -Double.MAX_VALUE, maxY = -Double.MAX_VALUE;

		for (Vertex v : getVertices()) {
			if (v.getX() < minX) {
				minX = v.getX();
			}
			if (v.getY() < minY) {
				minY = v.getY();
			}
			if (v.getX() > maxX) {
				maxX = v.getX();
			}
			if (v.getY() > maxY) {
				maxY = v.getY();
			}
		}

		return new Rectangle2D.Double(minX, minY, maxX - minX, maxY - minY);
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) {
			return true;
		}
		if (!(o instanceof Polygon)) {
			return false;
		}
		if (o instanceof Polygon) {
			boolean equal = true;
			boolean contained = false;
			for (Edge e1 : this.getEdges()) {
				for (Edge e2 : ((Polygon) o).getEdges()) {
					if (e2.equals(e1))
						contained = true;
				}
				equal = equal & contained;
			}
			for (Edge e2 : ((Polygon) o).getEdges()) {
				contained = false;
				for (Edge e1 : this.getEdges()) {
					if (e1.equals(e2))
						contained = true;
				}
				equal = equal & contained;
			}
		}
		return false;
	}

	@Override
	public Polygon clone() {
		Polygon clone = new Polygon();
		ArrayList<Vertex> cloneVertices = new ArrayList<Vertex>();
		for (Vertex v : this.getVertices()) {
			cloneVertices.add(v.clone());
		}
		clone.setVertices(cloneVertices);
		clone.setHoles(new ArrayList<Hole>(this.getHoles()));
		clone.generateEdges();
		return clone;
	}

}