package artgallery;

import java.util.*;

import artgallery.dataStructures.AVLTree;
import artgallery.geometricalElements.Edge;
import artgallery.geometricalElements.Hole;
import artgallery.geometricalElements.Polygon;
import artgallery.geometricalElements.Vertex;

import javax.activation.MimeTypeParameterList;

public class GeometricAlgorithms {

	private ArrayList<Edge> triangulationEdges = new ArrayList<Edge>();
	private ArrayList<Vertex> triangulationVertices = new ArrayList<Vertex>();

	public GeometricAlgorithms() {

	}

	private ArrayList<Polygon> trapezoidalDecomposition(Polygon fullPolygon) {
		// Not implemented - this simply adds a the whole vertices and attempts
		// to triangulate for debugging purposes
		ArrayList<Vertex> clonedVertices = new ArrayList<Vertex>(fullPolygon.getVertices());
		ArrayList<Edge> clonedEdges = new ArrayList<Edge>(fullPolygon.getEdges());
		ArrayList<Hole> clonedHoles = new ArrayList<Hole>(fullPolygon.getHoles());
		Polygon p = new Polygon(clonedVertices, clonedEdges, clonedHoles);

		// Decompose fullPolygon into array of trapezoids
		ArrayList<Polygon> trapezoids = new ArrayList<Polygon>();
		ArrayList<Polygon> monotonePolygons = new ArrayList<Polygon>();
		ArrayList<Polygon> triangulationMesh = new ArrayList<Polygon>();

		// Call monotone polygonization on said trapezoids
		trapezoids.add(p);
		for (Polygon trapezoid : trapezoids) {
			monotonePolygons.add(monotonePolygonization(trapezoid));
		}

		for (Polygon monotonePolygon : monotonePolygons) {
			triangulationMesh.addAll(triangulateMonotonePolygon(monotonePolygon));
		}

		return triangulationMesh;
	}

	private Polygon monotonePolygonization(Polygon trapezoid) {
		return trapezoid;
	}

	/**
	 * Triangulate using the method: {@link #makeMonotonePolygons(Polygon)} then
	 * {@link #triangulateMonotonePolygon(Polygon)}. Input polygon can have holes.
	 *
	 * @param polygon The simple polygon to triangulate
	 * @return Triangle mesh
	 */
	public ArrayList<Polygon> triangulateSimplePolygon(final Polygon polygon) {
		final ArrayList<Polygon> monotones = makeMonotonePolygons(polygon);
		final ArrayList<Polygon> triangles = new ArrayList<>();

		for (final Polygon poly : monotones) {
			final ArrayList<Polygon> tris = triangulateMonotonePolygon(poly);

			triangles.addAll(tris);
		}

		return triangles;
	}

	/**
	 * Break down a polygon into one or more monotone polygons.
	 *
	 * @param polygon The polygon to break down
	 * @return A list of monotone polygons
	 */
	public ArrayList<Polygon> makeMonotonePolygons(final Polygon polygon) {
		PriorityQueue<Vertex> q = new PriorityQueue<>(GeometricAlgorithms::vertexComparator);

		q.addAll(polygon.getVertices());

		for (Hole hole : polygon.getHoles()) {
			q.addAll(hole.getVertices());
		}

		AVLTree<Edge> t = new AVLTree<>();

		helper = new HashMap<>();

		ArrayList<Polygon> result = new ArrayList<>();

		while (!q.isEmpty()) {
			Vertex vi = q.poll();
			VertexType vertexType = judgeVertexType(vi);

			switch (vertexType) {
				case REGULAR: {
					// TODO: compute whether the interior is to the right of vi
					boolean interiorIsToTheRight = Math.random() > 0.5;

					if (interiorIsToTheRight) {
						Edge ei_m1 = vi.getPreviousAssociatedEdge();

						assert helper.containsKey(ei_m1);

						Vertex v = helper.get(ei_m1);
						VertexType vt = judgeVertexType(v);

						if (vt == VertexType.MERGE) {
							// TODO: insert diagonal from ei_m1 to vi
						}

						t.remove(ei_m1);

						Edge ei = vi.getAssociatedEdge();

						t.insert(ei);
						helper.put(ei, vi);
					} else {
						Edge ej = null;
						// TODO: search for ej

						Vertex v = helper.get(ej);
						VertexType vt = judgeVertexType(v);

						if (vt == VertexType.MERGE) {
							// TODO: insert diagonal from vi to v
						}

						helper.put(ej, vi);
					}
				}
				break;
				case START: {
					Edge ei = vi.getAssociatedEdge();
					t.insert(ei);
					helper.put(ei, vi);
				}
				break;
				case END: {
					Edge ei_m1 = vi.getPreviousAssociatedEdge();

					assert helper.containsKey(ei_m1);

					Vertex v = helper.get(ei_m1);
					VertexType vt = judgeVertexType(v);

					if (vt == VertexType.MERGE) {
						// TODO: insert diagonal from ei_m1 to vi
					}

					t.remove(ei_m1);
				}
				break;
				case MERGE: {
					Edge ei_m1 = vi.getPreviousAssociatedEdge();

					assert helper.containsKey(ei_m1);

					Vertex v = helper.get(ei_m1);
					VertexType vt = judgeVertexType(v);

					if (vt == VertexType.MERGE) {
						// TODO: insert diagonal from ei_m1 to vi
					}

					t.remove(ei_m1);

					Edge ej = null;
					// TODO: search for ej

					v = helper.get(ej);
					vt = judgeVertexType(v);

					if (vt == VertexType.MERGE) {
						// TODO: insert diagonal from vi to v
					}

					helper.put(ej, vi);
				}
				break;
				case SPLIT: {
					Edge ej = null;
					// TODO: search for ej

					Vertex dest = helper.get(ej);
					// TODO: insert diagonal from vi to dest

					helper.put(ej, vi);

					Edge ei = vi.getAssociatedEdge();

					t.insert(ei);
					helper.put(ei, vi);
				}
				break;
			}
		}

		return result;
	}

	private Map<Edge, Vertex> helper;

	private static VertexType judgeVertexType(Vertex vi) {
		ArrayList<Edge> edges = vi.getInEdges();

		assert edges.size() == 2;

		// Note: Here it is assumed that ref edges are stored in CCW order
		Edge e1 = edges.get(0);
		Edge e2 = edges.get(1);

		Vertex v1 = e1.getTheOtherVertexOf(vi);
		Vertex v2 = e2.getTheOtherVertexOf(vi);

		boolean allAbove = v1.getY() > vi.getY() && v2.getY() > vi.getY();
		boolean allBelow = v1.getY() < vi.getY() && v2.getY() < vi.getY();

		if (!allAbove && !allBelow) {
			return VertexType.REGULAR;
		}

		int x1 = v1.getX() - vi.getX();
		int y1 = v1.getY() - vi.getY();
		int x2 = v2.getX() - vi.getX();
		int y2 = v2.getY() - vi.getY();

		double a1 = Math.atan2(y1, x1);
		double a2 = Math.atan2(y2, x2);

		assert a1 != a2;

		if (allAbove) {
			if (a2 - a1 > 0) {
				return VertexType.MERGE;
			} else if (a2 - a1 < 0) {
				return VertexType.END;
			}
		} else if (allBelow) {
			if (a2 - a1 > 0) {
				return VertexType.SPLIT;
			} else if (a2 - a1 < 0) {
				return VertexType.START;
			}
		}

		throw new RuntimeException("Oops unexpected");
	}

	private enum VertexType {
		REGULAR, START, END, MERGE, SPLIT
	}

	private static int vertexComparator(final Vertex v1, final Vertex v2) {
		int result = v1.getY() - v2.getY();

		if (result != 0) {
			return result;
		} else {
			return v1.getX() - v2.getX();
		}
	}

	/*
	 * Implementation of "TriangulateMonotonePolygon(polygon P as DCEL)". From
	 * the course slides "3 - Art Gallery Triangulation", page 158. Input:
	 * Y-Monotone polygon "p". Output: ArrayList of triangular polygons.
	 */
	public ArrayList<Polygon> triangulateMonotonePolygon(Polygon p) {
		ArrayList<Vertex> polygonVertices = new ArrayList<Vertex>(p.getVertices());

		// Sorts the vertices on Y-Descending order.
		Collections.sort(polygonVertices, (v1, v2) -> v2.getY() - v1.getY());
		p.constructChains();
		triangulationVertices = polygonVertices;
		// Implementation translation
		Stack<Vertex> s = new Stack<Vertex>();
		s.push(polygonVertices.get(0));
		s.push(polygonVertices.get(1));
		for (int j = 2; j < polygonVertices.size() - 1; ++j) {
			Vertex uj = polygonVertices.get(j);
			// System.out.println(polygonVertices.indexOf(s.peek())+ " is
			// neighbor of "+ polygonVertices.indexOf(uj) + " : " +
			// s.peek().isNeighbor(uj));
			if (!s.peek().isNeighbor(uj)) {
				while (s.size() > 1) {
					Vertex v = s.pop();
					triangulationEdges.add(new Edge(uj, v));
				}
				s.pop();
				s.push(polygonVertices.get(j - 1));
				s.push(uj);
			} else {
				Vertex v = s.pop();
				// Must implement "sees"
				while (s.size() > 0) {
					v = s.pop();
					triangulationEdges.add(new Edge(uj, v));
				}
				s.push(v);
				s.push(uj);
			}
		}
		s.pop();
		while (s.size() > 1) {
			Vertex v = s.pop();
			triangulationEdges.add(new Edge(polygonVertices.get(polygonVertices.size() - 1), v));
		}

		ArrayList<Polygon> triangles = new ArrayList<Polygon>();
		triangulationEdges.addAll(p.getEdges());
		constructTriangles(triangles, triangulationVertices, triangulationEdges);

		// System.out.print(triangles.size());
		return triangles;
	}

	private void constructTriangles(ArrayList<Polygon> mesh, ArrayList<Vertex> vertices, ArrayList<Edge> edges) {
		for (Edge edge : edges) {
			for (Vertex vertex : vertices) {
				if (!vertex.equals(edge.getFirstVertex()) && !vertex.equals(edge.getSecondVertex())) {
					if (vertex.isNeighbor(edge.getFirstVertex()) && vertex.isNeighbor(edge.getSecondVertex())) {

						Edge edge1 = edges.stream()
								.filter(e -> e.containsVertex(vertex) && e.containsVertex(edge.getFirstVertex()))
								.findAny().orElse(null);
						Edge edge2 = edges.stream()
								.filter(e -> e.containsVertex(vertex) && e.containsVertex(edge.getSecondVertex()))
								.findAny().orElse(null);

						edge.getFirstVertex().removeNeighbor(vertex);
						edge.getSecondVertex().removeNeighbor(vertex);

						ArrayList<Edge> triangleEdges = new ArrayList<Edge>();
						triangleEdges.add(edge);
						triangleEdges.add(edge1);
						triangleEdges.add(edge2);
						Polygon newTriangle = new Polygon(triangleEdges);
						if (!mesh.stream().anyMatch(t -> t.edgeMatch(newTriangle))) {
							mesh.add(newTriangle);
						}
					}
				}
			}
		}
	}

	public ArrayList<Polygon> computeVisibilityPolygon(Vertex v, Polygon p) {
		//Array of triangles comprising the visibility polygon
		ArrayList<Polygon> visibilityPolygon = new ArrayList<Polygon>();

		ArrayList<Vertex> visibleVertices = visibleVertices(v, p);

		visibilityPolygon.add(p);
		//Gotta construct the triangles somehow...
		for (Vertex w : visibleVertices) {

		}

		return visibilityPolygon;
	}

	/*
	 * Implementation of "VisibleVertices(P, S)". From the book
	 * "Computational Geometry: Algorithms and Applications - Third Edition" by
	 * M. Berg , page 328. Input: Point P, set of polygonal obstacles S. Output:
	 * The set of visible vertices from P.
	 */
	private ArrayList<Vertex> visibleVertices(Vertex v, Polygon p) {
		ArrayList<Vertex> obstacleVertices = p.getVerticesAndHoles();

		obstacleVertices.sort((v1, v2) -> compareAngleAndProximity(v, v1, v2));

		AVLTree<Edge> tree = new AVLTree<Edge>();

		Edge sweepLine = new Edge(v, new Vertex((int) (Math.ceil(p.getMaxDistance())), v.getY()));

		ArrayList<Edge> intersectedEdges = new ArrayList<Edge>();
		for (Edge obstacleEdge : p.getEdges()) {
			if (edgesIntersect(sweepLine, obstacleEdge)) {
				intersectedEdges.add(obstacleEdge);
			}
		}
		intersectedEdges.sort((e1, e2) -> compareIntersectionDistance(v, sweepLine, e1, e2));
		for (Edge edge : intersectedEdges) {
			tree.insert(edge);
		}

		ArrayList<Vertex> visibleVertices = new ArrayList<Vertex>();

		for (int i = 0; i < obstacleVertices.size(); ++i) {
			Vertex vi = obstacleVertices.get(i);
			if (isVisible(obstacleVertices, i)) {
				visibleVertices.add(vi);
				ArrayList<Edge> incidentEdges = vi.getInEdges();
				for (Edge e : incidentEdges) {
					//Get the other vertex (not "vi") of the incident edge
					Vertex otherVertex = e.getFirstVertex().equals(vi) ? e.getSecondVertex() : e.getFirstVertex();
					//If the orientation if clockwise w.r.t. the sweepline, add the edge to the tree.
					if (orientation(v, vi, otherVertex) == 1) {
						tree.insert(e);
					}
					//If the orientation if counter-clockwise w.r.t. the sweepline, remove the edge to the tree.
					if (orientation(v, vi, otherVertex) == 2) {
						tree.remove(e);
					}
				}
			}
		}
		return visibleVertices;
	}

	/*
	 * Implementation of "Visible(W)" subroutine. From the book
	 * "Computational Geometry: Algorithms and Applications - Third Edition" by
	 * M. Berg , page 329.
	 * Input: Point W.
	 * Output: The set of visible vertices from P.
	 */
	private boolean isVisible(ArrayList<Vertex> w, int i) {
		/*
		 if() {
		 
			return false;
		}else if(i == 0 || w.get(i - 1)) {
			
		}
		*/
		return true;
	}

	private int compareAngleAndProximity(Vertex reference, Vertex v1, Vertex v2) {
		int result = Float.compare(v1.computeClockwiseAngle(reference), v2.computeClockwiseAngle(reference));
		if (result == 0) {
			Double d1 = Math.sqrt(Math.pow(v1.getY() - reference.getY(), 2) + Math.pow((v1.getX() - reference.getX()), 2));
			Double d2 = Math.sqrt(Math.pow(v2.getY() - reference.getY(), 2) + Math.pow((v2.getX() - reference.getX()), 2));
			return Double.compare(d1, d2);
		}
		return result;
	}

	private int compareIntersectionDistance(Vertex v, Edge e, Edge e1, Edge e2) {
		return Double.compare(getIntersectionDistance(v, e, e1), getIntersectionDistance(v, e, e2));
	}

	private double getIntersectionDistance(Vertex v, Edge e1, Edge e2) {
		int x1, x2, x3, x4, y1, y2, y3, y4;

		x1 = e1.getFirstVertex().getX();
		y1 = e1.getFirstVertex().getY();
		x2 = e1.getSecondVertex().getX();
		y2 = e1.getSecondVertex().getY();
		x3 = e2.getFirstVertex().getX();
		y3 = e2.getFirstVertex().getY();
		x4 = e2.getSecondVertex().getX();
		y4 = e2.getSecondVertex().getY();

		int x = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
		int y = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));

		return Math.sqrt(Math.pow(y - v.getY(), 2) + Math.pow((x - v.getX()), 2));
	}

	private boolean isVisible(Vertex viewer, ArrayList<Vertex> vertices) {
		return false;
	}

	private Polygon computeBoundingBox(Polygon p) {
		int minX, minY, maxX, maxY;
		minX = minY = Integer.MAX_VALUE;
		maxX = maxY = Integer.MIN_VALUE;

		for (Vertex v : p.getVertices()) {
			if (v.getX() > maxX) maxX = v.getX();
			if (v.getX() < minX) minX = v.getX();
			if (v.getY() > maxY) maxY = v.getY();
			if (v.getY() < minY) minY = v.getY();
		}
		ArrayList<Vertex> vertices = new ArrayList<Vertex>();
		vertices.add(new Vertex(maxX, maxY));
		vertices.add(new Vertex(maxX, minY));
		vertices.add(new Vertex(minX, minY));
		vertices.add(new Vertex(minX, maxY));
		Polygon boundingBox = new Polygon(vertices);
		return boundingBox;
	}

	private boolean insidePolygon(Vertex v, Polygon p) {
		return false;
	}

	private boolean edgesIntersect(Edge e1, Edge e2) {
		return intersect(e1.getFirstVertex(), e1.getSecondVertex(), e2.getFirstVertex(), e2.getSecondVertex());
	}

	/*
	 * Adapted from GeeksforGeeks
	 * "How to check if two given line segments intersect?"
	 * https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
	 */
	private boolean intersect(Vertex p1, Vertex q1, Vertex p2, Vertex q2) {
		// Find the four orientations needed for general and special cases
		int o1 = orientation(p1, q1, p2);
		int o2 = orientation(p1, q1, q2);
		int o3 = orientation(p2, q2, p1);
		int o4 = orientation(p2, q2, q1);

		// General case
		if (o1 != o2 && o3 != o4)
			return true;

		// Special Cases
		// p1, q1 and p2 are colinear and p2 lies on segment p1q1
		if (o1 == 0 && areColinear(p1, p2, q1))
			return true;

		// p1, q1 and q2 are colinear and q2 lies on segment p1q1
		if (o2 == 0 && areColinear(p1, q2, q1))
			return true;

		// p2, q2 and p1 are colinear and p1 lies on segment p2q2
		if (o3 == 0 && areColinear(p2, p1, q2))
			return true;

		// p2, q2 and q1 are colinear and q1 lies on segment p2q2
		if (o4 == 0 && areColinear(p2, q1, q2))
			return true;

		return false; // Doesn't fall in any of the above cases
	}

	private boolean areColinear(Vertex v1, Vertex v2, Vertex v3) {
		int area = v1.getX() * (v2.getY() - v3.getY()) + v2.getX() * (v3.getY() - v1.getY())
				+ v3.getX() * (v1.getY() - v2.getY());
		return (area == 0);
	}

	private int orientation(Vertex p, Vertex q, Vertex r) {
		int val = (q.getY() - p.getY()) * (r.getX() - q.getX()) - (q.getX() - p.getX()) * (r.getY() - q.getY());
		if (val == 0) {
			return 0; // colinear
		} else {
			return (val > 0) ? 1 : 2; // clock or counterclock wise
		}

	}
}
