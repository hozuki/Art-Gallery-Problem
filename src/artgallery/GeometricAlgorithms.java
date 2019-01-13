package artgallery;

import java.awt.geom.Rectangle2D;
import java.util.*;
import java.util.stream.Collectors;

import artgallery.dataStructures.AVLTree;
import artgallery.geometricalElements.Edge;
import artgallery.geometricalElements.Hole;
import artgallery.geometricalElements.Polygon;
import artgallery.geometricalElements.Vertex;

public final class GeometricAlgorithms {

	public ArrayList<Edge> trapezoidalEdges = new ArrayList<Edge>();

	public GeometricAlgorithms() {
	}

	// Public method to trigger the triangulation of a complex polygon through
	// the method proposed by Seidel here: http://gamma.cs.unc.edu/SEIDEL/#FIG
	// (Input) Complex Polygon -> Trapezoidal Decomposition -> Monotone Polygons
	// -> (Output) Monotonic Triangulation.
	public ArrayList<Polygon> computeTriangulation(Polygon polygon) {
		// Create a deep copy of the polygon object to avoid modifying the
		// original gallery.
		Polygon polygonCopy = polygon.clone();

		// Decompose the complex (copied) gallery polygon into trapezoids.
		ArrayList<Polygon> trapezoids = trapezoidalDecomposition(polygonCopy);
		// Obtain a set of monotone polygons out of said trapezoids.
		ArrayList<Polygon> monotonePolygons = monotonePolygonization(trapezoids);
		// Triangulate each of said monotone polygons.
		ArrayList<Polygon> triangulation = new ArrayList<Polygon>();
		for (Polygon mp : monotonePolygons) {
			triangulation.addAll(triangulateMonotonePolygon(mp));
		}

		// Finally return the list containing a triangle (Polygon) for each
		// triangle in the triangulation.
		return triangulation;
	}

	// Not implemented - for now simply returning the full polygon as a
	// "trapezoid".
	public ArrayList<Polygon> trapezoidalDecomposition(Polygon fullPolygon) {
		ArrayList<Polygon> trapezoids = new ArrayList<Polygon>();
		ArrayList<Vertex> polygonVertices = fullPolygon.getVerticesAndHoles();
		ArrayList<Edge> polygonEdges = fullPolygon.getEdgesAndHoles();
		// ArrayList<Edge> trapezoidalEdges = new ArrayList<Edge>();

		// Get all vertices and sort in descending Y-order.
		polygonVertices.sort((v1, v2) -> Double.compare(v2.getY(), v1.getY()));

		// Make a huge sweepline (edge) to test for intersections.
		for (Vertex vertex : polygonVertices) {
			Edge sweepLine = new Edge(fullPolygon,
					new Vertex(fullPolygon, (int) (-Math.floor(fullPolygon.getMaxDistance())), vertex.getY()),
					new Vertex(fullPolygon, (int) (Math.ceil(fullPolygon.getMaxDistance())), vertex.getY()));

			// Check for edge intersections
			ArrayList<Vertex> intersections = new ArrayList<>();
			for (Edge edge : polygonEdges) {
				Vertex intersection = getIntersectionPoint(sweepLine, edge);
				if (intersection != null && !intersections.contains(intersection)) {
					intersections.add(intersection);
				}
			}

			intersections.sort((v1, v2) -> Double.compare(v1.getX(), v2.getX()));

			for (int i = 0; i < intersections.size() - 1; ++i) {
				Vertex i1 = intersections.get(i);
				Vertex i2 = intersections.get(i + 1);
				Edge tempEdge = new Edge(fullPolygon, i1, i2);
				if (insidePolygon(tempEdge.getMidpoint(), fullPolygon)) {
					boolean isStructure1 = polygonVertices.contains(i1);
					boolean isStructure2 = polygonVertices.contains(i2);

					// Check if the midpoint is inside the polygon
					if (isStructure1 && isStructure2) {
						Vertex s1 = fullPolygon.getMatchingVertex(i1);
						Vertex s2 = fullPolygon.getMatchingVertex(i2);
						if (!s1.isNeighbor(s2)) {
							trapezoidalEdges.add(tempEdge);
						}
					} else if ((isStructure1 || isStructure2)) {
						trapezoidalEdges.add(tempEdge);
					} else {
						trapezoidalEdges.add(tempEdge);
					}
				}

			}
		}

		trapezoids.add(fullPolygon);
		return trapezoids;
	}

	// Not implemented - for now simply returning the full polygon as a
	// "trapezoid".
	private ArrayList<Polygon> monotonePolygonization(ArrayList<Polygon> trapezoid) {
		ArrayList<Polygon> monotonePolygons = new ArrayList<Polygon>();
		monotonePolygons.addAll(trapezoid);
		return monotonePolygons;
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
			final ArrayList<Polygon> tris = triangulateMonotonePolygon2(poly);

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
		polygon.tiltHorizontals();

		ArrayList<Vertex> vertexArray = new ArrayList<>(polygon.getVertices());

		for (Hole hole : polygon.getHoles()) {
			vertexArray.addAll(hole.getVertices());
		}

		vertexArray.sort(Vertex.TOP_TO_BOTTOM_LEFT_TO_RIGHT);

		AVLTree<Edge> t = new AVLTree<>();

		helper = new HashMap<>();

		ArrayList<Polygon> result = new ArrayList<>();

		result.add(polygon);

		LinkedList<Vertex> q = new LinkedList<>(vertexArray);

		while (!q.isEmpty()) {
			Vertex vi = q.poll();
			VertexType vertexType = judgeVertexType(vi);

			switch (vertexType) {
				case REGULAR: {
					final double tiltX = 0.01;

					// Compute whether the interior is to the right of vi
					Vertex tempVertex = new Vertex(polygon, vi.getX() + tiltX, vi.getY());
					boolean interiorIsToTheRight = insidePolygon(tempVertex, polygon);

					if (interiorIsToTheRight) {
						Edge ei_m1 = vi.getInEdge();

						assert helper.containsKey(ei_m1);

						Vertex v = helper.get(ei_m1);

						VertexType vt = judgeVertexType(v);

						if (vt == VertexType.MERGE) {
							breakPolygon(result, v, vi);
						}

						t.remove(ei_m1);

						Edge ei = vi.getOutEdge();

						t.insert(ei);
						helper.put(ei, vi);
					} else {
						Edge ej = getFirstLeftEdge(polygon, vi);

						if (ej == null) {
							throw new RuntimeException("Oops");
						}

						Vertex v = helper.get(ej);

						VertexType vt = judgeVertexType(v);

						if (vt == VertexType.MERGE) {
							breakPolygon(result, v, vi);
						}

						helper.put(ej, vi);
					}
				}
				break;
				case START: {
					Edge ei = vi.getOutEdge();
					t.insert(ei);
					helper.put(ei, vi);
				}
				break;
				case END: {
					Edge ei_m1 = vi.getInEdge();

					assert helper.containsKey(ei_m1);

					Vertex v = helper.get(ei_m1);

					VertexType vt = judgeVertexType(v);

					if (vt == VertexType.MERGE) {
						breakPolygon(result, v, vi);
					}

					t.remove(ei_m1);
				}
				break;
				case MERGE: {
					Edge ei_m1 = vi.getInEdge();

					assert helper.containsKey(ei_m1);

					Vertex v = helper.get(ei_m1);

					VertexType vt = judgeVertexType(v);

					if (vt == VertexType.MERGE) {
						breakPolygon(result, v, vi);
					}

					t.remove(ei_m1);

					Edge ej = getFirstLeftEdge(polygon, vi);

					if (ej != null) {
						v = helper.get(ej);

						vt = judgeVertexType(v);

						if (vt == VertexType.MERGE) {
							breakPolygon(result, v, vi);
						}

						helper.put(ej, vi);
					}
				}
				break;
				case SPLIT: {
					Edge ej = getFirstLeftEdge(polygon, vi);

					if (ej != null) {
						Vertex v = helper.get(ej);

						breakPolygon(result, v, vi);

						helper.put(ej, vi);
					}

					Edge ei = vi.getOutEdge();

					t.insert(ei);
					helper.put(ei, vi);
				}
				break;
			}
		}

		return result;
	}

	private static void breakPolygon(ArrayList<Polygon> newPolygons, Vertex v1, Vertex v2) {
		Optional<Polygon> toBreakOpt = newPolygons.stream().filter(poly -> {
			ArrayList<Vertex> vertices = poly.getVerticesAndHoles();
			return vertices.contains(v1) && vertices.contains(v2);
		}).findFirst();

		if (!toBreakOpt.isPresent()) {
			throw new IllegalArgumentException("There isn't a polygon containing these two vertices");
		}

		final Polygon toBreak = toBreakOpt.get();

		final boolean v1InBoundary = toBreak.getVertices().contains(v1);
		final boolean v2InBoundary = toBreak.getVertices().contains(v2);

		// We can't break holes
		assert v1InBoundary || v2InBoundary;

		final boolean shouldCreateNewPolygon = v1InBoundary == v2InBoundary;

		Polygon p1 = null, p2 = null;
		Hole breakingHole = null;

		if (!shouldCreateNewPolygon) {
			Vertex v;

			if (v1InBoundary) {
				v = v2;
			} else if (v2InBoundary) {
				v = v1;
			} else {
				throw new RuntimeException("Not possible");
			}

			breakingHole = toBreak.getHoles()
					.stream()
					.filter(h -> h.getVertices().contains(v))
					.findFirst()
					.get();
		}

		if (shouldCreateNewPolygon) {
			List<Vertex> vertices = toBreak.getVertices();
			int vertexCount = vertices.size();
			ArrayList<Vertex> boundaryVertices1 = new ArrayList<>();
			ArrayList<Vertex> boundaryVertices2 = new ArrayList<>();

			int v1Index = vertices.indexOf(v1);
			int v2Index = vertices.indexOf(v2);

			assert v1Index != v2Index;

			int iter1End, iter2End;

			if (v1Index > v2Index) {
				iter1End = v2Index + vertexCount;
				iter2End = v1Index;
			} else {
				iter1End = v2Index;
				iter2End = v1Index + vertexCount;
			}

			for (int i = v1Index; i <= iter1End; i = i + 1) {
				int k = i % vertexCount;
				boundaryVertices1.add(vertices.get(k).clone());
			}

			for (int i = v2Index; i <= iter2End; i = i + 1) {
				int k = i % vertexCount;
				boundaryVertices2.add(vertices.get(k).clone());
			}

			p1 = new Polygon();
			p1.setVertices(boundaryVertices1);
			p2 = new Polygon();
			p2.setVertices(boundaryVertices2);
		} else {
			// Connect the boundary and the hole
			Vertex vertexOnTheBoundary, vertexOnTheHole;

			if (toBreak.getVertices().contains(v1)) {
				vertexOnTheBoundary = v1;
				vertexOnTheHole = v2;
			} else {
				vertexOnTheBoundary = v2;
				vertexOnTheHole = v1;
			}

			ArrayList<Vertex> vertices = new ArrayList<>();

			int boundaryIndex = toBreak.getVertices().indexOf(vertexOnTheBoundary);
			int holeIndex = breakingHole.getVertices().indexOf(vertexOnTheHole);

			assert boundaryIndex >= 0;
			assert holeIndex >= 0;

			for (int i = boundaryIndex; i <= boundaryIndex + toBreak.getVertices().size(); ++i) {
				vertices.add(toBreak.getVertices().get(i % toBreak.getVertices().size()).clone());
			}

			for (int i = holeIndex; i <= holeIndex + breakingHole.getVertices().size(); ++i) {
				vertices.add(breakingHole.getVertices().get(i % breakingHole.getVertices().size()).clone());
			}

			p1 = new Polygon();
			p1.setVertices(vertices);
		}

		// Generate edges for insidePolygon test
		// Edges will be regenerated after adding holes
		p1.generateEdges();

		if (p2 != null) {
			p2.generateEdges();
		}

		// Handle holes if any; put them inside newly created polygons.
		if (toBreak.getHoles().size() > 0) {
			ArrayList<Hole> p1Holes = new ArrayList<>();
			ArrayList<Hole> p2Holes = new ArrayList<>();

			for (Hole hole : toBreak.getHoles()) {
				if (hole == breakingHole) {
					continue;
				}

				Vertex v = hole.getVertices().get(0);

				if (insidePolygon(v, p1)) {
					p1Holes.add(hole.clone(p1));
				} else {
					p2Holes.add(hole.clone(p2));
				}
			}

			p1.setHoles(p1Holes);

			if (p2 != null) {
				p2.setHoles(p2Holes);
			}
		}

		p1.generateEdges();

		if (p2 != null) {
			p2.generateEdges();
		}

		newPolygons.remove(toBreak);

		newPolygons.add(p1);

		if (p2 != null) {
			newPolygons.add(p2);
		}
	}

	private Edge getFirstLeftEdge(Polygon p, Vertex v) {
		if (v.getInEdge().isHorizontal() || v.getOutEdge().isHorizontal()) {
			return null;
		}

		Rectangle2D boundingBox = p.getBoundingBox();

		Edge line = new Edge(p, v, new Vertex(p, boundingBox.getMinX(), v.getY()));
		Vertex lastIntersection = null;
		Edge result = null;

		for (Edge edge : p.getEdges()) {
			Vertex intersectionPoint = getIntersectionPoint(line, edge);

			if (intersectionPoint != null && intersectionPoint.getX() < v.getX()) {
				if (result == null || computeDistance(v, intersectionPoint) < computeDistance(v, lastIntersection)) {
					result = edge;
					lastIntersection = intersectionPoint;
				}
			}
		}

		return result;
	}

	private Map<Edge, Vertex> helper;

	private static VertexType judgeVertexType(Vertex vi) {
		ArrayList<Vertex> neighbors = vi.getNeighbors();

		assert neighbors.size() == 2;

		// Note: Here it is assumed that boundary edges are stored in CCW order,
		// and holes in CW order
		Vertex v1 = neighbors.get(0);
		Vertex v2 = neighbors.get(1);

		boolean flat = v1.getY() == vi.getY() && v2.getY() == vi.getY();
		boolean allAbove = !flat && v1.getY() >= vi.getY() && v2.getY() >= vi.getY();
		boolean allBelow = !flat && v1.getY() <= vi.getY() && v2.getY() <= vi.getY();

		if (!allAbove && !allBelow) {
			return VertexType.REGULAR;
		}

		double x1 = v1.getX() - vi.getX();
		double y1 = v1.getY() - vi.getY();
		double x2 = v2.getX() - vi.getX();
		double y2 = v2.getY() - vi.getY();

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
			if (a1 == 0) {
				a1 = Math.PI * 2;
			} else {
				a1 = (a1 + Math.PI * 2) % (Math.PI * 2);
			}

			if (a2 == 0) {
				a2 = Math.PI * 2;
			} else {
				a2 = (a2 + Math.PI * 2) % (Math.PI * 2);
			}

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

	/**
	 * @param p A monotone polygon
	 * @return Triangles
	 */
	private static ArrayList<Polygon> triangulateMonotonePolygon2(Polygon p) {
		ArrayList<Vertex> u = new ArrayList<>(p.getVertices());
		final int n = u.size();

		u.sort(Vertex.TOP_TO_BOTTOM_LEFT_TO_RIGHT);

		ArrayList<Vertex> leftChain = new ArrayList<>();
		ArrayList<Vertex> rightChain = new ArrayList<>();

		findMonotoneChains(u, leftChain, rightChain);

		Stack<Vertex> s = new Stack<>();

		s.push(u.get(0));
		s.push(u.get(1));

		ArrayList<Polygon> d = new ArrayList<>();

		d.add(p);

		for (int j = 2; j < n - 1; ++j) {
			Vertex uj = u.get(j);
			Vertex top = s.peek();

			if ((leftChain.contains(uj) && rightChain.contains(top)) ||
					(rightChain.contains(uj) && leftChain.contains(top))) {
				while (!s.isEmpty()) {
					Vertex v = s.pop();
					boolean isLastOne = s.isEmpty();

					if (!isLastOne) {
						breakPolygon(d, uj, v);
					}
				}

				Vertex uj_m1 = u.get(j - 1);

				s.push(uj_m1);
				s.push(uj);
			} else {
				Vertex lastPopped = s.pop();

				while (!s.isEmpty()) {
					Vertex v2 = s.peek();

					Edge diagonal = new Edge(p, uj, v2);
					List<Vertex> intersections = edgeIntersectPolygon(diagonal, p);

					if (intersections.size() < 2) {
						throw new RuntimeException("Not expected");
					}

					if (intersections.size() != 2) {
						break;
					}

					lastPopped = s.pop();

					breakPolygon(d, uj, lastPopped);
				}

				s.push(lastPopped);

				s.push(uj);
			}
		}

		assert s.size() >= 2;

		Vertex $_ = s.pop();
		Vertex un = u.get(n - 1);

		while (!s.isEmpty()) {
			Vertex v = s.pop();
			boolean isLastOne = s.isEmpty();

			if (!isLastOne) {
				breakPolygon(d, un, v);
			}
		}

		return d;
	}

	private static void findMonotoneChains(ArrayList<Vertex> sorted, ArrayList<Vertex> left, ArrayList<Vertex> right) {
		final Vertex top = sorted.get(0);
		final Vertex bottom = sorted.get(sorted.size() - 1);

		Vertex nextLeft = top.getOutEdge().getOtherVertex(top);

		while (nextLeft != bottom) {
			left.add(nextLeft);
			nextLeft = nextLeft.getOutEdge().getOtherVertex(nextLeft);
		}

		Vertex nextRight = top.getInEdge().getOtherVertex(top);

		while (nextRight != bottom) {
			right.add(nextRight);
			nextRight = nextRight.getInEdge().getOtherVertex(nextRight);
		}

		boolean leftOrRight = left.size() < right.size();

		if (leftOrRight) {
			left.add(0, top);
		} else {
			right.add(0, top);
		}

		leftOrRight = left.size() < right.size();

		if (leftOrRight) {
			left.add(bottom);
		} else {
			right.add(bottom);
		}
	}

	/*
	 * Implementation of "TriangulateMonotonePolygon(polygon P as DCEL)". From
	 * the course slides "3 - Art Gallery Triangulation", page 158. Input:
	 * Y-Monotone polygon "p". Output: ArrayList of triangular polygons.
	 */
	private ArrayList<Polygon> triangulateMonotonePolygon(Polygon monotonePolygon) {
		ArrayList<Vertex> triangulationVertices = monotonePolygon.getVertices();
		ArrayList<Edge> triangulationEdges = monotonePolygon.getEdges();

		// Sorts the vertices on Y-Descending order.
		triangulationVertices.sort((v1, v2) -> Double.compare(v2.getY(), v1.getY()));
		monotonePolygon.constructChains();

		// Implementation translation from the book.
		Stack<Vertex> s = new Stack<Vertex>();
		s.push(triangulationVertices.get(0));
		s.push(triangulationVertices.get(1));

		// Begin generating the internal edges
		for (int j = 2; j < triangulationVertices.size() - 1; ++j) {
			Vertex uj = triangulationVertices.get(j);
			// System.out.println(triangulationVertices.indexOf(s.peek())+ " is
			// neighbor of "+ triangulationVertices.indexOf(uj) + " : " +
			// s.peek().isNeighbor(uj));
			if (!s.peek().isNeighbor(uj)) {
				while (s.size() > 1) {
					Vertex v = s.pop();
					triangulationEdges.add(new Edge(monotonePolygon, uj, v));
				}
				s.pop();
				s.push(triangulationVertices.get(j - 1));
				s.push(uj);
			} else {
				Vertex v = s.pop();
				// Must implement "sees"
				while (s.size() > 0) {
					v = s.pop();
					triangulationEdges.add(new Edge(monotonePolygon, uj, v));
				}
				s.push(v);
				s.push(uj);
			}
		}

		s.pop();
		while (s.size() > 1) {
			Vertex v = s.pop();
			triangulationEdges.add(new Edge(monotonePolygon, triangulationVertices.get(triangulationVertices.size() - 1), v));
		}

		// Finally convert the loose list of edges into solid triangular
		// polygons and return.
		ArrayList<Polygon> triangles = constructTriangles(triangulationEdges);
		return triangles;
	}

	// Iterates through all of the loose edges in the list and attempts to put
	// together triangles by matching their vertices.
	// Unfortunately O(n^3) as all edges have to be checked in each of the
	// levels required (triangles = 3 levels).
	// Perhaps could be improved through filtering, a better algorithm, or
	// better usage of references later on.
	private ArrayList<Polygon> constructTriangles(ArrayList<Edge> edges) {
		ArrayList<Polygon> triangles = new ArrayList<Polygon>();

		for (Edge firstEdge : edges) {
			Vertex firstVertex = firstEdge.getFirstVertex();
			for (Edge secondEdge : edges) {
				if (!secondEdge.equals(firstEdge) && secondEdge.containsVertex(firstVertex)) {
					Vertex secondVertex = secondEdge.getOtherVertex(firstVertex);
					for (Edge thirdEdge : edges) {
						if (!thirdEdge.equals(secondEdge) && !thirdEdge.equals(firstEdge)
								&& thirdEdge.containsVertex(secondVertex)) {
							if (thirdEdge.containsVertex(firstEdge.getOtherVertex(firstVertex))) {
								ArrayList<Edge> triangleEdges = new ArrayList<Edge>();
								triangleEdges.add(firstEdge);
								triangleEdges.add(secondEdge);
								triangleEdges.add(thirdEdge);

								Polygon triangle = new Polygon(triangleEdges);

								if (!triangles.stream().anyMatch(t -> t.edgeMatch(triangle))) {
									System.out.println("Adding triangle: " + firstEdge.getId() + " - "
											+ secondEdge.getId() + " - " + thirdEdge.getId());
									triangles.add(triangle);
								}
							}
						}
					}
				}
			}
		}
		return triangles;
	}

	public ArrayList<Polygon> computeVisibilityPolygon(Vertex viewPoint, Polygon p) {
		// Array of triangles comprising the visibility polygon
		ArrayList<Polygon> visiblePolygons = new ArrayList<Polygon>();

		// Array containing all points in the scene but the view point, and
		// sorted on angle+proximity.
		ArrayList<Vertex> visibleVertices = p.getVerticesAndHoles();
		visibleVertices.remove(viewPoint);
		visibleVertices.sort((v1, v2) -> compareAngleAndProximity(viewPoint, v1, v2));

		// Extend sightlines (edges) from the viewpoint into each of the scene
		// vertices.
		ArrayList<Edge> sightLines = new ArrayList<Edge>();
		for (Vertex w : visibleVertices) {
			double angle = computeCCWAngle(w, viewPoint);
			int endX = (int) (viewPoint.getX() + (Math.cos(Math.toRadians(angle)) * p.getMaxDistance()));
			int endY = (int) (viewPoint.getY() + (Math.sin(Math.toRadians(angle)) * p.getMaxDistance()));
			Vertex endOfLine = new Vertex(p, endX, endY);
			sightLines.add(new Edge(p, viewPoint, endOfLine));
			System.out.println(w.getId() + " - " + angle);
		}

		ArrayList<Edge> obstacles = p.getEdgesAndHoles();
		ArrayList<Vertex> visibilityPolygonVertices = new ArrayList<Vertex>();
		for (Edge sightLine : sightLines) {
			PriorityQueue<Vertex> queue = new PriorityQueue<Vertex>(
					(v1, v2) -> Double.compare(computeDistance(viewPoint, v1), computeDistance(viewPoint, v2)));
			for (Edge obstacle : obstacles) {
				Vertex intersection = getIntersectionPoint(sightLine, obstacle);
				if (intersection != null) {
					visibilityPolygonVertices.add(intersection);
				}
			}
		}

		for (int i = 0; i < visibilityPolygonVertices.size() - 1; ++i) {
			ArrayList<Vertex> visTriangleVertices = new ArrayList<Vertex>();
			visTriangleVertices.add(viewPoint);
			visTriangleVertices.add(visibilityPolygonVertices.get(i));
			visTriangleVertices.add(visibilityPolygonVertices.get((i + 1)));
			Polygon visibleTriangle = new Polygon(visTriangleVertices);
			visiblePolygons.add(visibleTriangle);
		}

		return visiblePolygons;
	}

	/*
	 * Implementation of "VisibleVertices(P, S)". From the book
	 * "Computational Geometry: Algorithms and Applications - Third Edition" by
	 * M. Berg , page 328. Input: Point P, set of polygonal obstacles S. Output:
	 * The set of visible vertices from P.
	 */
	private ArrayList<Vertex> visibleVertices(Vertex v, Polygon p) {
		// 1) Sort the obstacle vertices according to the clockwise angle that
		// the halfline from p to each vertex makes with the positive x-axis.
		ArrayList<Vertex> obstacleVertices = p.getVerticesAndHoles();
		obstacleVertices.sort((v1, v2) -> compareAngleAndProximity(v, v1, v2));

		// Find the obstacle edges that are properly intersected by p parallel
		// to the x-axis and store them in a
		// balanced search tree in the order they intersect p.
		AVLTree<Edge> tree = new AVLTree<Edge>();
		Edge sweepLine = new Edge(p, v, new Vertex(p, (int) (Math.ceil(p.getMaxDistance())), v.getY()));

		ArrayList<Edge> intersectedEdges = new ArrayList<Edge>();
		for (Edge obstacleEdge : p.getEdges()) {
			if (getIntersectionPoint(sweepLine, obstacleEdge) != null) {
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

				List<Edge> incidentEdges = (List<Edge>) intersectedEdges.stream().filter(e -> e.containsVertex(vi))
						.collect(Collectors.toList());
				for (Edge e : incidentEdges) {
					// Get the other vertex (not "vi") of the incident edge
					Vertex otherVertex = e.getOtherVertex(vi);
					// If the orientation if clockwise w.r.t. the sweepline, add
					// the edge to the tree.
					if (orientation(v, vi, otherVertex) == 1) {
						tree.insert(e);
					}
					// If the orientation if counter-clockwise w.r.t. the
					// sweepline, remove the edge to the tree.
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
	 * M. Berg , page 329. Input: Point W. Output: The set of visible vertices
	 * from P.
	 */
	private boolean isVisible(ArrayList<Vertex> w, int i) {
		/*
		 * if() {
		 *
		 * return false; }else if(i == 0 || w.get(i - 1)) {
		 *
		 * }
		 */
		return true;
	}

	private double computeDistance(Vertex v1, Vertex v2) {
		return Math.sqrt(Math.pow(v1.getY() - v2.getY(), 2) + Math.pow((v1.getX() - v2.getX()), 2));
	}

	private double computeCCWAngle(Vertex v, Vertex reference) {
		double deltaX = v.getX() - reference.getX();
		double deltaY = v.getY() - reference.getY();
		double angle = Math.toDegrees(Math.atan2(deltaY, deltaX));
		if (angle < 0)
			angle += 360;
		return angle;
	}

	// Comparison based on the angle
	// If the angles match (colinear) then prioritize the vertex closed to the
	// reference.
	public int compareAngleAndProximity(Vertex
			                                    reference, Vertex v1, Vertex v2) {
		int result = Double.compare(computeCCWAngle(v1, reference), computeCCWAngle(v2, reference));
		if (result == 0) {
			Double d1 = computeDistance(v1, reference);
			Double d2 = computeDistance(v2, reference);
			return d1 > d2 ? 1 : -1;
		}
		return result;
	}

	private static Vertex getIntersectionPoint(Edge e1, Edge e2) {
		assert e1.getPolygon() == e2.getPolygon();

		double x1 = e1.getFirstVertex().getX();
		double y1 = e1.getFirstVertex().getY();
		double x2 = e1.getSecondVertex().getX();
		double y2 = e1.getSecondVertex().getY();

		double x3 = e2.getFirstVertex().getX();
		double y3 = e2.getFirstVertex().getY();
		double x4 = e2.getSecondVertex().getX();
		double y4 = e2.getSecondVertex().getY();

		double ax = x2 - x1;
		double ay = y2 - y1;
		double bx = x4 - x3;
		double by = y4 - y3;

		double denominator = ax * by - ay * bx;

		if (denominator == 0)
			return null;

		double cx = x3 - x1;
		double cy = y3 - y1;

		double t = (cx * by - cy * bx) / (denominator);
		if (t < 0 || t > 1)
			return null;

		double u = (cx * ay - cy * ax) / (denominator);
		if (u < 0 || u > 1)
			return null;

		return new Vertex(e1.getPolygon(), (x1 + t * ax), (y1 + t * ay));
	}

	private double getIntersectionDistance(Vertex v, Edge e1, Edge e2) {
		double x1, x2, x3, x4, y1, y2, y3, y4;

		x1 = e1.getFirstVertex().getX();
		y1 = e1.getFirstVertex().getY();
		x2 = e1.getSecondVertex().getX();
		y2 = e1.getSecondVertex().getY();
		x3 = e2.getFirstVertex().getX();
		y3 = e2.getFirstVertex().getY();
		x4 = e2.getSecondVertex().getX();
		y4 = e2.getSecondVertex().getY();

		double x = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4))
				/ ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
		double y = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4))
				/ ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));

		return Math.sqrt(Math.pow(y - v.getY(), 2) + Math.pow((x - v.getX()), 2));
	}

	private int compareIntersectionDistance(Vertex
			                                        v, Edge e, Edge e1, Edge e2) {
		return Double.compare(getIntersectionDistance(v, e, e1), getIntersectionDistance(v, e, e2));
	}

	private Polygon computeBoundingBox(Polygon p) {
		double minX, minY, maxX, maxY;
		minX = minY = Integer.MAX_VALUE;
		maxX = maxY = Integer.MIN_VALUE;

		for (Vertex v : p.getVertices()) {
			if (v.getX() > maxX)
				maxX = v.getX();
			if (v.getX() < minX)
				minX = v.getX();
			if (v.getY() > maxY)
				maxY = v.getY();
			if (v.getY() < minY)
				minY = v.getY();
		}
		ArrayList<Vertex> vertices = new ArrayList<Vertex>();
		vertices.add(new Vertex(p, maxX, maxY));
		vertices.add(new Vertex(p, maxX, minY));
		vertices.add(new Vertex(p, minX, minY));
		vertices.add(new Vertex(p, minX, maxY));
		Polygon boundingBox = new Polygon(vertices);
		return boundingBox;
	}

	private static List<Vertex> edgeIntersectPolygon(Edge e, Polygon p) {
		ArrayList<Vertex> intersections = new ArrayList<>();
		ArrayList<Edge> edges = p.getEdgesAndHoles();

		for (Edge edge : edges) {
			Vertex intersectionPoint = getIntersectionPoint(e, edge);

			if (intersectionPoint != null && !intersections.contains(intersectionPoint)) {
				intersections.add(intersectionPoint);
			}
		}

		return intersections;
	}

	private static boolean insidePolygon(Vertex v, Polygon p) {
		ArrayList<Vertex> intersections = new ArrayList<>();
		Rectangle2D boundingBox = p.getBoundingBox();
//		Vertex v1 = new Vertex(p, boundingBox.getMinX(), v.getY());
		Vertex v2 = new Vertex(p, boundingBox.getMaxX() + 10, v.getY());
		Edge line = new Edge(p, v, v2);

//		// HACK: if the vertex to test is on an edge, add it to intersections first.
//		// HACK: in our case, all vertices "on" edges that we test should be endpoints
//		//       of edges, so simply test if it has in/out edges
//		if (v.getInEdge() != null && v.getOutEdge() != null) {
//			intersections.add(v);
//		}

		ArrayList<Edge> edges = p.getEdgesAndHoles();
		for (Edge edge : edges) {
			Vertex intersectionPoint = getIntersectionPoint(line, edge);

			if (intersectionPoint != null && !intersections.contains(intersectionPoint)) {
				intersections.add(intersectionPoint);
			}
		}

//		System.out.println(intersections.size() + " at " + v.getY());
		return intersections.size() % 2 == 1;
	}

	private boolean areColinear(Vertex v1, Vertex v2, Vertex v3) {
		double area = v1.getX() * (v2.getY() - v3.getY()) + v2.getX() * (v3.getY() - v1.getY())
				+ v3.getX() * (v1.getY() - v2.getY());
		return (area == 0);
	}

	private int orientation(Vertex p, Vertex q, Vertex r) {
		double val = (q.getY() - p.getY()) * (r.getX() - q.getX()) - (q.getX() - p.getX()) * (r.getY() - q.getY());
		if (val < 0) {
			return 0; // colinear
		} else {
			return (val > 0) ? 1 : 2; // clock or counterclock wise
		}

	}
}
