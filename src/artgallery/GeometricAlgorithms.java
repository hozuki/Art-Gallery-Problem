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

		for (int i = 0; i < monotones.size(); i++) {
			Polygon poly = monotones.get(i);
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
		polygon.untiltHorizontals();
		polygon.tiltHorizontals();

		ArrayList<Vertex> vertexArray = new ArrayList<>(polygon.getVertices());

		for (Hole hole : polygon.getHoles()) {
			vertexArray.addAll(hole.getVertices());
		}

		vertexArray.sort(Vertex.TOP_TO_BOTTOM_LEFT_TO_RIGHT);

		AVLTree<Edge> t = new AVLTree<>();

		final Map<Edge, Vertex> helper = new HashMap<>();

		ArrayList<Polygon> result = new ArrayList<>();

		result.add(polygon);

		LinkedList<Vertex> q = new LinkedList<>(vertexArray);

		while (!q.isEmpty()) {
			Vertex vi = q.poll();
			VertexType vertexType = judgeVertexType(vi);

			System.out.println(String.format("MakeMonotone: visiting vertex %s (%s)", vi.toString(), vertexType.toString()));

			try {
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
							Edge ej = getFirstLeftEdge(polygon, vi, t);

							if (ej == null) {
								throw new RuntimeException("Oops");
							}

							// Should return (-25,-124.9)-(-25,-200) @ vertex :-5(100,-124.9)
							// instead, returned (-100,-125)-(-25,-124.9)
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

						Edge ej = getFirstLeftEdge(polygon, vi, t);

						if (ej != null) {
							v = helper.get(ej);

							vt = judgeVertexType(v);

							if (vt == VertexType.MERGE) {
								breakPolygon(result, v, vi);
							} else {
//							if (vi.getOutEdge().getOtherVertex(vi).equals(v)) {
//								// Degenerated case (remove this condition and check AGS5)
//								Vertex otherV = ej.getOtherVertex(v);
//
//								if (otherV != null) {
//									// AGS5
//									vt = judgeVertexType(otherV);
//
//									if (vt == VertexType.MERGE) {
//										breakPolygon(result, otherV, vi);
//									}
//								} else {
//									// AGS4
//									for (Vertex vv : new Vertex[]{ej.getUpperVertex(), ej.getLowerVertex()}) {
//										vt = judgeVertexType(vv);
//
//										if (vt == VertexType.MERGE) {
//											breakPolygon(result, vv, vi);
//										}
//									}
//								}
//							}
							}

							helper.put(ej, vi);
						}
					}
					break;
					case SPLIT: {
						Edge ej = getFirstLeftEdge(polygon, vi, t);

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
			} catch (Exception e) {
				e.printStackTrace();
				break;
			}
		}

		return result;
	}

	private static boolean inSamePolygon(ArrayList<Polygon> polygons, Vertex v1, Vertex v2) {
		Optional<Polygon> resultOpt = polygons.stream().filter(poly -> {
			ArrayList<Vertex> vertices = poly.getVerticesAndHoles();
			return vertices.contains(v1) && vertices.contains(v2);
		}).findFirst();

		return resultOpt.isPresent();
	}

	private static void breakPolygon(ArrayList<Polygon> newPolygons, Vertex v1, Vertex v2) {
		breakPolygon(newPolygons, v1, v2, false);
	}

	private static void breakPolygon(ArrayList<Polygon> newPolygons, Vertex v1, Vertex v2, boolean debug) {
		System.out.println(String.format("Connecting %s to %s", v1.toString(), v2.toString()));

		Optional<Polygon> toBreakOpt = newPolygons.stream().filter(poly -> {
			ArrayList<Vertex> vertices = poly.getVerticesAndHoles();

			if (debug) {
				System.out.println("-------");

				for (Vertex v : vertices) {
					System.out.println(v.toString());
				}
			}

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

			Optional<Hole> breakingHoleOpt = toBreak.getHoles()
				.stream()
				.filter(h -> h.getVertices().contains(v))
				.findFirst();

			if (!breakingHoleOpt.isPresent()) {
				throw new RuntimeException("Query vertex is not on the hole");
			}

			breakingHole = breakingHoleOpt.get();
		}

		if (shouldCreateNewPolygon) {
			List<Vertex> vertices = toBreak.getVertices();
			int vertexCount = vertices.size();
			ArrayList<Vertex> boundaryVertices1 = new ArrayList<>();
			ArrayList<Vertex> boundaryVertices2 = new ArrayList<>();

			final Tuple<Integer, Integer> indices = findProperV1V2Indices(toBreak, vertices, v1, v2);

			final int v1Index = indices.item1;
			final int v2Index = indices.item2;

			assert v1Index != v2Index;

			int iter1End, iter2End;

			if (v1Index > v2Index) {
				iter1End = v2Index + vertexCount;
				iter2End = v1Index;
			} else {
				iter1End = v2Index;
				iter2End = v1Index + vertexCount;
			}

			p1 = new Polygon();
			p2 = new Polygon();

			for (int i = v1Index; i <= iter1End; i = i + 1) {
				int k = i % vertexCount;
				boundaryVertices1.add(vertices.get(k).clone(p1));
			}

			for (int i = v2Index; i <= iter2End; i = i + 1) {
				int k = i % vertexCount;
				boundaryVertices2.add(vertices.get(k).clone(p2));
			}

			p1.setVertices(boundaryVertices1);
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

			p1 = new Polygon();

			for (int i = boundaryIndex; i <= boundaryIndex + toBreak.getVertices().size(); ++i) {
				vertices.add(toBreak.getVertices().get(i % toBreak.getVertices().size()).clone(p1));
			}

			for (int i = holeIndex; i <= holeIndex + breakingHole.getVertices().size(); ++i) {
				vertices.add(breakingHole.getVertices().get(i % breakingHole.getVertices().size()).clone(p1));
			}

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

	// Problematic on AGS14, OK on other cases
	private static Tuple<Integer, Integer> findProperV1V2IndicesNaive(List<Vertex> vertices, Vertex v1, Vertex v2) {
		int v1Index = vertices.indexOf(v1);
		int v2Index = vertices.indexOf(v2);

		assert v1Index >= 0;
		assert v2Index >= 0;

		assert v1Index != v2Index;

		return new Tuple<>(v1Index, v2Index);
	}

	private static Tuple<Integer, Integer> findProperV1V2Indices(Polygon p, List<Vertex> vertices, Vertex v1, Vertex v2) {
		long v1Count = vertices.parallelStream().filter(v -> v.equals(v1)).count();
		long v2Count = vertices.parallelStream().filter(v -> v.equals(v2)).count();
		int v1Index, v2Index;

		assert (v1Count == 1 || v1Count == 2) && (v2Count == 1 || v2Count == 2);

		if (v1Count == 1 && v2Count == 1) {
			// Breaking a raw polygon (simple case)
			v1Index = vertices.indexOf(v1);
			assert v1Index >= 0;
			v2Index = vertices.indexOf(v2);
			assert v2Index >= 0;
		} else if ((v1Count == 1 && v2Count == 2) || (v2Count == 1 && v1Count == 2)) {
			// At least one of the vertices is broken by connecting boundary and a hole
			// Case: AGS14 (385,231)-(297,240) then error on (385,231)-(309,212)
			final Vertex vOne, vTwoPivot;

			if (v1Count == 1) {
				vOne = v1;
				vTwoPivot = v2;
			} else {
				vOne = v2;
				vTwoPivot = v1;
			}

			// Connect to the vTwo with lower edges (because we are sweeping from top to bottom)
			final Vertex[] vTwos = vertices.stream().filter(v -> v.equals(vTwoPivot)).toArray(Vertex[]::new);

			assert vTwos.length == 2;

			final double vTwo1Y = vTwos[0].getNeighbors().get(0).getY() + vTwos[0].getNeighbors().get(1).getY();
			final double vTwo2Y = vTwos[1].getNeighbors().get(0).getY() + vTwos[1].getNeighbors().get(1).getY();
			final double vTwo1X = vTwos[0].getNeighbors().get(0).getX() + vTwos[0].getNeighbors().get(1).getX();
			final double vTwo2X = vTwos[1].getNeighbors().get(0).getX() + vTwos[1].getNeighbors().get(1).getX();

			Vertex vToTilt;

			// Tilt the top-left one
			if (Util.notEquals(vTwo1Y, vTwo2Y)) {
				vToTilt = vTwo1Y > vTwo2Y ? vTwos[0] : vTwos[1];
			} else {
				assert Util.notEquals(vTwo1X, vTwo2X);
				vToTilt = vTwo1X < vTwo2X ? vTwos[0] : vTwos[1];
			}

			// TODO: What about (e1 is horizontal) || (e1 is vertical) || (e2 is horizontal) || (e2 is vertical)?
			//       Will require X-and-Y tilting on a non-intersect direction!
			final double tiltY = 0.01;
			Edge tempEdge = new Edge(p, vOne, new Vertex(p, vToTilt.getX(), vToTilt.getY() + tiltY));

			List<Vertex> intersections = edgeIntersectPolygon(tempEdge, p);

			final Vertex vTwo;

			if (intersections.isEmpty()) {
				vTwo = vToTilt;
			} else {
				// If there is at least one intersection, if we choose it as vTwo,
				// it will cause self intersection after breaking this polygon.
				// Here we need reference comparison.
				vTwo = vToTilt == vTwos[0] ? vTwos[1] : vTwos[0];
			}

			if (Objects.equals(vOne, v1)) {
				v1Index = vertices.indexOf(vOne);
				v2Index = Util.indexOfReference(vertices, vTwo);
			} else {
				v1Index = Util.indexOfReference(vertices, vTwo);
				v2Index = vertices.indexOf(vOne);
			}
		} else {
			throw new RuntimeException("Not expected");
		}

		return new Tuple<>(v1Index, v2Index);
	}

	private static final class Tuple<E1, E2> {

		public final E1 item1;

		public final E2 item2;

		public Tuple(final E1 element1, final E2 e2) {
			this.item1 = element1;
			this.item2 = e2;
		}

	}

	private static Edge getFirstLeftEdge(Polygon p, Vertex v, AVLTree<Edge> tree) {
		if (v.getInEdge().isHorizontal() || v.getOutEdge().isHorizontal()) {
			return null;
		}

		Rectangle2D boundingBox = p.getBoundingBox();

		Edge line = new Edge(p, v, new Vertex(p, boundingBox.getMinX(), v.getY()));
		Vertex lastIntersection = null;
		Edge result = null;

		ArrayList<Edge> edges = p.getEdgesAndHoles();

		for (Edge edge : edges) {
			if (edge.containsVertex(v)) {
				continue;
			}

			if (!tree.contains(edge)) {
				continue;
			}

			Vertex intersectionPoint = getIntersectionPoint(line, edge);

			if (intersectionPoint != null && intersectionPoint.getX() < v.getX()) {
				boolean shouldUpdate = false;

				if (result == null) {
					shouldUpdate = true;
				} else {
					double newDist = computeDistance(v, intersectionPoint);
					double oldDist = computeDistance(v, lastIntersection);

					if (newDist < oldDist) {
						shouldUpdate = true;
					} else if (Util.equals(newDist, oldDist)) {
						// Rare case (see notes).
						// That's why we can't filter out intersected points.
						// Return the one with the other vertex more to the right.
						Vertex newOther = edge.getOtherVertex(intersectionPoint);
						Vertex oldOther = result.getOtherVertex(lastIntersection);

						assert Util.notEquals(newOther.getY(), oldOther.getY());

						if ((newOther.getY() >= intersectionPoint.getY() && oldOther.getY() >= intersectionPoint.getY()) ||
							(newOther.getY() <= intersectionPoint.getY() && oldOther.getY() <= intersectionPoint.getY())) {
							assert Util.notEquals(newOther.getX(), oldOther.getX());

							// Degen case 1: v-shape
							// Select the one on the right
							shouldUpdate = newOther.getX() > oldOther.getX();
						} else {
							//Degen case 2: I-shape
							// Select the bottom one
							shouldUpdate = newOther.getY() < oldOther.getY();
						}
					}
				}

				if (shouldUpdate) {
					result = edge;
					lastIntersection = intersectionPoint;
				}
			}
		}

		return result;
	}

	private static VertexType judgeVertexType(Vertex vi) {
		ArrayList<Vertex> neighbors = vi.getNeighbors();

		assert neighbors.size() == 2;

		// Note: Here it is assumed that boundary edges are stored in CCW order,
		// and holes in CW order
		Vertex v1 = neighbors.get(0);
		Vertex v2 = neighbors.get(1);

		boolean flat = Util.equals(v1.getY(), vi.getY()) && Util.equals(v2.getY(), vi.getY());
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
			if (Util.equals(a1, 0)) {
				a1 = Math.PI * 2;
			} else {
				a1 = (a1 + Math.PI * 2) % (Math.PI * 2);
			}

			if (Util.equals(a2, 0)) {
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
					if (v.getX() == 682 && v.getY() == 234) {
						int i = 32;
					}
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

					if (v2.getX() == 682 && v2.getY() == 234) {
						int i = 32;
					}

					Edge diagonal = new Edge(p, uj, v2);
					List<Vertex> intersections = edgeIntersectPolygon(diagonal, p);

					if (intersections.size() < 2) {
						throw new RuntimeException("Not expected");
					}

					if (intersections.size() != 2) {
						// Intersects with other edge
						break;
					}

					{
						// Test whether the edge is out of the polygon (AGS14)
						Vertex tilted;
						final double tiltLength = 0.1;

						if (diagonal.isVertical()) {
							if (diagonal.getFirstVertex().getY() < diagonal.getSecondVertex().getY()) {
								tilted = new Vertex(p, diagonal.getFirstVertex().getX(), diagonal.getFirstVertex().getY() + tiltLength);
							} else {
								tilted = new Vertex(p, diagonal.getFirstVertex().getX(), diagonal.getFirstVertex().getY() - tiltLength);
							}
						} else {
							final double slope = diagonal.getSlope();
							final double cos = 1 / Math.sqrt(1 + slope * slope);
							final double sin = slope * cos;

							final double tiltX = tiltLength * cos;
							final double tiltY = tiltLength * sin;

							tilted = new Vertex(p, diagonal.getFirstVertex().getX() + tiltX, diagonal.getFirstVertex().getY() + tiltY);
						}

						if (!insidePolygon(tilted, p)) {
							break;
						}
					}

					lastPopped = s.pop();

					breakPolygon(d, uj, v2);
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

	// Simple visibility algorithm based on testing each vertex around the viewpoint.
	// Not finished and somewhat buggy on edge-cases.
	public ArrayList<Polygon> computeVisibilityPolygon(Vertex viewPoint, Polygon p) {
		// Array of triangles comprising the visibility polygon
		ArrayList<Polygon> visiblePolygons = new ArrayList<Polygon>();

		// Array containing all points in the scene but the view point, and
		// sorted on angle+proximity.
		ArrayList<Vertex> visibleVertices = p.getVerticesAndHoles();
		visibleVertices.remove(viewPoint);
		visibleVertices.sort((v1, v2) -> compareAngleAndProximity(viewPoint, v1, v2));

		visibleVertices.add(visibleVertices.get(0)); // Make it a full circle

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
					queue.offer(intersection);
				}
			}
			if (!queue.isEmpty()) {
				visibilityPolygonVertices.add(queue.poll());
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
	// Not currently working and the subroutine wasn't implemented.
	private ArrayList<Vertex> visibleVertices(Vertex v, Polygon p) {
		// 1) Sort the obstacle vertices according to the clockwise angle that
		// the half-line from p to each vertex makes with the positive x-axis.
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

	private static double computeDistance(Vertex v1, Vertex v2) {
		return Math.sqrt(Math.pow(v1.getY() - v2.getY(), 2) + Math.pow((v1.getX() - v2.getX()), 2));
	}

	private static double computeCCWAngle(Vertex v, Vertex reference) {
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
	public static int compareAngleAndProximity(Vertex reference, Vertex v1, Vertex v2) {
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

		if (Util.equals(denominator, 0))
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

	private static double getIntersectionDistance(Vertex v, Edge e1, Edge e2) {
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

	private static int compareIntersectionDistance(Vertex v, Edge e, Edge e1, Edge e2) {
		return Double.compare(getIntersectionDistance(v, e, e1), getIntersectionDistance(v, e, e2));
	}

	private static Polygon computeBoundingBox(Polygon p) {
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

		ArrayList<Vertex> vertices = p.getVerticesAndHoles();
		ArrayList<Edge> edges = p.getEdgesAndHoles();

		Set<Vertex> silentIntersections = new HashSet<>();

		edgeLoop:
		for (Edge edge : edges) {
			Vertex intersectionPoint = getIntersectionPoint(line, edge);

			if (silentIntersections.contains(intersectionPoint)) {
				continue;
			}

			if (intersections.contains(intersectionPoint)) {
				continue;
			}

			if (intersectionPoint != null) {
				for (final Vertex vi : vertices) {
					if (!vi.equals(intersectionPoint)) {
						continue;
					}

					List<Vertex> neighbors = vi.getNeighbors();

					assert neighbors.size() == 2;

					Vertex vn1 = neighbors.get(0);
					Vertex vn2 = neighbors.get(1);

					// Both above or below
					if ((vn1.getY() > vi.getY() && vn2.getY() > vi.getY()) ||
						(vn1.getY() < vi.getY() && vn2.getY() < vi.getY())) {
						// Ignore this intersection
						silentIntersections.add(v);
						continue edgeLoop;
					}

					// One is horizontal, only count intersection once
					if (Util.equals(vn1.getY(), vi.getY()) && Util.notEquals(vn2.getY(), vi.getY())) {
						silentIntersections.add(vn1);
					} else if (Util.equals(vn2.getY(), vi.getY()) && Util.notEquals(vn1.getY(), vi.getY())) {
						silentIntersections.add(vn2);
					}
				}

				intersections.add(intersectionPoint);
			}
		}

//		System.out.println(intersections.size() + " at " + v.getY());
		return intersections.size() % 2 == 1;
	}

	private static boolean areColinear(Vertex v1, Vertex v2, Vertex v3) {
		double area = v1.getX() * (v2.getY() - v3.getY()) + v2.getX() * (v3.getY() - v1.getY())
			+ v3.getX() * (v1.getY() - v2.getY());
		return (Util.equals(area, 0));
	}

	private static int orientation(Vertex p, Vertex q, Vertex r) {
		double val = (q.getY() - p.getY()) * (r.getX() - q.getX()) - (q.getX() - p.getX()) * (r.getY() - q.getY());
		if (val < 0) {
			return 0; // colinear
		} else {
			return (val > 0) ? 1 : 2; // clock or counterclock wise
		}

	}
}
