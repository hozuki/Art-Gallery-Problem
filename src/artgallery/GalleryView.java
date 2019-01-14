package artgallery;

import java.awt.*;
import java.awt.Polygon;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JPanel;

import artgallery.Actors.*;
import artgallery.geometricalElements.*;

@SuppressWarnings("serial")
public class GalleryView extends JPanel implements MouseListener, MouseMotionListener, MouseWheelListener {

	private GalleryModel gallery;
	private BufferedImage guardIcon;
	private BufferedImage thiefIcon;
	private double scale = 1;
	private Point startDrag, endDrag;
	private int offsetX = 0;
	private int offsetY = 0;

	GalleryView() {
		gallery = null;
		loadResources();
		addMouseListener(this);
		addMouseWheelListener(this);
	}

	GalleryView(GalleryModel gallery) {
		this.gallery = gallery;
		loadResources();
		addMouseWheelListener(this);
	}

	private void loadResources() {
		// Fetch the resources for the drawing of the library actors.
		try {
			guardIcon = ImageIO.read(getClass().getResource("/resources/guard.png"));
			thiefIcon = ImageIO.read(getClass().getResource("/resources/thief.png"));
		} catch (IOException e1) {
			e1.printStackTrace();
		}
	}

	public void updateGallery(GalleryModel updatedGallery) {
		this.gallery = updatedGallery;
		this.repaint();
	}

	@Override
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		Graphics2D g2d = (Graphics2D) g;

		// Hints for smooth rendering of shapes using anti-aliasing.
		RenderingHints rh = new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		rh.put(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		g2d.setRenderingHints(rh);

		// Get dimensions of the canvas for scaling of the components rendering.
		Dimension size = getSize();
		double w = size.getWidth();
		double h = size.getHeight();

		// Transform the Java coordinate system to the Cartesian system.
		g2d.translate(w / 2 + offsetX, h / 2 + offsetY);
		g2d.scale(scale, -scale);

//		g2d.setFont(new Font("sans-serif", Font.PLAIN, -10));

		if (this.gallery != null) {
			// Draw base surface polygon for the gallery.
			drawArea(g2d, gallery.getGallery(), null);

			// If selected, display the triangles from the gallery polygon
			// triangulation.
			if (gallery.showTriangulation()) {
				gallery.getGallery().computeTriangulation();
				drawTriangulation(g2d, gallery.getGallery().getTriangulation());
			}

			// If selected, display the visibility polygon of the guards.
			// Currently not fully implemented - only showing the polygon for
			// vertex 1.
			// Perhaps tweak in the future to work on a "clicked guard" to make
			// interactive & avoid overlapping.


			// Draw the gallery bounding edges.
			for (Edge edge : gallery.getGallery().getEdges()) {
				drawEdge(g2d, edge, 4, null);
			}
			// Draw the gallery vertices.
			for (Vertex vertex : gallery.getGallery().getVertices()) {
				drawVertex(g2d, vertex, 12, null);
			}
			// Draw the gallery holes.
			for (Hole hole : gallery.getGallery().getHoles()) {
				drawHole(g2d, hole);
			}

			// Could be modified with the uninplemented abstract class for
			// actors
			if (gallery.showActors()) {
				for (Guard guard : gallery.getGuards()) {
					if (gallery.showGuardsPaths()) {
						Color copBlue = new Color(100, 200, 250);
						for (Edge edge : guard.getRouteEdges(gallery.getGallery())) {
							drawEdge(g2d, edge, 2, copBlue);
						}
						for (RoutePoint point : guard.getRoutePoints()) {
							drawVertex(g2d, point.toVertex(gallery.getGallery()), 6, copBlue);
						}
					}
					drawGuard(g2d, guard);
				}

				for (Thief thief : gallery.getThieves()) {
					if (gallery.showThievesPaths()) {
						Color thiefOrange = new Color(250, 150, 0);
						for (Edge edge : thief.getRouteEdges(gallery.getGallery())) {
							drawEdge(g2d, edge, 2, thiefOrange);
						}
						for (RoutePoint point : thief.getRoutePoints()) {
							drawVertex(g2d, point.toVertex(gallery.getGallery()), 6, thiefOrange);
						}
					}
					drawThief(g2d, thief);
				}
			}
			if (gallery.showVisibility()) {
				//Vertex viewer = gallery.getGallery().getVertices().get(5);
				int numGuards = gallery.getGuards().size();

				for (int i = 0; i < numGuards; ++i) {
					ArrayList<artgallery.geometricalElements.Polygon> guardView = gallery.computeGuardVisibility(i);
					Color copBlue = new Color(255, 200, 10, 120);

					for (artgallery.geometricalElements.Polygon visibleTriangle : guardView) {
						drawArea(g2d, visibleTriangle, copBlue);
						visibleTriangle.getVertices().forEach(v -> drawVertex(g2d, v, 8, Color.ORANGE));
					}
				}
			}

//			for (Edge e : this.gallery.getTrapezoidalEdges()) {
//				drawVertex(g2d, e.getFirstVertex(), 4, Color.MAGENTA);
//				drawVertex(g2d, e.getSecondVertex(), 4, Color.MAGENTA);
//				drawEdge(g2d, e, 1, Color.MAGENTA);
//				drawEdgeLabel(g2d, e.getMidpoint(), Color.MAGENTA);
//			}

			if (gallery.showLegend()) {
				drawLegend(g2d, w, h);
				drawVertexPositions(g2d, w, h);
			}
		}
	}

	private void drawVertex(Graphics2D g2d, Vertex v, int size, Color c) {
		size = (int) (Math.round(size / scale));
		Ellipse2D e = new Ellipse2D.Double(v.getX() - size / 2, v.getY() - size / 2, size, size);
		if (c == null) {
			if (v.isArt()) {
				g2d.setColor(Color.YELLOW);
				g2d.fill(e);
				g2d.setColor(g2d.getColor().brighter());
			} else if (v.isExit()) {
				g2d.setColor(Color.RED);
				g2d.fill(e);
				g2d.setColor(g2d.getColor().brighter());
			} else {
				g2d.setColor(Color.LIGHT_GRAY);
				g2d.fill(e);
				g2d.setColor(g2d.getColor().brighter());
			}
		} else {
			g2d.setColor(c);
			g2d.fill(e);
			g2d.setColor(c.brighter());
		}
		g2d.setStroke(new BasicStroke(2));
		g2d.draw(e);
	}

	private void drawEdge(Graphics2D g2d, Edge e, int t, Color c) {
		Vertex v1 = e.getFirstVertex();
		Vertex v2 = e.getSecondVertex();
		c = c != null ? c : new Color(220, 220, 220, 127);
		g2d.setColor(c);
		g2d.setStroke(new BasicStroke(t));
		g2d.draw(new Line2D.Double(v1.getX(), v1.getY(), v2.getX(), v2.getY()));
	}

	private void drawEdgeLabel(Graphics2D g2d, Vertex v, Color c) {
		c = c != null ? c : new Color(220, 220, 220, 127);
		g2d.setColor(c);
		g2d.scale(1, -1);
		g2d.drawString(v.getId() + "", (float) (v.getX() - 8), (float) (-v.getY() - 6));
		g2d.scale(1, -1);
	}

	private void drawHole(Graphics2D g2d, Hole hole) {
		for (Edge e : hole.getEdges()) {
			drawEdge(g2d, e, 2, null);
		}
		for (Vertex v : hole.getVertices()) {
			drawVertex(g2d, v, 8, null);
		}
	}

	private void drawArea(Graphics2D g2d, artgallery.geometricalElements.Polygon galleryPolygon, Color c) {
		// Create a Java polygon and load the vertices contained in the gallery
		// (non-Java) Polygon.
		Polygon p = new Polygon();
		for (Vertex v : galleryPolygon.getVertices()) {
			p.addPoint((int) (v.getX()), (int) (v.getY()));
		}
		Area surface = new Area(p);

		// Subtract all the holes from the gallery surface
		for (Hole hole : galleryPolygon.getHoles()) {
			Polygon h = new Polygon();
			for (Vertex v : hole.getVertices()) {
				h.addPoint((int) (v.getX()), (int) (v.getY()));
			}
			Area currentHole = new Area(h);
			surface.subtract(currentHole);
		}

		if (c == null) {
			Color surfaceGray = new Color(45, 45, 45);
			c = surfaceGray;
		}
		g2d.setColor(c);
		g2d.fill(surface);
	}

	private void drawTriangulation(Graphics2D g2d, ArrayList<artgallery.geometricalElements.Polygon> triangulation) {
		for (artgallery.geometricalElements.Polygon p : triangulation) {
			int vertexSize = 4;
			Vertex v = new Vertex(gallery.getGallery(), p.getCentroid()[0], p.getCentroid()[1]);
			drawVertex(g2d, v, vertexSize, Color.GREEN);
			for (Edge e : p.getEdges()) {
				drawEdge(g2d, e, 2, new Color(255, 0, 0));
			}
		}

		/*
		 * //Delete for (Vertex vertex: gallery.getGallery().getVertices()){
		 * drawVertex(g2d, vertex, 12, null); g2d.scale(1, -1);
		 * g2d.drawString("" +
		 * gallery.getTriangulationVertices().indexOf(vertex),
		 * (float)vertex.getX() - 4, -(float)vertex.getY() - 15); g2d.scale(1,
		 * -1); }
		 */

	}

	private void drawGuard(Graphics2D g2d, Guard guard) {
		g2d.scale(0.2, -0.2);
		g2d.drawImage(guardIcon, guard.getX() * 5 - 50, -guard.getY() * 5 - 50, null);
		g2d.scale(5, -5);
	}

	private void drawThief(Graphics2D g2d, Thief thief) {
		g2d.scale(0.2, -0.2);
		g2d.drawImage(thiefIcon, thief.getX() * 5 - 50, -thief.getY() * 5 - 50, null);
		g2d.scale(5, -5);
	}

	private void drawLegend(Graphics2D g2d, double w, double h) {
		// Draw background rectangle
		Rectangle2D r = new Rectangle2D.Double(w / 2 - 185, h / 2 - 55, 180, 50);
		Color translucidGray = new Color(150, 150, 150, 120);
		g2d.translate((-offsetX) / scale, offsetY / scale);
		g2d.scale(1 / scale, -1 / scale);

		g2d.setColor(translucidGray);
		g2d.fill(r);
		g2d.draw(r);

		// Draw legend symbols
		int vertexSize = (int) (Math.round(10 * scale));
		Vertex symbol;
		symbol = new Vertex(gallery.getGallery(), (int) (w / 2 - 100), (int) ((h / 2 - 44)));
		drawVertex(g2d, symbol, vertexSize, null);

		symbol = new Vertex(gallery.getGallery(), (int) (w / 2 - 100), (int) ((h / 2 - 29)));
		symbol.setExit(true);
		drawVertex(g2d, symbol, vertexSize, null);

		symbol = new Vertex(gallery.getGallery(), (int) (w / 2 - 100), (int) ((h / 2 - 14)));
		symbol.setArt(true);
		drawVertex(g2d, symbol, vertexSize, null);

		g2d.scale(0.2, 0.2);
		g2d.drawImage(guardIcon, (int) (w / 2 - 175) * 5, (int) ((h / 2) - 26) * 5, null);
		g2d.drawImage(thiefIcon, (int) (w / 2 - 175) * 5, (int) ((h / 2) - 50) * 5, null);
		g2d.scale(5, 5);

		// Write Labels
		g2d.setColor(Color.WHITE);
		g2d.drawString("Normal Vertex", (float) (w / 2 - 90), (float) (h / 2 - 40));
		g2d.drawString("Exit Point", (float) (w / 2 - 90), (float) (h / 2 - 25));
		g2d.drawString("Art Piece", (float) (w / 2 - 90), (float) (h / 2 - 10));
		g2d.drawString("Guard", (float) (w / 2 - 150), (float) (h / 2 - 38));
		g2d.drawString("Thief", (float) (w / 2 - 150), (float) (h / 2 - 15));
		g2d.scale(1, -1);
	}

	private void drawVertexPositions(Graphics2D g2d, double w, double h) {
		if (true) {
			Paint paint = g2d.getPaint();

			g2d.scale(1, -1);
			g2d.translate(offsetX, offsetY);
			g2d.scale(scale, scale);

			g2d.setColor(Color.WHITE);

			for (final Vertex v : gallery.getGallery().getVerticesAndHoles()) {
				final int x = (int) v.getX();
				final int y = (int) (-v.getY());

				g2d.drawString(String.format("(%.1f,%.1f)", v.getX(), v.getY()), x, y);
			}

//			g2d.scale(1, -1);
			g2d.setPaint(paint);
		}
	}

	private void panView(MouseEvent e) {
		endDrag = e.getPoint();
		int dx = (int) (endDrag.getX() - startDrag.getX());
		int dy = (int) (endDrag.getY() - startDrag.getY());

//		if (offsetX + dx > getSize().getWidth() / 3) {
//			offsetX = (int) (getSize().getWidth() / 3);
//		} else if (offsetX + dx < -getSize().getWidth() / 3) {
//			offsetX = (int) (-getSize().getWidth() / 3);
//		} else {
//			offsetX += dx;
//		}
//
//		if (offsetY + dy > getSize().getHeight() / 3) {
//			offsetY = (int) (getSize().getHeight() / 3);
//		} else if (offsetY + dy < -getSize().getHeight() / 3) {
//			offsetY = (int) (-getSize().getHeight() / 3);
//		} else {
//			offsetY += dy;
//		}

		offsetX += dx;
		offsetY += dy;

		repaint();
	}

	public void resetView() {
		scale = 1;
		offsetX = 0;
		offsetY = 0;
	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent e) {
		scale -= (double) (e.getWheelRotation()) / 10;
		if (scale > 4) scale = 4;
		if (scale < 0.25) scale = 0.25;
		repaint();
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		Point p = e.getPoint();
		System.out.println((int) (p.getX() - this.getWidth() / 2) + ", " + (int) (-p.getY() + this.getHeight() / 2) + ", 0, 0");
	}

	@Override
	public void mouseEntered(MouseEvent e) {
	}

	@Override
	public void mouseExited(MouseEvent e) {
	}

	@Override
	public void mousePressed(MouseEvent e) {
		startDrag = e.getPoint();
		endDrag = null;
	}

	@Override
	public void mouseDragged(MouseEvent e) {
	}

	@Override
	public void mouseReleased(MouseEvent e) {
		panView(e);
	}

	@Override
	public void mouseMoved(MouseEvent e) {
	}

}
