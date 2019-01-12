package artgallery;

import java.util.ArrayList;

import artgallery.Actors.Guard;
import artgallery.Actors.Thief;
import artgallery.geometricalElements.Edge;
import artgallery.geometricalElements.Polygon;
import artgallery.geometricalElements.Vertex;

public class GalleryModel {
	private ArrayList<Guard> guards = new ArrayList<Guard>();
	private ArrayList<Thief> thieves = new ArrayList<Thief>();
	private Polygon galleryPolygon;
	private boolean actors = false;
	private boolean legend = true;
	private boolean guardsPaths = false;
	private boolean thievesPaths = false;
	private boolean triangulation = false;
	private boolean visibility = false;
	private int globalTime;

	public GalleryModel() {
		this.setGallery(new Polygon());
	}

	public GalleryModel(Polygon polygon, int globalTime) {
		this.setGallery(polygon);
		this.setGlobalTime(globalTime);
	}

	public Polygon getGallery() {
		return galleryPolygon;
	}

	public void setGallery(Polygon galleryPolygon) {
		this.galleryPolygon = galleryPolygon;
	}

	public int getGlobalTime() {
		return globalTime;
	}

	public void setGlobalTime(int globalTime) {
		this.globalTime = globalTime;
	}

	public ArrayList<Edge> getTrapezoidalEdges() {
		GeometricAlgorithms GA = new GeometricAlgorithms();
		GA.trapezoidalDecomposition(galleryPolygon);
		return GA.trapezoidalEdges;
	}

	public ArrayList<Polygon> computeGuardVisibility(int i) {
		GeometricAlgorithms GA = new GeometricAlgorithms();
		Guard guard = this.guards.get(i);
		Vertex guardPosition = new Vertex(galleryPolygon, guard.getX(), guard.getY());
		guard.setVisibilityPolygon(GA.computeVisibilityPolygon(guardPosition, galleryPolygon));
		return guard.getVisibilityPolygon();
	}

	public void toggleActors() {
		this.actors = !this.actors;
	}

	public void toggleLegend() {
		this.legend = !this.legend;
	}

	public void toggleGuardsPaths() {
		this.guardsPaths = !this.guardsPaths;
	}

	public void toggleThievesPaths() {
		this.thievesPaths = !this.thievesPaths;
	}

	public void toggleTriangulation() {
		this.triangulation = !this.triangulation;
	}

	public void toggleVisibility() {
		this.visibility = !this.visibility;
	}

	public boolean showActors() {
		return actors;
	}

	public boolean showLegend() {
		return legend;
	}

	public boolean showGuardsPaths() {
		return guardsPaths;
	}

	public boolean showThievesPaths() {
		return thievesPaths;
	}

	public boolean showTriangulation() {
		return triangulation;
	}

	public boolean showVisibility() {
		return visibility;
	}

	public ArrayList<Guard> getGuards() {
		return guards;
	}

	public void setGuards(ArrayList<Guard> guards) {
		this.guards = guards;
	}

	public ArrayList<Thief> getThieves() {
		return thieves;
	}

	public void setThieves(ArrayList<Thief> thieves) {
		this.thieves = thieves;
	}
}
