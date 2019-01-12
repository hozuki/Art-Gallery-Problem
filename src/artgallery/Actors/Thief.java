package artgallery.Actors;

import java.util.ArrayList;

import artgallery.geometricalElements.Edge;
import artgallery.geometricalElements.Polygon;
import artgallery.geometricalElements.RoutePoint;

public class Thief extends Actor {
	private int x;
	private int y;
	private int currentPointIndex;
	private int orientation;
	private int maxSpeed;
	private int observationTime;
	private ArrayList<RoutePoint> routePoints = new ArrayList<RoutePoint>();
	private ArrayList<Edge> routeEdges = new ArrayList<Edge>();
	private int routeLength;

	public Thief() {
		routePoints = new ArrayList<RoutePoint>();
	}

	public Thief(RoutePoint r, int ms, int ot) {
		this.maxSpeed = ms;
		this.observationTime = ot;
		routePoints = new ArrayList<RoutePoint>();
		routePoints.add(r);
		currentPointIndex = 0;
		this.setPosition(routePoints.get(currentPointIndex));
	}

	public void addPointToRoute(RoutePoint r) {
		routePoints.add(r);
		if (routePoints.size() == 1) {
			this.setX(routePoints.get(0).getX());
			this.setY(routePoints.get(0).getY());
		}
	}

	public int getRouteLength() {
		return routeLength;
	}

	public void setRouteLength(int routeLength) {
		this.routeLength = routeLength;
	}

	public ArrayList<RoutePoint> getRoutePoints() {
		return routePoints;
	}

	public ArrayList<Edge> getRouteEdges(Polygon polygon) {
		if (routePoints.size() > 1) {
			for (int i = 0; i < routePoints.size() - 1; ++i) {
				RoutePoint r1 = routePoints.get(i);
				RoutePoint r2 = routePoints.get((i + 1));
				routeEdges.add(new Edge(polygon, r1.toVertex(polygon), r2.toVertex(polygon)));
			}
		}
		return routeEdges;
	}

	public void advance(int time) {
		if (currentPointIndex < routePoints.size() - 1) {
			RoutePoint nextPoint = routePoints.get(currentPointIndex + 1);
			if (nextPoint.getTimeStamp() == time) {
				setPosition(nextPoint);
				currentPointIndex = (currentPointIndex + 1);
			}
		}
	}

	public void advanceWithInterpolation(int time) {
		RoutePoint nextPoint = getNext();
		if (nextPoint != null) {
			RoutePoint currentPoint = routePoints.get(currentPointIndex);

			int step = time - currentPoint.getTimeStamp() * 100;
			int interpolatedX = currentPoint.getX() + (nextPoint.getX() - currentPoint.getX()) / 100 * step;
			int interpolatedY = currentPoint.getY() + (nextPoint.getY() - currentPoint.getY()) / 100 * step;
			setPosition(interpolatedX, interpolatedY);

			if (time == nextPoint.getTimeStamp() * 100 && time > 0) {
				currentPointIndex++;
			}
		}
	}

	private RoutePoint getNext() {
		return (currentPointIndex < routePoints.size() - 1 ? routePoints.get(currentPointIndex + 1) : null);
	}

	private void setPosition(RoutePoint r) {
		this.setX(r.getX());
		this.setY(r.getY());
	}

	private void setPosition(int x, int y) {
		this.setX(x);
		this.setY(y);
	}

	public int getX() {
		return x;
	}

	public void setX(int x) {
		this.x = x;
	}

	public int getY() {
		return y;
	}

	public void setY(int y) {
		this.y = y;
	}
}
