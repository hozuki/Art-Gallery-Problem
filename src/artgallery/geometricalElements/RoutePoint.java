package artgallery.geometricalElements;

public class RoutePoint {	
	private int X;
	private int Y;
	private int timeStamp;
	private int observing;
	
	public RoutePoint(){
	}
	public RoutePoint(int x, int y){
		this.X = x;
		this.Y = y;
	}
	public RoutePoint(int x, int y, int timeStamp, int observing){
		this.X = x;
		this.Y = y;
		this.setTimeStamp(timeStamp);
		this.setObserving(observing);
	}
	
	public Vertex toVertex(){
		return new Vertex(X, Y, -1, 0, 0);
	}
		
	public int getX() {
		return X;
	}
	public void setX(int x) {
		X = x;
	}
	
	public int getY() {
		return Y;
	}
	public void setY(int y) {
		Y = y;
	}
	public int getTimeStamp() {
		return timeStamp;
	}
	public void setTimeStamp(int timeStamp) {
		this.timeStamp = timeStamp;
	}
	public int isObserving() {
		return observing;
	}
	public void setObserving(int observing) {
		this.observing = observing;
	}	
}