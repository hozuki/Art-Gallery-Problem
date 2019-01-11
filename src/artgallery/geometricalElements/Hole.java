package artgallery.geometricalElements;
import java.util.ArrayList;

public class Hole {
	private ArrayList<Vertex> vertices;
	private ArrayList<Edge> edges;
	private int size;
	
	public Hole(){
		this.vertices = new ArrayList<Vertex>();
		edges = new ArrayList<Edge>();
	}
	
	public Hole(ArrayList<Vertex> vertices){
		this.vertices = new ArrayList<Vertex>();
		this.setVertices(vertices);
		edges = new ArrayList<Edge>();
		for(int i = 0; i < vertices.size(); ++i){
		    Vertex v1 = vertices.get(i);
		    Vertex v2 = vertices.get((i+1) % vertices.size());
		    edges.add(new Edge(v1, v2));
		}
	}

	public ArrayList<Vertex> getVertices() {
		return vertices;
	}

	public void setVertices(ArrayList<Vertex> vertices) {
		this.vertices = vertices;
	}
	
	public void addVertex(Vertex v){
		if(this.vertices.size() < this.size){
			this.vertices.add(v);
		}			
		if(this.vertices.size() == this.size){
			for(int i = 0; i < vertices.size(); ++i){
			    Vertex v1 = vertices.get(i);
			    Vertex v2 = vertices.get((i+1) % vertices.size());
			    Edge edge = new Edge(v1, v2);
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
}
