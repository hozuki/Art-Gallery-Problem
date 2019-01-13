package artgallery;

import java.awt.Color;
import java.awt.Component;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.border.TitledBorder;

import artgallery.Actors.Guard;
import artgallery.Actors.Thief;
import artgallery.geometricalElements.Hole;
import artgallery.geometricalElements.Polygon;
import artgallery.geometricalElements.RoutePoint;
import artgallery.geometricalElements.Vertex;

public class ApplicationGUI implements Runnable {

	private JFrame frame;
	private JTextArea txtAreaStatus;
	private GalleryView canvas;
	private GalleryModel gallery;
	private JButton btnExecuteSimulation;
	private JButton btnRestart;
	private JTextField tfVertices;
	private JTextField tfEdges;
	private JTextField tfArt;
	private JTextField tfExits;
	private JTextField tfHoles;
	private JTextField tfGuards;
	private JTextField tfThieves;
	private JCheckBox cbDispActors;
	private JCheckBox cbDispGuardPt;
	private JCheckBox cbLoopPaths;
	private JCheckBox cbDispThiefPt;
	private JCheckBox cbDispTriangulation;
	private JCheckBox cbDispVisibility;
	private Thread animator;
	private boolean running = false;
	private int seconder = 0;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					ApplicationGUI window = new ApplicationGUI();
					window.frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the application.
	 */
	public ApplicationGUI() {
		initialize();
		txtAreaStatus.append("Simulator started\n");
		loadInput();
	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		frame = new JFrame();
		frame.setTitle("Group 3 - Art Gallery Simulator");
		frame.setBounds(100, 100, 1000, 700);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.getContentPane().setLayout(new BoxLayout(frame.getContentPane(), BoxLayout.X_AXIS));

		JSplitPane splitPane = new JSplitPane();
		splitPane.setEnabled(false);
		splitPane.setResizeWeight(1.0);
		splitPane.setDividerLocation(750);
		frame.getContentPane().add(splitPane);

		JSplitPane splitPane_1 = new JSplitPane();
		splitPane_1.setResizeWeight(1.0);
		splitPane_1.setOrientation(JSplitPane.VERTICAL_SPLIT);
		splitPane_1.setDividerLocation(500);
		splitPane.setLeftComponent(splitPane_1);

		canvas = new GalleryView();
		canvas.setBackground(Color.BLACK);
		splitPane_1.setLeftComponent(canvas);

		txtAreaStatus = new JTextArea();
		txtAreaStatus.setBackground(Color.WHITE);
		txtAreaStatus.setFont(new Font("Verdana", Font.PLAIN, 15));
		txtAreaStatus.setWrapStyleWord(true);
		txtAreaStatus.setText("Status:");
		txtAreaStatus.setEditable(false);
		JScrollPane scroll = new JScrollPane(txtAreaStatus);
		scroll.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);

		splitPane_1.setRightComponent(scroll);

		JPanel controlPanel = new JPanel();
		splitPane.setRightComponent(controlPanel);
		controlPanel.setLayout(null);

		JPanel informationPanel = new JPanel();
		informationPanel.setBounds(5, 0, 215, 220);
		informationPanel
				.setBorder(new TitledBorder(null, "Information", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		informationPanel.setLayout(null);

		JPanel geometryInformation = new JPanel();
		geometryInformation.setBounds(10, 21, 198, 123);
		geometryInformation
				.setBorder(new TitledBorder(null, "Geometry", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		informationPanel.add(geometryInformation);
		geometryInformation.setLayout(null);

		Box informationRow1 = Box.createHorizontalBox();
		informationRow1.setBounds(6, 16, 182, 20);
		geometryInformation.add(informationRow1);

		JLabel lblVertices = new JLabel("Vertices:");
		informationRow1.add(lblVertices);

		Component horizontalStrut = Box.createHorizontalStrut(16);
		informationRow1.add(horizontalStrut);

		tfVertices = new JTextField();
		tfVertices.setEditable(false);
		informationRow1.add(tfVertices);
		tfVertices.setColumns(10);

		Box informationRow2 = Box.createHorizontalBox();
		informationRow2.setBounds(6, 36, 182, 20);
		geometryInformation.add(informationRow2);

		JLabel lblEdges = new JLabel("Edges:");
		informationRow2.add(lblEdges);

		Component horizontalStrut_1 = Box.createHorizontalStrut(29);
		informationRow2.add(horizontalStrut_1);

		tfEdges = new JTextField();
		tfEdges.setEditable(false);
		informationRow2.add(tfEdges);
		tfEdges.setColumns(10);

		Box informationRow3 = Box.createHorizontalBox();
		informationRow3.setBounds(6, 56, 182, 20);
		geometryInformation.add(informationRow3);

		JLabel lblArt = new JLabel("Art:");
		informationRow3.add(lblArt);

		Component horizontalStrut_2 = Box.createHorizontalStrut(47);
		informationRow3.add(horizontalStrut_2);

		tfArt = new JTextField();
		tfArt.setAlignmentX(Component.LEFT_ALIGNMENT);
		tfArt.setEditable(false);
		informationRow3.add(tfArt);
		tfArt.setColumns(1);

		Box informationRow4 = Box.createHorizontalBox();
		informationRow4.setBounds(6, 76, 182, 20);
		geometryInformation.add(informationRow4);

		JLabel lblExits = new JLabel("Exits:");
		informationRow4.add(lblExits);

		Component horizontalStrut_3 = Box.createHorizontalStrut(36);
		informationRow4.add(horizontalStrut_3);

		tfExits = new JTextField();
		tfExits.setEditable(false);
		informationRow4.add(tfExits);
		tfExits.setColumns(10);

		Box InformationRow5 = Box.createHorizontalBox();
		InformationRow5.setBounds(6, 96, 182, 20);
		geometryInformation.add(InformationRow5);

		JLabel lblHoles = new JLabel("Holes:");
		InformationRow5.add(lblHoles);

		Component horizontalStrut_4 = Box.createHorizontalStrut(32);
		InformationRow5.add(horizontalStrut_4);

		tfHoles = new JTextField();
		tfHoles.setEditable(false);
		tfHoles.setColumns(10);
		InformationRow5.add(tfHoles);

		JPanel actorInformation = new JPanel();
		actorInformation.setBounds(10, 145, 198, 63);
		actorInformation
				.setBorder(new TitledBorder(null, "Actors", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		informationPanel.add(actorInformation);
		actorInformation.setLayout(null);

		Box informationRow6 = Box.createHorizontalBox();
		informationRow6.setBounds(6, 16, 182, 20);
		informationRow6.setAlignmentY(Component.CENTER_ALIGNMENT);
		actorInformation.add(informationRow6);

		JLabel lblGuards = new JLabel("Guards:");
		informationRow6.add(lblGuards);

		Component horizontalStrut_5 = Box.createHorizontalStrut(24);
		informationRow6.add(horizontalStrut_5);

		tfGuards = new JTextField();
		tfGuards.setEditable(false);
		tfGuards.setColumns(10);
		informationRow6.add(tfGuards);

		Box informationRow7 = Box.createHorizontalBox();
		informationRow7.setBounds(6, 36, 182, 20);
		actorInformation.add(informationRow7);

		JLabel lblThieves = new JLabel("Thieves:");
		informationRow7.add(lblThieves);

		Component horizontalStrut_6 = Box.createHorizontalStrut(21);
		informationRow7.add(horizontalStrut_6);

		tfThieves = new JTextField();
		tfThieves.setEditable(false);
		tfThieves.setColumns(10);
		informationRow7.add(tfThieves);
		controlPanel.add(informationPanel);

		JPanel optionsPanel = new JPanel();
		optionsPanel.setBounds(5, 220, 215, 212);
		controlPanel.add(optionsPanel);
		optionsPanel.setLayout(null);
		optionsPanel.setBorder(new TitledBorder(null, "Options", TitledBorder.LEADING, TitledBorder.TOP, null, null));

		JPanel panel_1 = new JPanel();
		panel_1.setLayout(null);
		panel_1.setBorder(new TitledBorder(null, "Route Options", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		panel_1.setBounds(6, 16, 202, 90);
		optionsPanel.add(panel_1);

		cbDispThiefPt = new JCheckBox("Display Thieves Paths");
		cbDispThiefPt.setEnabled(false);
		cbDispThiefPt.setBounds(6, 16, 155, 23);
		cbDispThiefPt.addActionListener(new ToggleDisplay("thief"));
		panel_1.add(cbDispThiefPt);

		cbDispGuardPt = new JCheckBox("Display Guards Paths");
		cbDispGuardPt.setEnabled(false);
		cbDispGuardPt.setBounds(6, 39, 155, 23);
		cbDispGuardPt.addActionListener(new ToggleDisplay("guard"));
		panel_1.add(cbDispGuardPt);

		cbLoopPaths = new JCheckBox("Loop Actor Routes");
		cbLoopPaths.setEnabled(false);
		cbLoopPaths.setBounds(6, 62, 155, 23);
		panel_1.add(cbLoopPaths);

		JPanel panel_2 = new JPanel();
		panel_2.setLayout(null);
		panel_2.setBorder(
				new TitledBorder(null, "Display Options", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		panel_2.setBounds(6, 110, 202, 90);
		optionsPanel.add(panel_2);

		cbDispTriangulation = new JCheckBox("Show Gallery Triangulation");
		cbDispTriangulation.setEnabled(false);
		cbDispTriangulation.setBounds(6, 35, 190, 23);
		cbDispTriangulation.addActionListener(new ToggleDisplay("triangulation"));
		panel_2.add(cbDispTriangulation);

		cbDispVisibility = new JCheckBox("Show Visibility Polygons");
		cbDispVisibility.setEnabled(false);
		cbDispVisibility.setBounds(6, 55, 190, 23);
		cbDispVisibility.addActionListener(new ToggleDisplay("visibility"));
		panel_2.add(cbDispVisibility);

		cbDispActors = new JCheckBox("Show Actors");
		cbDispActors.setEnabled(false);
		cbDispActors.setBounds(6, 15, 153, 23);
		cbDispActors.addActionListener(new ToggleDisplay("actors"));
		panel_2.add(cbDispActors);

		JPanel executionBtnsPanel = new JPanel();
		executionBtnsPanel.setBounds(0, 443, 228, 116);
		controlPanel.add(executionBtnsPanel);
		executionBtnsPanel.setLayout(null);

		JButton btnLoadFile = new JButton("Load File");
		btnLoadFile.setBounds(54, 5, 120, 30);
		executionBtnsPanel.add(btnLoadFile);
		btnLoadFile.addActionListener(new LoadInputFile());
		btnLoadFile.setAlignmentX(Component.CENTER_ALIGNMENT);

		btnExecuteSimulation = new JButton("Execute");
		btnExecuteSimulation.setBounds(54, 35, 120, 30);
		executionBtnsPanel.add(btnExecuteSimulation);
		btnExecuteSimulation.setEnabled(false);
		btnExecuteSimulation.addActionListener(new ExecuteAnimation());
		btnExecuteSimulation.setAlignmentX(Component.CENTER_ALIGNMENT);

		btnRestart = new JButton("Restart");
		btnRestart.setBounds(54, 65, 120, 30);
		btnRestart.addActionListener(new RestartGallery());
		btnRestart.setEnabled(false);
		btnRestart.setAlignmentX(0.5f);
		executionBtnsPanel.add(btnRestart);
	}

	private static String paths(final String first, final String... more) {
		return Paths.get(first, more).toString();
	}

	private void loadInput() {
		String directory = paths(System.getProperty("user.dir"), "data/sampleInputs/");
		JFileChooser fileChooser = new JFileChooser(directory);
		BufferedReader br = null;
		String line = "";
		String separator = ",";
		int lastHole = 0;
		int lastIndex = 0;

		//Array to add up objects while reading and later initialize the gallery.
		ArrayList<Vertex> vertices = new ArrayList<Vertex>();
		ArrayList<Hole> holes = new ArrayList<Hole>();
		ArrayList<Guard> guards = new ArrayList<Guard>();
		ArrayList<Thief> thieves = new ArrayList<Thief>();

		// Load Gallery Structure File
		//if (fileChooser.showOpenDialog(fileChooser) == JFileChooser.APPROVE_OPTION) {
		if (true) {
			//String fileName = fileChooser.getSelectedFile().getName();
			//fileNumber = fileName.charAt(fileName.length() - 1);
			Polygon galleryPolygon = new Polygon();

			char fileNumber = '3';
			String fileName = "AGS";
			try {
				txtAreaStatus.append("Attempting to load gallery file [" + paths(directory, fileName + fileNumber) + "]\n");
				br = new BufferedReader(new FileReader(paths(directory, fileName + fileNumber)));
				int count = 1;
				int n, h, e, a;
				int g, v, t, d;
				n = h = e = a = g = v = t = d = 0;

				while ((line = br.readLine()) != null) {
					String[] values = line.replaceAll(" ", "").split(separator);
					if (count == 1) {
						if (values.length != 4) {
							throw new IllegalArgumentException();
						}
						n = Integer.parseInt(values[0]);
						h = Integer.parseInt(values[1]);
						e = Integer.parseInt(values[2]);
						a = Integer.parseInt(values[3]);
						if (n < 3 || h < 0 || e < 1 || a < 1) {
							throw new IllegalArgumentException();
						}
					} else if (count == 2) {
						if (h > 0) {
							for (int i = 0; i < values.length; ++i) {
								Hole tempHole = new Hole(galleryPolygon);
								tempHole.setSize(Integer.parseInt(values[i]));
								holes.add(tempHole);
							}
						}
					} else if (count == 3) {
						g = Integer.parseInt(values[0]);
						v = Integer.parseInt(values[1]);
						t = Integer.parseInt(values[2]);
						d = Integer.parseInt(values[3]);
						if (g < 1 || v < 1 || t < 1 || d < 0) {
							throw new IllegalArgumentException();
						}
					} else if (count > 3 && count <= 3 + n) {
						Vertex tempVertex = new Vertex(galleryPolygon, Integer.parseInt(values[0]), Integer.parseInt(values[1]),
								count - 3, Integer.parseInt(values[2]), Integer.parseInt(values[3]));
						vertices.add(tempVertex);
					} else {
						Hole hole = holes.get(lastHole);
						if (count > 3 + n + lastIndex + hole.getSize()) {
							lastIndex += hole.getSize();
							lastHole++;
							hole = holes.get(lastHole);
						}
						Vertex tempVertex = new Vertex(galleryPolygon, Integer.parseInt(values[0]), Integer.parseInt(values[1]),
								count - 3, Integer.parseInt(values[2]), 0);
						hole.addVertex(tempVertex);

					}
					count++;
				}

				galleryPolygon.setVertices(vertices);
				galleryPolygon.setHoles(holes);
				galleryPolygon.generateEdges();

				for (Hole hole : holes) {
					// CCW to CW
					hole.reverseVerticesAndEdges();
				}

				galleryPolygon.fixVertexNeighbors();

				txtAreaStatus.append("File succesfully loaded!\n");
				this.gallery = new GalleryModel(galleryPolygon, t);

			} catch (FileNotFoundException e) {
				txtAreaStatus.append("No gallery structure file could be found with that name :(\n");
			} catch (IOException e) {
				e.printStackTrace();
			} catch (IllegalArgumentException e) {
				System.out.println("Incorrect input file format");
			} finally {
				if (br != null) {
					try {
						br.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			}

			// Load Guard Schedule File
			try {
				txtAreaStatus.append("Attempting to load guard schedule specification file...\n");
				br = new BufferedReader(new FileReader(paths(directory, "GSS" + fileNumber)));
				int count = 1;
				int g = 0;
				int currentGuard = -1;
				int routeLength = 0;
				int lastDef = 0;

				while ((line = br.readLine()) != null) {
					String[] values = line.replaceAll(" ", "").split(separator);
					if (count == 1) {
						if (values.length != 1) {
							throw new IllegalArgumentException();
						}
						g = Integer.parseInt(values[0]);
						for (int i = 0; i < g; ++i) {
							Guard tempGuard = new Guard();
							guards.add(tempGuard);
						}
					} else if (count == 2 || count > lastDef + guards.get(currentGuard).getRouteLength() + 1) {
						if (values.length != 1) {
							throw new IllegalArgumentException();
						}
						currentGuard++;
						lastDef = count;
						routeLength = Integer.parseInt(values[0]);
						guards.get(currentGuard).setRouteLength(routeLength);
					} else {
						RoutePoint tempPoint = new RoutePoint(Integer.parseInt(values[0]), Integer.parseInt(values[1]),
								Integer.parseInt(values[2]), Integer.parseInt(values[3]));
						guards.get(currentGuard).addPointToRoute(tempPoint);
					}
					count++;
				}
				gallery.setGuards(guards);
				txtAreaStatus.append("Guard schedule succesfully loaded!\n");
			} catch (FileNotFoundException e) {
				txtAreaStatus.append("No guard schedule file could be found with that name :(\n");
			} catch (IOException e) {
				e.printStackTrace();
			} catch (IllegalArgumentException e) {
				System.out.println("Incorrect input file format");
			} finally {
				if (br != null) {
					try {
						br.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			}

			// Load Thief Path File
			try {
				txtAreaStatus.append("Attempting to load thief path specification file...\n");
				br = new BufferedReader(new FileReader(paths(directory, "RPS" + fileNumber)));
				int count = 1;

				while ((line = br.readLine()) != null) {
					String[] values = line.replaceAll(" ", "").split(separator);
					if (count == 1) {
						if (values.length != 1) {
							throw new IllegalArgumentException();
						}
						thieves.add(new Thief());
					} else {
						RoutePoint tempPoint = new RoutePoint(Integer.parseInt(values[0]), Integer.parseInt(values[1]),
								Integer.parseInt(values[2]), 1);
						thieves.get(0).addPointToRoute(tempPoint);
					}
					count++;
				}
				gallery.setThieves(thieves);
				txtAreaStatus.append("Filthy thieves succesfully loaded!\n");
			} catch (FileNotFoundException e) {
				txtAreaStatus.append("No robber path file could be found with that name :(\n");
			} catch (IOException e) {
				e.printStackTrace();
			} catch (IllegalArgumentException e) {
				System.out.println("Incorrect input fileNumber format");
			} finally {
				if (br != null) {
					try {
						br.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			}
		}

		refreshCanvas();
		resetControls();
		enableControls();
	}

	private void refreshCanvas() {
		canvas.updateGallery(gallery);
		canvas.resetView();
		tfVertices.setText("" + gallery.getGallery().getVertices().size());
		tfEdges.setText("" + gallery.getGallery().getEdges().size());
		tfArt.setText("" + gallery.getGallery().getVertices().stream().filter(v -> v.isArt()).count());
		tfExits.setText("" + gallery.getGallery().getVertices().stream().filter(v -> v.isExit()).count());
		tfHoles.setText("" + gallery.getGallery().getHoles().size());
		tfGuards.setText("" + gallery.getGuards().size());
		tfThieves.setText("" + gallery.getThieves().size());
	}

	private void resetControls() {
		cbDispTriangulation.setSelected(false);
		cbDispActors.setSelected(false);
		cbDispGuardPt.setSelected(false);
		cbLoopPaths.setSelected(false);
		cbDispThiefPt.setSelected(false);
		cbDispTriangulation.setSelected(false);
		cbDispVisibility.setSelected(false);
	}

	public void toggleAnimation() {
		if (!running) {
			running = true;
			animator = new Thread(this);
			animator.start();
		} else {
			running = false;
			animator.interrupt();
		}
	}

	@Override
	public void run() {
		long beforeTime, timeDiff, sleep;
		beforeTime = System.currentTimeMillis();
		while (running) {
			if (seconder <= gallery.getGlobalTime() * 100) {
				advanceActors();
				canvas.repaint();

				timeDiff = System.currentTimeMillis() - beforeTime;
				sleep = 10 - timeDiff;

				if (sleep < 0) {
					sleep = 2;
				}
				try {
					Thread.sleep(sleep);
				} catch (InterruptedException e) {
				}

				beforeTime = System.currentTimeMillis();
				if (seconder % 100 == 0) {
					txtAreaStatus.append("Executing [" + seconder / 100 + "/" + gallery.getGlobalTime() + "]\n");
				}
				txtAreaStatus.setCaretPosition(txtAreaStatus.getDocument().getLength());
				seconder++;
			} else {
				txtAreaStatus.append("Global time ended\n");
				btnExecuteSimulation.setText("Execute");
				toggleAnimation();
				seconder = 0;
			}
		}
	}

	public boolean isRunning() {
		return running;
	}

	private void advanceActors() {
		for (Guard guard : gallery.getGuards()) {
			guard.advanceWithInterpolation(seconder);
		}

		for (Thief thief : gallery.getThieves()) {
			thief.advanceWithInterpolation(seconder);
		}
	}

	private class LoadInputFile implements ActionListener {
		@Override
		public void actionPerformed(ActionEvent e) {
			loadInput();
			enableControls();
		}
	}

	private void enableControls() {
		cbDispGuardPt.setEnabled(true);
		cbDispThiefPt.setEnabled(true);
		cbDispTriangulation.setEnabled(true);
		cbLoopPaths.setEnabled(true);
		cbDispActors.setEnabled(true);
		cbDispVisibility.setEnabled(true);
		btnExecuteSimulation.setEnabled(true);
		btnRestart.setEnabled(true);
	}

	private class ToggleDisplay implements ActionListener {
		private String type;

		public ToggleDisplay(String t) {
			this.type = t;
		}

		@Override
		public void actionPerformed(ActionEvent e) {
			switch (type) {
				case "actors":
					gallery.toggleActors();
					break;
				case "guard":
					gallery.toggleGuardsPaths();
					break;
				case "legend":
					break;
				case "thief":
					gallery.toggleThievesPaths();
					break;
				case "triangulation":
					gallery.toggleTriangulation();
					break;
				case "visibility":
					gallery.toggleVisibility();
					break;
			}
			canvas.repaint();
		}
	}

	private class ExecuteAnimation implements ActionListener {
		@Override
		public void actionPerformed(ActionEvent e) {
			if (!running) {
				btnExecuteSimulation.setText("Pause");
				txtAreaStatus.append("Simulation started\n");
			} else {
				btnExecuteSimulation.setText("Execute");
				txtAreaStatus.append("Simulation paused\n");
			}
			toggleAnimation();
		}
	}

	private class RestartGallery implements ActionListener {
		@Override
		public void actionPerformed(ActionEvent e) {
			seconder = 0;
			animator = new Thread();
			loadInput();
		}
	}
}
