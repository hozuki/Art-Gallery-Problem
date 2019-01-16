package artgallery.gui;

import javax.swing.filechooser.FileSystemView;
import java.io.File;
import java.io.IOException;

public class SingleRootFileSystemView extends FileSystemView {

	public SingleRootFileSystemView(final String path) {
		this(new File(path));
	}

	public SingleRootFileSystemView(final File path) {
		super();

		try {
			root = path.getCanonicalFile();
			roots[0] = root;
		} catch (IOException e) {
			throw new IllegalArgumentException(e);
		}

		if (!root.isDirectory()) {
			String message = root + " is not a directory";
			throw new IllegalArgumentException(message);
		}
	}

	@Override
	public File createNewFolder(File containingDir) {
		return null;
	}

	@Override
	public File getDefaultDirectory() {
		return root;
	}

	@Override
	public File getHomeDirectory() {
		return root;
	}

	@Override
	public File[] getRoots() {
		return roots.clone();
	}

	private final File root;
	private final File[] roots = new File[1];

}
