using Eto.Drawing;
using Eto.Forms;
using Eto.Veldrid;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using Veldrid;
using Veldrid.SPIRV;
using KDTree;
using System.Threading.Tasks;
using Eto;

namespace VeldridEto;

public struct VertexPositionColor
{
	private static uint _SizeInBytes = (uint)Marshal.SizeOf(typeof(VertexPositionColor));

	public static uint SizeInBytes => _SizeInBytes;

	private Vector3 Position;
	private RgbaFloat Color;

	public VertexPositionColor(Vector3 position, RgbaFloat color)
	{
		Position = position;
		Color = color;
	}
}

/// <summary>
/// A class that controls rendering to a VeldridSurface.
/// </summary>
/// <remarks>
/// VeldridSurface is only a basic control that lets you render to the screen
/// using Veldrid. How exactly to do that is up to you; this driver class is
/// only one possible approach, and in all likelihood not the most efficient.
/// </remarks>
public class VeldridDriver
{
	public string ExecutableDirectory { get; set; }
	public string ShaderSubdirectory { get; set; }

	public VeldridSurface Surface { get; set; }

	public UITimer Clock = new();

	public delegate void updateHost();
	public updateHost updateHostFunc { get; set; }

	public delegate void updateHostSelection(int index);
	public updateHostSelection updateHostSelectionFunc { get; set; }

	public bool ok { get; set; }
	public bool savedLocation_valid { get; set; }
	private PointF savedLocation;

	private uint[] polyFirst;
	private uint[] polyVertexCount;

	private uint[] tessFirst;
	private uint[] tessVertexCount;

	private uint[] lineFirst;
	private uint[] lineVertexCount;

	private uint[] pointsFirst;

	private uint[] gridIndices;

	private uint[] axesIndices;

	public OVPSettings ovpSettings;

	private CommandList CommandList;
	private DeviceBuffer GridVertexBuffer;
	private DeviceBuffer GridIndexBuffer;
	private DeviceBuffer AxesVertexBuffer;
	private DeviceBuffer AxesIndexBuffer;

	private DeviceBuffer LinesVertexBuffer;
	private DeviceBuffer PointsVertexBuffer;
	private DeviceBuffer PolysVertexBuffer;
	private DeviceBuffer TessVertexBuffer;

	private Pipeline PointsPipeline;
	private Pipeline LinePipeline;
	private Pipeline LinesPipeline;
	private Pipeline FilledPipeline;

	private Matrix4x4 ModelMatrix = Matrix4x4.Identity;
	private DeviceBuffer ModelBuffer;
	private ResourceSet ModelMatrixSet;

	private Matrix4x4 ViewMatrix;
	private DeviceBuffer ViewBuffer;
	private ResourceSet ViewMatrixSet;

	private bool Ready;
	private const float pointWidth = 0.50f;
	private bool hasFocus;
	private bool keyHandlerApplied;

	public VeldridDriver(ref OVPSettings svpSettings, ref VeldridSurface surface)
	{
		try
		{
			ovpSettings = svpSettings;
			Surface = surface;
			Surface.MouseDown += downHandler;
			Surface.MouseMove += dragHandler;
			Surface.MouseUp += upHandler;
			Surface.MouseWheel += zoomHandler;
			Surface.GotFocus += addKeyHandler;
			Surface.MouseEnter += setFocus;
			Surface.LostFocus += removeKeyHandler;
		}
		catch (Exception)
		{

		}
		Clock.Interval = 1.0f / 60.0f;
		Clock.Elapsed += Clock_Elapsed;

		Surface.Resize += (sender, e) =>
		{
			if (Surface.Visible && (e.Width > 0) && (e.Height > 0))
			{
				pUpdateViewport();
			}
		};
	}
		
	public void setFocus(object sender, EventArgs e)
	{
		switch (hasFocus)
		{
			case true:
				return;
			default:
				Surface.Focus();
				hasFocus = true;
				break;
		}
	}

	private void Clock_Elapsed(object sender, EventArgs e)
	{
		switch (ovpSettings.changed)
		{
			case false:
				return;
		}
		ovpSettings.changed = false;

		Draw();

		Surface.Invalidate();
	}

	private DateTime CurrentTime;
	private DateTime PreviousTime = DateTime.Now;

	private float axisZ;
	private float gridZ;

	// Use for drag handling.
	public bool dragging { get; set; }
	private float x_orig;
	private float y_orig;

	private ContextMenu menu;

	public void setContextMenu(ref ContextMenu menu_)
	{
		menu = menu_;
	}

	public void changeSettingsRef(ref OVPSettings newSettings)
	{
		ovpSettings = newSettings;
		updateViewport();
	}

	private Point WorldToScreen(float x, float y)
	{
		// int oX = (int)((x - ovpSettings.getCameraX() / (ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom())) + Surface.RenderWidth / 2);

		double oX_2 = (double)Surface.RenderWidth / 2;
		double oX_3 = ovpSettings.getCameraX() / (ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom());
		double oX_4 = x;

		int oXC = (int)(oX_4 - oX_3 + oX_2);

		// int oY = (int)((y - ovpSettings.getCameraY() / (ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom())) + Surface.RenderHeight / 2);

		double oY_2 = (double)Surface.RenderHeight / 2;
		double oY_3 = ovpSettings.getCameraY() / (ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom());
		double oY_4 = y;

		int oYC = (int)(oY_4 - oY_3 + oY_2);

		return new Point(oXC, oYC);
	}

	private Size WorldToScreen(SizeF pt)
	{
		Point pt1 = WorldToScreen(0, 0);
		Point pt2 = WorldToScreen(pt.Width, pt.Height);
		return new Size(pt2.X - pt1.X, pt2.Y - pt1.Y);
	}

	private SizeF ScreenToWorld(SizeF pt)
	{
		PointF pt1 = ScreenToWorld(0, 0);
		PointF pt2 = ScreenToWorld(pt.Width, pt.Height);
		return new SizeF(pt2.X - pt1.X, pt2.Y - pt1.Y);
	}

	private PointF ScreenToWorld(float x, float y)
	{
		double oX_2 = (double)Surface.RenderWidth / 2;
		double oX_3 = ovpSettings.getCameraX() / (ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom());
		double oX_4 = x;

		double oXC = oX_4 - oX_2 + oX_3;


		return new PointF((x - (float)Surface.RenderWidth / 2) * (ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom()) + ovpSettings.getCameraX(),
			((float)Surface.RenderHeight / 2 - y) * (ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom()) + ovpSettings.getCameraY());
	}

	private void downHandler(object sender, MouseEventArgs e)
	{
		switch (e.Buttons)
		{
			case MouseButtons.Primary:
				setDown(e.Location.X, e.Location.Y);
				break;
		}
		if (e.Buttons == MouseButtons.Middle || e.Modifiers == Keys.Control && e.Buttons == MouseButtons.Primary)
		{
			selectByClick(e.Location.X, e.Location.Y);
		}
		e.Handled = true;
	}

	private void setDown(float x, float y)
	{
		switch (dragging)
		{
			// might not be needed, but seemed like a safe approach to avoid re-setting these in a drag event.
			case false when !ovpSettings.isLocked():
				x_orig = x;
				y_orig = y;
				dragging = true;
				break;
		}
	}

	public void saveLocation()
	{
		savedLocation = new PointF(ovpSettings.getCameraX(), ovpSettings.getCameraY());
		savedLocation_valid = true;
	}

	public void zoomExtents(int index)
	{
		getExtents(index);

		if (ovpSettings.polyList.Count == 0 && ovpSettings.lineList.Count == 0 ||
		    ovpSettings.minX == 0 && ovpSettings.maxX == 0 ||
		    ovpSettings.minY == 0 && ovpSettings.maxY == 0)
		{
			reset();
			return;
		}

		// Locate camera at center of the polygon field.
		float dX = ovpSettings.maxX - ovpSettings.minX;
		float dY = ovpSettings.maxY - ovpSettings.minY;
		float cX = dX / 2.0f + ovpSettings.minX;
		float cY = dY / 2.0f + ovpSettings.minY;

		// Now need to get the zoom level organized.
		float zoomLevel_x = dX / Surface.Width;
		float zoomLevel_y = dY / Surface.Height;

		if (zoomLevel_x > zoomLevel_y)
		{
			ovpSettings.setZoomFactor(zoomLevel_x / ovpSettings.getBaseZoom());
		}
		else
		{
			ovpSettings.setZoomFactor(zoomLevel_y / ovpSettings.getBaseZoom());
		}

		goToLocation(cX, cY);
	}

	public void loadLocation()
	{
		switch (savedLocation_valid)
		{
			case true:
				ovpSettings.setCameraPos(savedLocation.X, savedLocation.Y);
				updateViewport();
				break;
		}
	}

	public void goToLocation(float x, float y)
	{
		ovpSettings.setCameraPos(x, y);
		updateViewport();
	}

	private void dragHandler(object sender, MouseEventArgs e)
	{
		if (!ovpSettings.isLocked())
		{
			switch (e.Buttons)
			{
				case MouseButtons.Primary:
				{
					PointF scaledLocation = e.Location * Surface.ParentWindow.LogicalPixelSize;

					switch (dragging)
					{
						case false:
							setDown(scaledLocation.X, scaledLocation.Y);
							break;
					}
					object locking = new();
					lock (locking)
					{
						float new_X = ovpSettings.getCameraX() - (scaledLocation.X - x_orig) * ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom();
						float new_Y = ovpSettings.getCameraY() + (scaledLocation.Y - y_orig) * ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom();
						ovpSettings.setCameraPos(new_X, new_Y);
						x_orig = scaledLocation.X;
						y_orig = scaledLocation.Y;
					}

					break;
				}
			}

			updateViewport();
		}
		e.Handled = true;
	}

	public void freeze_thaw()
	{
		ovpSettings.lockVP(!ovpSettings.isLocked());
		updateHostFunc?.Invoke();
	}
	
	private void selectByClick(float x, float y)
	{
		// Where did we click?
		PointF scaledLocation = new(x, y);
		scaledLocation = ScreenToWorld(scaledLocation.X * Surface.ParentWindow.LogicalPixelSize, scaledLocation.Y * Surface.ParentWindow.LogicalPixelSize);

		PointF cPos = ovpSettings.getCameraPos();

		// Populate our tree.
		int polyCount = ovpSettings.polyList.Count;
		switch (polyCount)
		{
			case > 0:
			{
				double[] distances = new double[polyCount];
				int[] indices = new int[polyCount];
				ParallelOptions po = new();
				//Parallel.For(0, polyCount, po, (poly, loopstate) =>
				for (int poly = 0; poly < ovpSettings.polyList.Count; poly++)
				{
					KDTree<PointF> pTree = new(2, ovpSettings.polyListPtCount[poly] + 1); // add one for the midpoint.
					foreach (PointF t1 in ovpSettings.polyList[poly].poly)
					{
						PointF t = new(t1.X, t1.Y);
						pTree.AddPoint(new double[] { t.X, t.Y }, t);
					}

					double maxX = ovpSettings.polyList[poly].poly.Max(p => p.X);
					double minX = ovpSettings.polyList[poly].poly.Min(p => p.X);
					double maxY = ovpSettings.polyList[poly].poly.Max(p => p.Y);
					double minY = ovpSettings.polyList[poly].poly.Min(p => p.Y);

					double deltaX = (maxX - minX) * 0.5f;
					double deltaY = (maxY - minY) * 0.5f;

					PointF midPoint = new((float)(minX + deltaX), (float)(minY + deltaY));
					pTree.AddPoint(new double[] { midPoint.X, midPoint.Y }, midPoint);

					// '1' forces a single nearest neighbor to be returned.
					NearestNeighbour<PointF> pIter = pTree.NearestNeighbors(new double[] { scaledLocation.X, scaledLocation.Y }, 1);
					while (pIter.MoveNext())
					{
						distances[poly] = Math.Abs(pIter.CurrentDistance);
						indices[poly] = ovpSettings.polySourceIndex[poly];
					}
				}
				//);

				int selIndex = indices[Array.IndexOf(distances, distances.Min())];

				updateHostSelectionFunc?.Invoke(selIndex);
				break;
			}
			default:
			{
				// Populate our tree.
				int lineCount = ovpSettings.lineList.Count;
				switch (lineCount)
				{
					case > 0:
					{
						double[] distances = new double[lineCount];
						int[] indices = new int[lineCount];
						ParallelOptions po = new();
						//Parallel.For(0, lineCount, po, (line, loopstate) =>
						for (int line = 0; line < ovpSettings.lineList.Count; line++)
						{
							KDTree<PointF> pTree = new(2, ovpSettings.lineListPtCount[line] + 1); // add one for the midpoint.
							foreach (PointF t1 in ovpSettings.lineList[line].poly)
							{
								PointF t = new(t1.X, t1.Y);
								pTree.AddPoint(new double[] { t.X, t.Y }, t);
							}

							double maxX = ovpSettings.lineList[line].poly.Max(p => p.X);
							double minX = ovpSettings.lineList[line].poly.Min(p => p.X);
							double maxY = ovpSettings.lineList[line].poly.Max(p => p.Y);
							double minY = ovpSettings.lineList[line].poly.Min(p => p.Y);

							double deltaX = (maxX - minX) * 0.5f;
							double deltaY = (maxY - minY) * 0.5f;

							PointF midPoint = new((float)(minX + deltaX), (float)(minY + deltaY));
							pTree.AddPoint(new double[] { midPoint.X, midPoint.Y }, midPoint);

							// '1' forces a single nearest neighbor to be returned.
							NearestNeighbour<PointF> pIter = pTree.NearestNeighbors(new double[] { scaledLocation.X, scaledLocation.Y }, 1);
							while (pIter.MoveNext())
							{
								distances[line] = Math.Abs(pIter.CurrentDistance);
								indices[line] = ovpSettings.lineSourceIndex[line];
							}
						}
						//);

						int selIndex = indices[Array.IndexOf(distances, distances.Min())];

						updateHostSelectionFunc?.Invoke(selIndex);
						break;
					}
				}

				break;
			}
		}
	}

	private void upHandler(object sender, MouseEventArgs e)
	{
		switch (e.Buttons)
		{
			case MouseButtons.Alternate:
			{
				menu?.Show(Surface);

				break;
			}
		}
		if (ovpSettings.isLocked())
		{
			return;
		}

		dragging = e.Buttons switch
		{
			MouseButtons.Primary => false,
			_ => dragging
		};
		e.Handled = true;
	}

	public void zoomIn(float delta)
	{
		ovpSettings.setZoomFactor(ovpSettings.getZoomFactor() + ovpSettings.getZoomStep() * 0.01f * delta);
		updateHostFunc?.Invoke();
	}

	public void zoomOut(float delta)
	{
		ovpSettings.setZoomFactor(ovpSettings.getZoomFactor() - ovpSettings.getZoomStep() * 0.01f * delta);
		updateHostFunc?.Invoke();
	}

	public void fastZoomIn(float delta)
	{
		ovpSettings.setZoomFactor(ovpSettings.getZoomFactor() * calcZoom(delta));
		updateHostFunc?.Invoke();
	}

	public void fastZoomOut(float delta)
	{
		ovpSettings.setZoomFactor(ovpSettings.getZoomFactor() / calcZoom(delta));
		updateHostFunc?.Invoke();
	}

	private float calcZoom (float delta)
	{
		float f = Math.Abs(delta) * 0.1f;
		f = delta switch
		{
			< 0 => 1.0f / f,
			_ => f
		};

		return f;
	}

	private void panVertical(float delta)
	{
		ovpSettings.setCameraY(ovpSettings.getCameraY() + delta / 10);
	}

	private void panHorizontal(float delta)
	{
		ovpSettings.setCameraX(ovpSettings.getCameraX() + delta / 10);
	}

	private void addKeyHandler(object sender, EventArgs e)
	{
		switch (keyHandlerApplied)
		{
			case true:
				return;
			default:
				Surface.KeyDown += keyHandler;
				keyHandlerApplied = true;
				break;
		}
	}

	private void removeKeyHandler(object sender, EventArgs e)
	{
		switch (keyHandlerApplied)
		{
			case false:
				return;
		}
		hasFocus = false;
		Surface.KeyDown -= keyHandler;
		keyHandlerApplied = false;
	}

	public void reset()
	{
		ovpSettings.resetCamera();
	}

	private void keyHandler(object sender, KeyEventArgs e)
	{
		//e.Handled = true;

		if (ovpSettings.isLocked())
		{
			if (e.Key != Keys.F)
			{
				return;
			}
		}

		switch (e.Key)
		{
			case Keys.F:
				ovpSettings.lockVP(!ovpSettings.isLocked());
				break;
			case Keys.R:
				reset();
				break;
		}

		float stepping = 10.0f * ovpSettings.getZoomFactor();

		bool doUpdate = true;
		switch (e.Key)
		{
			case Keys.A:
				panHorizontal(-stepping);
				break;
			case Keys.D:
				panHorizontal(stepping);
				break;
			case Keys.W:
				panVertical(stepping);
				break;
			case Keys.S:
				panVertical(-stepping);
				break;
			case Keys.N:
				zoomOut(-1);
				break;
			case Keys.M:
				zoomIn(-1);
				break;
			case Keys.X:
				zoomExtents(-1);
				doUpdate = false; // update performed in extents
				break;
			case Keys.Z:
				zoomExtents(ovpSettings.selectedIndex);
				doUpdate = false; // update performed in extents
				break;
		}

		switch (doUpdate)
		{
			case true when Platform.Instance.IsGtk:
				updateHostFunc?.Invoke();
				break;
			case true:
				updateViewport();
				break;
		}
	}

	private void zoomHandler(object sender, MouseEventArgs e)
	{
		if (!ovpSettings.isLocked())
		{
			float wheelZoom = e.Delta.Height; // SystemInformation.MouseWheelScrollLines;
			switch (wheelZoom)
			{
				case > 0:
					zoomIn(wheelZoom);
					break;
				case < 0:
					zoomOut(-wheelZoom);
					break;
			}

			updateViewport();
		}
		e.Handled = true;
	}

	private void getExtents(int index)
	{
		float minX = 0;
		float maxX = 0;
		float minY = 0, maxY = 0;

		bool set = false;

		switch (ovpSettings.polyList.Count)
		{
			case 0 when ovpSettings.lineList.Count == 0:
				ovpSettings.minX = 0;
				ovpSettings.maxX = 0;
				ovpSettings.minY = 0;
				ovpSettings.maxY = 0;
				return;
		}

		if (ovpSettings.polyList.Count != 0)
		{
			for (int poly = 0; poly < ovpSettings.polyList.Count; poly++)
			{
				if (index != -1 && (ovpSettings.polySourceIndex[poly] != index || !ovpSettings.polyMask[poly]))
				{
					continue;
				}

				switch (set)
				{
					case false:
						minX = ovpSettings.polyList[poly].poly[0].X;
						maxX = ovpSettings.polyList[poly].poly[0].X;
						minY = ovpSettings.polyList[poly].poly[0].Y;
						maxY = ovpSettings.polyList[poly].poly[0].Y;
						set = true;
						break;
				}
				float tMinX = ovpSettings.polyList[poly].poly.Min(p => p.X);
				float tMaxX = ovpSettings.polyList[poly].poly.Max(p => p.X);
				float tMinY = ovpSettings.polyList[poly].poly.Min(p => p.Y);
				float tMaxY = ovpSettings.polyList[poly].poly.Max(p => p.Y);
				minX = Math.Min(minX, tMinX);
				maxX = Math.Max(maxX, tMaxX);
				minY = Math.Min(minY, tMinY);
				maxY = Math.Max(maxY, tMaxY);
			}
		}

		if (ovpSettings.lineList.Count != 0)
		{
			for (int line = 0; line < ovpSettings.lineList.Count; line++)
			{
				if (index != -1 && (ovpSettings.lineSourceIndex[line] != index || !ovpSettings.lineMask[line]))
				{
					continue;
				}

				switch (set)
				{
					case false:
						minX = ovpSettings.lineList[line].poly[0].X;
						maxX = ovpSettings.lineList[line].poly[0].X;
						minY = ovpSettings.lineList[line].poly[0].Y;
						maxY = ovpSettings.lineList[line].poly[0].Y;
						set = true;
						break;
				}
				float tMinX = ovpSettings.lineList[line].poly.Min(p => p.X);
				float tMaxX = ovpSettings.lineList[line].poly.Max(p => p.X);
				float tMinY = ovpSettings.lineList[line].poly.Min(p => p.Y);
				float tMaxY = ovpSettings.lineList[line].poly.Max(p => p.Y);
				minX = Math.Min(minX, tMinX);
				maxX = Math.Max(maxX, tMaxX);
				minY = Math.Min(minY, tMinY);
				maxY = Math.Max(maxY, tMaxY);
			}
		}

		ovpSettings.minX = minX;
		ovpSettings.maxX = maxX;
		ovpSettings.minY = minY;
		ovpSettings.maxY = maxY;

		ovpSettings.changed = true;
	}

	private void drawPolygons()
	{
		int polyListCount = ovpSettings.polyList.Count;
		int bgPolyListCount = ovpSettings.bgPolyList.Count;
		int tessPolyListCount = ovpSettings.tessPolyList.Count;

		List<VertexPositionColor> polyList = new();

		List<VertexPositionColor> pointsList = new();

		List<VertexPositionColor> tessPolyList = new();

		try
		{

			// Carve our Z-space up to stack polygons
			int numPolys = 1;

			numPolys = polyListCount + bgPolyListCount;
			// Create our first and count arrays for the vertex indices, to enable polygon separation when rendering.
			polyFirst = new uint[numPolys];
			polyVertexCount = new uint[numPolys];

			tessFirst = new uint[tessPolyListCount];
			tessVertexCount = new uint[tessPolyListCount];

			List<uint> tFirst = new();

			uint tCounter = 0;

			if (ovpSettings.drawFilled())
			{
				numPolys += tessPolyListCount;
			}

			float polyZStep = 1.0f / Math.Max(1, numPolys + 1); // avoid a div by zero risk; pad the poly number also to reduce risk of adding a poly beyond the clipping range

			int counter = 0; // vertex count that will be used to define 'first' index for each polygon.
			int previouscounter = 0; // will be used to derive the number of vertices in each polygon.

			float polyZ = 0;

			if (ovpSettings.drawFilled())
			{
				for (int poly = 0; poly < tessPolyListCount; poly++)
				{
					tessFirst[poly] = (uint)(poly * 3);
					float alpha = ovpSettings.tessPolyList[poly].alpha;
					polyZ += polyZStep;
					for (int pt = 0; pt < 3; pt++)
					{
						tessPolyList.Add(new VertexPositionColor(new Vector3(ovpSettings.tessPolyList[poly].poly[pt].X, ovpSettings.tessPolyList[poly].poly[pt].Y, polyZ),
							new RgbaFloat(ovpSettings.tessPolyList[poly].color.R, ovpSettings.tessPolyList[poly].color.G, ovpSettings.tessPolyList[poly].color.B, alpha)));
					}
					tessVertexCount[poly] = 3;
				}
			}

			// Pondering options here - this would make a nice border construct around the filled geometry, amongst other things.
			for (int poly = 0; poly < polyListCount; poly++)
			{
				float alpha = ovpSettings.polyList[poly].alpha;
				if (ovpSettings.drawFilled())
				{
					alpha = 1.0f;
				}
				polyZ += polyZStep;
				polyFirst[poly] = (uint)counter;
				previouscounter = counter;
				int polyLength = ovpSettings.polyList[poly].poly.Length - 1;
				for (int pt = 0; pt < polyLength; pt++)
				{
					polyList.Add(new VertexPositionColor(new Vector3(ovpSettings.polyList[poly].poly[pt].X, ovpSettings.polyList[poly].poly[pt].Y, polyZ),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G, ovpSettings.polyList[poly].color.B, alpha)));
					counter++;
					polyList.Add(new VertexPositionColor(new Vector3(ovpSettings.polyList[poly].poly[pt + 1].X, ovpSettings.polyList[poly].poly[pt + 1].Y, polyZ),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G, ovpSettings.polyList[poly].color.B, alpha)));
					counter++;

					if (!ovpSettings.drawPoints())
					{
						continue;
					}

					tFirst.Add(tCounter);
					pointsList.Add(new VertexPositionColor(new Vector3(ovpSettings.polyList[poly].poly[pt].X - pointWidth / 2.0f, ovpSettings.polyList[poly].poly[pt].Y - pointWidth / 2.0f, 1.0f), new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G, ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
					pointsList.Add(new VertexPositionColor(new Vector3(ovpSettings.polyList[poly].poly[pt].X - pointWidth / 2.0f, ovpSettings.polyList[poly].poly[pt].Y + pointWidth / 2.0f, 1.0f), new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G, ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
					pointsList.Add(new VertexPositionColor(new Vector3(ovpSettings.polyList[poly].poly[pt].X + pointWidth / 2.0f, ovpSettings.polyList[poly].poly[pt].Y - pointWidth / 2.0f, 1.0f), new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G, ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;

					tFirst.Add(tCounter);
					pointsList.Add(new VertexPositionColor(new Vector3(ovpSettings.polyList[poly].poly[pt].X + pointWidth / 2.0f, ovpSettings.polyList[poly].poly[pt].Y - pointWidth / 2.0f, 1.0f), new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G, ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
					pointsList.Add(new VertexPositionColor(new Vector3(ovpSettings.polyList[poly].poly[pt].X - pointWidth / 2.0f, ovpSettings.polyList[poly].poly[pt].Y + pointWidth / 2.0f, 1.0f), new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G, ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
					pointsList.Add(new VertexPositionColor(new Vector3(ovpSettings.polyList[poly].poly[pt].X + pointWidth / 2.0f, ovpSettings.polyList[poly].poly[pt].Y + pointWidth / 2.0f, 1.0f), new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G, ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
				}
				polyVertexCount[poly] = (uint)(counter - previouscounter); // set our vertex count for the polygon.
			}

			polyZ = 0;
			for (int poly = 0; poly < bgPolyListCount; poly++)
			{
				float alpha = ovpSettings.bgPolyList[poly].alpha;
				polyZ += polyZStep;
				polyFirst[poly + polyListCount] = (uint)counter;
				previouscounter = counter;

				int bgPolyLength = ovpSettings.bgPolyList[poly].poly.Length - 1;
				for (int pt = 0; pt < bgPolyLength; pt++)
				{
					polyList.Add(new VertexPositionColor(new Vector3(ovpSettings.bgPolyList[poly].poly[pt].X, ovpSettings.bgPolyList[poly].poly[pt].Y, polyZ),
						new RgbaFloat(ovpSettings.bgPolyList[poly].color.R, ovpSettings.bgPolyList[poly].color.G, ovpSettings.bgPolyList[poly].color.B, alpha)));
					counter++;
					polyList.Add(new VertexPositionColor(new Vector3(ovpSettings.bgPolyList[poly].poly[pt + 1].X, ovpSettings.bgPolyList[poly].poly[pt + 1].Y, polyZ),
						new RgbaFloat(ovpSettings.bgPolyList[poly].color.R, ovpSettings.bgPolyList[poly].color.G, ovpSettings.bgPolyList[poly].color.B, alpha)));
					counter++;
				}
				polyVertexCount[poly + polyListCount] = (uint)(counter - previouscounter); // set our vertex count for the polygon.
			}

			pointsFirst = tFirst.ToArray();
		}
		catch (Exception)
		{
			// Can ignore - not critical.
		}

		if (polyListCount > 0 || bgPolyListCount > 0)
		{
			updateBuffer(ref PolysVertexBuffer, polyList.ToArray(), VertexPositionColor.SizeInBytes, BufferUsage.VertexBuffer);
		}
		if (ovpSettings.drawPoints() && polyListCount > 0)
		{
			updateBuffer(ref PointsVertexBuffer, pointsList.ToArray(), VertexPositionColor.SizeInBytes, BufferUsage.VertexBuffer);
		}
		if (ovpSettings.drawFilled() && tessPolyListCount > 0)
		{
			updateBuffer(ref TessVertexBuffer, tessPolyList.ToArray(), VertexPositionColor.SizeInBytes, BufferUsage.VertexBuffer);
		}
	}

	private void drawLines()
	{
		int tmp = ovpSettings.lineList.Count;

		switch (tmp)
		{
			// Create our first and count arrays for the vertex indices, to enable polygon separation when rendering.
			case > 0:
			{
				List<VertexPositionColor> lineList = new();

				// Carve our Z-space up to stack polygons
				float polyZStep = 1.0f / ovpSettings.lineList.Count;

				lineFirst = new uint[tmp];
				lineVertexCount = new uint[tmp];

				for (int poly = 0; poly < tmp; poly++)
				{
					float alpha = ovpSettings.lineList[poly].alpha;
					float polyZ = poly * polyZStep;
					lineFirst[poly] = (uint)lineList.Count;
					lineList.AddRange(ovpSettings.lineList[poly].poly.Select(t => new VertexPositionColor(new Vector3(t.X, t.Y, polyZ), new RgbaFloat(ovpSettings.lineList[poly].color.R, ovpSettings.lineList[poly].color.G, ovpSettings.lineList[poly].color.B, alpha))));
					lineVertexCount[poly] = (uint)ovpSettings.lineList[poly].poly.Length; // set our vertex count for the polygon.
				}

				updateBuffer(ref LinesVertexBuffer, lineList.ToArray(), VertexPositionColor.SizeInBytes, BufferUsage.VertexBuffer);
				break;
			}
			default:
				LinesVertexBuffer = null;
				break;
		}
	}

	private void drawGrid()
	{
		if (!ovpSettings.drawGrid())
		{
			return;
		}

		float spacing = ovpSettings.gridSpacing();
		if (ovpSettings.isGridDynamic())
		{
			while (WorldToScreen(new SizeF(spacing, 0.0f)).Width > 12.0f)
			{
				spacing /= 10.0f;
			}

			while (WorldToScreen(new SizeF(spacing, 0.0f)).Width < 4.0f)
			{
				spacing *= 10.0f;
			}
		}

		float zoom = ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom();
		float x = ovpSettings.getCameraX();
		float y = ovpSettings.getCameraY();

		List<VertexPositionColor> grid = new();

		if (WorldToScreen(new SizeF(spacing, 0.0f)).Width >= 4.0f)
		{
			int k = 0;
			for (float i = 0; i > -(Surface.RenderWidth * zoom) + x; i -= spacing)
			{
				float r = 0.0f;
				float g = 0.0f;
				float b = 0.0f;
				switch (k)
				{
					case <= 9:
						r = ovpSettings.minorGridColor.R;
						g = ovpSettings.minorGridColor.G;
						b = ovpSettings.minorGridColor.B;
						break;
					case 10:
						r = ovpSettings.majorGridColor.R;
						g = ovpSettings.majorGridColor.G;
						b = ovpSettings.majorGridColor.B;
						k = 0;
						break;
				}

				k++;
				grid.Add(new VertexPositionColor(new Vector3(i, y + zoom * Surface.RenderHeight, gridZ), new RgbaFloat(r, g, b, 1.0f)));
				grid.Add(new VertexPositionColor(new Vector3(i, y + zoom * -Surface.RenderHeight, gridZ), new RgbaFloat(r, g, b, 1.0f)));
			}
			k = 0;
			for (float i = 0; i < Surface.RenderWidth * zoom + x; i += spacing)
			{
				float r = 0.0f;
				float g = 0.0f;
				float b = 0.0f;
				switch (k)
				{
					case <= 9:
						r = ovpSettings.minorGridColor.R;
						g = ovpSettings.minorGridColor.G;
						b = ovpSettings.minorGridColor.B;
						break;
					case 10:
						r = ovpSettings.majorGridColor.R;
						g = ovpSettings.majorGridColor.G;
						b = ovpSettings.majorGridColor.B;
						k = 0;
						break;
				}

				k++;
				grid.Add(new VertexPositionColor(new Vector3(i, y + zoom * Surface.RenderHeight, gridZ), new RgbaFloat(r, g, b, 1.0f)));
				grid.Add(new VertexPositionColor(new Vector3(i, y + zoom * -Surface.RenderHeight, gridZ), new RgbaFloat(r, g, b, 1.0f)));
			}
			k = 0;
			for (float i = 0; i > -(Surface.RenderHeight * zoom) + y; i -= spacing)
			{
				float r = 0.0f;
				float g = 0.0f;
				float b = 0.0f;
				switch (k)
				{
					case <= 9:
						r = ovpSettings.minorGridColor.R;
						g = ovpSettings.minorGridColor.G;
						b = ovpSettings.minorGridColor.B;
						break;
					case 10:
						r = ovpSettings.majorGridColor.R;
						g = ovpSettings.majorGridColor.G;
						b = ovpSettings.majorGridColor.B;
						k = 0;
						break;
				}

				k++;
				grid.Add(new VertexPositionColor(new Vector3(x + zoom * Surface.RenderWidth, i, gridZ), new RgbaFloat(r, g, b, 1.0f)));
				grid.Add(new VertexPositionColor(new Vector3(x + zoom * -Surface.RenderWidth, i, gridZ), new RgbaFloat(r, g, b, 1.0f)));
			}
			k = 0;
			for (float i = 0; i < Surface.RenderHeight * zoom + y; i += spacing)
			{
				float r = 0.0f;
				float g = 0.0f;
				float b = 0.0f;
				switch (k)
				{
					case <= 9:
						r = ovpSettings.minorGridColor.R;
						g = ovpSettings.minorGridColor.G;
						b = ovpSettings.minorGridColor.B;
						break;
					case 10:
						r = ovpSettings.majorGridColor.R;
						g = ovpSettings.majorGridColor.G;
						b = ovpSettings.majorGridColor.B;
						k = 0;
						break;
				}

				k++;
				grid.Add(new VertexPositionColor(new Vector3(x + zoom * Surface.RenderWidth, i, gridZ), new RgbaFloat(r, g, b, 1.0f)));
				grid.Add(new VertexPositionColor(new Vector3(x + zoom * -Surface.RenderWidth, i, gridZ), new RgbaFloat(r, g, b, 1.0f)));
			}
		}

		uint gridCount = (uint)grid.Count;

		switch (gridCount)
		{
			case > 0:
			{
				gridIndices = new uint[gridCount];
				for (uint i = 0; i < gridIndices.Length; i++)
				{
					gridIndices[i] = i;
				}

				updateBuffer(ref GridVertexBuffer, grid.ToArray(), VertexPositionColor.SizeInBytes, BufferUsage.VertexBuffer);
				updateBuffer(ref GridIndexBuffer, gridIndices, sizeof(uint), BufferUsage.IndexBuffer);
				break;
			}
			default:
				GridVertexBuffer = null;
				GridIndexBuffer = null;
				break;
		}
	}

	private void drawAxes()
	{
		if (!ovpSettings.drawAxes())
		{
			return;
		}

		float zoom = ovpSettings.getBaseZoom() * ovpSettings.getZoomFactor();
		VertexPositionColor[] axesArray = new VertexPositionColor[4];
		axesArray[0] = new VertexPositionColor(new Vector3(0.0f, ovpSettings.getCameraY() + Surface.RenderHeight * zoom, axisZ), new RgbaFloat(ovpSettings.axisColor.R, ovpSettings.axisColor.G, ovpSettings.axisColor.B, 1.0f));
		axesArray[1] = new VertexPositionColor(new Vector3(0.0f, ovpSettings.getCameraY() - Surface.RenderHeight * zoom, axisZ), new RgbaFloat(ovpSettings.axisColor.R, ovpSettings.axisColor.G, ovpSettings.axisColor.B, 1.0f));
		axesArray[2] = new VertexPositionColor(new Vector3(ovpSettings.getCameraX() + Surface.RenderWidth * zoom, 0.0f, axisZ), new RgbaFloat(ovpSettings.axisColor.R, ovpSettings.axisColor.G, ovpSettings.axisColor.B, 1.0f));
		axesArray[3] = new VertexPositionColor(new Vector3(ovpSettings.getCameraX() - Surface.RenderWidth * zoom, 0.0f, axisZ), new RgbaFloat(ovpSettings.axisColor.R, ovpSettings.axisColor.G, ovpSettings.axisColor.B, 1.0f));

		axesIndices = new uint[4] { 0, 1, 2, 3 };

		updateBuffer(ref AxesVertexBuffer, axesArray, VertexPositionColor.SizeInBytes, BufferUsage.VertexBuffer);
		updateBuffer(ref AxesIndexBuffer, axesIndices, sizeof(uint), BufferUsage.IndexBuffer);
	}

	/// <summary>
	/// Fills the given buffer with the contents of 'data', creating or
	/// resizing it as necessary.
	/// </summary>
	/// <param name="buffer">The Veldrid.DeviceBuffer to fill.</param>
	/// <param name="data">The array of elements to put in the buffer.</param>
	/// <param name="elementSize">The size in bytes of each element.</param>
	/// <param name="usage">The Veldrid.BufferUsage type of 'buffer'.</param>
	public void updateBuffer<T>(ref DeviceBuffer buffer, T[] data, uint elementSize, BufferUsage usage)
		where T : unmanaged
	{
		switch (data.Length)
		{
			case > 0:
			{
				buffer?.Dispose();

				ResourceFactory factory = Surface.GraphicsDevice.ResourceFactory;

				buffer = factory.CreateBuffer(new BufferDescription(elementSize * (uint)data.Length, usage));

				Surface.GraphicsDevice.UpdateBuffer(buffer, 0, data);
				break;
			}
		}
	}

	public void updateViewport()
	{
		switch (ovpSettings.changed)
		{
			case false:
				return;
			default:
				pUpdateViewport();
				break;
		}
	}

	private void pUpdateViewport()
	{
		switch (Surface.GraphicsDevice)
		{
			case null:
				return;
		}

		if (Surface.Visible && (Surface.Width > 0) && (Surface.Height > 0))
		{
			drawAxes();
			drawGrid();
			drawLines();
			drawPolygons();
			updateHostFunc?.Invoke();
			Draw();
		}

	}

	public void Draw()
	{
		switch (Ready)
		{
			case false:
				return;
		}

		CommandList.Begin();

		ModelMatrix *= Matrix4x4.CreateFromAxisAngle(
			new Vector3(0, 0, 1), 0);
		CommandList.UpdateBuffer(ModelBuffer, 0, ModelMatrix);

		float zoom = ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom();

		float left = ovpSettings.getCameraX() - (float)Surface.RenderWidth / 2 * zoom;
		float right = ovpSettings.getCameraX() + (float)Surface.RenderWidth / 2 * zoom;
		float bottom = ovpSettings.getCameraY() + (float)Surface.RenderHeight / 2 * zoom;
		float top = ovpSettings.getCameraY() - (float)Surface.RenderHeight / 2 * zoom;

		ViewMatrix = Matrix4x4.CreateOrthographicOffCenter(left, right, bottom, top, 0.0f, 1.0f);
		CommandList.UpdateBuffer(ViewBuffer, 0, ViewMatrix);

		CommandList.SetFramebuffer(Surface.Swapchain.Framebuffer);

		// These commands differ from the stock Veldrid "Getting Started"
		// tutorial in two ways. First, the viewport is cleared to pink
		// instead of black so as to more easily distinguish between errors
		// in creating a graphics context and errors drawing vertices within
		// said context. Second, this project creates its swapchain with a
		// depth buffer, and that buffer needs to be reset at the start of
		// each frame.

		RgbaFloat bgColor = new(ovpSettings.backColor.R, ovpSettings.backColor.G, ovpSettings.backColor.B, 1.0f);

		CommandList.ClearColorTarget(0, bgColor);
		CommandList.ClearDepthStencil(1.0f);

		if (GridVertexBuffer != null)
		{
			lock (GridVertexBuffer)
			{
				try
				{
					CommandList.SetVertexBuffer(0, GridVertexBuffer);
					CommandList.SetIndexBuffer(GridIndexBuffer, IndexFormat.UInt32);
					CommandList.SetPipeline(LinePipeline);
					CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
					CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

					CommandList.DrawIndexed(
						indexCount: (uint)gridIndices.Length,
						instanceCount: 1,
						indexStart: 0,
						vertexOffset: 0,
						instanceStart: 0);
				}
				catch (Exception)
				{

				}
			}
		}

		if (AxesVertexBuffer != null)
		{
			lock (AxesVertexBuffer)
			{
				try
				{
					CommandList.SetVertexBuffer(0, AxesVertexBuffer);
					CommandList.SetIndexBuffer(AxesIndexBuffer, IndexFormat.UInt32);
					CommandList.SetPipeline(LinePipeline);
					CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
					CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

					CommandList.DrawIndexed(
						indexCount: (uint)axesIndices.Length,
						instanceCount: 1,
						indexStart: 0,
						vertexOffset: 0,
						instanceStart: 0);
				}
				catch (Exception)
				{

				}
			}
		}

		if (ovpSettings.drawFilled())
		{
			if (TessVertexBuffer != null)
			{
				lock (TessVertexBuffer)
				{
					try
					{
						CommandList.SetVertexBuffer(0, TessVertexBuffer);
						CommandList.SetPipeline(FilledPipeline);
						CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
						CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

						for (int l = 0; l < tessVertexCount.Length; l++)
						{
							CommandList.Draw(tessVertexCount[l], 1, tessFirst[l], 0);
						}
					}
					catch (Exception)
					{

					}
				}
			}
		}

		if (PolysVertexBuffer != null)
		{
			lock (PolysVertexBuffer)
			{
				try
				{
					CommandList.SetVertexBuffer(0, PolysVertexBuffer);
					CommandList.SetPipeline(LinesPipeline);
					CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
					CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

					for (int l = 0; l < polyVertexCount.Length; l++)
					{
						CommandList.Draw(polyVertexCount[l], 1, polyFirst[l], 0);
					}
				}
				catch (Exception)
				{

				}
			}
		}

		if (LinesVertexBuffer != null && ovpSettings.drawDrawn())
		{
			lock (LinesVertexBuffer)
			{
				try
				{
					CommandList.SetVertexBuffer(0, LinesVertexBuffer);
					CommandList.SetPipeline(LinesPipeline);
					CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
					CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

					for (int l = 0; l < lineVertexCount.Length; l++)
					{
						CommandList.Draw(lineVertexCount[l], 1, lineFirst[l], 0);
					}
				}
				catch (Exception)
				{

				}
			}
		}

		if (ovpSettings.drawPoints())
		{
			if (PointsVertexBuffer != null)
			{
				lock (PointsVertexBuffer)
				{
					try
					{
						CommandList.SetVertexBuffer(0, PointsVertexBuffer);
						CommandList.SetPipeline(FilledPipeline);
						CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
						CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

						foreach (uint t in pointsFirst)
						{
							CommandList.Draw(3, 1, t, 0);
						}
					}
					catch (Exception)
					{

					}
				}
			}
		}
		CommandList.End();

		try
		{
			lock (CommandList)
			{
				Surface.GraphicsDevice.SubmitCommands(CommandList);
			}
			Surface.GraphicsDevice.SwapBuffers(Surface.Swapchain);
		}
		catch (Exception)
		{

		}
	}

	public void SetUpVeldrid()
	{
		CreateResources();

		Ready = true;
	}

	private void CreateResources()
	{
		ResourceFactory factory = Surface.GraphicsDevice.ResourceFactory;

		ResourceLayout viewMatrixLayout = factory.CreateResourceLayout(
			new ResourceLayoutDescription(
				new ResourceLayoutElementDescription(
					"ViewMatrix",
					ResourceKind.UniformBuffer,
					ShaderStages.Vertex)));

		ViewBuffer = factory.CreateBuffer(
			new BufferDescription(64, BufferUsage.UniformBuffer));

		ViewMatrixSet = factory.CreateResourceSet(new ResourceSetDescription(
			viewMatrixLayout, ViewBuffer));

		ResourceLayout modelMatrixLayout = factory.CreateResourceLayout(
			new ResourceLayoutDescription(
				new ResourceLayoutElementDescription(
					"ModelMatrix",
					ResourceKind.UniformBuffer,
					ShaderStages.Vertex)));

		ModelBuffer = factory.CreateBuffer(
			new BufferDescription(64, BufferUsage.UniformBuffer));

		ModelMatrixSet = factory.CreateResourceSet(new ResourceSetDescription(
			modelMatrixLayout, ModelBuffer));

		// drawGrid();

		// Veldrid.SPIRV, when cross-compiling to HLSL, will always produce
		// TEXCOORD semantics; VertexElementSemantic.TextureCoordinate thus
		// becomes necessary to let D3D11 work alongside Vulkan and OpenGL.
		//
		//   https://github.com/mellinoe/veldrid/issues/121
		//
		VertexLayoutDescription vertexLayout = new(
			new VertexElementDescription("Position", VertexElementSemantic.TextureCoordinate, VertexElementFormat.Float3),
			new VertexElementDescription("Color", VertexElementSemantic.TextureCoordinate, VertexElementFormat.Float4));

		// Veldrid.SPIRV is an additional library that complements Veldrid
		// by simplifying the development of cross-backend shaders, and is
		// currently the recommended approach to doing so:
		//
		//   https://veldrid.dev/articles/portable-shaders.html
		//
		// If you decide against using it, you can try out Veldrid developer
		// mellinoe's other project, ShaderGen, or drive yourself crazy by
		// writing and maintaining custom shader code for each platform.
		byte[] vertexShaderSpirvBytes = LoadSpirvBytes(ShaderStages.Vertex);
		byte[] fragmentShaderSpirvBytes = LoadSpirvBytes(ShaderStages.Fragment);

		CrossCompileOptions options = new();
		switch (Surface.GraphicsDevice.BackendType)
		{
			// InvertVertexOutputY and FixClipSpaceZ address two major
			// differences between Veldrid's various graphics APIs, as
			// discussed here:
			//
			//   https://veldrid.dev/articles/backend-differences.html
			//
			// Note that the only reason those options are useful in this
			// example project is that the vertices being drawn are stored
			// the way Vulkan stores vertex data. The options will therefore
			// properly convert from the Vulkan style to whatever's used by
			// the destination backend. If you store vertices in a different
			// coordinate system, these may not do anything for you, and
			// you'll need to handle the difference in your shader code.
			case GraphicsBackend.Metal:
				options.InvertVertexOutputY = true;
				break;
			case GraphicsBackend.Direct3D11:
				options.InvertVertexOutputY = true;
				break;
			case GraphicsBackend.OpenGL:
				options.FixClipSpaceZ = true;
				options.InvertVertexOutputY = true;
				break;
		}

		ShaderDescription vertex = new(ShaderStages.Vertex, vertexShaderSpirvBytes, "main", true);
		ShaderDescription fragment = new(ShaderStages.Fragment, fragmentShaderSpirvBytes, "main", true);
		Shader[] shaders = factory.CreateFromSpirv(vertex, fragment, options);

		PointsPipeline = factory.CreateGraphicsPipeline(new GraphicsPipelineDescription
		{
			BlendState = BlendStateDescription.SingleOverrideBlend,
			DepthStencilState = new DepthStencilStateDescription(
				depthTestEnabled: false,
				depthWriteEnabled: false,
				comparisonKind: ComparisonKind.LessEqual),
			RasterizerState = new RasterizerStateDescription(
				cullMode: FaceCullMode.None,
				fillMode: PolygonFillMode.Solid,
				frontFace: FrontFace.Clockwise,
				depthClipEnabled: false,
				scissorTestEnabled: false),
			PrimitiveTopology = PrimitiveTopology.LineStrip,
			ResourceLayouts = new[] { viewMatrixLayout, modelMatrixLayout },
			ShaderSet = new ShaderSetDescription(
				vertexLayouts: new[] { vertexLayout },
				shaders: shaders),
			Outputs = Surface.Swapchain.Framebuffer.OutputDescription
		});

		LinePipeline = factory.CreateGraphicsPipeline(new GraphicsPipelineDescription
		{
			BlendState = BlendStateDescription.SingleOverrideBlend,
			DepthStencilState = new DepthStencilStateDescription(
				depthTestEnabled: true,
				depthWriteEnabled: true,
				comparisonKind: ComparisonKind.LessEqual),
			RasterizerState = new RasterizerStateDescription(
				cullMode: FaceCullMode.Back,
				fillMode: PolygonFillMode.Solid,
				frontFace: FrontFace.Clockwise,
				depthClipEnabled: true,
				scissorTestEnabled: false),
			PrimitiveTopology = PrimitiveTopology.LineList,
			ResourceLayouts = new[] { viewMatrixLayout, modelMatrixLayout },
			ShaderSet = new ShaderSetDescription(
				vertexLayouts: new[] { vertexLayout },
				shaders: shaders),
			Outputs = Surface.Swapchain.Framebuffer.OutputDescription
		});

		LinesPipeline = factory.CreateGraphicsPipeline(new GraphicsPipelineDescription
		{
			BlendState = BlendStateDescription.SingleAlphaBlend,
			DepthStencilState = new DepthStencilStateDescription(
				depthTestEnabled: false,
				depthWriteEnabled: false,
				comparisonKind: ComparisonKind.LessEqual),
			RasterizerState = new RasterizerStateDescription(
				cullMode: FaceCullMode.Back,
				fillMode: PolygonFillMode.Solid,
				frontFace: FrontFace.Clockwise,
				depthClipEnabled: false,
				scissorTestEnabled: false),
			PrimitiveTopology = PrimitiveTopology.LineStrip,
			ResourceLayouts = new[] { viewMatrixLayout, modelMatrixLayout },
			ShaderSet = new ShaderSetDescription(
				vertexLayouts: new[] { vertexLayout },
				shaders: shaders),
			Outputs = Surface.Swapchain.Framebuffer.OutputDescription
		});

		FilledPipeline = factory.CreateGraphicsPipeline(new GraphicsPipelineDescription
		{
			BlendState = BlendStateDescription.SingleAlphaBlend,
			DepthStencilState = new DepthStencilStateDescription(
				depthTestEnabled: false,
				depthWriteEnabled: false,
				comparisonKind: ComparisonKind.LessEqual),
			RasterizerState = new RasterizerStateDescription(
				cullMode: FaceCullMode.None,
				fillMode: PolygonFillMode.Solid,
				frontFace: FrontFace.CounterClockwise,
				depthClipEnabled: false,
				scissorTestEnabled: false),
			PrimitiveTopology = PrimitiveTopology.TriangleStrip,
			ResourceLayouts = new[] { viewMatrixLayout, modelMatrixLayout },
			ShaderSet = new ShaderSetDescription(
				vertexLayouts: new[] { vertexLayout },
				shaders: shaders),
			Outputs = Surface.Swapchain.Framebuffer.OutputDescription
		});

		CommandList = factory.CreateCommandList();
	}

	private byte[] LoadSpirvBytes(ShaderStages stage)
	{
		byte[] bytes;

		string name = $"VertexColor-{stage.ToString().ToLower()}.450.glsl";
		string full = Path.Combine(ExecutableDirectory, ShaderSubdirectory, name);

		// Precompiled SPIR-V bytecode can speed up program start by saving
		// the need to load text files and compile them before converting
		// the result to the final backend shader format. If they're not
		// available, though, the plain .glsl files will do just fine. Look
		// up glslangValidator to learn how to compile SPIR-V binary files.
		try
		{
			bytes = File.ReadAllBytes($"{full}.spv");
		}
		catch (FileNotFoundException)
		{
			bytes = File.ReadAllBytes(full);
		}

		return bytes;
	}
}