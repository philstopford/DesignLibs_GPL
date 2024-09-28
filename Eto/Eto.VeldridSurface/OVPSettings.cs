using Eto.Drawing;
using Eto.Forms;
using System;
using System.Collections.Generic;
using System.Linq;
using Clipper2Lib;
using LibTessDotNet.Double;

namespace VeldridEto;

public static class errorReporter
{
	public static void showMessage_OK(string stringToDisplay, string caption)
	{
		Application.Instance.Invoke(() =>
		{
			MessageBox.Show(stringToDisplay, caption, MessageBoxButtons.OK);
		});
	}
}

public class OVPSettings
{
	// This is true if something changed in the settings (set internally for query). The viewport itself should set this false when the changes are handled.
	public bool changed { get; set; }

	public int selectedIndex { get; set; }

	public float minX { get; set; }
	public float maxX { get; set; }
	public float minY { get; set; }
	public float maxY { get; set; }
	public Color minorGridColor { get; set; }
	public Color majorGridColor { get; set; }
	public Color axisColor { get; set; }
	public Color backColor { get; set; }
	public Color selectionColor { get; set; }
	public Color inverSelectionColor { get; set; }

	public List<ovp_Poly>? polyList { get; set; }
	public List<int>? polyListPtCount { get; set; }
	public List<int>? polySourceIndex { get; set; } // will eventually track source of polygon, allowing for layer generating, etc. in output.

	public List<bool>? polyMask { get; set; } // Masking boolean for overriding handling of polygons as-desired.
	public List<ovp_Poly>? bgPolyList { get; set; }
	public List<int>? bgPolyListPtCount { get; set; }
	public List<int>? bgPolySourceIndex { get; set; } // will eventually track source of polygon, allowing for layer generating, etc. in output.
	public List<ovp_Poly>? lineList { get; set; } // purely for lines.
	public List<int>? lineListPtCount { get; set; }
	public List<int>? lineSourceIndex { get; set; } // will eventually track source of polygon, allowing for layer generating, etc. in output.
	public List<bool>? lineMask { get; set; } // Masking boolean for overriding handling of polygons as-desired.

	public List<ovp_Poly>? tessPolyList { get; set; } // triangles, but also need to track color. This is decoupled to allow boundary extraction without triangles getting in the way.
	public List<bool>? drawnPoly { get; set; } // tracks whether the polygon corresponds to an enabled configuration or not.

	private bool enableFilledPolys;
	private bool showPoints;
	private float base_zoom;
	private float zoomFactor;
	private int zoomStep;
	private bool allowZoomAndPan;
	private bool dynamicGrid;
	private bool panning;
	private bool selecting;
	private bool showGrid;
	private bool showAxes;
	private bool showDrawn;
	private int grid_spacing;
	private PointF cameraPosition;
	private PointF default_cameraPosition;
	private bool antiAlias;
	private bool lockedViewport;

	public bool aA()
	{
		return antiAlias;
	}

	public void aA(bool val)
	{
		if (antiAlias == val)
		{
			return;
		}
		antiAlias = val;
		changed = true;
	}

	public bool drawGrid()
	{
		return showGrid;
	}

	public void drawGrid(bool val)
	{
		if (showGrid == val)
		{
			return;
		}
		showGrid = val;
		changed = true;
	}

	public bool drawPoints()
	{
		return showPoints;
	}

	public void drawPoints (bool val)
	{
		if (showPoints == val)
		{
			return;
		}
		showPoints = val;
		changed = true;
	}

	public bool drawDrawn()
	{
		return showDrawn;
	}

	public void drawDrawn(bool val)
	{
		if (showDrawn == val)
		{
			return;
		}
		showDrawn = val;
		changed = true;
	}

	public bool drawFilled()
	{
		return enableFilledPolys;
	}

	public void drawFilled(bool val)
	{
		if (enableFilledPolys == val)
		{
			return;
		}
		enableFilledPolys = val;
		changed = true;
	}

	public bool drawAxes()
	{
		return showAxes;
	}

	public void drawAxes(bool val)
	{
		if (showAxes == val)
		{
			return;
		}
		showAxes = val;
		changed = true;
	}

	public bool isGridDynamic()
	{
		return dynamicGrid;
	}

	public void isGridDynamic(bool val)
	{
		if (dynamicGrid == val)
		{
			return;
		}
		dynamicGrid = val;
		changed = true;
	}

	public float gridSpacing()
	{
		return grid_spacing;
	}

	public void resetCamera()
	{
		switch (lockedViewport)
		{
			case true:
				return;
			default:
				setCameraPos(default_cameraPosition.X, default_cameraPosition.Y);
				setZoomFactor(1.0f);
				break;
		}
	}

	public void setCameraPos(float x, float y)
	{
		setCameraX(x);
		setCameraY(y);
	}

	public PointF getCameraPos()
	{
		return cameraPosition;
	}

	public void setCameraX(float x)
	{
		switch (lockedViewport)
		{
			case true:
				return;
			default:
				cameraPosition.X = x;
				changed = true;
				break;
		}
	}

	public void setCameraY(float y)
	{
		switch (lockedViewport)
		{
			case true:
				return;
			default:
				cameraPosition.Y = y;
				changed = true;
				break;
		}
	}

	public float getCameraX()
	{
		return cameraPosition.X;
	}

	public float getCameraY()
	{
		return cameraPosition.Y;
	}

	public float getBaseZoom()
	{
		return base_zoom;
	}

	public void setBaseZoom(float val)
	{
		base_zoom = val;
	}
	public void setZoomFactor(float val)
	{
		switch (lockedViewport)
		{
			case true:
				return;
		}
		if (val < 0.0001)
		{
			val = 0.0001f; // avoid any chance of getting to zero.
		}
		zoomFactor = val;
		changed = true;
	}

	public float getZoomFactor()
	{
		return zoomFactor;
	}

	public int getZoomStep()
	{
		return zoomStep;
	}

	public void setZoomStep(int val)
	{
		zoomStep = val;
	}

	public void lockVP(bool val)
	{
		if (val == lockedViewport)
		{
			return;
		}
		lockedViewport = val;
		changed = true;
	}

	public bool isLocked()
	{
		return lockedViewport;
	}

	public float zoom()
	{
		return base_zoom * zoomFactor;
	}

	public void updateColors(Color newColor)
	{
		pUpdateColors(newColor);
	}

	private void pUpdateColors(Color newColor)
	{
		foreach (ovp_Poly t in polyList!)
		{
			t.color = newColor;
		}
		foreach (ovp_Poly t in tessPolyList!)
		{
			t.color = newColor;
		}
		changed = true;
	}

	public void reset(bool clearBG = true)
	{
		pReset(clearBG);
	}

	private void pReset(bool clearBG)
	{
		minX = 0;
		maxX = 0;
		minY = 0;
		maxY = 0;
		clear(clearBG);
		drawnPoly!.Clear();
		changed = true;
	}

	public void clear(bool clearBG = true)
	{
		pClear(clearBG);
	}

	private void pClear(bool clearBG)
	{
		polyList!.Clear();
		polySourceIndex!.Clear();
		polyMask!.Clear();
		polyListPtCount!.Clear();
		switch (clearBG)
		{
			case true:
				bgPolyList!.Clear();
				bgPolySourceIndex!.Clear();
				bgPolyListPtCount!.Clear();
				break;
		}
		lineList!.Clear();
		lineSourceIndex!.Clear();
		lineMask!.Clear();
		lineListPtCount!.Clear();
		tessPolyList!.Clear();
		changed = true;
	}

	public void addLine(PointF[] line, Color lineColor, float alpha, int layerIndex, bool mask = true)
	{
		pAddLine(line, lineColor, alpha, layerIndex, mask);
	}

	private void pAddLine(PointF[] line, Color lineColor, float alpha, int layerIndex, bool mask)
	{
		lineList!.Add(new ovp_Poly(line, lineColor, alpha));
		lineSourceIndex!.Add(layerIndex);
		lineListPtCount!.Add((line.Length - 1) * 2);
		lineMask!.Add(true);
		changed = true;
	}

	public void addPolygon(PointF[] poly, Color polyColor, float alpha, bool drawn, int layerIndex, bool mask = true)
	{
		switch (drawn)
		{
			case true:
				// Drawn polygons are to be treated as lines : they don't get filled.
				addLine(poly, polyColor, alpha, layerIndex, mask);
				break;
			default:
				pAddPolygon(poly, polyColor, alpha, drawn, layerIndex, mask);
				break;
		}
	}

	private void pAddPolygon(PointF[] poly, Color polyColor, float alpha, bool drawn, int layerIndex, bool mask)
	{
		switch (drawn)
		{
			// avoid tessellation unless really necessary.
			case false when enableFilledPolys:
				try
				{
					tessPoly(poly, polyColor, alpha);
				}
				catch (Exception)
				{

				}

				break;
		}

		polyList!.Add(new ovp_Poly(poly, polyColor, alpha));
		polySourceIndex!.Add(layerIndex);
		polyListPtCount!.Add(poly.Length);
		polyMask!.Add(mask);
		drawnPoly!.Add(drawn);
		changed = true;
	}

	public void addBGPolygon(PointF[] poly, Color polyColor, float alpha, int layerIndex)
	{
		pAddBGPolygon(poly, polyColor, alpha, layerIndex);
	}

	private void pAddBGPolygon(PointF[] poly, Color polyColor, float alpha, int layerIndex)
	{
		bgPolyList!.Add(new ovp_Poly(poly, polyColor, alpha));
		bgPolyListPtCount!.Add(poly.Length);
		bgPolySourceIndex!.Add(layerIndex);
		drawnPoly!.Add(false);
		changed = true;
	}

	public OVPSettings(float defX = 0.0f, float defY = 0.0f)
	{
		init(defX, defY);
	}

	private void init(float defX, float defY)
	{
		base_zoom = 1.0f;
		minorGridColor = new Color(1.0f, 0.8f, 0.3f);
		majorGridColor = new Color(0.2f, 0.8f, 0.3f);
		axisColor = new Color(0.1f, 0.1f, 0.1f);
		backColor = new Color(1.0f, 1.0f, 1.0f);
		selectionColor = SystemColors.Highlight;
		inverSelectionColor = SystemColors.Highlight;
		allowZoomAndPan = true;
		enableFilledPolys = false;
		showPoints = true;
		dynamicGrid = true;
		panning = false;
		selecting = false;
		showGrid = true;
		showAxes = true;
		grid_spacing = 10;
		antiAlias = true;
		zoomStep = 1;
		fullReset(defX, defY);
	}

	public void fullReset(float defX = 0.0f, float defY = 0.0f)
	{
		pFullReset(defX, defY);
	}

	private void pFullReset(float defX = 0.0f, float defY = 0.0f)
	{
		default_cameraPosition = new PointF(defX, defY);
		bgPolyList = new List<ovp_Poly>();
		bgPolySourceIndex = new List<int>();
		bgPolyListPtCount = new List<int>();
		polyList = new List<ovp_Poly>();
		polySourceIndex = new List<int>();
		polyMask = new List<bool>();
		polyListPtCount = new List<int>();
		lineList = new List<ovp_Poly>();
		lineSourceIndex = new List<int>();
		lineMask = new List<bool>();
		lineListPtCount = new List<int>();
		tessPolyList = new List<ovp_Poly>();
		drawnPoly = new List<bool>();
		zoomFactor = 1.0f;
		cameraPosition = new PointF(default_cameraPosition.X, default_cameraPosition.Y);
		changed = true;
	}

	private static PointF[] clockwiseOrder(PointF[] iPoints)
	{
		// Based on stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
		// Shoelace formula.

		double delta = 0;

		for (int pt = 0; pt < iPoints.Length; pt++)
		{
			double deltaX = 0;
			double deltaY = 0;
			if (pt == iPoints.Length - 1)
			{
				deltaX = iPoints[0].X - iPoints[pt].X;
				deltaY = iPoints[0].Y + iPoints[pt].Y;
			}
			else
			{
				deltaX = iPoints[pt + 1].X - iPoints[pt].X;
				deltaY = iPoints[pt + 1].Y + iPoints[pt].Y;
			}

			delta += deltaX * deltaY;
		}

		switch (delta)
		{
			case > 0:
				// clockwise
				break;
			default:
				// counter-clockwise.
				Array.Reverse(iPoints);
				break;
		}

		return iPoints;
	}

	public static PointF[] convertToClosedPointF(PathD poly)
	{
		bool forceClosed = false;
		int polyCount = poly.Count;
		if (!(Math.Abs(poly[0].x - poly[^1].x) > double.Epsilon) ||
		    !(Math.Abs(poly[0].y - poly[^1].y) > double.Epsilon))
		{
			forceClosed = true;
			polyCount++;
		}

		var tempPoly = new PointF[polyCount];
		Parallel.For(0, forceClosed ? polyCount : polyCount - 1, (pt) => // (int pt = 0; pt < poly.Count; pt++)
		{
			tempPoly[pt] = new PointF((float)poly[pt].x, (float)poly[pt].y);
		});

		if (forceClosed)
		{
			tempPoly[^1] = new PointF(tempPoly[0].X, tempPoly[0].Y);
		}

		return tempPoly;
	}

	public static PointF[] closePoly(PointF[] poly)
	{
		if (!(Math.Abs(poly[0].X - poly[^1].X) > double.Epsilon) ||
		    !(Math.Abs(poly[0].Y - poly[^1].Y) > double.Epsilon))
		{
			return poly;
		}

		int polyLength = poly.Length;
		PointF[] tempPoly = new PointF[poly.Length + 1];
		Parallel.For(0, polyLength, (pt) =>
		{
			tempPoly[pt] = new PointF(poly[pt].X, poly[pt].Y);
		});
		tempPoly[^1] = new PointF(tempPoly[0].X, tempPoly[0].Y);
		return tempPoly;
	}

	private void tessPoly(PointF[] source, Color polyColor, float alpha)
	{
		Tess tess = new();

		int sourceLength = source.Length;
		ContourVertex[] contour = new ContourVertex[sourceLength];
		Parallel.For(0, sourceLength, (pt) =>
		{
			contour[pt].Position = new Vec3 { X = source[pt].X, Y = source[pt].Y, Z = 0 };
		});
		tess.AddContour(contour, ContourOrientation.Clockwise); // keep our orientation to allow holes to be handled.

		// Triangulate.
		tess.Tessellate(WindingRule.NonZero); // We don't have any hole polygons here.

		// Iterate triangles and create output geometry
		for (int i = 0; i < tess.ElementCount; i++)
		{
			PointF[] tempPoly = new PointF[3]; // 3 points.
			tempPoly[0] = new PointF((float)tess.Vertices[tess.Elements[i * 3]].Position.X, (float)tess.Vertices[tess.Elements[i * 3]].Position.Y);
			tempPoly[1] = new PointF((float)tess.Vertices[tess.Elements[i * 3 + 1]].Position.X, (float)tess.Vertices[tess.Elements[i * 3 + 1]].Position.Y);
			tempPoly[2] = new PointF((float)tess.Vertices[tess.Elements[i * 3 + 2]].Position.X, (float)tess.Vertices[tess.Elements[i * 3 + 2]].Position.Y);

			tessPolyList!.Add(new ovp_Poly(clockwiseOrder(tempPoly).ToArray(), polyColor, alpha));
		}
	}
}
