using Eto.Drawing;
using Eto.Forms;
using LibTessDotNet.Double;
using System;
using System.Collections.Generic;
using System.Linq;

namespace VeldridEto
{
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

        public List<ovp_Poly> polyList { get; set; }
        public List<int> polyListPtCount { get; set; }
        public List<int> polySourceIndex { get; set; } // will eventually track source of polygon, allowing for layer generating, etc. in output.
        public List<ovp_Poly> bgPolyList { get; set; }
        public List<int> bgPolyListPtCount { get; set; }
        public List<int> bgPolySourceIndex { get; set; } // will eventually track source of polygon, allowing for layer generating, etc. in output.
        public List<ovp_Poly> lineList { get; set; } // purely for lines.

        public List<float> lineZBiasList { get; set; }
        public List<float> polyZBiasList { get; set; }
        public List<float> tessPolyZBiasList { get; set; }

        public List<int> lineListPtCount { get; set; }
        public List<int> lineSourceIndex { get; set; } // will eventually track source of polygon, allowing for layer generating, etc. in output.
        public List<ovp_Poly> tessPolyList { get; set; } // triangles, but also need to track color. This is decoupled to allow boundary extraction without triangles getting in the way.
        public List<bool> drawnPoly { get; set; } // tracks whether the polygon corresponds to an enabled configuration or not.

        bool enableFilledPolys;
		bool showPoints;
		float base_zoom;
		float zoomFactor;
		Int32 zoomStep;
		bool allowZoomAndPan;
		bool dynamicGrid;
		bool panning;
		bool selecting;
		bool showGrid;
		bool showAxes;
		bool showDrawn;
		int grid_spacing;
		PointF cameraPosition;
		PointF default_cameraPosition;
		bool antiAlias;
		bool lockedViewport;

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
            if (lockedViewport)
            {
                return;
            }
            setCameraPos(default_cameraPosition.X, default_cameraPosition.Y);
            setZoomFactor(1.0f);
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
            if (lockedViewport)
            {
                return;
            }
            cameraPosition.X = x;
            changed = true;
        }

        public void setCameraY(float y)
        {
            if (lockedViewport)
            {
                return;
            }
            cameraPosition.Y = y;
            changed = true;
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
            if (base_zoom == val)
            {
                return;
            }
            base_zoom = val;
        }
        public void setZoomFactor(float val)
        {
            if (lockedViewport)
            {
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
            if (zoomStep == val)
            {
                return;
            }
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

		void pUpdateColors(Color newColor)
		{
			for (int poly = 0; poly < polyList.Count(); poly++)
			{
				polyList[poly].color = newColor;
			}
			for (int poly = 0; poly < tessPolyList.Count(); poly++)
			{
				tessPolyList[poly].color = newColor;
			}
            changed = true;
        }

        public void reset(bool clearBG = true)
		{
			pReset(clearBG);
		}

		void pReset(bool clearBG)
		{
			minX = 0;
			maxX = 0;
			minY = 0;
			maxY = 0;
			clear(clearBG);
			drawnPoly.Clear();
            changed = true;
		}

		public void clear(bool clearBG = true)
		{
			pClear(clearBG);
		}

		void pClear(bool clearBG)
		{
			polyList.Clear();
            polyZBiasList.Clear();
			polySourceIndex.Clear();
			polyListPtCount.Clear();
			if (clearBG)
			{
				bgPolyList.Clear();
				bgPolySourceIndex.Clear();
				bgPolyListPtCount.Clear();
			}
			lineList.Clear();
            lineZBiasList.Clear();
			lineSourceIndex.Clear();
			lineListPtCount.Clear();
			tessPolyList.Clear();
            changed = true;
        }

        public void addLine(PointF[] line, Color lineColor, float alpha, int layerIndex, float zBias = 0.0f)
		{
			pAddLine(line, lineColor, alpha, layerIndex, zBias);
		}

		void pAddLine(PointF[] line, Color lineColor, float alpha, int layerIndex, float zBias)
		{
			lineList.Add(new ovp_Poly(line, lineColor, alpha));
            lineZBiasList.Add(zBias);
			lineSourceIndex.Add(layerIndex);
			lineListPtCount.Add((line.Length - 1) * 2);
            changed = true;
        }

        public void addPolygon(PointF[] poly, Color polyColor, float alpha, bool drawn, int layerIndex, float zBias = 0.0f)
		{
			if (drawn)
			{
				// Drawn polygons are to be treated as lines : they don't get filled.
				addLine(poly, polyColor, alpha, layerIndex, zBias);
			}
			else
			{
				pAddPolygon(poly, polyColor, alpha, drawn, layerIndex, zBias);
			}
		}

		void pAddPolygon(PointF[] poly, Color polyColor, float alpha, bool drawn, int layerIndex, float zBias)
		{
			if (!drawn && enableFilledPolys) // avoid tessellation unless really necessary.
			{
				try
				{
					tessPoly(poly, polyColor, alpha, zBias);
				}
				catch (Exception)
				{

				}
			}

			PointF[] polys = checkPoly(poly);

			polyList.Add(new ovp_Poly(poly, polyColor, alpha));
			polySourceIndex.Add(layerIndex);
            polyZBiasList.Add(zBias);
			polyListPtCount.Add(poly.Length);
			drawnPoly.Add(drawn);
            changed = true;
        }

        public void addBGPolygon(PointF[] poly, Color polyColor, float alpha, int layerIndex)
		{
			pAddBGPolygon(poly, polyColor, alpha, layerIndex);
		}

		void pAddBGPolygon(PointF[] poly, Color polyColor, float alpha, int layerIndex)
		{
			PointF[] polys = checkPoly(poly);

			bgPolyList.Add(new ovp_Poly(poly, polyColor, alpha));
			bgPolyListPtCount.Add(poly.Length);
			bgPolySourceIndex.Add(layerIndex);
			drawnPoly.Add(false);
            changed = true;
        }

        public OVPSettings(float defX = 0.0f, float defY = 0.0f)
		{
			init(defX, defY);
		}

		void init(float defX, float defY)
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

        void pFullReset(float defX = 0.0f, float defY = 0.0f)
        {
            default_cameraPosition = new PointF(defX, defY);
			bgPolyList = new List<ovp_Poly>();
			bgPolySourceIndex = new List<int>();
			bgPolyListPtCount = new List<int>();
			polyList = new List<ovp_Poly>();
			polySourceIndex = new List<int>();
			polyListPtCount = new List<int>();
			lineList = new List<ovp_Poly>();
			lineSourceIndex = new List<int>();
			lineListPtCount = new List<int>();
			tessPolyList = new List<ovp_Poly>();
            lineZBiasList = new List<float>();
            polyZBiasList = new List<float>();
			drawnPoly = new List<bool>();
			zoomFactor = 1.0f;
			cameraPosition = new PointF(default_cameraPosition.X, default_cameraPosition.Y);
            changed = true;
		}

		PointF[] clockwiseOrder(PointF[] iPoints)
		{
			// Based on stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
			// Shoelace formula.

			double delta = 0;

			for (Int32 pt = 0; pt < iPoints.Length; pt++)
			{
				double deltaX = 0;
				double deltaY = 0;
				if (pt == iPoints.Length - 1)
				{
					deltaX = (iPoints[0].X - iPoints[pt].X);
					deltaY = (iPoints[0].Y + iPoints[pt].Y);
				}
				else
				{
					deltaX = (iPoints[pt + 1].X - iPoints[pt].X);
					deltaY = (iPoints[pt + 1].Y + iPoints[pt].Y);
				}

				delta += deltaX * deltaY;
			}

			if (delta > 0)
			{
				// clockwise
			}
			else
			{
				// counter-clockwise.
				Array.Reverse(iPoints);
			}

			return iPoints;
		}

		PointF[] checkPoly(PointF[] poly)
		{
			PointF[] source = poly.ToArray();

			if ((poly[0].X != poly[poly.Length - 1].X) && (poly[0].Y != poly[poly.Length - 1].Y))
			{
				PointF[] tempPoly = new PointF[poly.Length + 1];
				for (int pt = 0; pt < poly.Length; pt++)
				{
					tempPoly[pt] = new PointF(poly[pt].X, poly[pt].Y);
				}
				tempPoly[tempPoly.Length - 1] = new PointF(tempPoly[0].X, tempPoly[0].Y);
				source = tempPoly.ToArray();
			}

			return source;
		}

		void tessPoly(PointF[] source, Color polyColor, float alpha, float zBias)
		{
			// Now we need to check for polyfill, and triangulate the polygon if needed.
			//if (enableFilledPolys)
			{
				var tess = new Tess();

				ContourVertex[] contour = new ContourVertex[source.Length];
				for (int pt = 0; pt < contour.Length; pt++)
				{
					contour[pt].Position = new Vec3 { X = source[pt].X, Y = source[pt].Y, Z = 0 };
				}
				tess.AddContour(contour, ContourOrientation.Clockwise); // keep our orientation to allow holes to be handled.

				// Triangulate.
				tess.Tessellate(WindingRule.NonZero, ElementType.Polygons, 3); // We don't have any hole polygons here.

				// Iterate triangles and create output geometry
				for (int i = 0; i < tess.ElementCount; i++)
				{
					PointF[] tempPoly = new PointF[3]; // 3 points.
					tempPoly[0] = new PointF((float)tess.Vertices[tess.Elements[i * 3]].Position.X, (float)tess.Vertices[tess.Elements[i * 3]].Position.Y);
					tempPoly[1] = new PointF((float)tess.Vertices[tess.Elements[(i * 3) + 1]].Position.X, (float)tess.Vertices[tess.Elements[(i * 3) + 1]].Position.Y);
					tempPoly[2] = new PointF((float)tess.Vertices[tess.Elements[(i * 3) + 2]].Position.X, (float)tess.Vertices[tess.Elements[(i * 3) + 2]].Position.Y);

                    tessPolyZBiasList.Add(zBias);
					tessPolyList.Add(new ovp_Poly(clockwiseOrder(tempPoly).ToArray(), polyColor, alpha));
				}
			}
		}
	}
}
