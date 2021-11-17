using Eto.Drawing;
using LibTessDotNet.Double;
using System;
using System.Collections.Generic;
using System.Linq;

namespace etoViewport;

public class OVPSettings
{
    public float minX { get; set; }
    public float maxX { get; set; }
    public float minY { get; set; }
    public float maxY { get; set; }
    public bool enableFilledPolys { get; set; }
    public bool immediateMode { get; set; } // if true, don't use VBOs.
    public bool drawPoints { get; set; }
    public RectangleF bounds { get; set; }
    public float base_zoom { get; set; }
    public float zoomFactor { get; set; }
    public int zoomStep { get; set; }
    public bool allowZoomAndPan { get; set; }
    public bool dynamicGrid { get; set; }
    public bool panning { get; set; }
    public bool selecting { get; set; }
    public bool showGrid { get; set; }
    public bool showAxes { get; set; }
    public bool showDrawn { get; set; }
    public int gridSpacing { get; set; }
    public Color minorGridColor { get; set; }
    public Color majorGridColor { get; set; }
    public Color axisColor { get; set; }
    public Color backColor { get; set; }
    public Color selectionColor { get; set; }
    public Color inverSelectionColor { get; set; }
    public PointF cameraPosition { get; set; }
    public PointF default_cameraPosition { get; set; }
    public bool antiAlias { get; set; }
    public List<ovp_Poly> polyList { get; set; }
    public List<int> polyListPtCount { get; set; }
    public List<int> polySourceIndex { get; set; } // will eventually track source of polygon, allowing for layer generating, etc. in output.
    public List<ovp_Poly> bgPolyList { get; set; }
    public List<int> bgPolyListPtCount { get; set; }
    public List<int> bgPolySourceIndex { get; set; } // will eventually track source of polygon, allowing for layer generating, etc. in output.
    public List<ovp_Poly> lineList { get; set; } // purely for lines.
    public List<int> lineListPtCount { get; set; }
    public List<int> lineSourceIndex { get; set; } // will eventually track source of polygon, allowing for layer generating, etc. in output.
    public List<ovp_Poly> tessPolyList { get; set; } // triangles, but also need to track color. This is decoupled to allow boundary extraction without triangles getting in the way.
    public List<bool> drawnPoly { get; set; } // tracks whether the polygon corresponds to an enabled configuration or not.
    public bool lockedViewport { get; set; }

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
        for (int poly = 0; poly < polyList.Count(); poly++)
        {
            polyList[poly].color = newColor;
        }
        for (int poly = 0; poly < tessPolyList.Count(); poly++)
        {
            tessPolyList[poly].color = newColor;
        }
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
        drawnPoly.Clear();
    }

    public void clear(bool clearBG = true)
    {
        pClear(clearBG);
    }

    private void pClear(bool clearBG)
    {
        polyList.Clear();
        polySourceIndex.Clear();
        polyListPtCount.Clear();
        switch (clearBG)
        {
            case true:
                bgPolyList.Clear();
                bgPolySourceIndex.Clear();
                bgPolyListPtCount.Clear();
                break;
        }
        lineList.Clear();
        lineSourceIndex.Clear();
        lineListPtCount.Clear();
        tessPolyList.Clear();
    }

    public void addLine(PointF[] line, Color lineColor, float alpha, int layerIndex)
    {
        pAddLine(line, lineColor, alpha, layerIndex);
    }

    private void pAddLine(PointF[] line, Color lineColor, float alpha, int layerIndex)
    {
        lineList.Add(new ovp_Poly(line, lineColor, alpha));
        lineSourceIndex.Add(layerIndex);
        lineListPtCount.Add((line.Length - 1) * 2);
    }

    public void addPolygon(PointF[] poly, Color polyColor, float alpha, bool drawn, int layerIndex)
    {
        switch (drawn)
        {
            case true:
                // Drawn polygons are to be treated as lines : they don't get filled.
                addLine(poly, polyColor, alpha, layerIndex);
                break;
            default:
                pAddPolygon(poly, polyColor, alpha, drawn, layerIndex);
                break;
        }
    }

    private void pAddPolygon(PointF[] poly, Color polyColor, float alpha, bool drawn, int layerIndex)
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

        PointF[] polys = checkPoly(poly);

        polyList.Add(new ovp_Poly(poly, polyColor, alpha));
        polySourceIndex.Add(layerIndex);
        polyListPtCount.Add(poly.Length);
        drawnPoly.Add(drawn);
    }

    public void addBGPolygon(PointF[] poly, Color polyColor, float alpha, int layerIndex)
    {
        pAddBGPolygon(poly, polyColor, alpha, layerIndex);
    }

    private void pAddBGPolygon(PointF[] poly, Color polyColor, float alpha, int layerIndex)
    {
        PointF[] polys = checkPoly(poly);

        bgPolyList.Add(new ovp_Poly(poly, polyColor, alpha));
        bgPolyListPtCount.Add(poly.Length);
        bgPolySourceIndex.Add(layerIndex);
        drawnPoly.Add(false);
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
        drawPoints = true;
        dynamicGrid = true;
        panning = false;
        selecting = false;
        showGrid = true;
        showAxes = true;
        gridSpacing = 10;
        antiAlias = true;
        zoomStep = 1;
        immediateMode = false;
        fullReset(defX, defY);
    }

    public void fullReset(float defX = 0.0f, float defY = 0.0f)
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
        drawnPoly = new List<bool>();
        zoomFactor = 1.0f;
        cameraPosition = new PointF(default_cameraPosition.X, default_cameraPosition.Y);
    }

    private PointF[] clockwiseOrder(PointF[] iPoints)
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

    private PointF[] checkPoly(PointF[] poly)
    {
        PointF[] source = poly.ToArray();

        if (poly[0].X != poly[poly.Length - 1].X && poly[0].Y != poly[poly.Length - 1].Y)
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

    private void tessPoly(PointF[] source, Color polyColor, float alpha)
    {
        // Now we need to check for polyfill, and triangulate the polygon if needed.
        //if (enableFilledPolys)
        {
            Tess tess = new();

            ContourVertex[] contour = new ContourVertex[source.Length];
            for (int pt = 0; pt < contour.Length; pt++)
            {
                contour[pt].Position = new Vec3 { X = source[pt].X, Y = source[pt].Y, Z = 0 };
            }
            tess.AddContour(contour, ContourOrientation.Clockwise); // keep our orientation to allow holes to be handled.

            // Triangulate.
            tess.Tessellate(WindingRule.Positive, ElementType.Polygons, 3); // We don't have any hole polygons here.

            // Iterate triangles and create output geometry
            for (int i = 0; i < tess.ElementCount; i++)
            {
                PointF[] tempPoly = new PointF[3]; // 3 points.
                tempPoly[0] = new PointF((float)tess.Vertices[tess.Elements[i * 3]].Position.X, (float)tess.Vertices[tess.Elements[i * 3]].Position.Y);
                tempPoly[1] = new PointF((float)tess.Vertices[tess.Elements[i * 3 + 1]].Position.X, (float)tess.Vertices[tess.Elements[i * 3 + 1]].Position.Y);
                tempPoly[2] = new PointF((float)tess.Vertices[tess.Elements[i * 3 + 2]].Position.X, (float)tess.Vertices[tess.Elements[i * 3 + 2]].Position.Y);

                tessPolyList.Add(new ovp_Poly(clockwiseOrder(tempPoly).ToArray(), polyColor, alpha));
            }
        }
    }
}