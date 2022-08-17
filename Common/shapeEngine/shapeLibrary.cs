using System.Diagnostics;
using geoLib;
using geoWrangler;
using utility;

namespace shapeEngine;

public class ShapeLibrary
{
    private static readonly List<string> availableShapes_all = new()
    {
        "(None)", "Rectangle/Square", "L-shape", "T-shape", "X-shape", "U-shape", "S-shape", "GDS/Oasis", "Boolean",
        "Text", "Bounding", "Layout"
    };

    public enum shapeNames_all
    {
        none,
        rect,
        Lshape,
        Tshape,
        Xshape,
        Ushape,
        Sshape,
        GEOCORE,
        BOOLEAN,
        text,
        bounding,
        complex
    }

    // Client sends an array that has mapping of the client shape index to the ShapeLibrary index.
    private int[]? shapeMapping_fromClient;

    public void shapesForClient(int[] clientShapeDefinition)
    {
        pShapesForClient(clientShapeDefinition);
    }

    void pShapesForClient(int[] clientShapeDefinition)
    {
        if (clientShapeDefinition.Length > availableShapes_all.Count)
        {
            throw new Exception("More shapes requested than are supported");
        }

        // Initialize to no-shapes
        shapeMapping_fromClient = new int[availableShapes_all.Count];
        for (int i = 0; i < availableShapes_all.Count; i++)
        {
            shapeMapping_fromClient[i] = -1;
        }

        // Set the support flags from the client.
        for (int i = 0; i < clientShapeDefinition.Length; i++)
        {
            shapeMapping_fromClient[i] = clientShapeDefinition[i];
        }

    }

    public static List<string> getAvailableShapes(int[] clientShapeDefinition)
    {
        // Set the support flags from the client.
        return clientShapeDefinition.Select(t => availableShapes_all[t]).ToList();
    }

    public int shapeIndex;
    public bool shapeValid { get; private set; }
    public bool geoCoreShapeOrthogonal { get; private set; }
    public MyVertex[] Vertex { get; private set; }
    public MyRound[] round1 { get; private set; }
    public bool[] tips { get; private set; }

    private class BoundingBox
    {
        private GeoLibPointF midPoint;

        public GeoLibPointF getMidPoint()
        {
            return pGetMidPoint();
        }

        private GeoLibPointF pGetMidPoint()
        {
            return midPoint;
        }

        public BoundingBox(List<GeoLibPointF> incomingPoints)
        {
            midPoint = new GeoLibPointF(0.0f, 0.0f);
            pBoundingBox(incomingPoints);
        }

        private void pBoundingBox(List<GeoLibPointF> incomingPoints)
        {
            if (incomingPoints.Count > 0)
            {
                var minX = incomingPoints.Min(p => p.X);
                var minY = incomingPoints.Min(p => p.Y);
                var maxX = incomingPoints.Max(p => p.X);
                var maxY = incomingPoints.Max(p => p.Y);
                midPoint = new GeoLibPointF(minX + (maxX - minX) / 2.0f, minY + (maxY - minY) / 2.0f);
            }
        }
    }

    public GeoLibPointF getPivotPoint()
    {
        return pGetPivotPoint();
    }

    private GeoLibPointF pGetPivotPoint()
    {
        int limit = Vertex.Length - 1;
        GeoLibPointF[] t = new GeoLibPointF[limit]; // closed shape, we don't need the final point
#if !SHAPELIBSINGLETHREADED
        Parallel.For(0, limit, i =>
#else
            for (int i = 0; i < t.Length; i++)
#endif
            {
                t[i] = new GeoLibPointF(Vertex[i].X, Vertex[i].Y);
            }
#if !SHAPELIBSINGLETHREADED
        );
#endif
        GeoLibPointF pivot = GeoWrangler.midPoint(t);

        return pivot;
    }

    private ShapeSettings layerSettings = new();

    public ShapeLibrary(int[] shapes, ShapeSettings shapeSettings)
    {
        Vertex = new MyVertex[1];
        round1 = new MyRound[1];
        tips = new bool[1];
        pShapeLibrary(shapes, shapeSettings);
    }

    private void pShapeLibrary(int[] shapes, ShapeSettings shapeSettings)
    {
        shapeValid = false;
        layerSettings = shapeSettings;
        pShapesForClient(shapes);
    }

    public ShapeLibrary(int[] shapes, int shapeIndex_, ShapeSettings shapeSettings)
    {
        Vertex = new MyVertex[1];
        round1 = new MyRound[1];
        tips = new bool[1];
        pShapeLibrary(shapes, shapeIndex_, shapeSettings);
    }

    private void pShapeLibrary(int[] shapes, int shapeIndex_, ShapeSettings shapeSettings)
    {
        pShapesForClient(shapes);
        shapeIndex = shapeIndex_;
        shapeValid = false;
        layerSettings = shapeSettings;
        pSetShape(shapeIndex);
    }

    public void setShape(int shapeIndex_, GeoLibPointF[]? sourcePoly = null)
    {
        pSetShape(shapeIndex_, sourcePoly);
    }

    private void pSetShape(int shapeIndex_, GeoLibPointF[]? sourcePoly = null)
    {
        try
        {
            Debug.Assert(shapeMapping_fromClient != null, nameof(shapeMapping_fromClient) + " != null");
            shapeIndex = shapeMapping_fromClient[shapeIndex_];
            switch (shapeIndex)
            {
                case (int)shapeNames_all.rect:
                case (int)shapeNames_all.text:
                case (int)shapeNames_all.bounding:
                    rectangle();
                    break;
                case (int)shapeNames_all.Lshape:
                    Lshape();
                    break;
                case (int)shapeNames_all.Tshape:
                    Tshape();
                    break;
                case (int)shapeNames_all.Xshape:
                    crossShape();
                    break;
                case (int)shapeNames_all.Ushape:
                    Ushape();
                    break;
                case (int)shapeNames_all.Sshape:
                    Sshape();
                    break;
                case (int)shapeNames_all.GEOCORE:
                    if (layerSettings.getInt(ShapeSettings.properties_i.gCSEngine) == 1)
                    {
                        customShape(sourcePoly);
                    }

                    break;
                case (int)shapeNames_all.complex:
                    customShape(sourcePoly);
                    break;
                default:
                    throw new Exception();
            }
        }
        catch (Exception)
        {
            Vertex = new MyVertex[1];
            tips = new[] { false };
            layerSettings = new();
            round1 = new MyRound[1];
        }
    }

    public static int getSubShapeCount(int index)
    {
        return pGetSubShapeCount(index);
    }

    private static int pGetSubShapeCount(int index)
    {
        switch (index)
        {
            case (int)shapeNames_all.Lshape:
            case (int)shapeNames_all.Tshape:
            case (int)shapeNames_all.Xshape:
            case (int)shapeNames_all.Ushape:
                return 2;
            case (int)shapeNames_all.Sshape:
                return 3;
            default:
                return 1;
        }
    }

    private void configureArrays()
    {
        double vertexCount = 0;
        switch (shapeIndex)
        {
            case (int)shapeNames_all.rect: // rectangle
            case (int)shapeNames_all.text:
            case (int)shapeNames_all.bounding:
                vertexCount = 9;
                break;
            case (int)shapeNames_all.Lshape: // L
                vertexCount = 13;
                break;
            case (int)shapeNames_all.Tshape: // T
                vertexCount = 17;
                break;
            case (int)shapeNames_all.Xshape: // Cross
                vertexCount = 25;
                break;
            case (int)shapeNames_all.Ushape: // U
                vertexCount = 17;
                break;
            case (int)shapeNames_all.Sshape: // S
                vertexCount = 25;
                break;
        }

        int arrayLength = (int)vertexCount;
        Vertex = new MyVertex[arrayLength];
        tips = new bool[arrayLength];
#if !SHAPELIBSINGLETHREADED
        Parallel.For(0, arrayLength, i =>
#else
            for (Int32 i = 0; i < arrayLength; i++)
#endif
            {
                tips[i] = false;
            }
#if !SHAPELIBSINGLETHREADED
        );
#endif
        arrayLength = (int)Math.Floor(vertexCount / 2) + 1;
        round1 = new MyRound[arrayLength];
#if !SHAPELIBSINGLETHREADED
        Parallel.For(0, arrayLength, i =>
#else
            for (Int32 i = 0; i < arrayLength; i++)
#endif
            {
                round1[i] = new MyRound();
            }
#if !SHAPELIBSINGLETHREADED
        );
#endif
    }

    private void rectangle()
    {
        configureArrays();

        // Sort out the tips by setting 'true' to center vertex that defines tip in each case.
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[1] = true;
                tips[8] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[1] = true;
                tips[8] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[3] = true;
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[1] = true;
                tips[3] = true;
                tips[8] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[3] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[1] = true;
                tips[3] = true;
                tips[5] = true;
                tips[8] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[1] = true;
                tips[7] = true;
                tips[8] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[5] = true;
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[1] = true;
                tips[5] = true;
                tips[7] = true;
                tips[8] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[3] = true;
                tips[7] = true;
                tips[1] = true;
                tips[8] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[3] = true;
                tips[7] = true;
                tips[5] = true;
                break;
            default: // All
                tips[3] = true;
                tips[7] = true;
                tips[1] = true;
                tips[8] = true;
                tips[5] = true;
                break;
        }

        // NOTE:
        // Subshape offsets are applied later to simplify ellipse generation
        // Horizontal and vertical global offsets are applied in callsite.
        // Repositioning with respect to subshape reference is also applied in callsite.

        double tmpX = 0.0;
        double tmpY = 0.0;
        Vertex[0] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0)) / 2;
        Vertex[1] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0));
        Vertex[2] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0)) / 2;
        Vertex[3] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));
        Vertex[4] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0)) / 2;
        Vertex[5] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY = 0.0;
        Vertex[6] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0)) / 2;
        tmpY = 0.0;
        Vertex[7] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        processEdgesForRounding();

        shapeValid = true;
    }

    private void Tshape()
    {
        configureArrays();
        // Sort out the tips by setting 'true' to center vertex that defines tip in each case.
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[1] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[5] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[1] = true;
                tips[5] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[3] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[3] = true;
                tips[1] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[3] = true;
                tips[5] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[1] = true;
                tips[3] = true;
                tips[5] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[15] = true;
                tips[1] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[15] = true;
                tips[5] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[1] = true;
                tips[15] = true;
                tips[5] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[1] = true;
                tips[3] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[3] = true;
                tips[5] = true;
                tips[13] = true;
                tips[15] = true;
                break;
            default: // All
                tips[1] = true;
                tips[3] = true;
                tips[5] = true;
                tips[13] = true;
                tips[15] = true;
                break;
        }

        switch (layerSettings.getInt(ShapeSettings.properties_i.subShape2TipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[7] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[7] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[7] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[11] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[11] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[7] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[7] = true;
                tips[9] = true;
                tips[11] = true;
                break;
            default: // All
                tips[7] = true;
                tips[9] = true;
                tips[11] = true;
                break;
        }

        // NOTE:
        // Subshape offsets are applied later to simplify ellipse generation
        // Horizontal and vertical global offsets are applied in callsite.

        double tmpX = 0.0;
        double tmpY = 0.0;
        Vertex[0] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0)) / 2;
        Vertex[1] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0));
        Vertex[2] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0)) / 2;
        Vertex[3] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0)) / 2;
        Vertex[4] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY = (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0)) -
                (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1)) +
                 (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1)) -
                  Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 0))))) / 2;
        Vertex[5] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1)) +
               Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1));
        Vertex[6] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1) / 2);
        Vertex[7] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1) / 2);
        Vertex[8] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1) / 2);
        Vertex[9] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1) / 2);
        Vertex[10] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1) / 2);
        Vertex[11] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        tmpX -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1) / 2);
        Vertex[12] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1)) / 2;
        Vertex[13] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY = 0;
        Vertex[14] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0) / 2);
        Vertex[15] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        processEdgesForRounding();

        shapeValid = true;
    }

    private void Lshape()
    {
        configureArrays();
        // Sort out the tips by setting 'true' to center vertex that defines tip in each case.
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[1] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[1] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[11] = true;
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[1] = true;
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[3] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[1] = true;
                tips[3] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[11] = true;
                tips[1] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[11] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[1] = true;
                tips[5] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[1] = true;
                tips[11] = true;
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[11] = true;
                tips[3] = true;
                tips[5] = true;
                break;
            default: // All
                tips[1] = true;
                tips[5] = true;
                tips[11] = true;
                tips[3] = true;
                break;
        }

        switch (layerSettings.getInt(ShapeSettings.properties_i.subShape2TipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[7] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[7] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[7] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[11] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[11] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[7] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[7] = true;
                tips[11] = true;
                tips[9] = true;
                break;
            default: // All
                tips[9] = true;
                tips[7] = true;
                tips[11] = true;
                break;
        }

        // NOTE:
        // Subshape offsets are applied later to simplify ellipse generation
        // Horizontal and vertical global offsets are applied in callsite.
        // Repositioning with respect to subshape reference is also applied in callsite.
        double tmpX = 0.0;
        double tmpY = 0.0;
        Vertex[0] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0)) / 2;
        Vertex[1] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0)) / 2;
        Vertex[2] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0)) / 2;
        Vertex[3] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0)) / 2;
        Vertex[4] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0)) -
                 Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1))) / 2;
        Vertex[5] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY -= (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0)) -
                 Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1))) / 2;
        Vertex[6] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1)) / 2;
        Vertex[7] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1)) / 2;
        Vertex[8] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1)) / 2;
        Vertex[9] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1)) / 2;
        Vertex[10] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX -= (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0)) +
                 Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1))) / 2;
        Vertex[11] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        processEdgesForRounding();

        shapeValid = true;
    }

    private void Ushape()
    {
        configureArrays();
        // Sort out the tips by setting 'true' to center vertex that defines tip in each case.
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[1] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[1] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[3] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[3] = true;
                tips[11] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[1] = true;
                tips[3] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[3] = true;
                tips[11] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[1] = true;
                tips[3] = true;
                tips[11] = true;
                tips[13] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[1] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[13] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[1] = true;
                tips[13] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[1] = true;
                tips[3] = true;
                tips[11] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[3] = true;
                tips[11] = true;
                tips[13] = true;
                tips[15] = true;
                break;
            default: // All
                tips[1] = true;
                tips[3] = true;
                tips[11] = true;
                tips[13] = true;
                tips[15] = true;
                break;
        }

        switch (layerSettings.getInt(ShapeSettings.properties_i.subShape2TipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[5] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[5] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[5] = true;
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[7] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[5] = true;
                tips[7] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[7] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[7] = true;
                tips[9] = true;
                break;
            default: // All
                tips[5] = true;
                tips[7] = true;
                tips[9] = true;
                break;
        }

        // NOTE:
        // Subshape offsets are applied later to simplify ellipse generation
        // Horizontal and vertical global offsets are applied in callsite.
        // Repositioning with respect to subshape reference is also applied in callsite.
        double tmpX = 0.0;
        double tmpY = 0.0;
        Vertex[0] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0)) / 2;
        Vertex[1] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0));
        Vertex[2] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horOffset, 1)) / 2;
        Vertex[3] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horOffset, 1)) / 2;
        Vertex[4] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1)) / 2;
        Vertex[5] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1)) / 2;
        Vertex[6] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1)) / 2;
        Vertex[7] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1)) / 2;
        Vertex[8] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1)) / 2;
        Vertex[9] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0));
        Vertex[10] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0)) -
                 (Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horOffset, 1)) +
                  Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1)))) / 2;
        Vertex[11] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));
        Vertex[12] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY /= 2;
        Vertex[13] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY = 0;
        Vertex[14] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX /= 2;
        Vertex[15] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        processEdgesForRounding();

        shapeValid = true;
    }

    private void crossShape()
    {
        configureArrays();
        // Sort out the tips.
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[1] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[1] = true;
                tips[9] = true;
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[11] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[1] = true;
                tips[9] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[11] = true;
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[1] = true;
                tips[9] = true;
                tips[11] = true;
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[1] = true;
                tips[9] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[13] = true;
                tips[21] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[1] = true;
                tips[9] = true;
                tips[13] = true;
                tips[21] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[11] = true;
                tips[1] = true;
                tips[9] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[13] = true;
                tips[11] = true;
                tips[23] = true;
                tips[21] = true;
                break;
            default: // All
                tips[1] = true;
                tips[9] = true;
                tips[13] = true;
                tips[21] = true;
                tips[11] = true;
                tips[23] = true;
                break;
        }

        switch (layerSettings.getInt(ShapeSettings.properties_i.subShape2TipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[5] = true;
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[7] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[3] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[3] = true;
                tips[7] = true;
                tips[15] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[5] = true;
                tips[7] = true;
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[7] = true;
                tips[15] = true;
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[5] = true;
                tips[7] = true;
                tips[15] = true;
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[3] = true;
                tips[5] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[3] = true;
                tips[17] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[3] = true;
                tips[5] = true;
                tips[17] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[5] = true;
                tips[3] = true;
                tips[7] = true;
                tips[15] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[3] = true;
                tips[7] = true;
                tips[15] = true;
                tips[17] = true;
                tips[19] = true;
                break;
            default: // All
                tips[5] = true;
                tips[17] = true;
                tips[3] = true;
                tips[7] = true;
                tips[15] = true;
                tips[19] = true;
                break;
        }

        // NOTE:
        // Subshape offsets are applied later to simplify ellipse generation
        // Horizontal and vertical global offsets are applied in callsite.

        double tmpX = 0.0;
        double tmpY = 0.0;
        Vertex[0] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1)) / 2;
        Vertex[1] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1));
        Vertex[2] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horOffset, 1));
        Vertex[3] = new MyVertex(tmpX / 2, tmpY, typeDirection.down1, false, false, typeVertex.center);

        Vertex[4] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        double tmpY2 = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1)) / 2;
        tmpY += tmpY2;
        Vertex[5] = new MyVertex(tmpX, tmpY / 2, typeDirection.left1, true, false, typeVertex.center);

        tmpY += tmpY2;
        Vertex[6] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX /= 2;
        Vertex[7] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX = 0.0;
        Vertex[8] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY2 = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1)) +
                Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1));
        tmpY2 = tmpY - tmpY2;
        tmpY -= (tmpY2 * 0.5);
        Vertex[9] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0));
        Vertex[10] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0)) / 2;
        Vertex[11] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));
        Vertex[12] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY2 = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1)) +
                Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1));
        tmpY2 = tmpY - tmpY2;
        tmpY -= (tmpY2 * 0.5);
        Vertex[13] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1)) +
               Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1));
        Vertex[14] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        // Need midpoint of edge
        double tmpX2 = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1)) +
                Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horOffset, 1));
        tmpX2 = tmpX - tmpX2;
        tmpX -= (tmpX2 * 0.5);
        Vertex[15] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1));
        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horOffset, 1));
        Vertex[16] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1)) / 2;
        Vertex[17] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1)) / 2;
        Vertex[18] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        // Need midpoint of edge
        tmpX2 = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1)) +
                       Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horOffset, 1));
        tmpX2 = tmpX - tmpX2;
        tmpX -= (tmpX2 * 0.5);
        Vertex[19] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));
        Vertex[20] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1)) / 2;
        Vertex[21] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);

        tmpY = 0.0;
        Vertex[22] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0)) / 2;
        Vertex[23] = new MyVertex(tmpX, tmpY, typeDirection.down1, true, false, typeVertex.center);

        processEdgesForRounding();

        shapeValid = true;
    }

    private void Sshape()
    {
        configureArrays();
        // Sort out the tips by setting 'true' to center vertex that defines tip in each case.
        switch (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex))
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[1] = true;
                tips[9] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[1] = true;
                tips[9] = true;
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[11] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[1] = true;
                tips[9] = true;
                tips[11] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[11] = true;
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[1] = true;
                tips[9] = true;
                tips[11] = true;
                tips[13] = true;
                tips[21] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[1] = true;
                tips[9] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[13] = true;
                tips[21] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[1] = true;
                tips[9] = true;
                tips[13] = true;
                tips[21] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[1] = true;
                tips[9] = true;
                tips[11] = true;
                tips[23] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[11] = true;
                tips[13] = true;
                tips[21] = true;
                tips[23] = true;
                break;
            default: // All
                tips[1] = true;
                tips[9] = true;
                tips[11] = true;
                tips[13] = true;
                tips[21] = true;
                tips[23] = true;
                break;
        }

        switch (layerSettings.getInt(ShapeSettings.properties_i.subShape2TipLocIndex)) // Bottom notch
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[3] = true;
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[3] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[3] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[3] = true;
                tips[5] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[5] = true;
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[5] = true;
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[3] = true;
                tips[7] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[3] = true;
                tips[5] = true;
                tips[7] = true;
                break;
            default: // All
                tips[3] = true;
                tips[5] = true;
                tips[7] = true;
                break;
        }

        switch (layerSettings.getInt(ShapeSettings.properties_i.subShape3TipLocIndex)) // Top notch
        {
            case (int)ShapeSettings.tipLocations.none: // None
                break;
            case (int)ShapeSettings.tipLocations.L: // Left
                break;
            case (int)ShapeSettings.tipLocations.R: // Right
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.LR: // Left and right
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.T: // Top
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.B: // Bottom
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.TB: // Top and Bottom
                tips[15] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TL: // Top and Left
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TR: // Top and Right
                tips[17] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TLR: // Top and Left and Right
                tips[17] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.BL: // Bottom and Left
                tips[15] = true;
                break;
            case (int)ShapeSettings.tipLocations.BR: // Bottom and Right
                tips[15] = true;
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.BLR: // Bottom and Left and Right
                tips[15] = true;
                tips[17] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBL: // Top and Bottom and Left
                tips[17] = true;
                tips[19] = true;
                break;
            case (int)ShapeSettings.tipLocations.TBR: // Top and Bottom and Right
                tips[15] = true;
                tips[17] = true;
                tips[19] = true;
                break;
            default: // All
                tips[15] = true;
                tips[17] = true;
                tips[19] = true;
                break;
        }

        // NOTE:
        // Subshape offsets are applied later to simplify ellipse generation
        // Horizontal and vertical global offsets are applied in callsite.
        // Repositioning with respect to subshape reference is also applied in callsite.
        double tmpX = 0.0;
        double tmpY = 0.0;
        Vertex[0] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1)) / 2;
        Vertex[1] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1));
        Vertex[2] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1)) / 2;
        Vertex[3] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1));
        Vertex[4] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1)) / 2;
        Vertex[5] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1)) +
               Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1));
        Vertex[6] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 1)) / 2;
        Vertex[7] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        tmpX = 0;
        Vertex[8] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        double tmpY2 = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 1)) +
                Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 1));
        tmpY2 = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0)) - tmpY2;
        tmpY += (tmpY2 * 0.5);

        Vertex[9] = new MyVertex(tmpX, tmpY, typeDirection.left1, true, false, typeVertex.center);

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0));
        Vertex[10] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0)) / 2;
        Vertex[11] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);
        // Center so no rounding definition

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));
        Vertex[12] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 2)) / 2;
        Vertex[13] = new MyVertex(tmpX, tmpY, typeDirection.right1, true, false, typeVertex.center);
        // Center so no rounding definition

        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0)) -
               Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 2));
        Vertex[14] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 2) / 2);
        Vertex[15] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horOffset, 2));
        Vertex[16] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY -= Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 2) / 2);
        Vertex[17] = new MyVertex(tmpX, tmpY, typeDirection.right1, false, false, typeVertex.center);

        tmpY2 = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 2)) +
                       Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verOffset, 2));
        tmpY = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.verLength, 0)) - tmpY2;

        Vertex[18] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX += Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 2) / 2);
        Vertex[19] = new MyVertex(tmpX, tmpY, typeDirection.up1, false, false, typeVertex.center);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0));
        Vertex[20] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpY /= 2;
        Vertex[21] = new MyVertex(tmpX, tmpY, typeDirection.right1, false, false, typeVertex.center);

        tmpY = 0;
        Vertex[22] = new MyVertex(tmpX, tmpY, typeDirection.tilt1, true, false, typeVertex.corner);

        tmpX = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.horLength, 0)) / 2;
        Vertex[23] = new MyVertex(tmpX, tmpY, typeDirection.down1, false, false, typeVertex.center);

        processEdgesForRounding();

        shapeValid = true;
    }

    private void processEdgesForRounding()
    {
        int horEdge = Vertex.Length - 2; // deal with padding.
        int verEdge = 1;
        for (int r = 0; r < round1.Length - 1; r++)
        {
            try
            {
                round1[r].index = r * 2;
                round1[r].horFace = horEdge;
                round1[r].verFace = verEdge;

                // Figure out our corner type. First is a special case.
                if (r == 0)
                {
                    round1[r].direction = typeRound.exter;
                    round1[r].MaxRadius =
                        Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.oCR));
                }
                else
                {
                    if (
                        Vertex[round1[r].verFace].direction == typeDirection.right1 &&
                        Vertex[round1[r].horFace].direction == typeDirection.up1 && horEdge < verEdge ||
                        Vertex[round1[r].verFace].direction == typeDirection.left1 &&
                        Vertex[round1[r].horFace].direction == typeDirection.up1 && horEdge > verEdge ||
                        Vertex[round1[r].verFace].direction == typeDirection.right1 &&
                        Vertex[round1[r].horFace].direction == typeDirection.down1 && horEdge > verEdge ||
                        Vertex[round1[r].verFace].direction == typeDirection.left1 &&
                        Vertex[round1[r].horFace].direction == typeDirection.down1 && horEdge < verEdge
                    )
                    {
                        round1[r].direction = typeRound.exter;
                        round1[r].MaxRadius =
                            Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.oCR));
                    }
                    else
                    {
                        round1[r].direction = typeRound.inner;
                        round1[r].MaxRadius =
                            Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.iCR));
                    }
                }

                // Small fudge for the 0 case
                if (r == 0)
                {
                    horEdge = -1;
                }

                // Change our edge configuration for the next loop. We need to handle overflow references as well
                if (r % 2 == 0)
                {
                    horEdge += 4;
                    horEdge %= Vertex.Length;
                }
                else
                {
                    verEdge += 4;
                    verEdge %= Vertex.Length;
                }
            }
            catch
            {
                // ignored
            }
        }

        // First and last are the same.
        round1[^1] = round1[0];
    }

    // Intended to take geometry from an external source and map it into our shape engine.
    private void customShape(GeoLibPointF[]? sourcePoly)
    {
        if (sourcePoly == null)
        {
            shapeValid = false;
            return;
        }

        // Note that we assume the point order matches our general primitives; might need upstream review to ensure this is being
        // fed correctly.
        // Upstream should trim array to ensure end point is different from start point, but we'll force the issue here for robustness.
        sourcePoly = GeoWrangler.stripTerminators(sourcePoly, true);
        sourcePoly = GeoWrangler.stripColinear(sourcePoly);
        // Remove duplicate points in case definition was badly constructed.
        sourcePoly = GeoWrangler.removeDuplicates(sourcePoly);
        //  Strip the terminator again to meet the requirements below.
        sourcePoly = GeoWrangler.stripTerminators(sourcePoly, false);
        sourcePoly = GeoWrangler.clockwise(sourcePoly);

        // We need to look at our incoming shape to see whether it's orthogonal and suitable for contouring.
        geoCoreShapeOrthogonal = GeoWrangler.orthogonal(sourcePoly, angularTolerance: 0.003);

        if (!geoCoreShapeOrthogonal)
        {
            customShape_nonOrthogonal(sourcePoly);
        }
        else
        {
            customShape_orthogonal(sourcePoly);
        }

        shapeValid = true;
    }

    private void customShape_nonOrthogonal(GeoLibPointF[] sourcePoly)
    {
        int sCount = sourcePoly.Length;
        Vertex = new MyVertex[sCount + 1]; // add one to close.
        tips = new bool[sCount + 1];
        // Assign shape vertices to Vertex and move on. EntropyShape will know what to do.
#if !SHAPELIBSINGLETHREADED
        Parallel.For(0, sCount, pt =>
#else
            for (int pt = 0; pt < sCount; pt++)
#endif
            {
                Vertex[pt] = new MyVertex(sourcePoly[pt].X, sourcePoly[pt].Y, typeDirection.tilt1, false, false,
                    typeVertex.corner);
                tips[pt] = false;
            }
#if !SHAPELIBSINGLETHREADED
        );
#endif
        // Close the shape.
        Vertex[^1] = new MyVertex(Vertex[0]);
        tips[^1] = false;
    }

    private void customShape_orthogonal(GeoLibPointF[] sourcePoly)
    {
        int sCount = sourcePoly.Length;
        int vertexCount = sCount * 2 + 1; // assumes no point in midpoint of edges, and 1 to close.
        Vertex = new MyVertex[vertexCount];
        tips = new bool[vertexCount];
        int vertexCounter = 0; // set up our vertex counter.

#if !SHAPELIBSINGLETHREADED
        Parallel.For(0, vertexCount, i =>
#else
            for (Int32 i = 0; i < vertexCount; i++)
#endif
            {
                tips[i] = false;
            }
#if !SHAPELIBSINGLETHREADED
        );
#endif

        int roundCount = sourcePoly.Length + 1;
        round1 = new MyRound[roundCount];
#if !SHAPELIBSINGLETHREADED
        Parallel.For(0, roundCount, i =>
#else
            for (Int32 i = 0; i < roundCount; i++)
#endif
            {
                round1[i] = new MyRound();
            }
#if !SHAPELIBSINGLETHREADED
        );
#endif
        // Set up first rounding entry
        round1[0].direction = typeRound.exter;
        round1[0].MaxRadius = Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.oCR));
        round1[0].verFace = 1;
        round1[0].horFace = vertexCount - 2;
        round1[^1] = round1[0]; // close the loop

        // Set up first vertex.
        Vertex[0] = new MyVertex(sourcePoly[0].X, sourcePoly[0].Y, typeDirection.tilt1, false, false,
            typeVertex.corner);
        vertexCounter++;
        // Set up first midpoint.
        Vertex[1] = new MyVertex((sourcePoly[0].X + sourcePoly[1].X) / 2.0f, (sourcePoly[0].Y + sourcePoly[1].Y) / 2.0f,
            typeDirection.left1, true, false, typeVertex.center);
        if (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.L ||
            layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
            (int)ShapeSettings.tipLocations.LR ||
            layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
            (int)ShapeSettings.tipLocations.BL ||
            layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
            (int)ShapeSettings.tipLocations.TL ||
            layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
            (int)ShapeSettings.tipLocations.TBL ||
            layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
            (int)ShapeSettings.tipLocations.BLR ||
            layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
            (int)ShapeSettings.tipLocations.TLR ||
            layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) == (int)ShapeSettings.tipLocations.all)
        {
            tips[vertexCounter] = true;
        }

        vertexCounter++;

        // Also set our end points
        Vertex[vertexCount - 2] = new MyVertex((sourcePoly[0].X + sourcePoly[^1].X) / 2.0f,
            (sourcePoly[0].Y + sourcePoly[^1].Y) / 2.0f, typeDirection.down1, false, false, typeVertex.center);

        // Figure out our rounding characteristics.

        // First edge is always vertical, left facing.
        bool left = true;
        bool up = false;

        for (int pt = 1; pt < roundCount - 1; pt++)
        {
            // Link to our vertical/horizontal edges
            round1[pt].index = vertexCounter;
            if (pt % 2 == 1)
            {
                round1[pt].verFace = vertexCounter - 1;
                round1[pt].horFace = vertexCounter + 1;
            }
            else
            {
                round1[pt].verFace = vertexCounter + 1;
                round1[pt].horFace = vertexCounter - 1;
            }

            // Register our corner point into the vertex array.
            Vertex[vertexCounter] = new MyVertex(sourcePoly[pt].X, sourcePoly[pt].Y, typeDirection.tilt1, false, false,
                typeVertex.corner);
            vertexCounter++;

            // Now we have to wrangle the midpoint.

            int next = (pt + 1) % sourcePoly.Length; // wrap to polygon length

            // Find the normal for the edge to the next point.

            double dx = sourcePoly[next].X - sourcePoly[pt].X;
            double dy = sourcePoly[next].Y - sourcePoly[pt].Y;

            // Set up our midpoint for convenience.
            GeoLibPointF midPt = new(sourcePoly[pt].X + dx / 2.0f, sourcePoly[pt].Y + dy / 2.0f);

            // The normal, to match convention in the distance calculation is assessed from this point to the next point.

            // Get average angle for this vertex based on angles from line segments.
            // http://stackoverflow.com/questions/1243614/how-do-i-calculate-the-normal-vector-of-a-line-segmen
            GeoLibPointF normalPt = new(-dy, dx);

            // Vertical edge has a normal with an X value non-zero and Y value ~0.
            // treating a 0.01 difference as being ~0
            bool vertical = Math.Abs(normalPt.X) > 0.01;

            // Assess the normal to establish direction
            if (vertical)
            {
                // left facing vertical edge has normal with negative X value.
                left = normalPt.X < 0;
            }
            else
            {
                // down facing horizontal edge has normal with negative Y value.
                up = !(normalPt.Y < 0);
            }

            if (!vertical)
            {
                if (up)
                {
                    Vertex[vertexCounter] = new MyVertex(midPt.X, midPt.Y, typeDirection.up1, vertical, false,
                        typeVertex.center);
                    if (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.T ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TB ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TL ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TBL ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TBR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TLR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.all)
                    {
                        tips[vertexCounter] = true;
                    }
                }
                else
                {
                    Vertex[vertexCounter] = new MyVertex(midPt.X, midPt.Y, typeDirection.down1, vertical, false,
                        typeVertex.center);
                    if (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.B ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TB ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.BL ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.BR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TBL ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TBR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.BLR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.all)
                    {
                        tips[vertexCounter] = true;
                    }
                }
            }
            else
            {
                if (left)
                {
                    Vertex[vertexCounter] = new MyVertex(midPt.X, midPt.Y, typeDirection.left1, vertical, false,
                        typeVertex.center);
                    if (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.L ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.LR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.BL ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TL ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TBL ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.BLR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TLR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.all)
                    {
                        tips[vertexCounter] = true;
                    }
                }
                else
                {
                    Vertex[vertexCounter] = new MyVertex(midPt.X, midPt.Y, typeDirection.right1, vertical, false,
                        typeVertex.center);
                    if (layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.R ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.LR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.BR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TBR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.BLR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.TLR ||
                        layerSettings.getInt(ShapeSettings.properties_i.subShapeTipLocIndex) ==
                        (int)ShapeSettings.tipLocations.all)
                    {
                        tips[vertexCounter] = true;
                    }
                }
            }

            vertexCounter++;
        }

        // Reprocess our corners for inner/outer rounding based on horFace/verFace directions
#if !SHAPELIBSINGLETHREADED
        Parallel.For(0, roundCount, pt =>
#else
            for (int pt = 0; pt < roundCount; pt++)
#endif
            {
                // Only certain changes in direction correspond to an outer vertex, for a clockwise ordered series of points.
                bool outerVertex = pt == 0 || pt == round1.Length - 1 ||
                                   round1[pt].verFace < round1[pt].horFace &&
                                   Vertex[round1[pt].verFace].direction == typeDirection.left1 &&
                                   Vertex[round1[pt].horFace].direction == typeDirection.up1 ||
                                   round1[pt].verFace > round1[pt].horFace &&
                                   Vertex[round1[pt].horFace].direction == typeDirection.up1 &&
                                   Vertex[round1[pt].verFace].direction == typeDirection.right1 ||
                                   round1[pt].verFace < round1[pt].horFace &&
                                   Vertex[round1[pt].verFace].direction == typeDirection.right1 &&
                                   Vertex[round1[pt].horFace].direction == typeDirection.down1 ||
                                   round1[pt].verFace > round1[pt].horFace &&
                                   Vertex[round1[pt].horFace].direction == typeDirection.down1 &&
                                   Vertex[round1[pt].verFace].direction == typeDirection.left1;

                if (outerVertex)
                {
                    round1[pt].direction = typeRound.exter;
                    round1[pt].MaxRadius =
                        Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.oCR));
                }
                else
                {
                    round1[pt].direction = typeRound.inner;
                    round1[pt].MaxRadius =
                        Convert.ToDouble(layerSettings.getDecimal(ShapeSettings.properties_decimal.iCR));
                }

                Vertex[round1[pt].index].inner = !outerVertex;
            }
#if !SHAPELIBSINGLETHREADED
        );
#endif
    }

    public void computeTips(double gTipBias, double hTipBias, double hTipBiasType, double hTipBiasNegVar,
        double hTipBiasPosVar, double vTipBias, double vTipBiasType, double vTipBiasNegVar, double vTipBiasPosVar)
    {
        // Wrangle the tips.
        for (int cp = 0;
             cp < Vertex.Length - 1;
             cp++) // We don't drive the last point directly - we'll close our shape.
        {
            if (!tips[cp])
            {
                continue;
            }

            // Note that these are reversed due to top-left origin of drawing!
            // A positive bias shrinks the shape in this code.....
            // Values below are correlated in the MCControl system for simulations. In preview mode, these are all zero.
            if (Vertex[cp].direction == typeDirection.down1 && Vertex[cp].yBiasApplied == false)
            {
                Vertex[cp].Y -= vTipBias;
                Vertex[cp].Y -= Convert.ToDouble(gTipBias);
                if (vTipBiasType < 0.5) // Need to use our negative variation value
                {
                    Vertex[cp].Y -= vTipBiasNegVar;
                }
                else
                {
                    Vertex[cp].Y -= vTipBiasPosVar;
                }

                Vertex[cp].yBiasApplied = true;
            }

            if (Vertex[cp].direction == typeDirection.up1 && Vertex[cp].yBiasApplied == false)
            {
                Vertex[cp].Y += vTipBias;
                Vertex[cp].Y += gTipBias;
                if (vTipBiasType < 0.5) // Need to use our negative variation value
                {
                    Vertex[cp].Y += vTipBiasNegVar;
                }
                else
                {
                    Vertex[cp].Y += vTipBiasPosVar;
                }

                Vertex[cp].yBiasApplied = true;
            }

            if (Vertex[cp].direction == typeDirection.left1 && Vertex[cp].xBiasApplied == false)
            {
                Vertex[cp].X -= hTipBias;
                Vertex[cp].X -= gTipBias;
                if (hTipBiasType < 0.5) // Need to use our negative variation value
                {
                    Vertex[cp].X -= hTipBiasNegVar;
                }
                else
                {
                    Vertex[cp].X -= hTipBiasPosVar;
                }

                Vertex[cp].xBiasApplied = true;
            }

            if (Vertex[cp].direction != typeDirection.right1 || Vertex[cp].xBiasApplied)
            {
                continue;
            }

            Vertex[cp].X += hTipBias;
            Vertex[cp].X += gTipBias;
            if (hTipBiasType < 0.5) // Need to use our negative variation value
            {
                Vertex[cp].X += hTipBiasNegVar;
            }
            else
            {
                Vertex[cp].X += hTipBiasPosVar;
            }

            Vertex[cp].xBiasApplied = true;
        }
    }

    public void computeBias(double gSideBias)
    {
        // Global bias for anything that isn't a tip.
        for (int cp = 0;
             cp < Vertex.Length - 1;
             cp++) // We don't drive the last point directly - we'll close our shape.
        {
            if (!Vertex[cp].xBiasApplied && !tips[cp])
            {
                switch (Vertex[cp].direction)
                {
                    case typeDirection.left1:
                        Vertex[cp].X -= gSideBias;
                        break;
                    case typeDirection.right1:
                        Vertex[cp].X += gSideBias;
                        break;
                }
            }

            if (!Vertex[cp].yBiasApplied && !tips[cp])
            {
                switch (Vertex[cp].direction)
                {
                    case typeDirection.up1:
                        Vertex[cp].Y += gSideBias;
                        break;
                    case typeDirection.down1:
                        Vertex[cp].Y -= gSideBias;
                        break;
                }
            }
        }
    }

    public void biasCorners()
    {
        if (shapeIndex == (int)shapeNames_all.complex)
        {
            return;
        }
        // Iterate the corners to apply the bias from the edges.
        foreach (MyRound t in round1)
        {
            Vertex[t.index].X = Vertex[t.verFace].X;
            Vertex[t.index].Y = Vertex[t.horFace].Y;
        }

        Vertex[^1] = Vertex[0]; // close the shape.
        round1[^1] = round1[0];
    }

    public void edgeMidpoints(int edgeSlide, double eTension)
    {
        // Set the midpoints of the edges to the average between the two corners
        for (int corner = 0; corner < round1.Length; corner++)
        {
            double previousEdgeLength;

            if (corner == 0)
            {
                previousEdgeLength = Math.Abs(
                    GeoWrangler.distanceBetweenPoints(
                        new GeoLibPointF(Vertex[round1[corner].index].X, Vertex[round1[corner].index].Y),
                        new GeoLibPointF(Vertex[round1[^1].index].X, Vertex[round1[^1].index].Y))
                );
            }
            else
            {
                previousEdgeLength = Math.Abs(
                    GeoWrangler.distanceBetweenPoints(
                        new GeoLibPointF(Vertex[round1[corner].index].X, Vertex[round1[corner].index].Y),
                        new GeoLibPointF(Vertex[round1[corner - 1].index].X, Vertex[round1[corner - 1].index].Y))
                );
            }

            // Wrap around if we exceed the length
            double nextEdgeLength = Math.Abs(
                GeoWrangler.distanceBetweenPoints(
                    new GeoLibPointF(Vertex[round1[(corner + 1) % (round1.Length - 1)].index].X,
                        Vertex[round1[(corner + 1) % (round1.Length - 1)].index].Y),
                    new GeoLibPointF(Vertex[round1[(corner + 2) % (round1.Length - 1)].index].X,
                        Vertex[round1[(corner + 2) % (round1.Length - 1)].index].Y))
            );

            double currentEdgeLength = Math.Abs(
                GeoWrangler.distanceBetweenPoints(
                    new GeoLibPointF(Vertex[round1[corner].index].X, Vertex[round1[corner].index].Y),
                    new GeoLibPointF(Vertex[round1[(corner + 1) % (round1.Length - 1)].index].X,
                        Vertex[round1[(corner + 1) % (round1.Length - 1)].index].Y))
            );

            double offset = 0.5f * currentEdgeLength;
            bool
                reverseSlide =
                    true; // used in the linear mode to handle reversed case (where ratio is > 1), and the no-slide case.)

            if (edgeSlide == 1 && previousEdgeLength > 0 && nextEdgeLength > 0)
            {
                // Now we need to figure out the weighting.
                double ratio = Math.Abs(nextEdgeLength / previousEdgeLength);
                bool doLinearSlide = false;

                if (ratio < 1)
                {
                    reverseSlide = false;
                    if (ratio < 1E-2)
                    {
                        ratio = 1E-2; // clamp
                    }

                    ratio = 1 / ratio; // normalize into our expected range
                }

                if (doLinearSlide)
                {
                    // Linear
                    offset = Math.Pow(currentEdgeLength / 2.0f, 1.0f / ratio); // ratio * currentEdgeLength / 2.0f;
                }
                else
                {
                    // Sigmoid function to try and provide some upper and lower resistance to the slide.
                    // center is to force 0.5 value of the scaling factor for a ratio of 1
                    // tension controls the shape of the curve, and thus the sensitivity of the response..
                    double center = 1.0f;
                    offset = currentEdgeLength * (1 / (1 + Math.Exp(-eTension * (center - ratio))));
                }
            }

            if (corner % 2 == 0)
            {
                // Get our associated vertical edge Y position
                double yPoint1;
                double yPoint2 = Vertex[round1[(corner + 1) % (round1.Length - 1)].horFace].Y;
                switch (corner)
                {
                    case 0:
                        // Need to wrap around for bias look-up
                        yPoint1 = Vertex[round1[^1].horFace].Y;
                        break;
                    default:
                        yPoint1 = Vertex[round1[corner].horFace].Y;
                        break;
                }

                if (yPoint1 < yPoint2)
                {
                    if (reverseSlide)
                    {
                        Vertex[round1[corner].verFace].Y = yPoint2 - offset;
                    }
                    else
                    {
                        Vertex[round1[corner].verFace].Y = yPoint1 + offset;
                    }
                }
                else
                {
                    if (reverseSlide)
                    {
                        Vertex[round1[corner].verFace].Y = yPoint2 + offset;
                    }
                    else
                    {
                        Vertex[round1[corner].verFace].Y = yPoint1 - offset;
                    }
                }
            }
            else
            {
                // Tweak horizontal edge
                double xPoint1 = Vertex[round1[corner].verFace].X;
                double xPoint2 = Vertex[round1[(corner + 1) % (round1.Length - 1)].verFace].X;

                if (xPoint1 < xPoint2)
                {
                    if (reverseSlide)
                    {
                        Vertex[round1[corner].horFace].X = xPoint2 - offset;
                    }
                    else
                    {
                        Vertex[round1[corner].horFace].X = xPoint1 + offset;
                    }
                }
                else
                {
                    if (reverseSlide)
                    {
                        Vertex[round1[corner].horFace].X = xPoint2 + offset;
                    }
                    else
                    {
                        Vertex[round1[corner].horFace].X = xPoint1 - offset;
                    }
                }
            }
        }
    }

    public List<GeoLibPointF> processCorners(bool previewMode, bool cornerCheck, bool ignoreCV, double s0HO,
        double s0VO, double iCR, double iCV, double iCVariation, bool iCPA, double oCR, double oCV, double oCVariation,
        bool oCPA, int cornerSegments, int optimizeCorners, double resolution, int scaleFactorForOperation)
    {
        Fragmenter fragment = new Fragmenter(resolution, scaleFactorForOperation);
        List<GeoLibPointF> mcPoints = new();
        List<GeoLibPointF>
            mcHorEdgePoints = new(); // corner coordinates list, used as a temporary container for each iteration
        List<List<GeoLibPointF>>
            mcHorEdgePointsList =
                new(); // Hold our lists of doubles for each corner in the shape, in order. We cast these to Int in the mcPoints list.
        List<List<GeoLibPointF>>
            mcVerEdgePointsList =
                new(); // Hold our lists of doubles for each edge in the shape, in order. We cast these to Int in the mcPoints list.

        for (int round = 0; round < round1.Length - 1; round++)
        {
            // Derive our basic coordinates for the three vertices on the edge.
            double start_x = Vertex[round1[round].index].X;
            double start_y = Vertex[round1[round].index].Y;
            double currentHorEdge_mid_x = Vertex[round1[round].horFace].X;
            double currentVerEdge_mid_y = Vertex[round1[round].verFace].Y;
            double end_x = Vertex[round1[round + 1].index].X;
            double end_y = Vertex[round1[round + 1].index].Y;
            double nextVerEdge_mid_y = Vertex[round1[round + 1].verFace].Y;

            switch (Math.Abs(start_y - end_y))
            {
                // Test whether we have a vertical edge or not. We only process horizontal edges to avoid doubling up
                case < double.Epsilon:
                {
                    double mcPX;
                    double mcPY = 0.0f;
                    // Establish corner rounding sign at start and end points of edge. Default is to move outwards (inner CRR)
                    bool startInnerRounding = true;
                    bool endInnerRounding = true;
                    if (round1[round].direction == typeRound.exter)
                    {
                        startInnerRounding = false;
                    }

                    if (round1[round + 1].direction == typeRound.exter)
                    {
                        endInnerRounding = false;
                    }

                    // Now sort out the shift based on face orientation.
                    bool horFaceUp = true;
                    switch (Vertex[round1[round].horFace].direction)
                    {
                        case typeDirection.up1:
                            break;
                        case typeDirection.down1:
                            horFaceUp = false;
                            break;
                    }

                    bool verFaceLeft = true;
                    switch (Vertex[round1[round].verFace].direction)
                    {
                        case typeDirection.left1:
                            break;
                        case typeDirection.right1:
                            verFaceLeft = false;
                            break;
                    }

                    // Segment 1

                    // Clamp radius in each direction, if needed, to available distance
                    double hRadius = round1[round].MaxRadius;
                    double x_Distance = Math.Sqrt(Utils.myPow(currentHorEdge_mid_x - start_x, 2));
                    double vRadius = round1[round].MaxRadius;
                    double y_Distance = Math.Sqrt(Utils.myPow(currentVerEdge_mid_y - start_y, 2));

                    // Add our random variation based on rounding type :
                    bool paSearchSetsCornerRoundingForThisCorner = false;

                    if (ignoreCV)
                    {
                        // Are we on a corner that has a PA-defined rounding value?
                        switch (startInnerRounding)
                        {
                            case true:
                                paSearchSetsCornerRoundingForThisCorner = iCPA;
                                break;
                            default:
                                paSearchSetsCornerRoundingForThisCorner = oCPA;
                                break;
                        }
                    }

                    // PA search works by setting rounding value directly in the settings, no variation needs to be added.
                    if (paSearchSetsCornerRoundingForThisCorner)
                    {
                        if (startInnerRounding)
                        {
                            hRadius = iCVariation * iCR;
                            vRadius = iCVariation * iCR;
                        }
                        else
                        {
                            hRadius = oCVariation * oCR;
                            vRadius = oCVariation * oCR;
                        }
                    }
                    else
                    {
                        if (!previewMode)
                        {
                            if (startInnerRounding)
                            {
                                hRadius += iCVariation * iCV;
                                vRadius += iCVariation * iCV;
                            }
                            else
                            {
                                hRadius += oCVariation * oCV;
                                vRadius += oCVariation * oCV;
                            }
                        }
                    }

                    if (hRadius > x_Distance)
                    {
                        hRadius = x_Distance;
                    }

                    if (vRadius > y_Distance)
                    {
                        vRadius = y_Distance;
                    }

                    // Clamp for negative radius values that would make no sense
                    if (hRadius < 0)
                    {
                        hRadius = 0;
                    }

                    if (vRadius < 0)
                    {
                        vRadius = 0;
                    }

                    double angleIncrement = 90.0f / cornerSegments;

                    // Sweep our corner.
                    double angle = 0.0f;
                    while (angle <= 90.0f)
                    {
                        // Set start condition
                        mcPX = start_x; // X position for new point.
                        mcPY = start_y; // this will hold our Y position for the new point.

                        // Remove full contribution from rounding.
                        if (verFaceLeft)
                        {
                            if (startInnerRounding)
                            {
                                mcPX -= hRadius;
                            }
                            else
                            {
                                mcPX += hRadius;
                            }
                        }
                        else
                        {
                            if (startInnerRounding)
                            {
                                mcPX += hRadius;
                            }
                            else
                            {
                                mcPX -= hRadius;
                            }
                        }

                        if (horFaceUp)
                        {
                            if (startInnerRounding)
                            {
                                mcPY += vRadius;
                            }
                            else
                            {
                                mcPY -= vRadius;
                            }
                        }
                        else
                        {
                            if (startInnerRounding)
                            {
                                mcPY -= vRadius;
                            }
                            else
                            {
                                mcPY += vRadius;
                            }
                        }

                        // Now process corner, adding back contribution from rounding
                        if (verFaceLeft)
                        {
                            if (startInnerRounding)
                            {
                                mcPX += hRadius * Math.Cos(Utils.toRadians(angle));
                            }
                            else
                            {
                                mcPX -= hRadius * Math.Cos(Utils.toRadians(angle));
                            }
                        }
                        else
                        {
                            if (startInnerRounding)
                            {
                                mcPX -= hRadius * Math.Cos(Utils.toRadians(angle));
                            }
                            else
                            {
                                mcPX += hRadius * Math.Cos(Utils.toRadians(angle));
                            }
                        }

                        if (horFaceUp)
                        {
                            if (startInnerRounding)
                            {
                                mcPY -= vRadius * Math.Sin(Utils.toRadians(angle));
                            }
                            else
                            {
                                mcPY += vRadius * Math.Sin(Utils.toRadians(angle));
                            }
                        }
                        else
                        {
                            if (startInnerRounding)
                            {
                                mcPY += vRadius * Math.Sin(Utils.toRadians(angle));
                            }
                            else
                            {
                                mcPY -= vRadius * Math.Sin(Utils.toRadians(angle));
                            }
                        }

                        GeoLibPointF cPt = new(mcPX, mcPY);
                        if (angle == 0 || Math.Abs(angle - 90) < double.Epsilon || optimizeCorners == 0 ||
                            optimizeCorners == 1 &&
                            Math.Abs(
                                GeoWrangler.distanceBetweenPoints(mcHorEdgePoints[^1], cPt)
                            )
                            > resolution
                           )
                        {
                            mcHorEdgePoints.Add(cPt);
                        }

                        angle += angleIncrement;
                    }

                    // OK. We now need to add points along the edge based on the simulation settings resolution.
                    // We need to add points from here to just before the midpoint

                    double bridgeX = mcHorEdgePoints[^1].X;

                    // Fragmenter returns first and last points in the point array.
                    GeoLibPointF[] fragments = fragment.fragmentPath(new[]
                        { new GeoLibPointF(bridgeX, mcPY), new GeoLibPointF(currentHorEdge_mid_x, mcPY) });

                    for (int i = 1; i < fragments.Length - 1; i++)
                    {
                        mcHorEdgePoints.Add(fragments[i]);
                    }

                    // Add our midpoint.
                    mcHorEdgePoints.Add(new GeoLibPointF(currentHorEdge_mid_x, mcPY));

                    // Segment 2, plus bridging on first pass through.

                    bool
                        firstPass =
                            true; // With this set, we bridge from midpoint to our first point in the first pass through
                    // segment 2 of the edge.
                    verFaceLeft = true;
                    switch (Vertex[round1[round + 1].verFace].direction)
                    {
                        case typeDirection.left1:
                            break;
                        case typeDirection.right1:
                            verFaceLeft = false;
                            break;
                    }

                    // Clamp radius to available distance, if needed.
                    hRadius = round1[round + 1].MaxRadius;
                    x_Distance = Math.Sqrt(Utils.myPow(currentHorEdge_mid_x - end_x, 2));
                    vRadius = round1[round + 1].MaxRadius;
                    y_Distance = Math.Sqrt(Utils.myPow(nextVerEdge_mid_y - end_y, 2));

                    // Add our random variation based on rounding type :

                    paSearchSetsCornerRoundingForThisCorner = false;
                    if (ignoreCV)
                    {
                        switch (startInnerRounding)
                        {
                            case true:
                                paSearchSetsCornerRoundingForThisCorner = iCPA;
                                break;
                            default:
                                paSearchSetsCornerRoundingForThisCorner = oCPA;
                                break;
                        }
                    }

                    if (paSearchSetsCornerRoundingForThisCorner)
                    {
                        if (startInnerRounding)
                        {
                            hRadius = iCVariation * iCR;
                            vRadius = iCVariation * iCR;
                        }
                        else
                        {
                            hRadius = oCVariation * oCR;
                            vRadius = oCVariation * oCR;
                        }
                    }
                    else
                    {
                        // Add our random variation based on rounding type :
                        if (!previewMode)
                        {
                            if (endInnerRounding)
                            {
                                hRadius += iCVariation * iCV;
                                vRadius += iCVariation * iCV;
                            }
                            else
                            {
                                hRadius += oCVariation * oCV;
                                vRadius += oCVariation * oCV;
                            }
                        }
                    }

                    if (hRadius > x_Distance)
                    {
                        hRadius = x_Distance;
                    }

                    if (vRadius > y_Distance)
                    {
                        vRadius = y_Distance;
                    }

                    // Clamp for negative radius values that would make no sense
                    if (hRadius < 0)
                    {
                        hRadius = 0;
                    }

                    if (vRadius < 0)
                    {
                        vRadius = 0;
                    }

                    // Sweep our end corner. We need to run the sweep in the opposite direction.
                    angle = 90.0f;
                    while (angle >= 0.0f)
                    {
                        // Set start conditions
                        mcPX = end_x;
                        mcPY = end_y;

                        // Remove full extent of rounding in each direction, based on face orientation
                        if (verFaceLeft)
                        {
                            if (endInnerRounding)
                            {
                                mcPX -= hRadius;
                            }
                            else
                            {
                                mcPX += hRadius;
                            }
                        }
                        else
                        {
                            if (endInnerRounding)
                            {
                                mcPX += hRadius;
                            }
                            else
                            {
                                mcPX -= hRadius;
                            }
                        }

                        if (horFaceUp)
                        {
                            if (endInnerRounding)
                            {
                                mcPY += vRadius;
                            }
                            else
                            {
                                mcPY -= vRadius;
                            }
                        }
                        else
                        {
                            if (endInnerRounding)
                            {
                                mcPY -= vRadius;
                            }
                            else
                            {
                                mcPY += vRadius;
                            }
                        }

                        // Process corners, adding back the contribution from the rounding based on the angle
                        if (verFaceLeft)
                        {
                            if (endInnerRounding)
                            {
                                mcPX += hRadius * Math.Cos(Utils.toRadians(angle));
                            }
                            else
                            {
                                mcPX -= hRadius * Math.Cos(Utils.toRadians(angle));
                            }
                        }
                        else
                        {
                            if (endInnerRounding)
                            {
                                mcPX -= hRadius * Math.Cos(Utils.toRadians(angle));
                            }
                            else
                            {
                                mcPX += hRadius * Math.Cos(Utils.toRadians(angle));
                            }
                        }

                        if (horFaceUp)
                        {
                            if (endInnerRounding)
                            {
                                mcPY -= vRadius * Math.Sin(Utils.toRadians(angle));
                            }
                            else
                            {
                                mcPY += vRadius * Math.Sin(Utils.toRadians(angle));
                            }
                        }
                        else
                        {
                            if (endInnerRounding)
                            {
                                mcPY += vRadius * Math.Sin(Utils.toRadians(angle));
                            }
                            else
                            {
                                mcPY -= vRadius * Math.Sin(Utils.toRadians(angle));
                            }
                        }

                        // If this is the first pass, we need to add points to the start of the rounding, from the midpoint.
                        if (firstPass)
                        {
                            bridgeX = currentHorEdge_mid_x;

                            // Fragmenter returns first and last points in the point array.
                            fragments = fragment.fragmentPath(new[]
                                { new GeoLibPointF(bridgeX, mcPY), new GeoLibPointF(mcPX, mcPY) });

                            for (int i = 1; i < fragments.Length - 1; i++)
                            {
                                mcHorEdgePoints.Add(fragments[i]);
                            }

                            firstPass = false;
                        }

                        GeoLibPointF cPt = new(mcPX, mcPY);
                        if (angle == 0 || Math.Abs(angle - 90) < double.Epsilon || optimizeCorners == 0 ||
                            optimizeCorners == 1 &&
                            Math.Abs(
                                GeoWrangler.distanceBetweenPoints(mcHorEdgePoints[^1], cPt)
                            )
                            > resolution
                           )
                        {
                            mcHorEdgePoints.Add(cPt);
                        }

                        angle -= angleIncrement;
                    }

                    mcHorEdgePointsList.Add(mcHorEdgePoints.ToList()); // make a deep copy of the points.
                    mcHorEdgePoints.Clear(); // clear our list of points to use on the next pass.
                    break;
                }
            }
        }

        if (cornerCheck)
        {
            mcPoints.Clear();
            foreach (List<GeoLibPointF> t in mcHorEdgePointsList)
            {
                foreach (GeoLibPointF t1 in t)
                {
                    mcPoints.Add(new GeoLibPointF(t1.X, t1.Y));
                }
            }

            return GeoWrangler.close(mcPoints);
        }

        // Now we have our corners, let's process the vertical edges. We need the corners in order to get our start/end on each vertical edge.
        for (int edge = 0; edge < mcHorEdgePointsList.Count; edge++)
        {
            // Get our start and end Y positions for our vertical edge.
            List<GeoLibPointF> startHorEdgePointList = mcHorEdgePointsList[edge];
            int endHorEdgePointListIndex;
            if (edge == 0)
            {
                endHorEdgePointListIndex = mcHorEdgePointsList.Count - 1; // need to wrap around.
            }
            else
            {
                endHorEdgePointListIndex = edge - 1;
            }

            List<GeoLibPointF> endHorEdgePointList = mcHorEdgePointsList[endHorEdgePointListIndex];
            double vert_x = endHorEdgePointList[^1].X;
            double startPoint_y = endHorEdgePointList[^1].Y;
            double endPoint_y = startHorEdgePointList[0].Y;

            // We get the start and end points here.
            List<GeoLibPointF> fragments = fragment.fragmentPath(new List<GeoLibPointF>
                { new GeoLibPointF(vert_x, startPoint_y), new GeoLibPointF(vert_x, endPoint_y) });
            mcVerEdgePointsList.Add(fragments);
        }

        // OK. We have our corners and edges. We need to walk them now. We'll apply the subshape 1 offset at the same time.
        for (int section = 0; section < mcVerEdgePointsList.Count; section++)
        {
            for (int point = 0; point < mcVerEdgePointsList[section].Count; point++)
            {
                double x = mcVerEdgePointsList[section][point].X + s0HO;
                double y = mcVerEdgePointsList[section][point].Y + s0VO;
                mcPoints.Add(new GeoLibPointF(x, y));
            }

            // Corner next.
            // Start and end points match those in the vertical edges, so we avoid them to eliminate duplicates.
            for (int point = 1; point < mcHorEdgePointsList[section].Count - 1; point++)
            {
                double x = mcHorEdgePointsList[section][point].X + s0HO;
                double y = mcHorEdgePointsList[section][point].Y + s0VO;
                mcPoints.Add(new GeoLibPointF(x, y));
            }
        }

        return mcPoints;

    }

    public List<GeoLibPointF> rotateShape(List<GeoLibPointF> input, ShapeSettings shapeSettings, double rotationVar,
        double rotationDirection, GeoLibPointF? pivot = null)
    {
        double rotationAngle = Convert.ToDouble(shapeSettings.getDecimal(ShapeSettings.properties_decimal.rot));
        if (rotationDirection <= 0.5)
        {
            rotationAngle -= rotationVar;
        }
        else
        {
            rotationAngle += rotationVar;
        }

        if (rotationAngle != 0 ||
            (shapeSettings.getInt(ShapeSettings.properties_i.flipH) == 1 ||
             shapeSettings.getInt(ShapeSettings.properties_i.flipV) == 1) &&
            (shapeSettings.getInt(ShapeSettings.properties_i.alignX) == 1 ||
             shapeSettings.getInt(ShapeSettings.properties_i.alignY) == 1))
        {
            // Get our bounding box.
            BoundingBox bb = new(input);

            switch (pivot)
            {
                case null:
                    pivot = new GeoLibPointF(bb.getMidPoint());
                    break;
            }

            // OK. Let's try some rotation and wobble.
            // Temporary separate container for our rotated points, just for now.
            return GeoWrangler.Rotate(pivot, input, rotationAngle);
        }

        return input;
    }
}
