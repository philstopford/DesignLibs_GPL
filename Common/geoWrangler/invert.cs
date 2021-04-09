using ClipperLib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;

namespace geoWrangler
{
    using Path = List<IntPoint>;
    using Paths = List<List<IntPoint>>;

    public static partial class GeoWrangler
    {
        // Bounds will force the negation to only work against the extents of the shape. Useful for capturing islands in negative tone.
        public static List<GeoLibPointF[]> invertTone(GeoLibPointF[] source, long scaleFactor, bool useTriangulation = false, bool useBounds = false)
        {
            return pPointFsFromPaths(pInvertTone(pPathsFromPointFs(new List<GeoLibPointF[]> { source }, scaleFactor), useTriangulation, useBounds), scaleFactor);
        }

        public static List<GeoLibPointF[]> invertTone(List<GeoLibPointF[]> source, long scaleFactor, bool useTriangulation = false, bool useBounds = false)
        {
            return pPointFsFromPaths(pInvertTone(pPathsFromPointFs(source, scaleFactor), useTriangulation, useBounds), scaleFactor);
        }

        public static Paths invertTone(Path sourcePath, bool useTriangulation = false, bool useBounds = false)
        {
            Paths t = new Paths();
            t.Add(sourcePath);
            return invertTone(t, useTriangulation, useBounds);
        }

        public static Paths invertTone(Paths sourcePaths, bool useTriangulation = false, bool useBounds = false)
        {
            return pInvertTone(sourcePaths, useTriangulation, useBounds);
        }

        static Paths pInvertTone(Paths sourcePaths, bool useTriangulation, bool useBounds)
        {
            if ((sourcePaths.Count() == 1) && !useBounds)
            {
                sourcePaths[0].Add(new IntPoint(-Int32.MaxValue, -Int32.MaxValue));
                sourcePaths[0].Add(new IntPoint(-Int32.MaxValue, Int32.MaxValue));
                sourcePaths[0].Add(new IntPoint(Int32.MaxValue, Int32.MaxValue));
                sourcePaths[0].Add(new IntPoint(Int32.MaxValue, -Int32.MaxValue));
                sourcePaths[0].Add(new IntPoint(-Int32.MaxValue, -Int32.MaxValue));

                return sourcePaths.ToList();
            }
            else
            {
                IntRect firstLayerBounds = ClipperBase.GetBounds(sourcePaths);
                Path firstLayerBP = new Path();
                if (!useBounds)
                {
                    firstLayerBP.Add(new IntPoint(-Int32.MaxValue, -Int32.MaxValue));
                    firstLayerBP.Add(new IntPoint(-Int32.MaxValue, Int32.MaxValue));
                    firstLayerBP.Add(new IntPoint(Int32.MaxValue, Int32.MaxValue));
                    firstLayerBP.Add(new IntPoint(Int32.MaxValue, -Int32.MaxValue));
                    firstLayerBP.Add(new IntPoint(-Int32.MaxValue, -Int32.MaxValue));
                }
                else
                {
                    IntRect bounds = ClipperBase.GetBounds(sourcePaths);
                    firstLayerBP.Add(new IntPoint(bounds.left, bounds.bottom));
                    firstLayerBP.Add(new IntPoint(bounds.left, bounds.top));
                    firstLayerBP.Add(new IntPoint(bounds.right, bounds.top));
                    firstLayerBP.Add(new IntPoint(bounds.right, bounds.bottom));
                    firstLayerBP.Add(new IntPoint(bounds.left, bounds.bottom));
                }

                Clipper c = new Clipper();
                c.PreserveCollinear = false;

                c.AddPath(firstLayerBP, PolyType.ptSubject, true);
                // Add hole polygons from our paths
                c.AddPaths(sourcePaths, PolyType.ptClip, true);

                Paths cutters = new Paths();
                c.Execute(ClipType.ctDifference, cutters);

                if (!useTriangulation)
                {
                    sourcePaths = cutters;
                    return sourcePaths;
                }
                else
                {
                    return pMakeKeyHole(sourcePaths);
                }
            }
        }
    }
}
