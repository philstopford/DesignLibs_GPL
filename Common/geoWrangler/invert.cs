using ClipperLib;
using System;
using System.Collections.Generic;
using System.Linq;

namespace geoWrangler
{
    using Path = List<IntPoint>;
    using Paths = List<List<IntPoint>>;

    public static partial class GeoWrangler
    {
        public static Paths invertTone(Paths sourcePaths, bool useTriangulation = false)
        {
            return pInvertTone(sourcePaths, useTriangulation);
        }

        static Paths pInvertTone(Paths sourcePaths, bool useTriangulation)
        {
            if (sourcePaths.Count() == 1)
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
                IntRect firstLayerBounds = Clipper.GetBounds(sourcePaths);
                IntPoint delta_to_midPoint = new IntPoint(((firstLayerBounds.right - firstLayerBounds.left) / 2), (firstLayerBounds.top - firstLayerBounds.bottom) / 2);
                Path firstLayerBP = new Path();
                firstLayerBP.Add(new IntPoint(-Int32.MaxValue, -Int32.MaxValue));
                firstLayerBP.Add(new IntPoint(-Int32.MaxValue, Int32.MaxValue));
                firstLayerBP.Add(new IntPoint(Int32.MaxValue, Int32.MaxValue));
                firstLayerBP.Add(new IntPoint(Int32.MaxValue, -Int32.MaxValue));
                firstLayerBP.Add(new IntPoint(-Int32.MaxValue, -Int32.MaxValue));

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
