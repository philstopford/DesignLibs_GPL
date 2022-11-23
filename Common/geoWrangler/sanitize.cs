using LibTessDotNet.Double;
using System;
using System.Collections.Generic;
using Clipper2Lib;
using utility;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PathsD fromSoup(PathsD source)
    {
        return pFromSoup(source, WindingRule.NonZero);
    }

    private static PathsD pFromSoup(PathsD source, WindingRule wr)
    {
        PathsD outers = new();
        PathsD cutters = new();

        foreach (PathD t in source)
        {
            if (Clipper.IsPositive(t) == Clipper.IsPositive(source[0]))
            {
                outers.Add(new (t));
            }
            else
            {
                cutters.Add(new (t));
            }
        }

        // Set up our contours. Since Clipper sets up the subject as the first item, we'll make that clockwise and force the rest to counterclockwise.	
        Tess tess = new();

        foreach (PathD t in outers)
        {
            // Skip tiny fragments. The tessellator has trouble with them.	
            /*
            if (Math.Abs(Clipper.Area(outers[poly])) < 1E-16)
            {
                continue;
            }
            */
            ContourVertex[] contour = new ContourVertex[t.Count];
            for (int i = 0; i < contour.Length; i++)
            {
                contour[i].Position = new Vec3 { X = t[i].x, Y = t[i].y, Z = 0 };
            }

            tess.AddContour(contour);
        }

        foreach (PathD t in cutters)
        {
            // Skip tiny fragments. The tessellator has trouble with them.
            /*
            if (Math.Abs(Clipper.Area(cutters[poly])) < 1E-16)
            {
                continue;
            }
            */
            ContourVertex[] contour = new ContourVertex[t.Count];
            for (int i = 0; i < contour.Length; i++)
            {
                contour[i].Position = new Vec3 { X = t[i].x, Y = t[i].y, Z = 0 };
            }

            tess.AddContour(contour, ContourOrientation.CounterClockwise);
        }

        try
        {
            // Triangulate. This gives us a triangle soup, which may not be entirely helpful for the case where we're punching multiple holes into a mesh. For now, this works, but the limitation will need to be reviewd.	
            // It might be that a PolyTree is needed as the input to the keyhole for these complex cases.... or clockwise paths may need to be evaluated piecewise using all counterclockwise paths for the tessellation.	
            const int polysize = 3;
            tess.Tessellate(wr, ElementType.Polygons, polysize);

            // Iterate triangles and create output geometry. We'll use clipper to simplify the output geometry.	
            ClipperD c = new() {PreserveCollinear = true};
            PathsD retPaths = new();

            PathsD cPaths = new();
            PathsD aPaths = new();

            for (int i = 0; i < tess.ElementCount; i++)
            {
                PathD trianglePath = new();
                for (int p = 0; p < polysize; p++)
                {
                    PointD tmpPt = new(tess.Vertices[tess.Elements[i * polysize + p]].Position.X, tess.Vertices[tess.Elements[i * polysize + p]].Position.Y);
                    trianglePath.Add(tmpPt);
                }

                if (Clipper.IsPositive(trianglePath))
                {
                    cPaths.Add(trianglePath);
                }
                else
                {
                    aPaths.Add(trianglePath);
                }
            }

            // Add paths to the clipper.	
            c.AddSubject(cPaths);
            c.AddClip(aPaths);

            c.Execute(ClipType.Union, FillRule.NonZero, retPaths);

            retPaths = pReorderXY(retPaths);

            retPaths = pClose(retPaths);

            return retPaths;
        }
        catch (Exception)
        {
            return new(source);
        }
    }

    public static PathsD clean_and_flatten(PathsD source, double customSizing = 0, double extension = 0)
    {
        return pClean_and_flatten(source, customSizing, extension);
    }

    private static PathsD pClean_and_flatten(PathsD source, double customSizing = 0, double extension = 0)
    {
        ClipperD c = new();
        c.AddSubject(source);
        PathsD solution = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, solution);

        solution = pReorderXY(solution);

        PathsD keyHoled = pMakeKeyHole(solution, reverseEval:false, biDirectionalEval:true, customSizing: customSizing, extension: extension);

        return pClockwiseAndReorderXY(keyHoled);
    }

    public static PathsD removeDuplicatePaths(PathsD source)
    {
        return pRemoveDuplicatePaths(source);
    }
    
    private static PathsD pRemoveDuplicatePaths(PathsD source)
    {
        PathsD ret = new();
        List<string> polyHashCodes = new();

        foreach (PathD p in source)
        {
            string polyHash = Utils.GetMD5Hash(p);
            if (polyHashCodes.IndexOf(polyHash) != -1)
            {
                continue;
            }
            polyHashCodes.Add(polyHash);
            ret.Add(new(p));
        }

        return ret;
    }
    
}