using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Clipper2Lib;
using geoLib;

namespace geoWrangler;

using Path = List<Point64>;
using Paths = List<List<Point64>>;

public static class ProcessOverlaps
{
    public static GeometryResult processOverlaps(List<GeoLibPointF[]> sourceData, List<bool> drawn, double extension, double resolution, Int64 scaleFactorForOperation, double customSizing = 0, bool forceOverride = false, FillRule pft = FillRule.NonZero)
    {
        // Filter drawn, process those, then do not-drawn. This allows for element counts to change.
        List<GeoLibPointF[]> drawnStuff = new();
        List<GeoLibPointF[]> notDrawnStuff = new();
        int sCount = sourceData.Count;
        for (int i = 0; i < sCount; i++)
        {
            if (drawn[i])
            {
                drawnStuff.Add(sourceData[i]);
            }
            else
            {
                notDrawnStuff.Add(sourceData[i]);
            }
        }
        
        List<GeoLibPointF[]> processed_Drawn = processOverlaps_core(drawnStuff, customSizing:customSizing, extension:extension, resolution:resolution, scaleFactorForOperation:scaleFactorForOperation, forceOverride, pft);

        List<GeoLibPointF[]> processed_NotDrawn = processOverlaps_core( notDrawnStuff, customSizing:customSizing, extension: extension, resolution:resolution, scaleFactorForOperation:scaleFactorForOperation, forceOverride, pft);


        GeometryResult ret = new();
        
        int pdCount = processed_Drawn.Count;
        int pndCount = processed_NotDrawn.Count;

        for (int i = 0; i < pdCount; i++)
        {
            ret.geometry.Add(processed_Drawn[i]);
            ret.drawn.Add(true);
        }

        for (int i = 0; i < pndCount; i++)
        {
            ret.geometry.Add(processed_NotDrawn[i]);
            ret.drawn.Add(false);
        }

        return ret;

    }

    private static List<GeoLibPointF[]> processOverlaps_core(List<GeoLibPointF[]> sourceData, double customSizing, double extension, double resolution, Int64 scaleFactorForOperation, bool forceOverride = false, FillRule pft = FillRule.NonZero)
    {
        if (sourceData.Count == 0)
        {
            return sourceData;
        }
        try
        {
            Clipper64 c = new() {PreserveCollinear = true};
            Paths sourcePolyData = GeoWrangler.pathsFromPointFs(sourceData, scaleFactorForOperation);
            Paths mergedPolyData = new();
            
            // Union isn't always robust, so get a bounding box and run an intersection boolean to rationalize the geometry.
            Rect64 bounds = Clipper.GetBounds(sourcePolyData);
            Path bounding = new()
            {
                new Point64(bounds.left, bounds.bottom),
                new Point64(bounds.left, bounds.top),
                new Point64(bounds.right, bounds.top),
                new Point64(bounds.right, bounds.bottom)
            };

            c.AddClip(sourcePolyData);
            c.AddSubject(bounding);
            
            c.Execute(ClipType.Intersection, pft, mergedPolyData);
            mergedPolyData = GeoWrangler.reOrderXY(mergedPolyData);
            mergedPolyData = GeoWrangler.close(mergedPolyData);

            double area = Clipper.Area(mergedPolyData);

            // Avoid the keyhole handling if we don't actually need to do it.
            bool noKeyHolesNeeded = true;
            // Decompose to outers and cutters
            Paths[] decomp = GeoWrangler.getDecomposed(mergedPolyData);

            Paths outers = decomp[(int)GeoWrangler.type.outer];
            Paths cutters = decomp[(int)GeoWrangler.type.cutter];

            int oCount = outers.Count;
            int cCount = cutters.Count;

            // Is any cutter fully enclosed within an outer?
            for (int outer = 0; outer < oCount; outer++)
            {
                double origArea = Math.Abs(Clipper.Area(outers[outer]));
                for (int cutter = 0; cutter < cCount; cutter++)
                {
                    c.Clear();
                    c.AddSubject(outers[outer]);
                    c.AddClip(cutters[cutter]);
                    Paths test = new();
                    c.Execute(ClipType.Union, FillRule.Positive, test);
                    test = GeoWrangler.reOrderXY(test);

                    double uArea = test.Sum(t => Math.Abs(Clipper.Area(t)));

                    if (!(Math.Abs(uArea - origArea) < double.Epsilon))
                    {
                        continue;
                    }

                    noKeyHolesNeeded = false;
                    break;
                }
                if (!noKeyHolesNeeded)
                {
                    break;
                }
            }

            if (noKeyHolesNeeded)
            {
                // Send back our merged data to the caller; no keyholes needed and overlaps are reconciled.
                return GeoWrangler.pointFsFromPaths(mergedPolyData, scaleFactorForOperation);
            }
            
            // So it turns out we need to worry about keyholes if we didn't return already. Let's get started.
            
            // Here, we can run into trouble....we might have a set of polygons which need to get keyholed. For example, where we have fully enclosed 'cutters' within an outer boundary.
            // Can geoWrangler help us out here?

            // We need to run the fragmenter here because the keyholer / raycaster pipeline needs points for emission.
            Fragmenter f = new(Convert.ToDouble(resolution) * scaleFactorForOperation);
            Paths toKeyHole = f.fragmentPaths(mergedPolyData);
            // Nudge the default keyhole size to avoid collapsing existing keyholes.
            Paths keyHoled = GeoWrangler.makeKeyHole(toKeyHole, true, RayCast.inversionMode.x, customSizing:customSizing, extension:Convert.ToDouble(extension));

            if (!keyHoled.Any())
            {
                return sourceData; // might need to be mergedPolyData here.
            }

            double newArea = Clipper.Area(keyHoled);
            
            if (newArea > area) // We gained area, which suggests we lost a hole. Let's try with the alternate inversion scheme.
            {
                keyHoled = GeoWrangler.makeKeyHole(toKeyHole, true, RayCast.inversionMode.y, customSizing:customSizing, extension:Convert.ToDouble(extension));
            }

            mergedPolyData = GeoWrangler.close(keyHoled);

            // We got some resulting geometry from our Boolean so let's process it to send back to the caller.
            List<GeoLibPointF[]> refinedData = new();
            
            int rpdCount = mergedPolyData.Count;

            // Switch our fragmenter to use a new configuration for the downsized geometry.
            f = new Fragmenter(Convert.ToDouble(resolution), scaleFactorForOperation);

            // Convert back our geometry.                
            for (int rPoly = 0; rPoly < rpdCount; rPoly++)
            {
                // We have to re-fragment as the overlap processing changed the geometry heavily.
                refinedData.Add(f.fragmentPath(GeoWrangler.pointFFromPath(mergedPolyData[rPoly], scaleFactorForOperation)));
            }

            return refinedData;
        }
        catch (Exception)
        {
            return sourceData;
        }
    }
    
}