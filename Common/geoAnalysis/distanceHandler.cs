using System.Globalization;
using Clipper2Lib;
using geoWrangler;
using KDTree;
using utility;
#pragma warning disable CS8618

namespace geoAnalysis;

public class DistanceHandler
{
    public enum spacingCalcModes { spacing, enclosure, spacingOld, enclosureOld } // exp triggers projection from shortest edge for overlap evaluation.

    public delegate void ErrorRep(string text, string header);
    public ErrorRep errorRep { get; set; }
    private const bool debug = false;
    public PathsD resultPaths { get; private set; }
    private double resultDistance;
    public string distanceString { get; private set; }
    private RayCast.inversionMode invert;

    private static void ZFillCallback(PointD bot1, PointD top1, PointD bot2, PointD top2, ref PointD pt)
    {
        pt.z = -1; // Tag our intersection points.
    }

    private class spaceResult
    {
        public bool done { get; set; } // to flag early return.
        public double distance { get; set; }
        public PathsD resultPaths { get; set; }

        public spaceResult()
        {
            resultPaths = [];
        }
    }

    public DistanceHandler(bool debugCalc, PathsD aPaths, PathsD bPaths, int mode, bool runThreaded)
    {
        distanceHandlerLogic(debugCalc, aPaths, bPaths, mode, runThreaded);
    }

    private void distanceHandlerLogic(bool debugCalc, PathsD aPaths_, PathsD bPaths_, int mode, bool runThreaded)
    {
        PathsD aPaths = new (aPaths_);
        PathsD bPaths = new (bPaths_);
        resultPaths = [];
        // Safety check for no active layers.
        if (aPaths is [{ Count: 1 }] ||
            bPaths is [{ Count: 1 }])
        {
            resultPaths.Add([new PointD(0, 0)]);
            distanceString = "N/A";
        }
        else
        {
            invert = RayCast.inversionMode.x;
            // Experimental check to see whether we can simplify this approach.
            // bool isEnclosed = true; // GeoWrangler.enclosed(aPaths, bPaths);

            bool anyOverlap = GeoWrangler.anyPartialOverlap(aPaths, bPaths);

            if (anyOverlap)
            {
                // Overlap method sets result fields.
                overlap(debugCalc, aPaths, bPaths, mode, runThreaded);
            }
            else
            {
                spaceResult result = fastKDTree(aPaths, bPaths, mode);

                resultDistance = result.distance;
                resultPaths = result.resultPaths;
            }

            distanceString = resultDistance.ToString(CultureInfo.InvariantCulture);
        }
    }

    private spaceResult fastKDTree(PathsD aPaths, PathsD bPaths, int mode)
    {
        int numberOfPoints = 0;
        double currentMinimum = 0;
        PathD minimumDistancePath = [];
        bool resultNeedsInversion = false;

        double refArea = 0;

        foreach (PathD t in bPaths)
        {
            numberOfPoints += t.Count;
            refArea += Clipper.Area(t);
        }
        // KDTree to store the points from our target shape(s)
        KDTree<PointD> pTree = new(2, numberOfPoints);

        for (int shapeA = 0; shapeA < aPaths.Count; shapeA++)
        {
            foreach (PathD t in bPaths)
            {
                for (int pointB = 0; pointB < t.Count; pointB++)
                {
                    try
                    {
                        pTree.AddPoint([t[pointB].x, t[pointB].y], new PointD(t[pointB].x, t[pointB].y));
                    }
                    catch (Exception ex)
                    {
                        errorRep?.Invoke("Oops", "jobEngine() KDTree error: " + ex);
                    }
                }
            }

            // Do we need to invert the result?
            ClipperD c = new(Constants.roundingDecimalPrecision) {PreserveCollinear = true};
            PathsD oCheck = [];
            c.AddClip(bPaths);
            c.AddSubject(aPaths[shapeA]);
            c.Execute(ClipType.Union, FillRule.EvenOdd, oCheck);
            oCheck = GeoWrangler.reOrderXY(oCheck);

            double oCheckArea = oCheck.Sum(Clipper.Area);

            resultNeedsInversion = mode switch
            {
                // negative value since we're fully inside a containing polygon.
                (int)spacingCalcModes.enclosure or (int)spacingCalcModes.enclosureOld => Math.Abs(Math.Abs(oCheckArea) -
                    Math.Abs(refArea)) > double.Epsilon,
                // negative value since we're fully outside a containing polygon.
                (int)spacingCalcModes.spacing or (int)spacingCalcModes.spacingOld => Math.Abs(Math.Abs(oCheckArea) -
                    Math.Abs(refArea)) < double.Epsilon,
                _ => resultNeedsInversion
            };

            // We can query here for the minimum distance for each shape combination.
            for (int pointA = 0; pointA < aPaths[shapeA].Count; pointA++)
            {
                // '1' forces a single nearest neighbor to be returned.
                NearestNeighbour<PointD> pIter = pTree.NearestNeighbors([aPaths[shapeA][pointA].x, aPaths[shapeA][pointA].y
                ], 1);
                while (pIter.MoveNext())
                {
                    double currentDistance = pIter.CurrentDistance;

                    if ((shapeA != 0 || pointA != 0) && !(currentDistance < currentMinimum))
                    {
                        continue;
                    }

                    minimumDistancePath.Clear();
                    minimumDistancePath.Add(new PointD(aPaths[shapeA][pointA]));
                    minimumDistancePath.Add(new PointD((pIter.Current.x + aPaths[shapeA][pointA].x) / 2.0f, (pIter.Current.y + aPaths[shapeA][pointA].y) / 2.0f));
                    minimumDistancePath.Add(new PointD(pIter.Current.x, pIter.Current.y));
                    currentMinimum = currentDistance;
                }
            }
        }

        spaceResult result = new()
        {
            resultPaths = [minimumDistancePath],
            // k-d tree distance is the squared distance. Need to scale and sqrt
            distance = Math.Sqrt(currentMinimum)
        };

        if (resultNeedsInversion)
        {
            result.distance *= -1;
        }

        return result;
    }

    private void overlap(bool debugCalc, PathsD aPaths, PathsD bPaths, int mode, bool runThreaded)
    {
        bool completeOverlap = false;
        foreach (PathD a in aPaths)
        {
            PathD layerAPath = [];
            for (int pa = 0; pa < a.Count; pa++)
            {
                layerAPath.Add(new PointD(a[pa].x, a[pa].y, 1));
            }
            foreach (PathD b in bPaths)
            {
                PathD layerBPath = [];
                for (int pb = 0; pb < b.Count; pb++)
                {
                    layerBPath.Add(new PointD(b[pb].x, b[pb].y, 2));
                }

                PathsD overlapShape = [];

                // Check for complete overlap
                ClipperD c = new(Constants.roundingDecimalPrecision) {PreserveCollinear = true, ZCallback = ZFillCallback};
                // Try a union and see whether the point count of the perimeter changes. This might break for re-entrant cases, but those are problematic anyway.
                PathsD fullOverlapCheck = [];
                c.AddSubject(layerAPath);
                c.AddClip(layerBPath);
                c.Execute(ClipType.Union, FillRule.EvenOdd, fullOverlapCheck);
                double aArea = Math.Abs(Clipper.Area(layerAPath));
                double bArea = Math.Abs(Clipper.Area(layerBPath));
                double uArea = fullOverlapCheck.Sum(Clipper.Area);
                uArea = Math.Abs(uArea);

                // If overlap area matches either of the input areas, we have a full overlap
                if (Math.Abs(aArea - uArea) < double.Epsilon || Math.Abs(bArea - uArea) < double.Epsilon)
                {
                    completeOverlap = true;
                }

                if (mode is (int)spacingCalcModes.spacing or (int)spacingCalcModes.spacingOld) // spacing
                {
                    // Perform an area check in case of overlap.
                    // Overlap means X/Y negative space needs to be reported.
                    AreaHandler aH = new([layerAPath], [layerBPath], maySimplify: false, perPoly: false);
                    overlapShape = aH.listOfOutputPoints;
                }

                if (!completeOverlap && mode is (int)spacingCalcModes.enclosure or (int)spacingCalcModes.enclosureOld) // enclosure
                {
                    // Need to find the region outside our enclosure shape. We use the modifier to handle this.
                    c.Clear();
                    c.ZCallback = ZFillCallback;
                    // Try cutting layerB from layerA
                    c.AddSubject(layerAPath);
                    c.AddClip(layerBPath);
                    c.Execute(ClipType.Difference, FillRule.EvenOdd, overlapShape);
                }

                if (completeOverlap || !overlapShape.Any())
                {
                    continue;
                }

                // This is needed to ensure the downstream evaluation works.
                overlapShape = GeoWrangler.reOrderYX(overlapShape);

                spaceResult result = doPartialOverlap(debugCalc, overlapShape, layerAPath, layerBPath, mode, runThreaded);
                if (!result.done || resultPaths.Any() && !(result.distance < resultDistance))
                {
                    continue;
                }

                resultPaths = result.resultPaths;
                resultDistance = result.distance;
            }
        }
    }

    private spaceResult doPartialOverlap(bool debugCalc, PathsD overlapShape, PathD aPath, PathD bPath, int mode, bool runThreaded)
    {
        spaceResult result = new();
        int oCount = overlapShape.Count;

        try
        {
            // Process the overlap shape polygon(s) to evaluate the overlap.
            for (int poly = 0; poly < oCount; poly++)
            {
                bool osOrient = Clipper.IsPositive(overlapShape[poly]);
                if (osOrient)
                {
                    // Reverse to accommodate raycaster needs for emission (normal computation). 
                    overlapShape[poly].Reverse();
                }

                // Use Z value to figure out which part of the overlap shape belongs to which original input source.
                int ptCount = overlapShape[poly].Count;
                bool[] residesOnAEdge = new bool[ptCount];
                bool[] residesOnBEdge = new bool[ptCount];

                for (int pt = 0; pt < ptCount; pt++)
                {
                    residesOnAEdge[pt] = false;
                    residesOnBEdge[pt] = false;
                    switch (overlapShape[poly][pt].z)
                    {
                        case -1:
                            // residesOnAEdge[pt] = true;
                            // residesOnBEdge[pt] = true;
                            break;
                        case 1:
                            residesOnAEdge[pt] = true;
                            break;
                        case 2:
                            residesOnBEdge[pt] = true;
                            break;
                    }
                }

                // Now we need to construct our Paths for the overlap edges, based on the true/false case for each array.
                PathsD aOverlapEdge = [];
                PathsD bOverlapEdge = [];

                PathD tempPath = [];
                for (int i = 0; i < ptCount; i++)
                {
                    if (residesOnAEdge[i])
                    {
                        tempPath.Add(new PointD(overlapShape[poly][i]));
                    }
                    else // not found on A edge, probably resides on B edge, but we'll check that separately to keep code readable.
                    {
                        // If we have a tempPath with more than a single point in it, we have been creating a segment.
                        // Since we haven't found this point in the A edge, commit the segment and clear the temporary path.
                        if (tempPath.Count < 1)
                        {
                            continue;
                        }

                        // We have some points, but now we're getting a new segment.
                        aOverlapEdge.Add(new PathD(tempPath));
                        tempPath.Clear();
                    }
                }
                if (tempPath.Count > 1)
                {
                    aOverlapEdge.Add(tempPath);
                }

                tempPath.Clear();
                for (int i = 0; i < ptCount; i++)
                {
                    if (residesOnBEdge[i])
                    {
                        tempPath.Add(new PointD(overlapShape[poly][i]));
                    }
                    else
                    {
                        if (tempPath.Count <= 1)
                        {
                            continue;
                        }

                        // We have some points, but now we're getting a new segment.
                        bOverlapEdge.Add(new PathD(tempPath));
                        tempPath.Clear();
                    }
                }
                if (tempPath.Count > 1)
                {
                    bOverlapEdge.Add(tempPath);
                }

                // Strip any single point path from the paths.
                bool changed = true;
                while (changed)
                {
                    changed = false;
                    for (int i = aOverlapEdge.Count; i == -1; i--)
                    {
                        if (aOverlapEdge[i].Count > 1)
                        {
                            continue;
                        }
                        aOverlapEdge.RemoveAt(i);
                        changed = true;
                    }
                }

                changed = true;
                while (changed)
                {
                    changed = false;
                    for (int i = bOverlapEdge.Count; i == -1; i--)
                    {
                        if (bOverlapEdge[i].Count > 1)
                        {
                            continue;
                        }
                        bOverlapEdge.RemoveAt(i);
                        changed = true;
                    }
                }

                // Walk our edges to figure out the overlap.
                var poly1 = poly;
                foreach (spaceResult tResult in aOverlapEdge.SelectMany(t => bOverlapEdge.Select(t1 => overlapAssess(debugCalc, overlapShape[poly1], t, t1, aPath, bPath, mode, runThreaded)).Where(tResult => result.resultPaths.Count == 0 || tResult.distance > result.distance)))
                {
                    result.distance = tResult.distance;
                    if (!debugCalc)
                    {
                        result.resultPaths = tResult.resultPaths;
                    }
                    else
                    {
                        result.resultPaths.AddRange(tResult.resultPaths);
                    }
                }
            }
        }
        catch (Exception e)
        {
            errorRep.Invoke(e.ToString(), "Oops");
        }
        
        result.done = true;
        if (!debugCalc)
        {
            result.distance = -result.distance;
        }
        return result;
    }

    private spaceResult overlapAssess(bool debugCalc, PathD overlapPoly, PathD aOverlapEdge, PathD bOverlapEdge, PathD aPath, PathD bPath, int mode, bool runThreaded)
    {
        spaceResult result = new();
        if (aOverlapEdge.Count == 0)
        {
            return result;
        }
        if (bOverlapEdge.Count == 0)
        {
            return result;
        }
        double maxDistance_orthogonalFallback = 0; // fallback for the overlap case of orthogonal geometry, where self-intersections occur.
        int maxDistance_fallbackIndex = 0;
        double maxDistance = 0; // for the overlap cases.

        double lengthA = 0;
        double lengthB = 0;

        aOverlapEdge = GeoWrangler.reOrderXY(aOverlapEdge);
        bOverlapEdge = GeoWrangler.reOrderXY(bOverlapEdge);
        PathD extractedPath = [];
        bool usingAEdge = false;

        /* Prior to 1.7, we used the edge from the layerB combination to raycast. This was not ideal.
            We should have used the shortest edge. 1.7 adds an option to make this the behavior, and the old 
            approach is retained for now as well.
        */

        PointD shortestPathBeforeStartPoint = new(0, 0);
        PointD shortestPathAfterEndPoint = new(0, 0);

        switch (mode)
        {
            case (int)spacingCalcModes.spacing:
            {
                // Find the shortest edge length and use that for the projection reference
                // Calculate lengths and check.
                for (int pt = 0; pt < aOverlapEdge.Count - 1; pt++)
                {
                    lengthA += Math.Sqrt(Utils.myPow(aOverlapEdge[pt + 1].x - aOverlapEdge[pt].x, 2) + Utils.myPow(aOverlapEdge[pt + 1].y - aOverlapEdge[pt].y, 2));
                }
                for (int pt = 0; pt < bOverlapEdge.Count - 1; pt++)
                {
                    lengthB += Math.Sqrt(Utils.myPow(bOverlapEdge[pt + 1].x - bOverlapEdge[pt].x, 2) + Utils.myPow(bOverlapEdge[pt + 1].y - bOverlapEdge[pt].y, 2));
                }

                extractedPath = bOverlapEdge; // need a default in case of equivalent lengths

                if (extractedPath.Count == 0)
                {
                    // No common overlap.
                    result.distance = 0;
                    return result;
                }

                // Here we need to go back to our input polygons to derive the edge segment normal for the start/end of each edge.
                if (lengthA < lengthB)
                {
                    usingAEdge = true;
                    extractedPath = aOverlapEdge;
                    int startIndex = aPath.FindIndex(p => p.Equals(extractedPath[0]));
                    if (startIndex == 0)
                    {
                        startIndex = aPath.Count - 1;
                    }
                    if (startIndex == -1)
                    {
                        // Failed to find it cheaply. Let's try a different approach.
                        double startDistanceCheck = GeoWrangler.distanceBetweenPoints(extractedPath[0], aPath[0]);
                        startIndex = 0;
                        for (int pt = 1; pt < aPath.Count; pt++)
                        {
                            double distanceCheck = GeoWrangler.distanceBetweenPoints(extractedPath[0], aPath[pt]);
                            if (!(distanceCheck < startDistanceCheck))
                            {
                                continue;
                            }

                            startDistanceCheck = distanceCheck;
                            startIndex = pt;
                        }
                    }

                    int endIndex = aPath.FindIndex(p => p.Equals(extractedPath[^1]));
                    if (endIndex == aPath.Count - 1)
                    {
                        endIndex = 0;
                    }
                    if (endIndex == -1)
                    {
                        // Failed to find it cheaply. Let's try a different approach.
                        double endDistanceCheck = GeoWrangler.distanceBetweenPoints(extractedPath[^1], aPath[0]);
                        endIndex = 0;
                        for (int pt = 1; pt < aPath.Count; pt++)
                        {
                            double distanceCheck = GeoWrangler.distanceBetweenPoints(extractedPath[^1], aPath[pt]);
                            if (!(distanceCheck < endDistanceCheck))
                            {
                                continue;
                            }

                            endDistanceCheck = distanceCheck;
                            endIndex = pt;
                        }
                    }

                    shortestPathBeforeStartPoint = new PointD(aPath[startIndex]);
                    shortestPathAfterEndPoint = new PointD(aPath[endIndex]);
                }
                else
                {
                    int startIndex = bPath.FindIndex(p => p.Equals(extractedPath[0]));
                    if (startIndex == 0)
                    {
                        startIndex = bPath.Count - 1;
                    }
                    if (startIndex == -1)
                    {
                        // Failed to find it cheaply. Let's try a different approach.
                        double startDistanceCheck = GeoWrangler.distanceBetweenPoints(extractedPath[0], bPath[0]);
                        startIndex = 0;

                        for (int pt = 1; pt < bPath.Count; pt++)
                        {
                            double distanceCheck = GeoWrangler.distanceBetweenPoints(extractedPath[0], bPath[pt]);
                            if (!(distanceCheck < startDistanceCheck))
                            {
                                continue;
                            }

                            startDistanceCheck = distanceCheck;
                            startIndex = pt;
                        }
                    }

                    int endIndex = bPath.FindIndex(p => p.Equals(extractedPath[^1]));
                    if (endIndex == bPath.Count - 1)
                    {
                        endIndex = 0;
                    }
                    if (endIndex == -1)
                    {
                        // Failed to find it cheaply. Let's try a different approach.
                        double endDistanceCheck = GeoWrangler.distanceBetweenPoints(extractedPath[^1], bPath[0]);
                        endIndex = 0;
                        for (int pt = 1; pt < bPath.Count; pt++)
                        {
                            double distanceCheck = GeoWrangler.distanceBetweenPoints(extractedPath[^1], bPath[pt]);
                            if (!(distanceCheck < endDistanceCheck))
                            {
                                continue;
                            }

                            endDistanceCheck = distanceCheck;
                            endIndex = pt;
                        }
                    }

                    shortestPathBeforeStartPoint = new PointD(bPath[startIndex]);
                    shortestPathAfterEndPoint = new PointD(bPath[endIndex]);
                }

                break;
            }
            case (int)spacingCalcModes.spacingOld or (int)spacingCalcModes.enclosure or (int)spacingCalcModes.enclosureOld:
            {
                extractedPath = bOverlapEdge;

                if (extractedPath.Count == 0) // No common edge.
                {
                    result.distance = 0;
                    return result;
                }

                int startIndex = bPath.FindIndex(p => p.Equals(extractedPath[0]));
                if (startIndex == 0)
                {
                    startIndex = bPath.Count - 1;
                }
                if (startIndex == -1)
                {
                    // Failed to find it cheaply. Let's try a different approach.
                    double startDistanceCheck = GeoWrangler.distanceBetweenPoints(extractedPath[0], bPath[0]);
                    startIndex = 0;
                    for (int pt = 1; pt < bPath.Count; pt++)
                    {
                        double distanceCheck = GeoWrangler.distanceBetweenPoints(extractedPath[0], bPath[pt]);
                        if (!(distanceCheck < startDistanceCheck))
                        {
                            continue;
                        }

                        startDistanceCheck = distanceCheck;
                        startIndex = pt;
                    }
                }

                int endIndex = bPath.FindIndex(p => p.Equals(extractedPath[^1]));
                if (endIndex == bPath.Count - 1)
                {
                    endIndex = 0;
                }
                if (endIndex == -1)
                {
                    // Failed to find it cheaply. Let's try a different approach.
                    double endDistanceCheck = GeoWrangler.distanceBetweenPoints(extractedPath[^1], bPath[0]);
                    for (int pt = 1; pt < bPath.Count; pt++)
                    {
                        double distanceCheck = GeoWrangler.distanceBetweenPoints(extractedPath[^1], bPath[pt]);
                        endIndex = 0;
                        if (!(distanceCheck < endDistanceCheck))
                        {
                            continue;
                        }

                        endDistanceCheck = distanceCheck;
                        endIndex = pt;
                    }
                }

                try
                {
                    shortestPathBeforeStartPoint = new PointD(bPath[startIndex]);
                }
                catch (Exception)
                {
                    errorRep.Invoke("distanceHandler: shortestPathBeforeStartPoint failed.", "Oops");
                }
                try
                {
                    shortestPathAfterEndPoint = new PointD(bPath[endIndex]);
                }
                catch (Exception)
                {
                    errorRep.Invoke("distanceHandler: shortestPathAfterEndPoint failed.", "Oops");
                }

                break;
            }
        }

        // No blurry rays, so no point running the inner loop threaded. We thread the outer loop (the emission edge raycast), though. Testing showed small performance improvement for this approach.
        RayCast rc = new(extractedPath, overlapPoly, 1000, runThreaded, invert, 0, true, false, shortestPathBeforeStartPoint, shortestPathAfterEndPoint);

        if (debug || debugCalc)
        {
            result.done = true;
            result.resultPaths = rc.getRays();
            result.distance = -1;
            return result;
        }

        PathsD clippedLines = rc.getClippedRays();

        bool validOverlapFound = false;

        // Need to scan for the maximum path length to record the overlap.
        for (int line = 0; line < clippedLines.Count; line++)
        {
            double lineLength = rc.getRayLength(line);
            // With the new LER implementation, normals can be fired in problematic directions.
            // Valid overlaps don't have start and end points on the same polygon.
            // Happily, we feed this polygon-wise, which makes the evaluation much easier.
            bool validOverlap = true;

            double startPointCheck_A_dist = GeoWrangler.distanceBetweenPoints(clippedLines[line][0], aOverlapEdge[0]);
            double endPointCheck_A_dist = GeoWrangler.distanceBetweenPoints(clippedLines[line][1], aOverlapEdge[0]);
            foreach (PointD t in aOverlapEdge)
            {
                double temp_startPointCheck_A_dist = GeoWrangler.distanceBetweenPoints(clippedLines[line][0], t);
                double temp_endPointCheck_A_dist = GeoWrangler.distanceBetweenPoints(clippedLines[line][1], t);
                if (temp_startPointCheck_A_dist < startPointCheck_A_dist)
                {
                    startPointCheck_A_dist = temp_startPointCheck_A_dist;
                }
                if (temp_endPointCheck_A_dist < endPointCheck_A_dist)
                {
                    endPointCheck_A_dist = temp_endPointCheck_A_dist;
                }
            }

            double startPointCheck_B_dist = GeoWrangler.distanceBetweenPoints(clippedLines[line][0], bOverlapEdge[0]);
            double endPointCheck_B_dist = GeoWrangler.distanceBetweenPoints(clippedLines[line][1], bOverlapEdge[0]);
            foreach (PointD t in bOverlapEdge)
            {
                double temp_startPointCheck_B_dist = GeoWrangler.distanceBetweenPoints(clippedLines[line][0], t);
                double temp_endPointCheck_B_dist = GeoWrangler.distanceBetweenPoints(clippedLines[line][1], t);
                if (temp_startPointCheck_B_dist < startPointCheck_B_dist)
                {
                    startPointCheck_B_dist = temp_startPointCheck_B_dist;
                }
                if (temp_endPointCheck_B_dist < endPointCheck_B_dist)
                {
                    endPointCheck_B_dist = temp_endPointCheck_B_dist;
                }
            }

            bool ortho = Math.Abs(clippedLines[line][0].x - clippedLines[line][1].x) < Constants.tolerance || Math.Abs(clippedLines[line][0].y - clippedLines[line][1].y) < Constants.tolerance;

            const double threshold = 0.01500; // arbitrary, dialed in by hand.

            if (!(ortho || usingAEdge && startPointCheck_A_dist < threshold && endPointCheck_A_dist < threshold || !usingAEdge && startPointCheck_B_dist < threshold && endPointCheck_B_dist < threshold
               ))
            {
                // This is a special situation, it turns out.
                // There is one specific scenario where this overlap case (start and end on the same geometry) is valid - orthogonal shapes with a bisection.
                // The orthogonal shapes cause the rays to hit the opposite side of the shape. We don't want to reject this.
                if (!(lineLength > maxDistance_orthogonalFallback))
                {
                    continue;
                }

                maxDistance_fallbackIndex = line;
                maxDistance_orthogonalFallback = lineLength;
            }
            else
            {
                if (usingAEdge && startPointCheck_A_dist < threshold && endPointCheck_B_dist > threshold ||
                    !usingAEdge && startPointCheck_B_dist < threshold && endPointCheck_A_dist > threshold)
                {
                    validOverlap = false;
                }

                if (!validOverlap || !(lineLength > maxDistance))
                {
                    continue;
                }

                validOverlapFound = true;
                maxDistance = lineLength;
                result.resultPaths.Clear();
                result.resultPaths.Add([
                    new PointD(clippedLines[line][0]),
                    new PointD(clippedLines[line][0]),
                    new PointD(clippedLines[line][clippedLines[line].Count - 1])
                ]);
            }
        }

        try
        {
            if (clippedLines.Any() && !validOverlapFound || maxDistance_orthogonalFallback > maxDistance)
            {
                // Couldn't find a valid overlap so assume the orthogonal fallback is needed.
                maxDistance = maxDistance_orthogonalFallback;
                result.resultPaths.Clear();
                result.resultPaths.Add([
                    new PointD(clippedLines[maxDistance_fallbackIndex][0]),
                    new PointD(clippedLines[maxDistance_fallbackIndex][0]),
                    new PointD(
                        clippedLines[maxDistance_fallbackIndex]
                            [clippedLines[maxDistance_fallbackIndex].Count - 1])
                ]);
            }
        }
        catch (Exception)
        {
            // Harmless - we'll reject the case and move on.
        }

        result.distance = maxDistance;
        return result;
    }
}