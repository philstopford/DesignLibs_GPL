using System.Globalization;
using Clipper2Lib;
using geoLib;
using geoWrangler;
using KDTree;
using utility;

namespace geoAnalysis;

using Path = List<Point64>;
using Paths = List<List<Point64>>;

public class DistanceHandler
{
    public delegate void ErrorRep(string text, string header);
    public ErrorRep errorRep { get; set; }

    public enum spacingCalcModes { spacing, enclosure, spacingOld, enclosureOld } // exp triggers projection from shortest edge for overlap evaluation.

    private const bool debug = false;
    public Paths resultPaths { get; private set; }
    private double resultDistance;
    public string distanceString { get; private set; }
    private RayCast.inversionMode invert;

    private void ZFillCallback(Point64 bot1, Point64 top1, Point64 bot2, Point64 top2, ref Point64 pt)
    {
        pt.Z = -1; // Tag our intersection points.
    }

    private class spaceResult
    {
        public bool done { get; set; } // to flag early return.
        public double distance { get; set; }
        public Paths resultPaths { get; set; }

        public spaceResult()
        {
            resultPaths = new Paths();
        }
    }

    public DistanceHandler(Paths aPaths, Paths bPaths, double resolution, int scaleFactorForOperation, int subMode, bool runThreaded, bool debug)
    {
        distanceHandlerLogic(aPaths, bPaths, resolution, scaleFactorForOperation, subMode, runThreaded, debug);
    }

    private void distanceHandlerLogic(Paths aPaths, Paths bPaths, double resolution, int scaleFactorForOperation, int subMode, bool runThreaded, bool debugCalc)
    {
        resultPaths = new Paths();
        // Safety check for no active layers.
        if (aPaths.Count == 1 && aPaths[0].Count == 1 ||
            bPaths.Count == 1 && bPaths[0].Count == 1)
        {
            resultPaths.Add(new() {new Point64(0, 0)});
            distanceString = "N/A";
        }
        else
        {
            invert = RayCast.inversionMode.x;
            // Experimental check to see whether we can simplify this approach.
            bool isEnclosed = GeoWrangler.enclosed(aPaths, bPaths);


            if (!isEnclosed)
            {
                // Overlap method sets result fields.
                overlap(aPaths, bPaths, resolution, scaleFactorForOperation, subMode, runThreaded, debugCalc);
            }
            else
            {
                spaceResult result = fastKDTree(aPaths, bPaths, scaleFactorForOperation, subMode);

                resultDistance = result.distance;
                resultPaths = result.resultPaths;
            }

            distanceString = resultDistance.ToString(CultureInfo.InvariantCulture);
        }
    }

    private spaceResult fastKDTree(Paths aPaths, Paths bPaths, int scaleFactorForOperation, int subMode)
    {
        int numberOfPoints = 0;
        double currentMinimum = 0;
        Path minimumDistancePath = new();
        bool resultNeedsInversion = false;

        double refArea = 0;

        foreach (Path t in bPaths)
        {
            numberOfPoints += t.Count;
            refArea += Clipper.Area(t);
        }
        // KDTree to store the points from our target shape(s)
        KDTree<GeoLibPointF> pTree = new(2, numberOfPoints);

        for (int shapeA = 0; shapeA < aPaths.Count; shapeA++)
        {
            foreach (Path t in bPaths)
            {
                for (int pointB = 0; pointB < t.Count; pointB++)
                {
                    try
                    {
                        pTree.AddPoint(new double[] { t[pointB].X, t[pointB].Y }, new GeoLibPointF(t[pointB].X, t[pointB].Y));
                    }
                    catch (Exception ex)
                    {
                        errorRep?.Invoke("Oops", "jobEngine() KDTree error: " + ex);
                    }
                }
            }

            // Do we need to invert the result?
            Clipper64 c = new() {PreserveCollinear = true};
            Paths oCheck = new();
            c.AddClip(bPaths);
            c.AddSubject(aPaths[shapeA]);
            c.Execute(ClipType.Union, FillRule.EvenOdd, oCheck);
            oCheck = GeoWrangler.reOrderXY(oCheck);

            double oCheckArea = oCheck.Sum(t => Clipper.Area(t));

            if (subMode == (int)spacingCalcModes.enclosure || subMode == (int)spacingCalcModes.enclosureOld) // negative value since we're fully inside a containing polygon.
            {
                resultNeedsInversion = Math.Abs(Math.Abs(oCheckArea) - Math.Abs(refArea)) > double.Epsilon;
            }
            if (subMode == (int)spacingCalcModes.spacing || subMode == (int)spacingCalcModes.spacingOld) // negative value since we're fully outside a containing polygon.
            {
                resultNeedsInversion = Math.Abs(Math.Abs(oCheckArea) - Math.Abs(refArea)) < double.Epsilon;
            }

            // We can query here for the minimum distance for each shape combination.
            for (int pointA = 0; pointA < aPaths[shapeA].Count; pointA++)
            {
                // '1' forces a single nearest neighbor to be returned.
                NearestNeighbour<GeoLibPointF> pIter = pTree.NearestNeighbors(new double[] { aPaths[shapeA][pointA].X, aPaths[shapeA][pointA].Y }, 1);
                while (pIter.MoveNext())
                {
                    double currentDistance = pIter.CurrentDistance;

                    if ((shapeA != 0 || pointA != 0) && !(currentDistance < currentMinimum))
                    {
                        continue;
                    }

                    minimumDistancePath.Clear();
                    minimumDistancePath.Add(new Point64(aPaths[shapeA][pointA]));
                    minimumDistancePath.Add(new Point64((pIter.Current.X + aPaths[shapeA][pointA].X) / 2.0f, (pIter.Current.Y + aPaths[shapeA][pointA].Y) / 2.0f));
                    minimumDistancePath.Add(new Point64(pIter.Current.X, pIter.Current.Y));
                    currentMinimum = currentDistance;
                }
            }
        }

        spaceResult result = new()
        {
            resultPaths = new Paths {minimumDistancePath},
            // k-d tree distance is the squared distance. Need to scale and sqrt
            distance = Math.Sqrt(currentMinimum / Utils.myPow(scaleFactorForOperation, 2))
        };

        if (resultNeedsInversion)
        {
            result.distance *= -1;
        }

        return result;
    }

    private void overlap(Paths aPaths, Paths bPaths, double resolution, int scaleFactorForOperation, int subMode, bool runThreaded, bool debugCalc)
    {
        bool completeOverlap = false;
        foreach (Path a in aPaths)
        {
            Path layerAPath = new();
            for (int pa = 0; pa < a.Count; pa++)
            {
                layerAPath.Add(new(a[pa].X, a[pa].Y, 1));
            }
            foreach (Path b in bPaths)
            {
                Path layerBPath = new();
                for (int pb = 0; pb < b.Count; pb++)
                {
                    layerBPath.Add(new(b[pb].X, b[pb].Y, 2));
                }

                Paths overlapShape = new();

                // Check for complete overlap
                Clipper64 c = new() {PreserveCollinear = true, ZFillFunc = ZFillCallback};
                // Try a union and see whether the point count of the perimeter changes. This might break for re-entrant cases, but those are problematic anyway.
                Paths fullOverlapCheck = new();
                c.AddSubject(layerAPath);
                c.AddClip(layerBPath);
                c.Execute(ClipType.Union, FillRule.EvenOdd, fullOverlapCheck);
                double aArea = Math.Abs(Clipper.Area(layerAPath));
                double bArea = Math.Abs(Clipper.Area(layerBPath));
                double uArea = fullOverlapCheck.Sum(t => Clipper.Area(t));
                uArea = Math.Abs(uArea);

                // If overlap area matches either of the input areas, we have a full overlap
                if (Math.Abs(aArea - uArea) < double.Epsilon || Math.Abs(bArea - uArea) < double.Epsilon)
                {
                    completeOverlap = true;
                }

                if (subMode == (int)spacingCalcModes.spacing || subMode == (int)spacingCalcModes.spacingOld) // spacing
                {
                    // Perform an area check in case of overlap.
                    // Overlap means X/Y negative space needs to be reported.
                    AreaHandler aH = new(new Paths { layerAPath }, new Paths { layerBPath }, maySimplify: false, perPoly: false, scaleFactorForOperation: 1.0);
                    overlapShape = aH.listOfOutputPoints;
                }

                if (!completeOverlap && (subMode == (int)spacingCalcModes.enclosure || subMode == (int)spacingCalcModes.enclosureOld)) // enclosure
                {
                    // Need to find the region outside our enclosure shape. We use the modifier to handle this.
                    c.Clear();
                    c.ZFillFunc = ZFillCallback;
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

                spaceResult result = doPartialOverlap(overlapShape, layerAPath, layerBPath, resolution, scaleFactorForOperation, subMode, runThreaded, debugCalc);
                if (!result.done || resultPaths.Any() && !(result.distance < resultDistance))
                {
                    continue;
                }

                resultPaths = result.resultPaths;
                resultDistance = result.distance;
            }
        }
    }

    private spaceResult doPartialOverlap(Paths overlapShape, Path aPath, Path bPath, double resolution, int scaleFactorForOperation, int subMode, bool runThreaded, bool debugCalc)
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
                    switch (overlapShape[poly][pt].Z)
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
                Paths aOverlapEdge = new();
                Paths bOverlapEdge = new();

                Path tempPath = new();
                for (int i = 0; i < ptCount; i++)
                {
                    if (residesOnAEdge[i])
                    {
                        tempPath.Add(new Point64(overlapShape[poly][i]));
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
                        aOverlapEdge.Add(tempPath.ToList());
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
                        tempPath.Add(new Point64(overlapShape[poly][i]));
                    }
                    else
                    {
                        if (tempPath.Count <= 1)
                        {
                            continue;
                        }

                        // We have some points, but now we're getting a new segment.
                        bOverlapEdge.Add(tempPath.ToList());
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
                        if (aOverlapEdge[i].Count <= 1)
                        {
                            aOverlapEdge.RemoveAt(i);
                            changed = true;
                        }
                    }
                }

                changed = true;
                while (changed)
                {
                    changed = false;
                    for (int i = bOverlapEdge.Count; i == -1; i--)
                    {
                        if (bOverlapEdge[i].Count <= 1)
                        {
                            bOverlapEdge.RemoveAt(i);
                            changed = true;
                        }
                    }
                }

                // Walk our edges to figure out the overlap.
                foreach (spaceResult tResult in aOverlapEdge.SelectMany(t => bOverlapEdge.Select(t1 => overlapAssess(overlapShape[poly], t, t1, aPath, bPath, scaleFactorForOperation, subMode, runThreaded, debugCalc)).Where(tResult => result.resultPaths.Count == 0 || tResult.distance > result.distance)))
                {
                    result.distance = tResult.distance;
                    result.resultPaths = tResult.resultPaths;
                }
            }
        }
        catch (Exception e)
        {
            errorRep?.Invoke(e.ToString(), "Oops");
        }
        
        result.done = true;
        if (debugCalc)
        {
            result.distance = -result.distance / scaleFactorForOperation;
        }
        return result;
    }

    private spaceResult overlapAssess(Path overlapPoly, Path aOverlapEdge, Path bOverlapEdge, Path aPath, Path bPath, int scaleFactorForOperation, int subMode, bool runThreaded, bool debugCalc)
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
        Path extractedPath = new();
        bool usingAEdge = false;

        /* Prior to 1.7, we used the edge from the layerB combination to raycast. This was not ideal.
            We should have used the shortest edge. 1.7 adds an option to make this the behavior, and the old 
            approach is retained for now as well.
        */

        Point64 shortestPathBeforeStartPoint = new(0, 0);
        Point64 shortestPathAfterEndPoint = new(0, 0);

        if (subMode == (int)spacingCalcModes.spacing)
        {
            // Find the shortest edge length and use that for the projection reference
            // Calculate lengths and check.
            for (int pt = 0; pt < aOverlapEdge.Count - 1; pt++)
            {
                lengthA += Math.Sqrt(Utils.myPow(aOverlapEdge[pt + 1].X - aOverlapEdge[pt].X, 2) + Utils.myPow(aOverlapEdge[pt + 1].Y - aOverlapEdge[pt].Y, 2));
            }
            for (int pt = 0; pt < bOverlapEdge.Count - 1; pt++)
            {
                lengthB += Math.Sqrt(Utils.myPow(bOverlapEdge[pt + 1].X - bOverlapEdge[pt].X, 2) + Utils.myPow(bOverlapEdge[pt + 1].Y - bOverlapEdge[pt].Y, 2));
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

                shortestPathBeforeStartPoint = new Point64(aPath[startIndex]);
                shortestPathAfterEndPoint = new Point64(aPath[endIndex]);
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

                shortestPathBeforeStartPoint = new Point64(bPath[startIndex]);
                shortestPathAfterEndPoint = new Point64(bPath[endIndex]);
            }
        }

        if (subMode == (int)spacingCalcModes.spacingOld || subMode == (int)spacingCalcModes.enclosure || subMode == (int)spacingCalcModes.enclosureOld)
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
                shortestPathBeforeStartPoint = new Point64(bPath[startIndex]);
            }
            catch (Exception)
            {
                errorRep?.Invoke("distanceHandler: shortestPathBeforeStartPoint failed.", "Oops");
            }
            try
            {
                shortestPathAfterEndPoint = new Point64(bPath[endIndex]);
            }
            catch (Exception)
            {
                errorRep?.Invoke("distanceHandler: shortestPathAfterEndPoint failed.", "Oops");
            }
        }

        bool extrOrient = Clipper.IsPositive(extractedPath);
        bool overOrient = Clipper.IsPositive(overlapPoly);

        // No blurry rays, so no point running the inner loop threaded. We thread the outer loop (the emission edge raycast), though. Testing showed small performance improvement for this approach.
        RayCast rc = new(extractedPath, overlapPoly, scaleFactorForOperation * scaleFactorForOperation, runThreaded, invert, 0, true, false, shortestPathBeforeStartPoint, shortestPathAfterEndPoint);

        if (debug || debugCalc)
        {
            result.done = true;
            result.resultPaths = rc.getRays();
            result.distance = -1;
            return result;
        }

        Paths clippedLines = rc.getClippedRays();

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
            foreach (Point64 t in aOverlapEdge)
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
            foreach (Point64 t in bOverlapEdge)
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

            bool ortho = clippedLines[line][0].X == clippedLines[line][1].X || clippedLines[line][0].Y == clippedLines[line][1].Y;

            double threshold = 1500; // arbitrary, dialed in by hand.

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
                result.resultPaths.Add(new()
                {
                    new Point64(clippedLines[line][0]),
                    new Point64(clippedLines[line][0]),
                    new Point64(clippedLines[line][clippedLines[line].Count - 1])
                });
            }
        }

        try
        {
            if (clippedLines.Any() && !validOverlapFound || maxDistance_orthogonalFallback > maxDistance)
            {
                // Couldn't find a valid overlap so assume the orthogonal fallback is needed.
                maxDistance = maxDistance_orthogonalFallback;
                result.resultPaths.Clear();
                result.resultPaths.Add(new()
                {
                    new Point64(clippedLines[maxDistance_fallbackIndex][0]),
                    new Point64(clippedLines[maxDistance_fallbackIndex][0]),
                    new Point64(
                        clippedLines[maxDistance_fallbackIndex]
                            [clippedLines[maxDistance_fallbackIndex].Count - 1])
                });
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