using Clipper2Lib;
using geoWrangler;
using utility;
#pragma warning disable CS8618

namespace geoAnalysis;

/// <summary>
/// Analyzes the minimum intersection angle between two sets of polygonal paths.
/// This class computes geometric intersections and determines the smallest angle
/// formed between intersecting polygon edges.
/// </summary>
public class angleHandler
{
    /// <summary>
    /// Gets the minimum intersection angle in degrees between the analyzed polygon sets.
    /// Returns 180.0 if no intersections are found.
    /// </summary>
    public double minimumIntersectionAngle { get; private set; }
    
    /// <summary>
    /// Internal storage for intersection polygon results from boolean operations.
    /// </summary>
    private PathsD listOfOutputPoints;
    
    /// <summary>
    /// Gets the result paths highlighting the location of the minimum intersection angle.
    /// Contains visualization geometry for the minimum angle intersection point.
    /// </summary>
    public PathsD resultPaths { get; private set; } // will only have one path, for minimum angle.

    /// <summary>
    /// Callback function for Z-fill operations during polygon intersection.
    /// Tags intersection points with a specific Z value for identification.
    /// </summary>
    /// <param name="bot1">Bottom point of first edge</param>
    /// <param name="top1">Top point of first edge</param>
    /// <param name="bot2">Bottom point of second edge</param>
    /// <param name="top2">Top point of second edge</param>
    /// <param name="pt">Reference to the intersection point to be tagged</param>
    private static void ZFillCallback(PointD bot1, PointD top1, PointD bot2, PointD top2, ref PointD pt)
    {
        pt.z = -1; // Tag our intersection points.
    }

    /// <summary>
    /// Minimum distance threshold for intersection marker scaling.
    /// Used to ensure visibility of intersection points in visualization.
    /// </summary>
    private const double minDistance = 10.0;

    /// <summary>
    /// Initializes a new instance of the angleHandler class and performs intersection angle analysis.
    /// </summary>
    /// <param name="layerAPath">First set of polygonal paths for intersection analysis</param>
    /// <param name="layerBPath">Second set of polygonal paths for intersection analysis</param>
    public angleHandler(PathsD layerAPath, PathsD layerBPath)
    {
        angleHandlerLogic(layerAPath, layerBPath);
    }

    /// <summary>
    /// Core logic for computing minimum intersection angles between two polygon sets.
    /// Performs boolean intersection, analyzes edge angles, and creates visualization geometry.
    /// </summary>
    /// <param name="layerAPath">First polygon set</param>
    /// <param name="layerBPath">Second polygon set</param>
    private void angleHandlerLogic(PathsD layerAPath, PathsD layerBPath)
    {
        listOfOutputPoints = [];
        resultPaths = [];
        PathD resultPath = [];
        ClipperD c = new(Constants.roundingDecimalPrecision) {ZCallback = ZFillCallback};
        c.AddSubject(layerAPath);
        c.AddClip(layerBPath);
        
        // Boolean AND of the two levels for the area operation.
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, listOfOutputPoints);

        listOfOutputPoints = GeoWrangler.reOrderXY(listOfOutputPoints);

        // Set initial output value for the case there are no intersections
        minimumIntersectionAngle = 180.0; // no intersection angle.

        double tmpVal = listOfOutputPoints.Sum(Clipper.Area);
        if (tmpVal == 0.0)
        {
            // No overlap
            // Set output path and avoid heavy lifting
            resultPath.Add(new PointD(0, 0));
            resultPaths.Add(resultPath);
        }
        else
        {
            double temporaryResult = 180.0;
            PathD temporaryPath = [new PointD(0, 0), new PointD(0, 0), new PointD(0, 0)];
            foreach (PathD t in listOfOutputPoints)
            {
                PathD overlapPath = GeoWrangler.clockwise(t);

                int pt = 0;
                while (pt < overlapPath.Count)
                {
                    if (overlapPath[pt].z == -1)
                    {
                        // intersection point found - let's get our three points to find the angle.
                        // http://en.wikipedia.org/wiki/Law_of_cosines
                        PointD interSection_B;
                        PointD interSection_C;
                        PointD interSection_A;
                        if (pt == 0)
                        {
                            // Find preceding not-identical point.
                            int refPt = overlapPath.Count - 1;
                            while (Math.Abs(GeoWrangler.distanceBetweenPoints(overlapPath[refPt], overlapPath[pt])) == 0)
                            {
                                refPt--;
                                if (refPt == 0)
                                {
                                    break;
                                }
                            }
                            interSection_B = overlapPath[refPt]; // map to last point
                            interSection_C = overlapPath[pt];
                            // Find following not-identical point.
                            refPt = 0;
                            while (Math.Abs(GeoWrangler.distanceBetweenPoints(overlapPath[refPt], overlapPath[pt])) == 0)
                            {
                                refPt++;
                                if (refPt == overlapPath.Count - 1)
                                {
                                    break;
                                }
                            }
                            interSection_A = overlapPath[refPt];
                        }
                        else if (pt == overlapPath.Count - 1) // last point in the list
                        {
                            // Find preceding not-identical point.
                            int refPt = pt;
                            while (Math.Abs(GeoWrangler.distanceBetweenPoints(overlapPath[refPt], overlapPath[pt])) == 0)
                            {
                                refPt--;
                                if (refPt == 0)
                                {
                                    break;
                                }
                            }
                            interSection_B = overlapPath[refPt];
                            interSection_C = overlapPath[pt];
                            // Find following not-identical point.
                            refPt = 0;
                            while (Math.Abs(GeoWrangler.distanceBetweenPoints(overlapPath[refPt], overlapPath[pt])) == 0)
                            {
                                refPt++;
                                if (refPt == overlapPath.Count - 1)
                                {
                                    break;
                                }
                            }
                            interSection_A = overlapPath[0]; // map to the first point
                        }
                        else
                        {
                            // Find preceding not-identical point.
                            int refPt = pt;
                            while (Math.Abs(GeoWrangler.distanceBetweenPoints(overlapPath[refPt], overlapPath[pt])) == 0)
                            {
                                refPt--;
                                if (refPt == 0)
                                {
                                    break;
                                }
                            }
                            interSection_B = overlapPath[refPt];
                            interSection_C = overlapPath[pt];
                            // Find following not-identical point.
                            refPt = pt;
                            while (Math.Abs(GeoWrangler.distanceBetweenPoints(overlapPath[refPt], overlapPath[pt])) == 0)
                            {
                                refPt++;
                                if (refPt == overlapPath.Count - 1)
                                {
                                    break;
                                }
                            }
                            interSection_A = overlapPath[refPt];
                        }

                        PointD cBVector = new(interSection_B.x - interSection_C.x, interSection_B.y - interSection_C.y);
                        PointD cAVector = new(interSection_A.x - interSection_C.x, interSection_A.y - interSection_C.y);

                        double xComponents = cBVector.x * cAVector.x;
                        double yComponents = cBVector.y * cAVector.y;

                        double scalarProduct = xComponents + yComponents;

                        double cBMagnitude = Math.Sqrt(Utils.myPow(cBVector.x, 2) + Utils.myPow(cBVector.y, 2));
                        double cAMagnitude = Math.Sqrt(Utils.myPow(cAVector.x, 2) + Utils.myPow(cAVector.y, 2));

                        double theta = Math.Abs(Utils.toDegrees(Math.Acos(scalarProduct / (cBMagnitude * cAMagnitude)))); // Avoid falling into a trap with negative angles.

                        if (theta < temporaryResult)
                        {
                            temporaryResult = theta;
                            temporaryPath.Clear();
                            temporaryPath.Add(new PointD(interSection_A.x, interSection_A.y));
                            temporaryPath.Add(new PointD(interSection_C.x, interSection_C.y));
                            temporaryPath.Add(new PointD(interSection_B.x, interSection_B.y));
                        }
                    }
                    pt++;
                }
            }
            minimumIntersectionAngle = temporaryResult;

            // Check our temporary path to see if we need to scale it up.
            double distance = GeoWrangler.distanceBetweenPoints(temporaryPath[0], temporaryPath[1]);
            PointD distancePoint64 = GeoWrangler.PointD_distanceBetweenPoints(temporaryPath[0], temporaryPath[1]); // A to C
            if (distance < minDistance)
            {
                double X = temporaryPath[0].x;
                double Y = temporaryPath[0].y;
                if (Math.Abs(temporaryPath[1].x - temporaryPath[0].x) > ConstantsGA.tolerance)
                {
                    X = temporaryPath[1].x + distancePoint64.x * (minDistance / distance);
                }
                if (Math.Abs(temporaryPath[1].y - temporaryPath[0].y) > ConstantsGA.tolerance)
                {
                    Y = temporaryPath[1].y + distancePoint64.y * (minDistance / distance);
                }
                temporaryPath[0] = new PointD(X, Y);
            }
            distance = GeoWrangler.distanceBetweenPoints(temporaryPath[2], temporaryPath[1]);
            distancePoint64 = GeoWrangler.PointD_distanceBetweenPoints(temporaryPath[2], temporaryPath[1]); // B to C
            if (distance < minDistance)
            {
                double X = temporaryPath[2].x;
                double Y = temporaryPath[2].y;
                if (Math.Abs(temporaryPath[1].y - temporaryPath[2].y) > ConstantsGA.tolerance)
                {
                    Y = temporaryPath[1].y + distancePoint64.y * (minDistance / distance);
                }
                if (Math.Abs(temporaryPath[1].x - temporaryPath[2].x) > ConstantsGA.tolerance)
                {
                    X = temporaryPath[1].x + distancePoint64.x * (minDistance / distance);
                }
                temporaryPath[2] = new PointD(X, Y);
            }
            resultPaths.Add(temporaryPath);
        }
    }
}