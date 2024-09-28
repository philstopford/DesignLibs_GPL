using Clipper2Lib;
using geoWrangler;
using utility;
#pragma warning disable CS8618

namespace geoAnalysis;

public class angleHandler
{
    public double minimumIntersectionAngle { get; private set; }
    private PathsD listOfOutputPoints;
    public PathsD resultPaths { get; private set; } // will only have one path, for minimum angle.

    private static void ZFillCallback(PointD bot1, PointD top1, PointD bot2, PointD top2, ref PointD pt)
    {
        pt.z = -1; // Tag our intersection points.
    }

    // Distance functions to drive scale-up of intersection marker if needed.
    private const double minDistance = 10.0;

    public angleHandler(PathsD layerAPath, PathsD layerBPath)
    {
        angleHandlerLogic(layerAPath, layerBPath);
    }

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