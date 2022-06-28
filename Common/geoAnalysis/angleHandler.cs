using Clipper2Lib;
using geoWrangler;
using utility;

namespace geoAnalysis;

using Path = List<Point64>;
using Paths = List<List<Point64>>;

public class angleHandler
{
    public double minimumIntersectionAngle { get; private set; }
    private Paths listOfOutputPoints;
    public Paths resultPaths { get; private set; } // will only have one path, for minimum angle.

    private void ZFillCallback(Point64 bot1, Point64 top1, Point64 bot2, Point64 top2, ref Point64 pt)
    {
        pt.Z = -1; // Tag our intersection points.
    }

    // Distance functions to drive scale-up of intersection marker if needed.
    private double minDistance = 10.0;

    public angleHandler(Paths layerAPath, Paths layerBPath, int scaleFactorForOperation)
    {
        angleHandlerLogic(layerAPath, layerBPath, scaleFactorForOperation);
    }

    private void angleHandlerLogic(Paths layerAPath, Paths layerBPath, int scaleFactorForOperation)
    {
        listOfOutputPoints = new Paths();
        resultPaths = new Paths();
        Path resultPath = new();
        Clipper64 c = new() {ZFillFunc = ZFillCallback};
        c.AddSubject(layerAPath);
        c.AddClip(layerBPath);
        
        // Boolean AND of the two levels for the area operation.
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, listOfOutputPoints);

        listOfOutputPoints = GeoWrangler.reOrderXY(listOfOutputPoints);

        // Set initial output value for the case there are no intersections
        minimumIntersectionAngle = 180.0; // no intersection angle.

        double tmpVal = listOfOutputPoints.Sum(t => Clipper.Area(t));
        if (tmpVal == 0.0)
        {
            // No overlap
            // Set output path and avoid heavy lifting
            resultPath.Add(new Point64(0, 0));
            resultPaths.Add(resultPath);
        }
        else
        {
            double temporaryResult = 180.0;
            Path temporaryPath = new() {new Point64(0, 0), new Point64(0, 0), new Point64(0, 0)};
            foreach (Path t in listOfOutputPoints)
            {
                Path overlapPath = GeoWrangler.clockwise(t);

                int pt = 0;
                while (pt < overlapPath.Count)
                {
                    if (overlapPath[pt].Z == -1)
                    {
                        // intersection point found - let's get our three points to find the angle.
                        // http://en.wikipedia.org/wiki/Law_of_cosines
                        Point64 interSection_B;
                        Point64 interSection_C;
                        Point64 interSection_A;
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

                        Point64 cBVector = new(interSection_B.X - interSection_C.X, interSection_B.Y - interSection_C.Y);
                        Point64 cAVector = new(interSection_A.X - interSection_C.X, interSection_A.Y - interSection_C.Y);

                        long xComponents = cBVector.X * cAVector.X;
                        long yComponents = cBVector.Y * cAVector.Y;

                        long scalarProduct = xComponents + yComponents;

                        double cBMagnitude = Math.Sqrt(Utils.myPow(cBVector.X, 2) + Utils.myPow(cBVector.Y, 2));
                        double cAMagnitude = Math.Sqrt(Utils.myPow(cAVector.X, 2) + Utils.myPow(cAVector.Y, 2));

                        double theta = Math.Abs(Utils.toDegrees(Math.Acos(scalarProduct / (cBMagnitude * cAMagnitude)))); // Avoid falling into a trap with negative angles.

                        if (theta < temporaryResult)
                        {
                            temporaryResult = theta;
                            temporaryPath.Clear();
                            temporaryPath.Add(new Point64(interSection_A.X, interSection_A.Y));
                            temporaryPath.Add(new Point64(interSection_C.X, interSection_C.Y));
                            temporaryPath.Add(new Point64(interSection_B.X, interSection_B.Y));
                        }
                    }
                    pt++;
                }
            }
            minimumIntersectionAngle = temporaryResult;

            // Check our temporary path to see if we need to scale it up.
            double distance = GeoWrangler.distanceBetweenPoints(temporaryPath[0], temporaryPath[1]) / scaleFactorForOperation;
            Point64 distancePoint64 = GeoWrangler.Point64_distanceBetweenPoints(temporaryPath[0], temporaryPath[1]); // A to C
            if (distance < minDistance)
            {
                double X = temporaryPath[0].X;
                double Y = temporaryPath[0].Y;
                if (temporaryPath[1].X != temporaryPath[0].X)
                {
                    X = temporaryPath[1].X + distancePoint64.X * (minDistance / distance);
                }
                if (temporaryPath[1].Y != temporaryPath[0].Y)
                {
                    Y = temporaryPath[1].Y + distancePoint64.Y * (minDistance / distance);
                }
                temporaryPath[0] = new Point64((long)X, (long)Y);
            }
            distance = GeoWrangler.distanceBetweenPoints(temporaryPath[2], temporaryPath[1]) / scaleFactorForOperation;
            distancePoint64 = GeoWrangler.Point64_distanceBetweenPoints(temporaryPath[2], temporaryPath[1]); // B to C
            if (distance < minDistance)
            {
                double X = temporaryPath[2].X;
                double Y = temporaryPath[2].Y;
                if (temporaryPath[1].Y != temporaryPath[2].Y)
                {
                    Y = temporaryPath[1].Y + distancePoint64.Y * (minDistance / distance);
                }
                if (temporaryPath[1].X != temporaryPath[2].X)
                {
                    X = temporaryPath[1].X + distancePoint64.X * (minDistance / distance);
                }
                temporaryPath[2] = new Point64((long)X, (long)Y);
            }
            resultPaths.Add(temporaryPath);
        }
    }
}