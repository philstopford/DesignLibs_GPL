using ClipperLib;
using System;
using System.Collections.Generic;
using System.Linq;
using utility;

namespace geoWrangler
{
    using Path = List<IntPoint>;
    using Paths = List<List<IntPoint>>;

    public static partial class GeoWrangler
    {
        // Sizing is used to define the keyhole width (default) and will be used for the sliver/gap removal.
        // Use of a custom value will cause headaches.
        static double sizing = 500;
        static double default_nudge = 1.03;

        public static Paths makeKeyHole(Paths outers, Paths cutters, double customSizing = 0, double extension = 0, double angularTolerance = 0)
        {
            return pMakeKeyHole(outers, cutters, customSizing: customSizing, extension: extension, angularTolerance: angularTolerance);
        }

        public static Paths makeKeyHole(Path source, double customSizing = 0, double extension = 0, double angularTolerance = 0)
        {
            return pMakeKeyHole(new Paths() { source }, customSizing, extension, angularTolerance);
        }

        public static Paths makeKeyHole(Paths source, double customSizing = 0, double extension = 0, double angularTolerance = 0)
        {
            return pMakeKeyHole(source, customSizing: customSizing, extension: extension, angularTolerance: angularTolerance);
        }

        static Paths pMakeKeyHole(Paths source, double customSizing = 0, double extension = 0, double angularTolerance = 0)
        {
            if (source.Count < 1)
            {
                return source;
            }

            if (customSizing == 0)
            {
                customSizing = sizing;
            }

            // Limit the offset used in the removal otherwise we can cause self-intersections that
            // result in lots of artifacts and trouble.
            Paths input = pRemoveFragments(source, customSizing, extension);

            if (input.Count < 2)
            {
                return pClose(input);
            }

            Paths[] decomp = pGetDecomposed(input);
            for (int i = 0; i < decomp.Length; i++)
            {
                decomp[i] = pClose(decomp[i]);
            }
            Paths ret = pMakeKeyHole(decomp[(int)type.outer], decomp[(int)type.cutter], customSizing, extension, angularTolerance);

            if (ret.Count == 0)
            {
                return source; // something blew up. Send back the original geometry.
            }

            ret = pClose(ret);

            return ret;
        }

        static Paths pMakeKeyHole(Paths outers, Paths cutters, double customSizing = 0, double extension = 0, double angularTolerance = 0)
        {
            if (customSizing == 0)
            {
                customSizing = sizing;
            }
            try
            {
                bool outerOrient = Clipper.Orientation(outers[0]);
                // Use raycaster to project from holes to outer, to try and find a keyhole path that is minimal length, and ideally orthogonal.
                for (int hole = 0; hole < cutters.Count; hole++)
                {
                    Paths extraCutters = new Paths();
                    Path projCheck = pStripTerminators(cutters[hole], true);
                    projCheck = pStripColinear(projCheck);
                    //  Strip the terminator again to meet the requirements below.
                    projCheck = pStripTerminators(projCheck, false);
                    projCheck = pClockwise(projCheck);

                    bool projectCorners = orthogonal(projCheck, angularTolerance);
                    RayCast rc = new RayCast(cutters[hole], outers, 1000000, invert: true, projectCorners: projectCorners);
                    Paths clipped = rc.getClippedRays();
                    // Need to find minimal length, ideally orthogonal.
                    double minLength = -1;
                    bool minLength_ortho = false;
                    int cutPathIndex = -1;
                    for (int r = 0; r < clipped.Count; r++)
                    {
                        bool ray_isOrtho = (clipped[r][0].X == clipped[r][1].X) || (clipped[r][0].Y == clipped[r][1].Y);

                        // Don't replace an orthogonal ray with a non-orthogonal ray.
                        if (minLength_ortho && !ray_isOrtho)
                        {
                            continue;
                        }

                        double ray_length = Math.Abs(distanceBetweenPoints(clipped[r][0], clipped[r][1]));

                        if (ray_length < double.Epsilon)
                        {
                            continue;
                        }

                        minLength_ortho = minLength_ortho || ray_isOrtho;

                        // First ray or a smaller distance causes us to make this the keyhole edge.
                        if ((r == 0) || (ray_length < minLength))
                        {
                            cutPathIndex = r;
                            minLength = ray_length;
                        }
                    }

                    if (cutPathIndex != -1)
                    {
                        // Offset our cutter and assign to the clipping scenario.
                        Paths sPaths = pInflateEdge(clipped[cutPathIndex], customSizing);

                        extraCutters.AddRange(new Paths(sPaths));
                    }
                    else
                    {
                        int xxx = 2;
                    }

                    Clipper c = new Clipper();
                    c.AddPaths(cutters, PolyType.ptSubject, true);
                    c.AddPaths(extraCutters, PolyType.ptClip, true);

                    Paths mergedCutters = new Paths();
                    c.Execute(ClipType.ctUnion, mergedCutters);

                    c.Clear();
                    c.AddPaths(outers, PolyType.ptSubject, true);
                    c.AddPaths(mergedCutters, PolyType.ptClip, true);

                    // Reduce our geometry back to the simplest form.
                    Paths new_outers = new Paths();
                    c.Execute(ClipType.ctDifference, new_outers);//, PolyFillType.pftNonZero, PolyFillType.pftNegative);

                    outers.Clear();
                    for (int p = 0; p < new_outers.Count; p++)
                    {
                        if (Clipper.Orientation(new_outers[p]) == outerOrient)
                        {
                            outers.Add(new_outers[p]);
                        }
                    }
                }

                return pClockwiseAndReorder(outers);
            }
            catch (Exception)
            {
                throw new Exception("pMakeKeyHole error");
            }
        }

        static Paths pInflateEdge(Path edge, double customSizing)
        {
            if (customSizing == 0)
            {
                customSizing = sizing;
            }
            // Force clockwise, which should get us something consistent to work with.

            double dTmp0 = pDistanceBetweenPoints(new IntPoint(0, 0), edge[0]);
            double dTmp1 = pDistanceBetweenPoints(new IntPoint(0, 0), edge[1]);

            if (dTmp1 < dTmp0)
            {
                edge.Reverse();
            }

            // Get sorted out for dx, dy and normalization.
            double dx = edge[0].X - edge[1].X;
            double dy = edge[0].Y - edge[1].Y;

            double length = Math.Sqrt(Utils.myPow(dx, 2) + Utils.myPow(dy, 2));

            dx /= length;
            dy /= length;

            // Extend the line slightly.
            edge[0] = new IntPoint((Int64)(edge[0].X - Math.Abs(dx * sizing)), (Int64)(edge[0].Y - Math.Abs(dy * sizing)));
            edge[1] = new IntPoint((Int64)(edge[1].X + Math.Abs(dx * sizing)), (Int64)(edge[1].Y + Math.Abs(dy * sizing)));

            ClipperOffset co = new ClipperOffset();
            co.PreserveCollinear = true;
            co.AddPath(edge, JoinType.jtMiter, EndType.etOpenSquare);
            PolyTree solution = new PolyTree();
            co.Execute(ref solution, customSizing);
            Paths sPaths = Clipper.ClosedPathsFromPolyTree(solution);

            return sPaths;
        }

        public static Paths sliverGapRemoval(Path source, double customSizing = 0, double extension = 0, bool maySimplify = false, bool doSomething = true)
        {
            return pSliverGapRemoval(new Paths() { source }, customSizing, extension, maySimplify: maySimplify, doSomething: doSomething);
        }

        public static Paths sliverGapRemoval(Paths source, double customSizing = 0, double extension = 0, bool maySimplify = false, bool doSomething = true)
        {
            return pSliverGapRemoval(source, customSizing, extension, maySimplify: maySimplify, doSomething: doSomething);
        }

        static Paths pSliverGapRemoval(Paths source, double customSizing, double extension, bool maySimplify, bool doSomething)
        {
            if (customSizing == 0)
            {
                customSizing = sizing;
            }
            // Remove gaps, then remove slivers. Same process, different direction for sizing.
            Paths ret = pRemoveFragments(pRemoveFragments(source, customSizing, extension, maySimplify), -customSizing, extension, maySimplify: maySimplify, doSomething: doSomething);

            if (ret.Count == 0)
            {
                // Probably not what we intended, so return the original as we assume the incoming geometry got crushed out of existence.
                return source;
            }

            return ret;
        }

        public static Paths gapRemoval(Path source, double customSizing = 0, double extension = 0, bool maySimplify = false, bool doSomething = true)
        {
            return gapRemoval(new Paths() { source }, customSizing, extension, maySimplify, doSomething);
        }

        public static Paths gapRemoval(Paths source, double customSizing = 0, double extension = 0, bool maySimplify = false, bool doSomething = true)
        {
            if (source.Count < 1)
            {
                return source;
            }

            bool orig_orient_gw = GeoWrangler.isClockwise(source[0]);
            bool orig_orient_c = Clipper.Orientation(source[0]);

            Paths ret = pRemoveFragments(source, customSizing, extension, maySimplify, doSomething: doSomething);

            if (ret.Count == 0)
            {
                return source; // something blew up. Send back the original geometry.
            }

            // Clean-up the geometry.
            ret = Clipper.SimplifyPolygons(ret);

            // Validate orientations.
            bool gR_orient_gw = GeoWrangler.isClockwise(ret[0]);
            bool gR_orient_c = Clipper.Orientation(ret[0]);

            bool reverseNeeded = (gR_orient_gw != orig_orient_gw) || (gR_orient_c != orig_orient_c);

            // Respin if needed. 1st path is always an outer, so we need to ensure the orientation is consistent.
            if (reverseNeeded)
            {
                for (int i = 0; i < ret.Count; i++)
                {
                    ret[i].Reverse();
                }
            }

            // Clean-up stripped colinear vertices; need to refragment.
            ret = pClose(ret);

            return ret;
        }

        public static Paths sliverRemoval(Paths source, double customSizing = 0, double extension = 0, bool maySimplify = false, bool doSomething = true)
        {
            if (customSizing == 0)
            {
                customSizing = sizing;
            }
            double oArea = 0;
            for (int i = 0; i < source.Count; i++)
            {
                oArea += Clipper.Area(source[i]);
            }
            Paths ret = pRemoveFragments(source, -customSizing, extension, maySimplify: maySimplify, doSomething: doSomething);
            double nArea = 0;
            for (int i = 0; i < ret.Count; i++)
            {
                nArea += Clipper.Area(ret[i]);
            }

            if ((Math.Abs(oArea) - Math.Abs(nArea)) < 10000)
            {
                return source;
            }

            return ret;
        }

        // Positive incoming value removes gaps (keyholes); negative incoming value will remove slivers.
        static Paths pRemoveFragments(Paths source, double customSizing, double extension, bool maySimplify = false, JoinType joinType = JoinType.jtMiter, bool doSomething = true)
        {
            if (customSizing == 0)
            {
                customSizing = sizing;
            }

            if (extension == 0)
            {
                extension = default_nudge;
            }

            // Used to try and avoid residual fragments; empirically derived.
            customSizing *= extension;

            Paths cGeometry = new Paths();

            ClipperOffset co = new ClipperOffset();
            co.PreserveCollinear = !maySimplify;
            co.AddPaths(source, joinType, EndType.etClosedPolygon);
            co.Execute(ref cGeometry, customSizing);
            co.Clear();
            co.AddPaths(cGeometry.ToList(), joinType, EndType.etClosedPolygon);
            cGeometry.Clear();
            co.Execute(ref cGeometry, -customSizing); // Size back to original dimensions

            double newArea = 0;
            for (int i = 0; i < cGeometry.Count; i++)
            {
                newArea += Clipper.Area(cGeometry[i]);
            }

            if (Math.Abs(newArea) < double.Epsilon)
            {
                // We crushed our geometry, it seems.
                // This was probably not the plan, so send back the original geometry instead.
                return source;
            }

            return cGeometry;
        }

        static Paths pRemoveFragments(Path source, double customSizing, bool maySimplify = false, JoinType joinType = JoinType.jtMiter, bool doSomething = true)
        {
            Paths cGeometry = new Paths();
            if (!doSomething)
            {
                cGeometry.Add(source);
                return cGeometry;
            }
            else
            {
                if (customSizing == 0)
                {
                    customSizing = sizing;
                }
                double sourceArea = Clipper.Area(source);

                ClipperOffset co = new ClipperOffset();
                co.PreserveCollinear = !maySimplify;
                co.AddPath(source, joinType, EndType.etClosedPolygon);
                co.Execute(ref cGeometry, customSizing);
                co.Clear();
                co.AddPaths(cGeometry.ToList(), joinType, EndType.etClosedPolygon);
                cGeometry.Clear();
                co.Execute(ref cGeometry, -customSizing); // Size back to original dimensions

                double newArea = 0;
                for (int i = 0; i < cGeometry.Count; i++)
                {
                    newArea += Clipper.Area(cGeometry[i]);
                }

                if (Math.Abs(newArea) < double.Epsilon)
                {
                    // We crushed our geometry, it seems.
                    // This was probably not the plan, so send back the original geometry instead.
                    cGeometry.Clear();
                    cGeometry.Add(source);
                }
                else
                {
                    // Do we need to flip the direction to match the original orientation?
                    if (((sourceArea < 0) && (newArea > 0)) || ((sourceArea > 0) && (newArea < 0)))
                    {
                        // Multi-path handling gets interesting. The first path is assumed to be the outer. Let's compare that with the original geometry. If the orientation is different, reverse the full set.
                        bool orientation = Clipper.Orientation(source);
                        bool origCG0_o = Clipper.Orientation(cGeometry[0]);
                        for (int i = 0; i < cGeometry.Count; i++)
                        {
                            if (origCG0_o != orientation)
                            {
                                cGeometry[i].Reverse();
                            }
                        }
                    }
                }

                return cGeometry;
            }
        }
    }
}
