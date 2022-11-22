using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Clipper2Lib;
using geoLib;
using utility;

namespace geoWrangler;

public static class Proximity
{
    // Drawn poly allows for geometry to be excluded from consideration in the input list
    public static GeometryResult proximityBias(PathsD input, List<bool> drawnPoly_, decimal pBias,
        decimal pBiasDist, int proxRays, int proxSideRaysFallOff, decimal proxSideRaysMultiplier, decimal rayExtension,
        double fragmenterResolution, bool doCleanUp, double cleanUpEpsilon)
    {
        // Proximity biasing - where isolated edges get bias based on distance to nearest supporting edge.
        bool proxBiasNeeded = (pBias != 0) && (pBiasDist != 0);

        if (!proxBiasNeeded)
        {
            return new() { geometry = new (input), drawn = drawnPoly_.ToList() };
        }

        bool debug = false;
        bool linear = false;

        PathsD preOverlapMergePolys = new();

        Paths64 dRays = new();

        // Scale up our geometry for processing. Force a clockwise point order here due to potential upstream point order changes (e.g. polygon merging)
        Paths64 sourceGeometry = GeoWrangler.paths64FromPathsD(input);

        List<bool> overlapDrawnList = new();

        // Ensure geometry meets our needs. This is important for the normal computation and calculations.
        // If this is not done, sawtooth profiles are seen due to problematic ray casts.
        for (int p1 = 0; p1 < sourceGeometry.Count; p1++)
        {
            if (Clipper.IsPositive(sourceGeometry[p1]))
            {
                sourceGeometry[p1].Reverse();
            }
        }
        sourceGeometry = GeoWrangler.close(sourceGeometry);

        int sCount = sourceGeometry.Count;

        Fragmenter f = new(fragmenterResolution);

        for (int poly = 0; poly < sCount; poly++)
        {
            if (sourceGeometry[poly].Count <= 1 || drawnPoly_[poly])
            {
                // Nothing to do with drawn or zero count entries.
                continue;
            }
            
            overlapDrawnList.Add(false);

            Path64 sourcePoly = new(sourceGeometry[poly]);
            Paths64 collisionGeometry = new(sourceGeometry);
            // collisionGeometry.RemoveAt(poly); // Don't actually want to remove the emission as self-aware proximity matters.
            Path64 deformedPoly = new();

            // Threading operation here gets more tricky than the distance handler. We have a less clear trade off of threading based on the emission edge (the polygon being biased) vs the multisampling emission.
            // In batch calculation mode, this tradeoff gets more awkward.
            // Threading both options also causes major performance degradation as far too many threads are spawned for the host system.
            bool multiSampleThread = false;
            bool emitThread = false;

            if (proxRays > 1)
            {
                multiSampleThread = true;
                // for multipolygon scenarios, avoid threading the multisampling and instead favor threading emitting edge.
                if (sourceGeometry.Count > 1)
                {
                    emitThread = true;
                    multiSampleThread = false;
                }
            }
            else
            {
                emitThread = true;
            }
            
            sourcePoly = GeoWrangler.stripColinear(sourcePoly);
            sourcePoly = GeoWrangler.removeDuplicates(sourcePoly);
            sourcePoly = GeoWrangler.close(sourcePoly);
            sourcePoly = f.fragmentPath(sourcePoly);

            collisionGeometry = f.fragmentPaths(collisionGeometry);

            RayCast rc = new(sourcePoly, collisionGeometry, Convert.ToInt32(pBiasDist), false, invert:0, proxRays, emitThread, multiSampleThread, sideRayFallOff: (RayCast.falloff)proxSideRaysFallOff, sideRayFallOffMultiplier: Convert.ToDouble(proxSideRaysMultiplier));

            Paths64 clippedLines = rc.getClippedRays();
            // ReSharper disable once ConditionIsAlwaysTrueOrFalse
            if (debug)
            {
                dRays.AddRange(clippedLines);
            }

            // We hope to get the same number of clipped lines back as the number of points that went in....
            int cLCount = clippedLines.Count;
            for (int line = 0; line < cLCount; line++)
            {
                long displacedX = sourcePoly[line].X;
                long displacedY = sourcePoly[line].Y;

                double lineLength = rc.getRayLength(line);

                switch (lineLength)
                {
                    // No biasing - ray never made it beyond the surface. Short-cut the 
                    case 0:
                        deformedPoly.Add(new Point64(clippedLines[line][0]));
                        continue;
                    case < 0:
                        lineLength *= -1;
                        break;
                }

                // Calculate our bias based on this distance and apply it.
                double biasScaling = lineLength / 1.0 / Convert.ToDouble(pBiasDist);

                if (biasScaling > 1)
                {
                    biasScaling = 1;
                }

                // Probably should be a sigmoid, but using this for now.
                double displacedAmount;

                // ReSharper disable once ConditionIsAlwaysTrueOrFalse
                if (linear)
                {
                    displacedAmount = biasScaling * Convert.ToDouble(pBias);
                }
                else
                {
                    // Using sine to make a ease-in/ease-out effect.
                    displacedAmount = Math.Sin(Utils.toRadians(biasScaling * 90.0f)) * Convert.ToDouble(pBias);
                }

                // Use our cast ray from rc to get a normalized average 
                Point64 averagedEdgeNormal = GeoWrangler.Point64_distanceBetweenPoints(clippedLines[line][clippedLines[line].Count - 1], clippedLines[line][0]);

                // Normalize our vector length.
                double aX = averagedEdgeNormal.X / lineLength;
                double aY = averagedEdgeNormal.Y / lineLength;

                displacedY += (long)(displacedAmount * aY);
                displacedX += (long)(displacedAmount * aX);

                deformedPoly.Add(new Point64(displacedX, displacedY));
            }
            
            // Experimental clean-up
            if (doCleanUp)
            {
                Path64 rdpPath = Clipper.RamerDouglasPeucker(GeoWrangler.close(deformedPoly), cleanUpEpsilon); // scale from int to make a double
                deformedPoly = f.fragmentPath(GeoWrangler.close(rdpPath));
            }

            // deformedPoly.Add(new Point64(deformedPoly[0]));
            preOverlapMergePolys.Add(GeoWrangler.pathDFromPath64(deformedPoly));
        }

        // Check for overlaps and process as needed post-biasing.
        GeometryResult ret = ProcessOverlaps.processOverlaps(preOverlapMergePolys, overlapDrawnList, extension:Convert.ToDouble(rayExtension), resolution:fragmenterResolution, forceOverride: false);

        // ReSharper disable once ConditionIsAlwaysTrueOrFalse
        if (debug)
        {
            foreach (Path64 t in dRays)
            {
                ret.geometry.Add(GeoWrangler.pathDFromPath64(t));
                ret.drawn.Add(true);
            }
        }

        return ret;
    }
    
}