using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Clipper2Lib;
using geoLib;
using utility;

namespace geoWrangler;

using Path = List<Point64>;
using Paths = List<List<Point64>>;

public static class Proximity
{
    // Drawn poly allows for geometry to be excluded from consideration in the input list
    public static GeometryResult proximityBias(List<GeoLibPointF[]> input, List<bool> drawnPoly_, decimal pBias, decimal pBiasDist, int proxRays, int proxSideRaysFallOff, decimal proxSideRaysMultiplier, decimal rayExtension, double fragmenterResolution, Int64 scaleFactorForOperation)
    {
        // Proximity biasing - where isolated edges get bias based on distance to nearest supporting edge.
        bool proxBiasNeeded = (pBias != 0) && (pBiasDist != 0);

        if (!proxBiasNeeded)
        {
            return new() {geometry = input.ToList(), drawn = drawnPoly_.ToList()};
        }
        
        bool debug = false;
        bool linear = false;

        List<GeoLibPointF[]> preOverlapMergePolys = new();

        Paths dRays = new();

        // Scale up our geometry for processing. Force a clockwise point order here due to potential upstream point order changes (e.g. polygon merging)
        Paths sourceGeometry = GeoWrangler.pathsFromPointFs(input, scaleFactorForOperation);

        List<bool> overlapDrawnList = new();
        
        int sCount = sourceGeometry.Count;
        for (int poly = 0; poly < sCount; poly++)
        {
            if (sourceGeometry[poly].Count <= 1 || drawnPoly_[poly])
            {
                // Nothing to do with drawn or zero count entries.
                continue;
            }
            
            overlapDrawnList.Add(false);

            Path sourcePoly = sourceGeometry[poly].ToList();
            Paths collisionGeometry = sourceGeometry.ToList();
            // collisionGeometry.RemoveAt(poly); // Don't actually want to remove the emission as self-aware proximity matters.
            Path deformedPoly = new();

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

            Fragmenter f = new(fragmenterResolution * scaleFactorForOperation);

            sourcePoly = f.fragmentPath(sourcePoly);

            collisionGeometry = f.fragmentPaths(collisionGeometry);

            RayCast rc = new(sourcePoly, collisionGeometry, Convert.ToInt32(pBiasDist * scaleFactorForOperation), false, invert:0, proxRays, emitThread, multiSampleThread, sideRayFallOff: (RayCast.falloff)proxSideRaysFallOff, sideRayFallOffMultiplier: Convert.ToDouble(proxSideRaysMultiplier));

            Paths clippedLines = rc.getClippedRays();
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
                double biasScaling = lineLength / scaleFactorForOperation / Convert.ToDouble(pBiasDist);

                if (biasScaling > 1)
                {
                    biasScaling = 1;
                }

                // Probably should be a sigmoid, but using this for now.
                double displacedAmount;

                // ReSharper disable once ConditionIsAlwaysTrueOrFalse
                if (linear)
                {
                    displacedAmount = biasScaling * Convert.ToDouble(pBias) * scaleFactorForOperation;
                }
                else
                {
                    // Using sine to make a ease-in/ease-out effect.
                    displacedAmount = Math.Sin(Utils.toRadians(biasScaling * 90.0f)) * Convert.ToDouble(pBias) * scaleFactorForOperation;
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
            preOverlapMergePolys.Add(GeoWrangler.pointFFromPath(deformedPoly, scaleFactorForOperation));
            deformedPoly.Add(new Point64(deformedPoly[0]));
        }

        // Check for overlaps and process as needed post-biasing.
        GeometryResult ret = ProcessOverlaps.processOverlaps(preOverlapMergePolys, overlapDrawnList, extension:Convert.ToDouble(rayExtension), resolution:fragmenterResolution, scaleFactorForOperation:scaleFactorForOperation, forceOverride: false);

        // ReSharper disable once ConditionIsAlwaysTrueOrFalse
        if (debug)
        {
            foreach (Path t in dRays)
            {
                ret.geometry.Add(GeoWrangler.pointFFromPath(t, scaleFactorForOperation));
                ret.drawn.Add(true);
            }
        }

        return ret;
    }
    
}