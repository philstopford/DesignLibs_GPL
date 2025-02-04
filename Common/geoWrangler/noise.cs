using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Clipper2Lib;
using Noise;
using utility;

namespace geoWrangler;

public static class NoiseC
{
    public static readonly List<string> noiseTypes = ["Perlin", "Simplex", "OpenSimplex"];

    public enum noiseIndex { perlin, simplex, opensimplex }

    public static PathsD doNoise(PathsD input, List<bool> drawnPoly, int noiseType, int seed, double freq, double jitterScale)
    {
        PathsD ret = new ();
        // Gets a -1 to +1 noise field. We get a seed from our RNG of choice unless the layer preview mode is set, where a fixed seed is used.
        // Random constants to mitigate continuity effects in noise that cause nodes in the noise across multiple layers, due to periodicity.
        const double x_const = 123489.1928734;
        const double y_const = 891243.0982134;

        object noiseSource = noiseType switch
        {
            (int)noiseIndex.opensimplex => new OpenSimplexNoise(seed),
            (int)noiseIndex.simplex => new SimplexNoise(seed),
            _ => new PerlinNoise(seed)
        };

        // Need to iterate our preview points.
        for (int poly = 0; poly < input.Count; poly++)
        {
            if (input[poly].Count <= 1 || drawnPoly[poly])
            {
                continue;
            }

            // Let's see if this is a keyholed polygon; if so, we need to be a bit clever.
            PathsD geometry = GeoWrangler.sliverGapRemoval(input[poly]);

            foreach (PathD mcPoints in geometry)
            {
                int ptCount = mcPoints.Count;

                // Create our jittered polygon in a new list to avoid breaking normal computation, etc. by modifying the source.
                PathD jitteredPoints = Helper.initedPathD(ptCount);

                // We could probably simply cast rays in the raycaster and use those, but for now reinvent the wheel here...
                PathD normals = Helper.initedPathD(ptCount);
                PathD previousNormals = Helper.initedPathD(ptCount);
                // Pre-calculate these for the threading to be an option.
                // This is a serial evaluation as we need both the previous and the current normal for each point.
#if !NOISESINGLETHREADED
                Parallel.For(0, ptCount - 1, pt =>
#else
                for (Int32 pt = 0; pt < ptCount - 1; pt++)
#endif
                    {
                        double dx;
                        double dy;
                        if (pt == 0)
                        {
                            // First vertex needs special care.
                            dx = mcPoints[0].x - mcPoints[ptCount - 2].x;
                            dy = mcPoints[0].y - mcPoints[ptCount - 2].y;
                        }
                        else
                        {
                            dx = mcPoints[pt + 1].x - mcPoints[pt].x;
                            dy = mcPoints[pt + 1].y - mcPoints[pt].y;
                        }

                        normals[pt] = new PointD(-dy, dx);
                    }
#if !NOISESINGLETHREADED
                );
#endif
                normals[^1] = new PointD(normals[0]);

                int nLength = normals.Count;
#if !NOISESINGLETHREADED
                Parallel.For(1, nLength, pt =>
#else
                for (int pt = 1; pt < nLength; pt++)
#endif
                    {
                        previousNormals[pt] = new PointD(normals[pt - 1]);
                    }
#if !NOISESINGLETHREADED
                );
#endif

                previousNormals[0] = new PointD(normals[^2]);

#if !NOISESINGLETHREADED
                Parallel.For(0, ptCount - 1, pt =>
#else
                for (int pt = 0; pt < ptCount - 1; pt++)
#endif
                    {
                        // We need to average the normals of two edge segments to get the vector we need to displace our point along.
                        // This ensures that we handle corners and fluctuations in a reasonable manner.
                        PointD averagedEdgeNormal = new((previousNormals[pt].x + normals[pt].x) / 2.0f,
                            (previousNormals[pt].y + normals[pt].y) / 2.0f);
                        // Normalize our vector length.
                        double length = Math.Sqrt(Utils.myPow(averagedEdgeNormal.x, 2) +
                                                  Utils.myPow(averagedEdgeNormal.y, 2));
                        const double normalTolerance = 1E-3;
                        if (length < normalTolerance)
                        {
                            length = normalTolerance;
                        }

                        averagedEdgeNormal = new PointD(averagedEdgeNormal.x / length, averagedEdgeNormal.y / length);

                        // Use a tolerance as we're handling floats; we don't expect a normalized absolute value generally above 1.0, ignoring the float error.
                        /*
                        if ((Math.Abs(averagedEdgeNormal.X) - 1 > normalTolerance) || (Math.Abs(averagedEdgeNormal.Y) - 1 > normalTolerance))
                        {
                            ErrorReporter.showMessage_OK("averageNormal exceeded limits: X:" + averagedEdgeNormal.X.ToString() + ",Y:" + averagedEdgeNormal.Y.ToString(), "oops");
                        }
                        */

                        // We can now modify the position of our point and stuff it into our jittered list.

                        double jitterAmount = noiseType switch
                        {
                            (int)noiseIndex.opensimplex => ((OpenSimplexNoise)noiseSource).Evaluate(
                                freq * (mcPoints[pt].x + x_const), freq * (mcPoints[pt].y + y_const)),
                            (int)noiseIndex.simplex => ((SimplexNoise)noiseSource).GetNoise(
                                freq * (mcPoints[pt].x + x_const), freq * (mcPoints[pt].y + y_const)),
                            _ => ((PerlinNoise)noiseSource).Noise(freq * (mcPoints[pt].x + x_const),
                                freq * (mcPoints[pt].y + y_const), 0)
                        };

                        jitterAmount *= jitterScale;

                        double jitteredX = mcPoints[pt].x;
                        jitteredX += jitterAmount * averagedEdgeNormal.x;

                        double jitteredY = mcPoints[pt].y;
                        jitteredY += jitterAmount * averagedEdgeNormal.y;

                        jitteredPoints[pt] = new PointD(jitteredX, jitteredY);
                    }
#if !NOISESINGLETHREADED
                );
#endif
                jitteredPoints[ptCount - 1] = new PointD(jitteredPoints[0]);

                // Push back to mcPoints for further processing.
                ret.Add(jitteredPoints);
            }
        }

        /*
        if (ret.Count != input.Count)
        {
            SvgWriter svgWriter = new SvgWriter();
            SvgUtils.AddSubject(svgWriter, input);
            SvgUtils.AddSolution(svgWriter, ret, false);
            SvgUtils.SaveToFile(svgWriter, "/d/development/debug_prekh.svg", FillRule.Positive);
            int xx = 2;
        }
        */

        ret = GeoWrangler.makeKeyHole(ret, false, false);

        /*
        if (ret.Count != input.Count)
        {
            SvgWriter svgWriter = new SvgWriter();
            SvgUtils.AddSubject(svgWriter, input);
            SvgUtils.AddSolution(svgWriter, ret, false);
            SvgUtils.SaveToFile(svgWriter, "/d/development/debug.svg", FillRule.Positive);
            int xx = 2;
        }
        */

        return ret;
    }    
}