using System;
using System.Collections.Generic;
using System.Numerics;

namespace RectilinearCenteredTurns
{
    class Program
    {
        static readonly List<Vector2> Orig = new List<Vector2> {
            new Vector2(-0.36600f, -0.39900f),
            new Vector2(-0.59900f,  0.36600f),
            new Vector2( 0.16600f,  0.59900f),
            new Vector2( 0.39900f, -0.16600f),
            new Vector2( 0.26800f, -0.20700f),
            new Vector2( 0.13400f,  0.06500f),
            new Vector2(-0.23600f, -0.06100f),
            new Vector2(-0.21300f, -0.34200f)
        };

        static void Main()
        {
            var origin = Orig[0];
            var (bestTheta, bestErr) = FindBestTheta(origin);
            var finalPts = RectilinearizeCentered(bestTheta, origin);

            Console.WriteLine($"Optimal θ = {bestTheta * 180.0 / Math.PI:F2}° (Error = {bestErr:F6})\n");

            Console.WriteLine("Idx   Orig X,Y           →  Final X,Y    Angle (°)");
            var angles = ComputeAngles(finalPts);
            for (int i = 0; i < Orig.Count; i++)
            {
                var o = Orig[i];
                var f = finalPts[i];
                Console.WriteLine(
                    $"{i,2}: ({o.X,7:F5},{o.Y,7:F05}) → ({f.X,7:F05},{f.Y,7:F05})    {angles[i],7:F2}"
                );
            }
        }

        static (double theta, double error) FindBestTheta(Vector2 origin)
        {
            double bestTheta = 0, bestErr = double.MaxValue;
            int steps = 900;
            double dθ = (Math.PI / 2) / steps;

            for (int i = 0; i <= steps; i++)
            {
                double θ = i * dθ;
                double err = ComputeErrorDP(θ, origin);
                if (err < bestErr)
                {
                    bestErr = err;
                    bestTheta = θ;
                }
            }
            return (bestTheta, bestErr);
        }

        static double ComputeErrorDP(double θ, Vector2 origin)
        {
            float c = (float)Math.Cos(θ), s = (float)Math.Sin(θ);
            int N = Orig.Count;

            // 1) Rotate into frame
            var P = new Vector2[N];
            for (int i = 0; i < N; i++)
            {
                var d = Orig[i] - origin;
                P[i] = new Vector2(c * d.X - s * d.Y, s * d.X + c * d.Y);
            }

            // 2) Deltas
            var D = new Vector2[N];
            for (int i = 1; i < N; i++)
                D[i] = P[i] - P[i - 1];

            // 3) DP cost: 0 = horizontal, 1 = vertical
            const double INF = 1e9;
            var cost = new double[N, 2];
            cost[0, 0] = cost[0, 1] = 0;

            for (int i = 1; i < N; i++)
            {
                for (int axis = 0; axis < 2; axis++)
                {
                    double projErr = axis == 0
                        ? Math.Abs(D[i].Y)
                        : Math.Abs(D[i].X);

                    // only alternate
                    cost[i, axis] = Math.Min(
                        cost[i - 1, 1 - axis] + projErr,
                        INF
                    );
                }
            }

            return Math.Min(cost[N - 1, 0], cost[N - 1, 1]);
        }

        static List<Vector2> RectilinearizeCentered(double θ, Vector2 origin)
        {
            float c = (float)Math.Cos(θ), s = (float)Math.Sin(θ);
            int N = Orig.Count;

            // 1) Rotate into frame
            var P = new Vector2[N];
            for (int i = 0; i < N; i++)
            {
                var d = Orig[i] - origin;
                P[i] = new Vector2(c * d.X - s * d.Y, s * d.X + c * d.Y);
            }

            // 2) Deltas
            var D = new Vector2[N];
            for (int i = 1; i < N; i++)
                D[i] = P[i] - P[i - 1];

            // 3) DP & backtrack
            const double INF = 1e9;
            var cost   = new double[N, 2];
            var parent = new int[N, 2];
            cost[0, 0] = cost[0, 1] = 0;
            parent[0, 0] = parent[0, 1] = -1;

            for (int i = 1; i < N; i++)
            {
                for (int axis = 0; axis < 2; axis++)
                {
                    double projErr = axis == 0
                        ? Math.Abs(D[i].Y)
                        : Math.Abs(D[i].X);

                    double best = INF;
                    int    par  = -1;

                    for (int prev = 0; prev < 2; prev++)
                    {
                        if (prev == axis) continue;
                        double cst = cost[i - 1, prev] + projErr;
                        if (cst < best)
                        {
                            best = cst;
                            par  = prev;
                        }
                    }
                    cost[i, axis]   = best;
                    parent[i, axis] = par;
                }
            }

            var axisChoice = new int[N];
            axisChoice[N - 1] = cost[N - 1, 0] < cost[N - 1, 1] ? 0 : 1;
            for (int i = N - 1; i > 0; i--)
                axisChoice[i - 1] = parent[i, axisChoice[i]];

            // 4) Reconstruct with centering
            var Q = new List<Vector2> { P[0] };
            for (int i = 1; i < N; i++)
            {
                var prev = Q[i - 1];
                var curr = P[i];
                if (axisChoice[i] == 0)
                {
                    // horizontal: y = midpoint of original Ys
                    float midY = (prev.Y + curr.Y) * 0.5f;
                    Q.Add(new Vector2(curr.X, midY));
                }
                else
                {
                    // vertical: x = midpoint of original Xs
                    float midX = (prev.X + curr.X) * 0.5f;
                    Q.Add(new Vector2(midX, curr.Y));
                }
            }

            // 5) Rotate back to world
            var result = new List<Vector2>(N);
            float ic =  c, is_ = -s;
            foreach (var q in Q)
            {
                result.Add(new Vector2(
                    ic * q.X - is_ * q.Y,
                    is_ * q.X + ic * q.Y
                ) + origin);
            }

            return result;
        }
        
        // Compute interior angles at each point in the polyline
        static double[] ComputeAngles(List<Vector2> pts)
        {
            int N = pts.Count;
            var angles = new double[N];

            for (int i = 0; i < N; i++)
            {
                if (i == 0 || i == N - 1)
                {
                    // no angle at endpoints in an open polyline
                    angles[i] = double.NaN;
                    continue;
                }

                Vector2 v1 = Vector2.Normalize(pts[i - 1] - pts[i]);
                Vector2 v2 = Vector2.Normalize(pts[i + 1] - pts[i]);
                double dot = Math.Clamp(Vector2.Dot(v1, v2), -1.0f, 1.0f);
                // interior angle in radians, convert to degrees
                angles[i] = Math.Acos(dot) * (180.0 / Math.PI);
            }

            return angles;
        }
        
    }
}
