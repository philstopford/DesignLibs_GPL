using System;
using System.Collections.Generic;
using System.Numerics;

namespace RectilinearCenteredTurns
{
    /// <summary>
    /// Demonstrates a rectilinear path approximation algorithm using dynamic programming.
    /// 
    /// This prototype converts an arbitrary polygon path into a rectilinear (Manhattan-style)
    /// path that minimizes projection error while ensuring all segments are either purely
    /// horizontal or vertical. The algorithm uses rotation optimization and dynamic programming
    /// to find the best axis-aligned approximation.
    /// 
    /// Algorithm Overview:
    /// 1. Test multiple rotation angles to find optimal coordinate frame
    /// 2. Transform input points to rotated coordinate system
    /// 3. Use dynamic programming to find optimal H/V segment assignment
    /// 4. Transform back to original coordinate system
    /// 
    /// Applications:
    /// - PCB routing and trace layout optimization
    /// - VLSI physical design and Manhattan routing
    /// - Architectural floor plan simplification
    /// - CNC machining path optimization for axis-parallel cuts
    /// - Game pathfinding on grid-based maps
    /// 
    /// Mathematical Foundation:
    /// - Rotation matrices for coordinate transformation
    /// - Dynamic programming for optimal substructure problems
    /// - L1 norm minimization for Manhattan distance metrics
    /// - Cross product calculations for angle measurement
    /// </summary>
    public static class Program
    {
        /// <summary>
        /// Test polygon vertices representing an irregular path to be rectilinearized.
        /// These points simulate a typical input that might come from hand-drawn sketches,
        /// sensor data, or other non-axis-aligned sources.
        /// </summary>
        private static readonly List<Vector2> TestVertices = new List<Vector2> {
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
            var origin = TestVertices[0];
            var (bestTheta, bestErr) = FindOptimalRotation(origin);
            var finalPts = RectilinearizeCentered(bestTheta, origin);

            Console.WriteLine($"Optimal θ = {bestTheta * 180.0 / Math.PI:F2}° (Error = {bestErr:F6})\n");
            Console.WriteLine("Idx   Orig X,Y           →  Final X,Y    Angle (°)");

            var angles = ComputeAngles(finalPts);
            for (int i = 0; i < TestVertices.Count; i++)
            {
                var o = TestVertices[i];
                var f = finalPts[i];
                Console.WriteLine(
                    $"{i,2}: ({o.X,7:F5},{o.Y,7:F05}) → ({f.X,7:F05},{f.Y,7:F05})    {angles[i],7:F2}"
                );
            }
        }

        /// <summary>
        /// Finds the optimal rotation angle that minimizes rectilinear approximation error.
        /// 
        /// Tests rotation angles from 0° to 90° in small increments to find the coordinate
        /// frame that best aligns with the dominant directions in the input path.
        /// </summary>
        /// <param name="origin">Reference point for rotation center</param>
        /// <returns>Tuple of optimal angle in radians and minimum error achieved</returns>
        public static (double theta, double error) FindOptimalRotation(Vector2 origin)
        {
            double bestTheta = 0, bestErr = double.MaxValue;
            int steps = 900;  // 0.1 degree increments for high precision
            double dθ = (Math.PI / 2) / steps;

            for (int i = 0; i <= steps; i++)
            {
                double θ = i * dθ;
                double err = ComputeRectilinearError(θ, origin);
                if (err < bestErr)
                {
                    bestErr = err;
                    bestTheta = θ;
                }
            }
            return (bestTheta, bestErr);
        }

        /// <summary>
        /// Computes the total approximation error for a given rotation angle using dynamic programming.
        /// 
        /// This function rotates the input points into a candidate coordinate frame and then
        /// uses dynamic programming to find the optimal assignment of horizontal/vertical
        /// segments that minimizes the total projection error.
        /// </summary>
        /// <param name="θ">Rotation angle in radians</param>
        /// <param name="origin">Center point for rotation</param>
        /// <returns>Total L1 projection error for this rotation</returns>
        public static double ComputeRectilinearError(double θ, Vector2 origin)
        {
            float c = (float)Math.Cos(θ), s = (float)Math.Sin(θ);
            int N = TestVertices.Count;

            // Step 1: Rotate points into the candidate coordinate frame
            var P = new Vector2[N];
            for (int i = 0; i < N; i++)
            {
                var d = TestVertices[i] - origin;
                P[i] = new Vector2(c * d.X - s * d.Y, s * d.X + c * d.Y);
            }

            // Step 2: Compute displacement vectors between consecutive points
            var D = new Vector2[N];
            for (int i = 1; i < N; i++)
                D[i] = P[i] - P[i - 1];

            // Step 3: Dynamic programming to find optimal axis assignment
            // State: cost[i, axis] = minimum cost to reach point i using 'axis' for segment i
            // where axis 0 = horizontal constraint, axis 1 = vertical constraint
            const double INF = 1e9;
            var cost = new double[N, 2];
            cost[0, 0] = cost[0, 1] = 0;  // No cost for starting point

            for (int i = 1; i < N; i++)
            {
                for (int axis = 0; axis < 2; axis++)
                {
                    // Project displacement onto chosen axis and compute error
                    double projectionError = axis == 0
                        ? Math.Abs(D[i].Y)  // Horizontal: penalize Y component
                        : Math.Abs(D[i].X); // Vertical: penalize X component

                    // Must alternate axes between consecutive segments
                    cost[i, axis] = cost[i - 1, 1 - axis] + projectionError;
                }
            }

            // Return minimum cost over both possible final axes
            return Math.Min(cost[N - 1, 0], cost[N - 1, 1]);
        }

        /// <summary>
        /// Generates the final rectilinear path using the optimal rotation and dynamic programming.
        /// 
        /// This is the main algorithm that produces the actual rectilinear approximation.
        /// It performs the same dynamic programming as the error computation, but also
        /// reconstructs the optimal path by backtracking through the DP table.
        /// </summary>
        /// <param name="θ">Optimal rotation angle</param>
        /// <param name="origin">Center point for rotation</param>
        /// <returns>List of points forming the rectilinear path</returns>
        public static List<Vector2> RectilinearizeCentered(double θ, Vector2 origin)
        {
            float c = (float)Math.Cos(θ), s = (float)Math.Sin(θ);
            int N = TestVertices.Count;

            // Step 1: Rotate points into optimal coordinate frame
            var P = new Vector2[N];
            for (int i = 0; i < N; i++)
            {
                var d = TestVertices[i] - origin;
                P[i] = new Vector2(c * d.X - s * d.Y, s * d.X + c * d.Y);
            }

            // Step 2: Compute displacement vectors
            var D = new Vector2[N];
            for (int i = 1; i < N; i++) D[i] = P[i] - P[i - 1];

            // Step 3: Dynamic programming with backtracking information
            const double INF = 1e9;
            var cost = new double[N, 2];
            var parent = new int[N, 2];  // Track optimal previous axis choice
            cost[0, 0] = cost[0, 1] = 0;
            parent[0, 0] = parent[0, 1] = -1;

            for (int i = 1; i < N; i++)
            {
                for (int axis = 0; axis < 2; axis++)
                {
                    double projErr = axis == 0
                        ? Math.Abs(D[i].Y)
                        : Math.Abs(D[i].X);

                    // Find best previous axis (must be different from current)
                    double best = INF;
                    int par = -1;
                    for (int prevAxis = 0; prevAxis < 2; prevAxis++)
                    {
                        if (prevAxis == axis) continue;  // Force alternation
                        double cst = cost[i - 1, prevAxis] + projErr;
                        if (cst < best)
                        {
                            best = cst;
                            par = prevAxis;
                        }
                    }

                    cost[i, axis] = best;
                    parent[i, axis] = par;
                }
            }

            // Step 4: Backtrack to find optimal axis sequence
            var axisChoice = new int[N];
            axisChoice[N - 1] = cost[N - 1, 0] < cost[N - 1, 1] ? 0 : 1;
            for (int i = N - 1; i > 0; i--)
                axisChoice[i - 1] = parent[i, axisChoice[i]];

            // Step 5: Reconstruct rectilinear path with axis constraints
            var Q = new List<Vector2>(N) { P[0] };
            for (int i = 1; i < N; i++)
            {
                var prev = Q[i - 1];
                var curr = P[i];

                // Apply axis constraint by zeroing the non-aligned component
                if (axisChoice[i] == 0)
                {
                    // Horizontal segment: maintain Y coordinate, use target X
                    Q.Add(new Vector2(curr.X, prev.Y));
                }
                else
                {
                    // Vertical segment: maintain X coordinate, use target Y
                    Q.Add(new Vector2(prev.X, curr.Y));
                }
            }

            if (Q.Count != N)
                throw new InvalidOperationException($"Expected {N} points, got {Q.Count}");

            // Step 6: Rotate back to original coordinate system
            var result = new List<Vector2>(N);
            float ic = c, is_ = -s;  // Inverse rotation matrix
            foreach (var q in Q)
                result.Add(new Vector2(ic * q.X - is_ * q.Y,
                    is_ * q.X + ic * q.Y) + origin);

            return result;
        }

        /// <summary>
        /// Computes interior angles at each vertex of the final rectilinear path.
        /// 
        /// This is useful for validating the output and understanding the geometric
        /// properties of the resulting path. Angles close to 90° or 270° indicate
        /// good rectilinear approximation.
        /// </summary>
        /// <param name="pts">List of path vertices</param>
        /// <returns>Array of interior angles in degrees (NaN for endpoints)</returns>
        public static double[] ComputeAngles(List<Vector2> pts)
        {
            int N = pts.Count;
            var angles = new double[N];

            for (int i = 0; i < N; i++)
            {
                // Endpoints don't have well-defined interior angles
                if (i == 0 || i == N - 1)
                {
                    angles[i] = double.NaN;
                    continue;
                }

                // Compute normalized vectors for incoming and outgoing segments
                Vector2 v1 = Vector2.Normalize(pts[i - 1] - pts[i]);
                Vector2 v2 = Vector2.Normalize(pts[i + 1] - pts[i]);
                
                // Compute angle using dot product (clamped for numerical stability)
                double dot = Math.Clamp(Vector2.Dot(v1, v2), -1.0f, 1.0f);
                angles[i] = Math.Acos(dot) * (180.0 / Math.PI);
            }

            return angles;
        }
    }
}
