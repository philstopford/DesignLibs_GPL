using System;
using System.Collections.Generic;

namespace PolygonClassification
{
    struct Point
    {
        public double X, Y;
        public Point(double x, double y) { X = x; Y = y; }
    }

    static class Program
    {
        static void Main()
        {
            // Threshold for short edges
            double shortEdgeLength = 0.011;
            // Tolerance for orthogonality (dot product ≈ 0)
            double orthoEpsilon = 1e-6;

            var vertices = new List<Point>
            {
                new Point(0.00000, 0.00000),
                new Point(0.00000, 0.10000),
                new Point(0.02000, 0.10000),
                new Point(0.02000, 0.09000),
                new Point(0.03000, 0.09000),
                new Point(0.03000, 0.08000),
                new Point(0.04000, 0.08000),
                new Point(0.04000, 0.07000),
                new Point(0.05000, 0.07000),
                new Point(0.05000, 0.06000),
                new Point(0.08000, 0.06000),
                new Point(0.08000, 0.04000),
                new Point(0.09000, 0.04000),
                new Point(0.09000, 0.03000),
                new Point(0.10000, 0.03000),
                new Point(0.10000, 0.00000)
            };

            // Determine polygon orientation: positive = CCW, negative = CW
            double area2 = 0;
            for (int i = 0; i < vertices.Count; i++)
            {
                var p1 = vertices[i];
                var p2 = vertices[(i + 1) % vertices.Count];
                area2 += p1.X * p2.Y - p2.X * p1.Y;
            }
            bool isCCW = area2 > 0;

            var status = new List<string>(new string[vertices.Count]);

            for (int i = 0; i < vertices.Count; i++)
            {
                var prev = vertices[(i - 1 + vertices.Count) % vertices.Count];
                var curr = vertices[i];
                var next = vertices[(i + 1) % vertices.Count];

                // Edge vectors
                double vx1 = curr.X - prev.X;
                double vy1 = curr.Y - prev.Y;
                double vx2 = next.X - curr.X;
                double vy2 = next.Y - curr.Y;

                // Edge lengths
                double len1 = Math.Sqrt(vx1 * vx1 + vy1 * vy1);
                double len2 = Math.Sqrt(vx2 * vx2 + vy2 * vy2);

                // Check for two short edges that are orthogonal
                if (len1 <= shortEdgeLength && len2 <= shortEdgeLength)
                {
                    double dot = vx1 * vx2 + vy1 * vy2;
                    if (Math.Abs(dot) <= orthoEpsilon)
                    {
                        status[i] = "ShortEdge";
                        continue;
                    }
                }

                // Compute cross product Z for convex/concave
                double crossZ = vx1 * vy2 - vy1 * vx2;
                bool isVertexConvex = isCCW ? (crossZ > 0) : (crossZ < 0);
                status[i] = isVertexConvex ? "Convex" : "Concave";
            }

            // Print statuses
            for (int i = 0; i < vertices.Count; i++)
            {
                var p = vertices[i];
                Console.WriteLine(
                    $"Vertex {i + 1}: ({p.X:F5}, {p.Y:F5}) → **{status[i]}**"
                );
            }
        }
    }
}
