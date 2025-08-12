using System;
using System.Collections.Generic;

namespace PolygonClassification
{
    struct Point
    {
        public double X, Y;
        public Point(double x, double y) { X = x; Y = y; }
    }

    class Program
    {
        static void Main()
        {
            var vertices = new List<Point>
            {
                new Point(0.41600, -0.35500),
                new Point(0.14100, -0.29100),
                new Point(0.06100,  0.22800),
                new Point(-0.28900, 0.15300),
                new Point(-0.39400,-0.10700),
                new Point(-0.69600, 0.01000),
                new Point(-0.40400, 0.50100),
                new Point(0.46000,  0.44200),
                new Point(0.46200,  0.35500),
                new Point(0.67400,  0.24200),
                new Point(0.50100,  0.01000)
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

            // Prepare status list
            var status = new List<string>(new string[vertices.Count]);

            for (int i = 0; i < vertices.Count; i++)
            {
                var prev = vertices[(i - 1 + vertices.Count) % vertices.Count];
                var curr = vertices[i];
                var next = vertices[(i + 1) % vertices.Count];

                // Vectors: prev→curr and curr→next
                double vx1 = curr.X - prev.X;
                double vy1 = curr.Y - prev.Y;
                double vx2 = next.X - curr.X;
                double vy2 = next.Y - curr.Y;

                // Z component of 3D cross product
                double crossZ = vx1 * vy2 - vy1 * vx2;

                // For CCW polygon, positive crossZ = convex. For CW, negative = convex.
                bool isVertexConvex = isCCW ? (crossZ > 0) : (crossZ < 0);
                status[i] = isVertexConvex ? "Convex" : "Concave";
            }

            // Print statuses in input order
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
