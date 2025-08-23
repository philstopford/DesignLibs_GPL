using System;
using System.Collections.Generic;

namespace PolygonClassification
{
    /// <summary>
    /// Represents a 2D point with double precision coordinates.
    /// </summary>
    public struct Point
    {
        public double X, Y;
        
        /// <summary>
        /// Initializes a new Point with the specified coordinates.
        /// </summary>
        /// <param name="x">The X coordinate</param>
        /// <param name="y">The Y coordinate</param>
        public Point(double x, double y) { X = x; Y = y; }
    }

    /// <summary>
    /// Demonstrates polygon vertex classification using cross product analysis.
    /// 
    /// This prototype implements a corner categorization algorithm that determines
    /// whether each vertex in a polygon is convex or concave by analyzing the
    /// turning angle at each vertex using vector cross products.
    /// 
    /// Algorithm Overview:
    /// 1. Determine polygon orientation (clockwise or counter-clockwise) using the shoelace formula
    /// 2. For each vertex, compute vectors from previous vertex to current and current to next
    /// 3. Calculate cross product to determine turning direction
    /// 4. Classify as convex/concave based on cross product sign and polygon orientation
    /// 
    /// Mathematical Foundation:
    /// - Uses the shoelace formula: Area = 1/2 * Σ(x_i * y_(i+1) - x_(i+1) * y_i)
    /// - Cross product z-component: v1 × v2 = v1.x * v2.y - v1.y * v2.x
    /// - For CCW polygons: positive cross product = convex corner, negative = concave
    /// - For CW polygons: negative cross product = convex corner, positive = concave
    /// </summary>
    public static class Program
    {
        static void Main()
        {
            // Test polygon vertices - represents an irregular polygon for demonstration
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

            var classifications = ClassifyVertices(vertices);
            
            // Display results
            for (int i = 0; i < vertices.Count; i++)
            {
                var p = vertices[i];
                Console.WriteLine(
                    $"Vertex {i + 1}: ({p.X:F5}, {p.Y:F5}) → **{classifications[i]}**"
                );
            }
        }

        /// <summary>
        /// Classifies each vertex in a polygon as convex or concave.
        /// </summary>
        /// <param name="vertices">List of polygon vertices in order</param>
        /// <returns>Array of strings indicating "Convex" or "Concave" for each vertex</returns>
        public static string[] ClassifyVertices(List<Point> vertices)
        {
            if (vertices == null || vertices.Count < 3)
                throw new ArgumentException("Polygon must have at least 3 vertices");

            bool isCCW = DetermineOrientation(vertices);
            var status = new string[vertices.Count];

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

                // Z component of 3D cross product (determines turning direction)
                double crossZ = vx1 * vy2 - vy1 * vx2;

                // For CCW polygon, positive crossZ = convex. For CW, negative = convex.
                bool isVertexConvex = isCCW ? (crossZ > 0) : (crossZ < 0);
                status[i] = isVertexConvex ? "Convex" : "Concave";
            }

            return status;
        }

        /// <summary>
        /// Determines if a polygon is oriented counter-clockwise using the shoelace formula.
        /// </summary>
        /// <param name="vertices">List of polygon vertices</param>
        /// <returns>True if counter-clockwise, false if clockwise</returns>
        public static bool DetermineOrientation(List<Point> vertices)
        {
            if (vertices == null || vertices.Count < 3)
                throw new ArgumentException("Polygon must have at least 3 vertices");

            // Calculate twice the signed area using the shoelace formula
            double area2 = 0;
            for (int i = 0; i < vertices.Count; i++)
            {
                var p1 = vertices[i];
                var p2 = vertices[(i + 1) % vertices.Count];
                area2 += p1.X * p2.Y - p2.X * p1.Y;
            }
            
            // Positive area indicates counter-clockwise orientation
            return area2 > 0;
        }
    }
}
