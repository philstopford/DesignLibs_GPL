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
    /// Advanced polygon vertex classification that handles special case of short orthogonal edges.
    /// 
    /// This prototype extends basic corner categorization by detecting and specially handling
    /// vertices where two short edges meet at approximately 90 degrees. This is particularly
    /// useful for detecting features like small notches or steps in polygon geometry.
    /// 
    /// Algorithm Enhancement:
    /// 1. Standard convex/concave classification using cross products
    /// 2. Additional check for vertices with two short adjacent edges
    /// 3. Orthogonality test using dot product ≈ 0
    /// 4. Special "ShortEdge" classification for qualifying vertices
    /// 
    /// Applications:
    /// - PCB trace routing where small steps are common
    /// - Architectural drawings with fine details
    /// - Any geometry with mixed large/small features
    /// 
    /// Technical Details:
    /// - Short edge threshold: configurable distance threshold
    /// - Orthogonality tolerance: dot product absolute value threshold
    /// - Preserves standard classification for normal vertices
    /// </summary>
    public static class Program
    {
        static void Main()
        {
            // Configuration parameters for short edge detection
            double shortEdgeLength = 0.011;  // Maximum length for "short" edges
            double orthoEpsilon = 1e-6;      // Tolerance for orthogonality check

            // Test polygon with mix of normal and short edges (staircase pattern)
            var vertices = new List<Point>
            {
                new Point(0.00000, 0.00000),
                new Point(0.00000, 0.10000),
                new Point(0.02000, 0.10000),
                new Point(0.02000, 0.09000),  // Short vertical edge
                new Point(0.03000, 0.09000),  // Short horizontal edge
                new Point(0.03000, 0.08000),  // Pattern continues...
                new Point(0.04000, 0.08000),
                new Point(0.04000, 0.07000),
                new Point(0.05000, 0.07000),
                new Point(0.05000, 0.06000),
                new Point(0.08000, 0.06000),  // Longer horizontal edge
                new Point(0.08000, 0.04000),  // Longer vertical edge
                new Point(0.09000, 0.04000),
                new Point(0.09000, 0.03000),
                new Point(0.10000, 0.03000),
                new Point(0.10000, 0.00000)
            };

            var classifications = ClassifyVerticesWithShortEdges(vertices, shortEdgeLength, orthoEpsilon);
            
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
        /// Classifies vertices with special handling for short orthogonal edges.
        /// </summary>
        /// <param name="vertices">List of polygon vertices in order</param>
        /// <param name="shortEdgeLength">Maximum length considered a "short" edge</param>
        /// <param name="orthoEpsilon">Tolerance for orthogonality check (dot product threshold)</param>
        /// <returns>Array of classification strings: "Convex", "Concave", or "ShortEdge"</returns>
        public static string[] ClassifyVerticesWithShortEdges(List<Point> vertices, 
            double shortEdgeLength, double orthoEpsilon)
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

                // Calculate edge vectors
                double vx1 = curr.X - prev.X;
                double vy1 = curr.Y - prev.Y;
                double vx2 = next.X - curr.X;
                double vy2 = next.Y - curr.Y;

                // Calculate edge lengths
                double len1 = Math.Sqrt(vx1 * vx1 + vy1 * vy1);
                double len2 = Math.Sqrt(vx2 * vx2 + vy2 * vy2);

                // Special case: two short edges that are orthogonal
                if (len1 <= shortEdgeLength && len2 <= shortEdgeLength)
                {
                    // Check orthogonality using dot product
                    double dot = vx1 * vx2 + vy1 * vy2;
                    if (Math.Abs(dot) <= orthoEpsilon)
                    {
                        status[i] = "ShortEdge";
                        continue;
                    }
                }

                // Standard convex/concave classification
                double crossZ = vx1 * vy2 - vy1 * vx2;
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
