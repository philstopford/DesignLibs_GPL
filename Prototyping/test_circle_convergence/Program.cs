using System;
using System.Collections.Generic;
using System.IO;
using System.Globalization;
using System.Linq;
using Clipper2Lib;
using shapeEngine;

namespace test_circle_convergence
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Testing ContourGen Circular Convergence Fix");
            Console.WriteLine("==========================================");

            // Test Case 1: Square with small radius (should produce rounded square)
            TestSquareSmallRadius();
            
            // Test Case 2: Square with large radius (should converge to circle)
            TestSquareLargeRadius();
            
            // Test Case 3: Rectangle with large radius 
            TestRectangleLargeRadius();
            
            // Test Case 4: Hexagon with large radius
            TestHexagonLargeRadius();

            Console.WriteLine("\nAll tests completed. Check generated SVG files for visual verification.");
        }

        static void TestSquareSmallRadius()
        {
            Console.WriteLine("\n1. Testing Square with Small Radius (10):");
            
            var square = new PathD
            {
                new PointD(0, 0),
                new PointD(100, 0),
                new PointD(100, 100),
                new PointD(0, 100),
                new PointD(0, 0)
            };

            var result = contourGen.makeContour(
                square, 
                concaveRadius: 0, 
                convexRadius: 10, 
                edgeResolution: 1.0, 
                angularResolution: 0.1, 
                shortEdgeLength: 5, 
                maxShortEdgeLength: 10, 
                optimizeCorners: 0);

            CreateSvgOutput("square_small_radius.svg", square, result, "Square with small radius (10)");
            AnalyzeCircularity(result, "Small radius");
        }

        static void TestSquareLargeRadius()
        {
            Console.WriteLine("\n2. Testing Square with Large Radius (60):");
            
            var square = new PathD
            {
                new PointD(0, 0),
                new PointD(100, 0),
                new PointD(100, 100),
                new PointD(0, 100),
                new PointD(0, 0)
            };

            // Debug: Check edge half lengths  
            double edgeLength = 100; // Each edge of the square
            double halfEdge = edgeLength / 2.0; // = 50
            double radius = 60;
            Console.WriteLine($"  Edge length: {edgeLength}, Half-edge: {halfEdge}, Radius: {radius}");
            Console.WriteLine($"  Radius > half-edge * 1.2: {radius} > {halfEdge * 1.2} = {radius > halfEdge * 1.2}");

            var result = contourGen.makeContour(
                square, 
                concaveRadius: 0, 
                convexRadius: 60, 
                edgeResolution: 1.0, 
                angularResolution: 0.1, 
                shortEdgeLength: 5, 
                maxShortEdgeLength: 10, 
                optimizeCorners: 0);

            CreateSvgOutput("square_large_radius.svg", square, result, "Square with large radius (60) - should be circular");
            AnalyzeCircularity(result, "Large radius (should be circular)");
        }

        static void TestRectangleLargeRadius()
        {
            Console.WriteLine("\n3. Testing Rectangle with Large Radius (40):");
            
            var rectangle = new PathD
            {
                new PointD(0, 0),
                new PointD(120, 0),
                new PointD(120, 60),
                new PointD(0, 60),
                new PointD(0, 0)
            };

            var result = contourGen.makeContour(
                rectangle, 
                concaveRadius: 0, 
                convexRadius: 40, // Greater than half of shorter edge (30)
                edgeResolution: 1.0, 
                angularResolution: 0.1, 
                shortEdgeLength: 5, 
                maxShortEdgeLength: 10, 
                optimizeCorners: 0);

            CreateSvgOutput("rectangle_large_radius.svg", rectangle, result, "Rectangle with large radius (40) - should be circular");
            AnalyzeCircularity(result, "Rectangle large radius");
        }

        static void TestHexagonLargeRadius()
        {
            Console.WriteLine("\n4. Testing Hexagon with Large Radius (50):");
            
            // Create a regular hexagon centered at origin, then translate
            var hexagon = new PathD();
            for (int i = 0; i < 6; i++)
            {
                double angle = i * Math.PI / 3; // 60 degrees apart
                hexagon.Add(new PointD(
                    50 + 40 * Math.Cos(angle),
                    50 + 40 * Math.Sin(angle)
                ));
            }
            hexagon.Add(hexagon[0]); // Close the path

            var result = contourGen.makeContour(
                hexagon, 
                concaveRadius: 0, 
                convexRadius: 50, 
                edgeResolution: 1.0, 
                angularResolution: 0.1, 
                shortEdgeLength: 5, 
                maxShortEdgeLength: 10, 
                optimizeCorners: 0);

            CreateSvgOutput("hexagon_large_radius.svg", hexagon, result, "Hexagon with large radius (50) - should be circular");
            AnalyzeCircularity(result, "Hexagon large radius");
        }

        static void CreateSvgOutput(string filename, PathD original, PathD result, string title)
        {
            var svg = new System.Text.StringBuilder();
            
            // Find bounds for both paths
            double minX = Math.Min(original.Min(p => p.x), result.Min(p => p.x)) - 20;
            double minY = Math.Min(original.Min(p => p.y), result.Min(p => p.y)) - 20;
            double maxX = Math.Max(original.Max(p => p.x), result.Max(p => p.x)) + 20;
            double maxY = Math.Max(original.Max(p => p.y), result.Max(p => p.y)) + 20;
            
            double width = maxX - minX;
            double height = maxY - minY;
            
            svg.AppendLine($"<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            svg.AppendLine($"<svg width=\"{width}\" height=\"{height}\" viewBox=\"{minX} {minY} {width} {height}\" xmlns=\"http://www.w3.org/2000/svg\">");
            svg.AppendLine($"<title>{title}</title>");
            
            // Add original shape in light gray
            svg.Append("<path d=\"M ");
            foreach (var pt in original)
            {
                svg.Append($"{pt.x.ToString(CultureInfo.InvariantCulture)},{pt.y.ToString(CultureInfo.InvariantCulture)} ");
            }
            svg.AppendLine("\" stroke=\"#cccccc\" stroke-width=\"1\" fill=\"none\" />");
            
            // Add result shape in blue
            svg.Append("<path d=\"M ");
            foreach (var pt in result)
            {
                svg.Append($"{pt.x.ToString(CultureInfo.InvariantCulture)},{pt.y.ToString(CultureInfo.InvariantCulture)} ");
            }
            svg.AppendLine("\" stroke=\"blue\" stroke-width=\"2\" fill=\"none\" />");
            
            // Add center point for analysis
            if (result.Count > 0)
            {
                double centerX = result.Average(p => p.x);
                double centerY = result.Average(p => p.y);
                svg.AppendLine($"<circle cx=\"{centerX.ToString(CultureInfo.InvariantCulture)}\" cy=\"{centerY.ToString(CultureInfo.InvariantCulture)}\" r=\"2\" fill=\"red\" />");
            }
            
            svg.AppendLine("</svg>");
            
            File.WriteAllText(filename, svg.ToString());
        }

        static void AnalyzeCircularity(PathD path, string description)
        {
            if (path.Count < 3)
            {
                Console.WriteLine($"  {description}: Too few points to analyze circularity");
                return;
            }

            // Find center by averaging all points
            double centerX = path.Average(p => p.x);
            double centerY = path.Average(p => p.y);
            var center = new PointD(centerX, centerY);

            // Calculate distances from center
            var distances = path.Select(p => Helper.Length(Helper.Minus(p, center))).ToArray();
            
            double avgRadius = distances.Average();
            double minRadius = distances.Min();
            double maxRadius = distances.Max();
            double radiusVariation = (maxRadius - minRadius) / avgRadius;

            Console.WriteLine($"  Points: {path.Count}");
            Console.WriteLine($"  Center: ({centerX:F2}, {centerY:F2})");
            Console.WriteLine($"  Average radius: {avgRadius:F2}");
            Console.WriteLine($"  Radius variation: {radiusVariation:F4} ({radiusVariation * 100:F2}%)");
            Console.WriteLine($"  Circularity score: {1 - radiusVariation:F4} (1.0 = perfect circle)");
        }
    }
}