using NUnit.Framework;
using Clipper2Lib;
using shapeEngine;
using System;
using System.Linq;

namespace UnitTests
{
    [TestFixture]
    public class CircularConvergenceTests
    {
        [Test]
        public void SquareWithLargeRadius_ConvergesToCircle()
        {
            // Arrange
            var square = new PathD
            {
                new PointD(0, 0),
                new PointD(100, 0),
                new PointD(100, 100),
                new PointD(0, 100),
                new PointD(0, 0)
            };

            // Act - Use radius larger than half-edge (50 * 1.1 = 55)
            var result = contourGen.makeContour(
                square, 
                concaveRadius: 0, 
                convexRadius: 60, 
                edgeResolution: 2.0, 
                angularResolution: 0.2, 
                shortEdgeLength: 5, 
                maxShortEdgeLength: 10, 
                optimizeCorners: 0);

            // Assert - Should produce a reasonably circular result
            AssertCircularity(result, maxVariation: 0.07, description: "Square with large radius");
            
            // Should produce significantly fewer points than the non-circular case
            Assert.That(result.Count, Is.LessThan(100), 
                "Circular convergence should produce fewer points");
        }

        [Test]
        public void SquareWithSmallRadius_DoesNotConvergeToCircle()
        {
            // Arrange
            var square = new PathD
            {
                new PointD(0, 0),
                new PointD(100, 0),
                new PointD(100, 100),
                new PointD(0, 100),
                new PointD(0, 0)
            };

            // Act - Use radius smaller than convergence threshold
            var result = contourGen.makeContour(
                square, 
                concaveRadius: 0, 
                convexRadius: 20, 
                edgeResolution: 2.0, 
                angularResolution: 0.2, 
                shortEdgeLength: 5, 
                maxShortEdgeLength: 10, 
                optimizeCorners: 0);

            // Assert - Should NOT be circular (should have more variation)
            var circularity = CalculateCircularityScore(result);
            Assert.That(circularity, Is.LessThan(0.98), 
                "Small radius should not produce perfect circular result");
                
            // Should produce more points than the circular case
            Assert.That(result.Count, Is.GreaterThan(200), 
                "Non-circular case should produce more points");
        }

        [Test]
        public void RectangleWithLargeRadius_ConvergesToCircle()
        {
            // Arrange
            var rectangle = new PathD
            {
                new PointD(0, 0),
                new PointD(120, 0),
                new PointD(120, 60),
                new PointD(0, 60),
                new PointD(0, 0)
            };

            // Act - Use radius larger than half of shortest edge (30 * 1.1 = 33)
            var result = contourGen.makeContour(
                rectangle, 
                concaveRadius: 0, 
                convexRadius: 40, 
                edgeResolution: 2.0, 
                angularResolution: 0.2, 
                shortEdgeLength: 5, 
                maxShortEdgeLength: 10, 
                optimizeCorners: 0);

            // Assert
            AssertCircularity(result, maxVariation: 0.07, description: "Rectangle with large radius");
            Assert.That(result.Count, Is.LessThan(100), 
                "Circular convergence should produce fewer points");
        }

        [Test]
        public void ConcavePolygon_DoesNotConvergeToCircle()
        {
            // Arrange - Create a concave L-shape
            var lShape = new PathD
            {
                new PointD(0, 0),
                new PointD(100, 0),
                new PointD(100, 50),
                new PointD(50, 50),
                new PointD(50, 100),
                new PointD(0, 100),
                new PointD(0, 0)
            };

            // Act - Use large radius 
            var result = contourGen.makeContour(
                lShape, 
                concaveRadius: 0, 
                convexRadius: 60, 
                edgeResolution: 2.0, 
                angularResolution: 0.2, 
                shortEdgeLength: 5, 
                maxShortEdgeLength: 10, 
                optimizeCorners: 0);

            // Assert - Should NOT be circular because it's concave
            var circularity = CalculateCircularityScore(result);
            Assert.That(circularity, Is.LessThan(0.95), 
                "Concave shapes should not converge to circles even with large radius");
        }

        [Test]
        public void ComplexPolygon_DoesNotConvergeToCircle()
        {
            // Arrange - Create a polygon with many corners (should exceed the 8-corner limit)
            var star = new PathD();
            int points = 10; // 10-pointed star = 20 corners, exceeds our 8-corner limit
            for (int i = 0; i < points * 2; i++)
            {
                double angle = i * Math.PI / points;
                double radius = (i % 2 == 0) ? 50 : 25; // Alternating outer/inner radius
                star.Add(new PointD(
                    50 + radius * Math.Cos(angle),
                    50 + radius * Math.Sin(angle)
                ));
            }
            star.Add(star[0]); // Close the path

            // Act - Use large radius 
            var result = contourGen.makeContour(
                star, 
                concaveRadius: 0, 
                convexRadius: 40, 
                edgeResolution: 2.0, 
                angularResolution: 0.2, 
                shortEdgeLength: 5, 
                maxShortEdgeLength: 10, 
                optimizeCorners: 0);

            // Assert - Should NOT be circular because it has too many corners
            var circularity = CalculateCircularityScore(result);
            Assert.That(circularity, Is.LessThan(0.95), 
                "Complex polygons with many corners should not converge to circles");
        }

        private void AssertCircularity(PathD path, double maxVariation, string description)
        {
            var circularity = CalculateCircularityScore(path);
            Assert.That(circularity, Is.GreaterThan(1.0 - maxVariation), 
                $"{description} should be circular (circularity score: {circularity:F4})");
        }

        private double CalculateCircularityScore(PathD path)
        {
            if (path.Count < 3)
                return 0.0;

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

            return 1.0 - radiusVariation; // 1.0 = perfect circle, lower values = less circular
        }
    }
}