using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;

namespace UnitTests
{
    /// <summary>
    /// Test data structures that mirror the prototype implementations for testing purposes.
    /// These allow us to test the algorithms without complex cross-project references.
    /// </summary>
    public static class PrototypeTestHelpers
    {
        /// <summary>
        /// Point structure for geometric calculations
        /// </summary>
        public struct Point
        {
            public double X, Y;
            public Point(double x, double y) { X = x; Y = y; }
            
            public static Point operator +(Point a, Point b) => new Point(a.X + b.X, a.Y + b.Y);
            public static Point operator -(Point a, Point b) => new Point(a.X - b.X, a.Y - b.Y);
            public static Point operator *(Point a, double s) => new Point(a.X * s, a.Y * s);
            
            public double Length() => Math.Sqrt(X * X + Y * Y);
            public Point Normalized()
            {
                double len = Length();
                return len > 0 ? new Point(X / len, Y / len) : new Point(0, 0);
            }
        }

        /// <summary>
        /// Polygon corner classification algorithms (mirror of prototype implementations)
        /// </summary>
        public static class CornerClassification
        {
            /// <summary>
            /// Classifies each vertex in a polygon as convex or concave.
            /// </summary>
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

                    // Z component of 3D cross product
                    double crossZ = vx1 * vy2 - vy1 * vx2;

                    // For CCW polygon, positive crossZ = convex. For CW, negative = convex.
                    bool isVertexConvex = isCCW ? (crossZ > 0) : (crossZ < 0);
                    status[i] = isVertexConvex ? "Convex" : "Concave";
                }

                return status;
            }

            /// <summary>
            /// Classifies vertices with special handling for short orthogonal edges.
            /// </summary>
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

        /// <summary>
        /// Quadratic Bézier curve sampling algorithms (mirror of prototype implementations)
        /// </summary>
        public static class BezierSampling
        {
            public enum SamplingMode
            {
                ByMaxSegmentLength,
                ByMaxAngle
            }

            /// <summary>
            /// Computes the control point for a quadratic Bézier curve given endpoints and tangent directions.
            /// </summary>
            public static Point ComputeControlPoint(Point start, Point startDir, Point end, Point endDir)
            {
                double dx = end.X - start.X;
                double dy = end.Y - start.Y;
                double det = startDir.X * endDir.Y - startDir.Y * endDir.X;
                
                if (Math.Abs(det) < 1e-6)
                    throw new Exception("Tangent lines are parallel; no unique control point.");

                double tParam = (dx * endDir.Y - dy * endDir.X) / det;
                return start + startDir * tParam;
            }

            /// <summary>
            /// Samples a quadratic Bézier curve using the specified sampling strategy.
            /// </summary>
            public static List<Point> SampleCurve(Point start, Point control, Point end, 
                SamplingMode mode, double scale)
            {
                return mode switch
                {
                    SamplingMode.ByMaxSegmentLength => SampleByMaxSegmentLength(start, control, end, 0.05 * scale),
                    SamplingMode.ByMaxAngle => SampleByMaxAngle(start, control, end, 5.0 * Math.PI / 180.0),
                    _ => throw new ArgumentOutOfRangeException(nameof(mode))
                };
            }

            /// <summary>
            /// Samples curve by recursively subdividing until all segments are below length threshold.
            /// </summary>
            public static List<Point> SampleByMaxSegmentLength(Point P0, Point P1, Point P2, double maxSegLen)
            {
                var pts = new List<Point> { P0 };
                SubdivideByLength(P0, P1, P2, maxSegLen, pts);
                pts.Add(P2);
                return pts;
            }

            /// <summary>
            /// Samples curve by subdividing based on angular changes between tangent vectors.
            /// </summary>
            public static List<Point> SampleByMaxAngle(Point P0, Point P1, Point P2, double maxAngle)
            {
                var pts = new List<Point> { P0 };
                SubdivideByAngle(P0, P1, P2, maxAngle, pts);
                pts.Add(P2);
                return pts;
            }

            private static void SubdivideByLength(Point p0, Point p1, Point p2, double maxSegLen, List<Point> outPts)
            {
                if ((p2 - p0).Length() <= maxSegLen)
                {
                    outPts.Add(p2);
                    return;
                }
                
                Point p01 = Mid(p0, p1), p12 = Mid(p1, p2), p012 = Mid(p01, p12);
                SubdivideByLength(p0,  p01,  p012, maxSegLen, outPts);
                SubdivideByLength(p012, p12,  p2,   maxSegLen, outPts);
            }

            private static void SubdivideByAngle(Point p0, Point p1, Point p2, double maxAngle, List<Point> outPts)
            {
                var tan0 = (p1 - p0).Normalized();
                var tan1 = (p2 - p1).Normalized();
                
                double dot = Math.Max(-1.0, Math.Min(1.0, tan0.X * tan1.X + tan0.Y * tan1.Y));
                double angle = Math.Acos(dot);
                
                if (angle <= maxAngle)
                {
                    outPts.Add(p2);
                    return;
                }
                
                Point p01 = Mid(p0, p1), p12 = Mid(p1, p2), p012 = Mid(p01, p12);
                SubdivideByAngle(p0,  p01,  p012, maxAngle, outPts);
                SubdivideByAngle(p012, p12,  p2,   maxAngle, outPts);
            }

            private static Point Mid(Point a, Point b) => new Point((a.X + b.X) / 2, (a.Y + b.Y) / 2);
        }
    }

    /// <summary>
    /// Comprehensive test suite for prototyping algorithms covering polygon corner
    /// categorization, Bézier curve sampling, and related geometric operations.
    /// </summary>
    [TestFixture]
    public class PrototypingTests
    {
        #region Corner Categorization Tests

        [Test]
        public void TestBasicSquareConvexClassification()
        {
            // Arrange: Simple square (all corners should be convex)
            var square = new List<PrototypeTestHelpers.Point>
            {
                new(0, 0),
                new(1, 0),
                new(1, 1),
                new(0, 1)
            };

            // Act
            var result = PrototypeTestHelpers.CornerClassification.ClassifyVertices(square);

            // Assert: All vertices should be convex for a simple square
            Assert.That(result.Length, Is.EqualTo(4));
            Assert.That(result.All(r => r == "Convex"), Is.True, 
                "All vertices of a square should be classified as convex");
        }

        [Test]
        public void TestConcavePolygonClassification()
        {
            // Arrange: L-shaped polygon with one concave corner
            var lShape = new List<PrototypeTestHelpers.Point>
            {
                new(0, 0),
                new(2, 0),
                new(2, 1),
                new(1, 1),  // This creates the concave corner
                new(1, 2),
                new(0, 2)
            };

            // Act
            var result = PrototypeTestHelpers.CornerClassification.ClassifyVertices(lShape);

            // Assert: Should have exactly one concave vertex
            int concaveCount = result.Count(r => r == "Concave");
            Assert.That(concaveCount, Is.EqualTo(1), "L-shape should have exactly one concave vertex");
        }

        [Test]
        public void TestPolygonOrientationCCW()
        {
            // Arrange: Counter-clockwise square
            var ccwSquare = new List<PrototypeTestHelpers.Point>
            {
                new(0, 0),
                new(1, 0),
                new(1, 1),
                new(0, 1)
            };

            // Act
            bool isCCW = PrototypeTestHelpers.CornerClassification.DetermineOrientation(ccwSquare);

            // Assert
            Assert.That(isCCW, Is.True, "Square with CCW vertex order should be detected as CCW");
        }

        [Test]
        public void TestPolygonOrientationCW()
        {
            // Arrange: Clockwise square (reverse order)
            var cwSquare = new List<PrototypeTestHelpers.Point>
            {
                new(0, 0),
                new(0, 1),
                new(1, 1),
                new(1, 0)
            };

            // Act
            bool isCCW = PrototypeTestHelpers.CornerClassification.DetermineOrientation(cwSquare);

            // Assert
            Assert.That(isCCW, Is.False, "Square with CW vertex order should be detected as CW");
        }

        [Test]
        public void TestClassifyVerticesEdgeCases()
        {
            // Test with minimum valid polygon
            var triangle = new List<PrototypeTestHelpers.Point>
            {
                new(0, 0),
                new(1, 0),
                new(0.5, 1)
            };

            var result = PrototypeTestHelpers.CornerClassification.ClassifyVertices(triangle);
            Assert.That(result.Length, Is.EqualTo(3));

            // Test with null input
            Assert.Throws<ArgumentException>(() => 
                PrototypeTestHelpers.CornerClassification.ClassifyVertices(null));

            // Test with insufficient vertices
            Assert.Throws<ArgumentException>(() => 
                PrototypeTestHelpers.CornerClassification.ClassifyVertices(new List<PrototypeTestHelpers.Point>
                {
                    new(0, 0),
                    new(1, 0)
                }));
        }

        #endregion

        #region Short Edge Corner Categorization Tests

        [Test]
        public void TestShortEdgeDetection()
        {
            // Arrange: Polygon with short orthogonal edges creating a step pattern
            var stepPattern = new List<PrototypeTestHelpers.Point>
            {
                new(0, 0),
                new(0, 0.1),
                new(0.005, 0.1),    // Short horizontal edge
                new(0.005, 0.105),  // Short vertical edge  
                new(0.05, 0.105),
                new(0.05, 0)
            };

            double shortEdgeThreshold = 0.01;
            double orthoTolerance = 1e-6;

            // Act
            var result = PrototypeTestHelpers.CornerClassification.ClassifyVerticesWithShortEdges(
                stepPattern, shortEdgeThreshold, orthoTolerance);

            // Assert: Should detect the short edge vertex
            Assert.That(result.Any(r => r == "ShortEdge"), Is.True,
                "Should detect at least one short edge vertex");
        }

        [Test]
        public void TestShortEdgeOrthogonalityRequirement()
        {
            // Arrange: Short edges that are NOT orthogonal
            var nonOrthoShort = new List<PrototypeTestHelpers.Point>
            {
                new(0, 0),
                new(0.005, 0),      // Short edge
                new(0.01, 0.005),   // Short edge but not orthogonal
                new(0.1, 0.005),
                new(0.1, 0.1),
                new(0, 0.1)
            };

            double shortEdgeThreshold = 0.01;
            double orthoTolerance = 1e-6;

            // Act
            var result = PrototypeTestHelpers.CornerClassification.ClassifyVerticesWithShortEdges(
                nonOrthoShort, shortEdgeThreshold, orthoTolerance);

            // Assert: Should NOT classify as short edge due to lack of orthogonality
            Assert.That(result.Any(r => r == "ShortEdge"), Is.False,
                "Non-orthogonal short edges should not be classified as ShortEdge");
        }

        #endregion

        #region Quadratic Bézier Curve Tests

        [Test]
        public void TestControlPointCalculation()
        {
            // Arrange: Known configuration where control point can be verified
            var start = new PrototypeTestHelpers.Point(0, 0);
            var startDir = new PrototypeTestHelpers.Point(1, 0); // Horizontal
            var end = new PrototypeTestHelpers.Point(2, 2);
            var endDir = new PrototypeTestHelpers.Point(0, 1);   // Vertical

            // Act
            var control = PrototypeTestHelpers.BezierSampling.ComputeControlPoint(
                start, startDir, end, endDir);

            // Assert: Control point should be at intersection (2, 0)
            Assert.That(control.X, Is.EqualTo(2).Within(1e-10));
            Assert.That(control.Y, Is.EqualTo(0).Within(1e-10));
        }

        [Test]
        public void TestParallelTangentLinesThrowException()
        {
            // Arrange: Parallel tangent directions
            var start = new PrototypeTestHelpers.Point(0, 0);
            var startDir = new PrototypeTestHelpers.Point(1, 0);
            var end = new PrototypeTestHelpers.Point(0, 1);
            var endDir = new PrototypeTestHelpers.Point(1, 0); // Same direction

            // Act & Assert
            Assert.Throws<Exception>(() =>
                PrototypeTestHelpers.BezierSampling.ComputeControlPoint(start, startDir, end, endDir));
        }

        [Test]
        public void TestSampleByMaxSegmentLength()
        {
            // Arrange: Simple quadratic curve
            var start = new PrototypeTestHelpers.Point(0, 0);
            var control = new PrototypeTestHelpers.Point(1, 2);
            var end = new PrototypeTestHelpers.Point(2, 0);
            double maxLength = 0.5;

            // Act
            var samples = PrototypeTestHelpers.BezierSampling.SampleByMaxSegmentLength(
                start, control, end, maxLength);

            // Assert
            Assert.That(samples.Count, Is.GreaterThan(2), "Should have more than just start and end points");
            Assert.That(samples[0].X, Is.EqualTo(start.X).Within(1e-10));
            Assert.That(samples[0].Y, Is.EqualTo(start.Y).Within(1e-10));
            Assert.That(samples[^1].X, Is.EqualTo(end.X).Within(1e-10));
            Assert.That(samples[^1].Y, Is.EqualTo(end.Y).Within(1e-10));

            // Verify segment lengths are within threshold
            for (int i = 1; i < samples.Count; i++)
            {
                var prev = samples[i - 1];
                var curr = samples[i];
                double segLength = Math.Sqrt(Math.Pow(curr.X - prev.X, 2) + Math.Pow(curr.Y - prev.Y, 2));
                Assert.That(segLength, Is.LessThanOrEqualTo(maxLength * 1.1), // Small tolerance for subdivision
                    $"Segment {i} length {segLength} exceeds threshold {maxLength}");
            }
        }

        [Test]
        public void TestSampleByMaxAngle()
        {
            // Arrange: Quadratic curve with significant curvature
            var start = new PrototypeTestHelpers.Point(0, 0);
            var control = new PrototypeTestHelpers.Point(0, 2);
            var end = new PrototypeTestHelpers.Point(2, 0);
            double maxAngle = Math.PI / 6; // 30 degrees

            // Act
            var samples = PrototypeTestHelpers.BezierSampling.SampleByMaxAngle(
                start, control, end, maxAngle);

            // Assert
            Assert.That(samples.Count, Is.GreaterThan(2));
            Assert.That(samples[0].X, Is.EqualTo(start.X).Within(1e-10));
            Assert.That(samples[0].Y, Is.EqualTo(start.Y).Within(1e-10));
            Assert.That(samples[^1].X, Is.EqualTo(end.X).Within(1e-10));
            Assert.That(samples[^1].Y, Is.EqualTo(end.Y).Within(1e-10));

            // Verify the curve is properly sampled (more points in high curvature areas)
            Assert.That(samples.Count, Is.GreaterThan(4), 
                "Curved path should require multiple subdivision steps");
        }

        [Test]
        public void TestSampleCurveWithDifferentModes()
        {
            // Arrange
            var start = new PrototypeTestHelpers.Point(0, 0);
            var control = new PrototypeTestHelpers.Point(1, 1);
            var end = new PrototypeTestHelpers.Point(2, 0);
            double scale = 1.0;

            // Act
            var lengthSamples = PrototypeTestHelpers.BezierSampling.SampleCurve(
                start, control, end, PrototypeTestHelpers.BezierSampling.SamplingMode.ByMaxSegmentLength, scale);
            var angleSamples = PrototypeTestHelpers.BezierSampling.SampleCurve(
                start, control, end, PrototypeTestHelpers.BezierSampling.SamplingMode.ByMaxAngle, scale);

            // Assert: Both should produce valid samples but potentially different counts
            Assert.That(lengthSamples.Count, Is.GreaterThan(1));
            Assert.That(angleSamples.Count, Is.GreaterThan(1));

            // Both should start and end at the same points
            Assert.That(lengthSamples[0].X, Is.EqualTo(angleSamples[0].X).Within(1e-10));
            Assert.That(lengthSamples[^1].X, Is.EqualTo(angleSamples[^1].X).Within(1e-10));
        }

        [Test]
        public void TestInvalidSamplingMode()
        {
            var start = new PrototypeTestHelpers.Point(0, 0);
            var control = new PrototypeTestHelpers.Point(1, 1);
            var end = new PrototypeTestHelpers.Point(2, 0);

            Assert.Throws<ArgumentOutOfRangeException>(() =>
                PrototypeTestHelpers.BezierSampling.SampleCurve(start, control, end, 
                    (PrototypeTestHelpers.BezierSampling.SamplingMode)999, 1.0));
        }

        #endregion

        #region Point and Vector Operation Tests

        [Test]
        public void TestPointArithmetic()
        {
            // Test Point operations
            var p1 = new PrototypeTestHelpers.Point(1, 2);
            var p2 = new PrototypeTestHelpers.Point(3, 4);

            var sum = p1 + p2;
            Assert.That(sum.X, Is.EqualTo(4));
            Assert.That(sum.Y, Is.EqualTo(6));

            var diff = p2 - p1;
            Assert.That(diff.X, Is.EqualTo(2));
            Assert.That(diff.Y, Is.EqualTo(2));

            var scaled = p1 * 2.0;
            Assert.That(scaled.X, Is.EqualTo(2));
            Assert.That(scaled.Y, Is.EqualTo(4));
        }

        [Test]
        public void TestPointLength()
        {
            var p = new PrototypeTestHelpers.Point(3, 4);
            Assert.That(p.Length(), Is.EqualTo(5.0).Within(1e-10));

            var zero = new PrototypeTestHelpers.Point(0, 0);
            Assert.That(zero.Length(), Is.EqualTo(0.0));
        }

        [Test]
        public void TestPointNormalization()
        {
            var p = new PrototypeTestHelpers.Point(3, 4);
            var normalized = p.Normalized();
            
            Assert.That(normalized.Length(), Is.EqualTo(1.0).Within(1e-10));
            Assert.That(normalized.X, Is.EqualTo(0.6).Within(1e-10));
            Assert.That(normalized.Y, Is.EqualTo(0.8).Within(1e-10));

            // Test zero vector normalization
            var zero = new PrototypeTestHelpers.Point(0, 0);
            var normalizedZero = zero.Normalized();
            Assert.That(normalizedZero.X, Is.EqualTo(0));
            Assert.That(normalizedZero.Y, Is.EqualTo(0));
        }

        #endregion
    }
}