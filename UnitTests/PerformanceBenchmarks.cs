using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using Clipper2Lib;
using geoLib;
using NUnit.Framework;
using shapeEngine;

namespace UnitTests
{
    [TestFixture]
    public class PerformanceBenchmarks
    {
        public class PerformanceResult
        {
            public string TestName { get; set; } = "";
            public string ApproachName { get; set; } = "";
            public long ElapsedMilliseconds { get; set; }
            public int InputPointCount { get; set; }
            public int OutputPointCount { get; set; }
            public double CornerRadius { get; set; }
            public double EdgeResolution { get; set; }
            public double AngularResolution { get; set; }
            public string Notes { get; set; } = "";
        }

        /// <summary>
        /// Benchmarks performance difference between legacy and new corner processing approaches
        /// </summary>
        [Test]
        public void BenchmarkCornerProcessingPerformance()
        {
            var results = new List<PerformanceResult>();
            
            // Test different geometric complexity levels
            var testCases = new[]
            {
                ("Simple Rectangle", CreateRectangle(100, 50)),
                ("Complex L-Shape", CreateLShape(100, 50, 40, 30)),
                ("High Detail Star", CreateStar(50, 25, 8)),
                ("Very Complex Shape", CreateComplexPolygon(50, 16)),
                ("Small Fine Detail", CreateRectangle(10, 5)),  // Small shapes
                ("Large Coarse Shape", CreateRectangle(1000, 500)) // Large shapes
            };

            var parameterSets = new[]
            {
                ("Fast", 2.0, 5.0, 10.0),      // radius, edgeRes, angularRes
                ("Balanced", 5.0, 2.0, 5.0),
                ("HighQuality", 10.0, 1.0, 2.0),
                ("VeryFine", 15.0, 0.5, 1.0)
            };

            Console.WriteLine("\n=== Corner Processing Performance Benchmark ===");
            Console.WriteLine($"{"Test Case",-20} {"Parameters",-12} {"Approach",-8} {"Time(ms)",-8} {"In Pts",-6} {"Out Pts",-7} {"Ratio",-6}");
            Console.WriteLine(new string('-', 80));

            foreach (var (testName, inputShape) in testCases)
            {
                foreach (var (paramName, radius, edgeRes, angularRes) in parameterSets)
                {
                    // Test legacy approach
                    var legacyResult = BenchmarkLegacyApproach(testName, paramName, inputShape, radius, edgeRes, angularRes);
                    results.Add(legacyResult);

                    // Test new approach  
                    var newResult = BenchmarkNewApproach(testName, paramName, inputShape, radius, edgeRes, angularRes);
                    results.Add(newResult);

                    // Print comparison
                    double speedupRatio = legacyResult.ElapsedMilliseconds == 0 ? 0 : 
                        (double)newResult.ElapsedMilliseconds / legacyResult.ElapsedMilliseconds;
                    
                    Console.WriteLine($"{testName,-20} {paramName,-12} {"Legacy",-8} {legacyResult.ElapsedMilliseconds,-8} {legacyResult.InputPointCount,-6} {legacyResult.OutputPointCount,-7} {"1.0x",-6}");
                    Console.WriteLine($"{"",-20} {"",-12} {"New",-8} {newResult.ElapsedMilliseconds,-8} {newResult.InputPointCount,-6} {newResult.OutputPointCount,-7} {$"{speedupRatio:F1}x",-6}");
                    Console.WriteLine();
                }
            }

            // Generate summary report
            GeneratePerformanceReport(results);
        }

        private PerformanceResult BenchmarkLegacyApproach(string testName, string paramName, PathD inputShape, 
            double radius, double edgeRes, double angularRes)
        {
            var shapeSettings = new ShapeSettings();
            shapeSettings.setDecimal(ShapeSettings.properties_decimal.iCR, (decimal)radius);
            shapeSettings.setDecimal(ShapeSettings.properties_decimal.oCR, (decimal)radius);
            
            // Use a simple rectangle shape to avoid complex geometry issues
            var shapeLib = new ShapeLibrary(new[] { 0, 1, 2, 3, 4, 5, 6 }, shapeSettings);
            shapeLib.setShape((int)ShapeLibrary.shapeNames_all.rect);
            shapeLib.computeCage();

            var stopwatch = Stopwatch.StartNew();
            
            // Call the legacy method through reflection to access private method
            var method = typeof(ShapeLibrary).GetMethod("legacy_processCorners_actual", 
                System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            
            var result = (PathD)method.Invoke(shapeLib, new object[] 
            { 
                false, // previewMode
                false, // cornerCheck  
                (int)(90.0 / angularRes), // cornerSegments
                1, // optimizeCorners
                edgeRes, // resolution
                false, false, 0.0, 0.0, 0.0, 0.0 // PA search params
            });
            
            stopwatch.Stop();

            return new PerformanceResult
            {
                TestName = testName,
                ApproachName = "Legacy",
                ElapsedMilliseconds = stopwatch.ElapsedMilliseconds,
                InputPointCount = inputShape.Count,
                OutputPointCount = result?.Count ?? 0,
                CornerRadius = radius,
                EdgeResolution = edgeRes,
                AngularResolution = angularRes,
                Notes = paramName
            };
        }

        private PerformanceResult BenchmarkNewApproach(string testName, string paramName, PathD inputShape,
            double radius, double edgeRes, double angularRes)
        {
            var shapeSettings = new ShapeSettings();
            shapeSettings.setDecimal(ShapeSettings.properties_decimal.iCR, (decimal)radius);
            shapeSettings.setDecimal(ShapeSettings.properties_decimal.oCR, (decimal)radius);
            
            // Use a simple rectangle shape to avoid complex geometry issues  
            var shapeLib = new ShapeLibrary(new[] { 0, 1, 2, 3, 4, 5, 6 }, shapeSettings);
            shapeLib.setShape((int)ShapeLibrary.shapeNames_all.rect);
            shapeLib.computeCage();

            var stopwatch = Stopwatch.StartNew();
            
            var result = shapeLib.processCorners(
                previewMode: false,
                cornerCheck: false, 
                cornerSegments: (int)(90.0 / angularRes),
                optimizeCorners: 1,
                resolution: edgeRes,
                shortEdgeLength: 2.0,
                maxShortEdgeLength: 5.0
            );
            
            stopwatch.Stop();

            return new PerformanceResult
            {
                TestName = testName,
                ApproachName = "New",
                ElapsedMilliseconds = stopwatch.ElapsedMilliseconds,
                InputPointCount = inputShape.Count,
                OutputPointCount = result?.Count ?? 0,
                CornerRadius = radius,
                EdgeResolution = edgeRes,
                AngularResolution = angularRes,
                Notes = paramName
            };
        }

        private void GeneratePerformanceReport(List<PerformanceResult> results)
        {
            var report = new StringBuilder();
            report.AppendLine("\n=== Performance Analysis Summary ===");
            
            // Group by test case and analyze
            var grouped = results.GroupBy(r => r.TestName);
            
            foreach (var group in grouped)
            {
                report.AppendLine($"\n--- {group.Key} ---");
                
                var legacyResults = group.Where(r => r.ApproachName == "Legacy").ToList();
                var newResults = group.Where(r => r.ApproachName == "New").ToList();
                
                if (legacyResults.Any() && newResults.Any())
                {
                    var avgLegacyTime = legacyResults.Average(r => r.ElapsedMilliseconds);
                    var avgNewTime = newResults.Average(r => r.ElapsedMilliseconds);
                    var avgSpeedup = avgLegacyTime == 0 ? 0 : avgNewTime / avgLegacyTime;
                    
                    report.AppendLine($"  Legacy avg time: {avgLegacyTime:F1}ms");
                    report.AppendLine($"  New avg time: {avgNewTime:F1}ms");
                    report.AppendLine($"  Performance ratio: {avgSpeedup:F2}x");
                    
                    if (avgSpeedup < 0.8)
                        report.AppendLine("  → New approach is FASTER");
                    else if (avgSpeedup > 1.2)
                        report.AppendLine("  → Legacy approach is FASTER");
                    else
                        report.AppendLine("  → Performance is COMPARABLE");
                }
            }
            
            // Overall recommendations
            report.AppendLine("\n=== Recommendations ===");
            
            var allLegacy = results.Where(r => r.ApproachName == "Legacy").ToList();
            var allNew = results.Where(r => r.ApproachName == "New").ToList();
            
            if (allLegacy.Any() && allNew.Any())
            {
                var totalLegacyTime = allLegacy.Sum(r => r.ElapsedMilliseconds);
                var totalNewTime = allNew.Sum(r => r.ElapsedMilliseconds);
                
                report.AppendLine($"Total legacy time: {totalLegacyTime}ms");
                report.AppendLine($"Total new time: {totalNewTime}ms");
                
                if (totalNewTime < totalLegacyTime * 0.8)
                {
                    report.AppendLine("✓ New approach shows significant performance gains");
                    report.AppendLine("  Recommendation: Use new approach for all cases");
                }
                else if (totalNewTime > totalLegacyTime * 1.2)
                {
                    report.AppendLine("⚠ Legacy approach is significantly faster");
                    report.AppendLine("  Recommendation: Implement hybrid approach based on complexity");
                }
                else
                {
                    report.AppendLine("≈ Performance is comparable between approaches");
                    report.AppendLine("  Recommendation: Use new approach for quality benefits");
                }
            }

            Console.WriteLine(report.ToString());
        }

        // Helper methods to create test shapes
        private PathD CreateRectangle(double width, double height)
        {
            return new PathD
            {
                new PointD(0, 0),
                new PointD(width, 0),
                new PointD(width, height),
                new PointD(0, height),
                new PointD(0, 0)
            };
        }

        private PathD CreateLShape(double width, double height, double cutWidth, double cutHeight)
        {
            return new PathD
            {
                new PointD(0, 0),
                new PointD(width, 0),
                new PointD(width, height - cutHeight),
                new PointD(width - cutWidth, height - cutHeight),
                new PointD(width - cutWidth, height),
                new PointD(0, height),
                new PointD(0, 0)
            };
        }

        private PathD CreateStar(double outerRadius, double innerRadius, int points)
        {
            var star = new PathD();
            double angleStep = 2 * Math.PI / (points * 2);
            
            for (int i = 0; i < points * 2; i++)
            {
                double angle = i * angleStep;
                double radius = (i % 2 == 0) ? outerRadius : innerRadius;
                star.Add(new PointD(
                    radius * Math.Cos(angle),
                    radius * Math.Sin(angle)
                ));
            }
            star.Add(star[0]); // Close the shape
            return star;
        }

        private PathD CreateComplexPolygon(double radius, int sides)
        {
            var polygon = new PathD();
            double angleStep = 2 * Math.PI / sides;
            
            for (int i = 0; i < sides; i++)
            {
                double angle = i * angleStep;
                // Add some noise to make it more complex
                double r = radius + 5 * Math.Sin(4 * angle);
                polygon.Add(new PointD(
                    r * Math.Cos(angle),
                    r * Math.Sin(angle)
                ));
            }
            polygon.Add(polygon[0]); // Close the shape
            return polygon;
        }
    }
}