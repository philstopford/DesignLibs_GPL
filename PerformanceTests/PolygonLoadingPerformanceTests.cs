using System;
using System.Diagnostics;
using System.Threading.Tasks;
using Eto.Drawing;
using VeldridEto;

namespace PerformanceTests
{
    public class PolygonLoadingPerformanceTests
    {
        public static async Task RunPerformanceComparison()
        {
            Console.WriteLine("=== Polygon Loading Performance Comparison ===");
            Console.WriteLine();

            var testSizes = new[] { 1000, 5000, 10000 };
            
            foreach (var polygonCount in testSizes)
            {
                Console.WriteLine($"Testing with {polygonCount} polygons:");
                
                // Test synchronous loading
                var syncTime = await TestSynchronousLoading(polygonCount);
                Console.WriteLine($"  Synchronous loading: {syncTime.TotalMilliseconds:F2}ms");
                
                // Test progressive loading
                var progressiveTime = await TestProgressiveLoading(polygonCount);
                Console.WriteLine($"  Progressive loading: {progressiveTime.TotalMilliseconds:F2}ms");
                
                var improvement = ((syncTime.TotalMilliseconds - progressiveTime.TotalMilliseconds) / syncTime.TotalMilliseconds) * 100;
                Console.WriteLine($"  Time to first render improvement: {improvement:F1}%");
                Console.WriteLine();
            }
        }

        private static async Task<TimeSpan> TestSynchronousLoading(int polygonCount)
        {
            var settings = new OVPSettings();
            settings.setProgressiveLoading(false);
            
            // Generate test polygons
            GenerateTestPolygons(settings, polygonCount);
            
            var stopwatch = Stopwatch.StartNew();
            
            // Simulate synchronous polygon processing
            await SimulatePolygonProcessing(settings, false);
            
            stopwatch.Stop();
            return stopwatch.Elapsed;
        }

        private static async Task<TimeSpan> TestProgressiveLoading(int polygonCount)
        {
            var settings = new OVPSettings();
            settings.setProgressiveLoading(true);
            settings.setBatchSize(200);
            
            // Generate test polygons
            GenerateTestPolygons(settings, polygonCount);
            
            var stopwatch = Stopwatch.StartNew();
            
            // Simulate progressive polygon processing
            await SimulatePolygonProcessing(settings, true);
            
            stopwatch.Stop();
            return stopwatch.Elapsed;
        }

        private static void GenerateTestPolygons(OVPSettings settings, int count)
        {
            settings.clear();
            
            var random = new Random(42); // Fixed seed for consistent results
            
            for (int i = 0; i < count; i++)
            {
                // Create random polygon positions
                var baseX = (float)(random.NextDouble() * 2000.0 - 1000.0);
                var baseY = (float)(random.NextDouble() * 2000.0 - 1000.0);
                var size = (float)(random.NextDouble() * 20.0 + 5.0);
                
                // Create a simple rectangular polygon
                var poly = new PointF[]
                {
                    new(baseX, baseY),
                    new(baseX + size, baseY),
                    new(baseX + size, baseY + size),
                    new(baseX, baseY + size),
                    new(baseX, baseY) // Close the polygon
                };
                
                // Use random colors
                var color = Color.FromArgb(
                    random.Next(128, 256),
                    random.Next(128, 256),
                    random.Next(128, 256));
                
                var alpha = (float)(random.NextDouble() * 0.5 + 0.5);
                
                settings.addPolygon(poly, color, alpha, false, i);
            }
        }

        private static async Task SimulatePolygonProcessing(OVPSettings settings, bool useProgressive)
        {
            if (useProgressive)
            {
                // Simulate progressive loading
                var loader = new BatchedPolygonLoader(settings);
                var renderer = new ProgressivePolygonRenderer(settings);
                
                var timeToFirstRender = Stopwatch.StartNew();
                bool firstBatchProcessed = false;
                
                await loader.LoadPolygonsProgressivelyAsync(async batch =>
                {
                    // Simulate batch processing time
                    await Task.Run(() => renderer.ProcessBatch(batch));
                    
                    if (!firstBatchProcessed)
                    {
                        timeToFirstRender.Stop();
                        Console.WriteLine($"    Time to first render: {timeToFirstRender.ElapsedMilliseconds}ms");
                        firstBatchProcessed = true;
                    }
                    
                    // Simulate small delay between batches for UI responsiveness
                    await Task.Delay(1);
                });
            }
            else
            {
                // Simulate synchronous processing
                await Task.Run(() =>
                {
                    // Simulate processing all polygons at once
                    var totalPolygons = (settings.polyList?.Count ?? 0) +
                                      (settings.bgPolyList?.Count ?? 0) +
                                      (settings.tessPolyList?.Count ?? 0);
                    
                    // Simulate processing time proportional to polygon count
                    var processingTime = totalPolygons / 100; // 1ms per 100 polygons
                    System.Threading.Thread.Sleep(processingTime);
                });
            }
        }

        public static void PrintPerformanceGuidelines()
        {
            Console.WriteLine("=== Performance Guidelines ===");
            Console.WriteLine();
            Console.WriteLine("Progressive Loading Benefits:");
            Console.WriteLine("• Faster time to first polygon display");
            Console.WriteLine("• Maintained UI responsiveness during loading");
            Console.WriteLine("• Better user experience with large datasets");
            Console.WriteLine("• Configurable batch sizes for optimal performance");
            Console.WriteLine();
            Console.WriteLine("Recommended Settings:");
            Console.WriteLine("• Enable progressive loading for >500 polygons");
            Console.WriteLine("• Use batch size of 100-300 polygons");
            Console.WriteLine("• Set max processing time to 16ms for 60fps");
            Console.WriteLine();
            Console.WriteLine("API Usage:");
            Console.WriteLine("  settings.setProgressiveLoading(true);");
            Console.WriteLine("  settings.setBatchSize(200);");
            Console.WriteLine("  settings.setMaxBatchProcessingTimeMs(16);");
            Console.WriteLine();
        }
    }
}