using System;
using System.Threading.Tasks;
using PerformanceTests;

namespace PerformanceTests
{
    class Program
    {
        static async Task Main(string[] args)
        {
            Console.WriteLine("Progressive Polygon Loading Performance Tests");
            Console.WriteLine("============================================");
            Console.WriteLine();

            try
            {
                // Print guidelines
                PolygonLoadingPerformanceTests.PrintPerformanceGuidelines();
                
                // Run performance comparison
                await PolygonLoadingPerformanceTests.RunPerformanceComparison();
                
                Console.WriteLine("Performance testing completed!");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error during performance testing: {ex.Message}");
                Console.WriteLine(ex.StackTrace);
            }
            
            Console.WriteLine();
            Console.WriteLine("Press any key to exit...");
            Console.ReadKey();
        }
    }
}