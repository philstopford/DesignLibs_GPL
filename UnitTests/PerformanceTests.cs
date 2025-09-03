using System;
using System.Diagnostics;
using System.Threading.Tasks;
using entropyRNG;
using utility;
using NUnit.Framework;

namespace UnitTests;

/// <summary>
/// Performance benchmarking tests to measure optimization improvements
/// </summary>
public class PerformanceTests
{
    private const int WarmupIterations = 1000;
    private const int BenchmarkIterations = 50000;
    
    [Test]
    public static void BenchmarkRNGPerformance()
    {
        Console.WriteLine("=== RNG Performance Benchmarks ===");
        
        // Warm up
        WarmupRNG();
        
        // Benchmark System RNG
        var systemTime = BenchmarkSystemRNG();
        Console.WriteLine($"System RNG: {systemTime:F2} ms ({BenchmarkIterations / systemTime * 1000:F0} ops/sec)");
        
        // Benchmark Mersenne Twister RNG
        var mtTime = BenchmarkMersenneTwisterRNG();
        Console.WriteLine($"Mersenne Twister RNG: {mtTime:F2} ms ({BenchmarkIterations / mtTime * 1000:F0} ops/sec)");
        
        // Benchmark Crypto RNG
        var cryptoTime = BenchmarkCryptoRNG();
        Console.WriteLine($"Crypto RNG: {cryptoTime:F2} ms ({BenchmarkIterations / cryptoTime * 1000:F0} ops/sec)");
        
        // All tests should complete without errors
        Assert.Pass("Performance benchmarks completed successfully");
    }
    
    [Test]
    public static void BenchmarkUtilityPerformance()
    {
        Console.WriteLine("=== Utility Performance Benchmarks ===");
        
        // Benchmark myPow function
        var powTime = BenchmarkMyPow();
        Console.WriteLine($"myPow function: {powTime:F2} ms ({BenchmarkIterations / powTime * 1000:F0} ops/sec)");
        
        // Benchmark hash functions
        var hashTime = BenchmarkHashFunctions();
        Console.WriteLine($"Hash functions: {hashTime:F2} ms ({BenchmarkIterations / hashTime * 1000:F0} ops/sec)");
        
        // Benchmark FrExp function
        var frexpTime = BenchmarkFrExp();
        Console.WriteLine($"FrExp function: {frexpTime:F2} ms ({BenchmarkIterations / frexpTime * 1000:F0} ops/sec)");
        
        Assert.Pass("Utility performance benchmarks completed successfully");
    }
    
    private static void WarmupRNG()
    {
        for (int i = 0; i < WarmupIterations; i++)
        {
            RNG.random_gauss();
            MersenneTwister_RNG.random_gauss();
            Crypto_RNG.random_gauss();
        }
    }
    
    private static double BenchmarkSystemRNG()
    {
        var sw = Stopwatch.StartNew();
        
        for (int i = 0; i < BenchmarkIterations; i++)
        {
            var result = RNG.random_gauss3();
            // Use result to prevent optimization
            if (result[0] > 1000.0) Console.Write("");
        }
        
        sw.Stop();
        return sw.Elapsed.TotalMilliseconds;
    }
    
    private static double BenchmarkMersenneTwisterRNG()
    {
        var sw = Stopwatch.StartNew();
        
        for (int i = 0; i < BenchmarkIterations; i++)
        {
            var result = MersenneTwister_RNG.random_gauss3();
            // Use result to prevent optimization
            if (result[0] > 1000.0) Console.Write("");
        }
        
        sw.Stop();
        return sw.Elapsed.TotalMilliseconds;
    }
    
    private static double BenchmarkCryptoRNG()
    {
        var sw = Stopwatch.StartNew();
        
        for (int i = 0; i < BenchmarkIterations; i++)
        {
            var result = Crypto_RNG.random_gauss3();
            // Use result to prevent optimization
            if (result[0] > 1000.0) Console.Write("");
        }
        
        sw.Stop();
        return sw.Elapsed.TotalMilliseconds;
    }
    
    private static double BenchmarkMyPow()
    {
        var sw = Stopwatch.StartNew();
        var random = new Random(42);
        
        for (int i = 0; i < BenchmarkIterations; i++)
        {
            var baseValue = random.NextDouble() * 10;
            var exponent = random.Next(0, 16);
            var result = Utils.myPow(baseValue, exponent);
            // Use result to prevent optimization
            if (result > 1000000.0) Console.Write("");
        }
        
        sw.Stop();
        return sw.Elapsed.TotalMilliseconds;
    }
    
    private static double BenchmarkHashFunctions()
    {
        var sw = Stopwatch.StartNew();
        var testData = "Performance test data for hashing";
        
        for (int i = 0; i < BenchmarkIterations / 3; i++)
        {
            var md5 = Utils.GetMD5Hash(testData + i);
            var sha1 = Utils.GetSHA1Hash(testData + i);
            var sha256 = Utils.GetSHA256Hash(testData + i);
            // Use results to prevent optimization
            if (md5.Length == 0 || sha1.Length == 0 || sha256.Length == 0) Console.Write("");
        }
        
        sw.Stop();
        return sw.Elapsed.TotalMilliseconds;
    }
    
    private static double BenchmarkFrExp()
    {
        var sw = Stopwatch.StartNew();
        var random = new Random(42);
        
        for (int i = 0; i < BenchmarkIterations; i++)
        {
            var value = random.NextDouble() * 1000;
            var result = FrExp.calculate(value);
            // Use result to prevent optimization
            if (result.exponent > 1000) Console.Write("");
        }
        
        sw.Stop();
        return sw.Elapsed.TotalMilliseconds;
    }
}