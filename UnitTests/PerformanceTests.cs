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

        // Run multiple iterations for more reliable results
        double systemTotal = 0, mtTotal = 0, cryptoTotal = 0;
        const int runs = 3;

        for (int run = 0; run < runs; run++)
        {
            // Benchmark System RNG
            systemTotal += BenchmarkSystemRNG();

            // Benchmark Mersenne Twister RNG
            mtTotal += BenchmarkMersenneTwisterRNG();

            // Benchmark Crypto RNG
            cryptoTotal += BenchmarkCryptoRNG();
        }

        double systemTime = systemTotal / runs;
        double mtTime = mtTotal / runs;
        double cryptoTime = cryptoTotal / runs;

        Console.WriteLine($"System RNG: {systemTime:F2} ms ({BenchmarkIterations / systemTime * 1000:F0} ops/sec)");
        Console.WriteLine($"Mersenne Twister RNG: {mtTime:F2} ms ({BenchmarkIterations / mtTime * 1000:F0} ops/sec)");
        Console.WriteLine($"Crypto RNG: {cryptoTime:F2} ms ({BenchmarkIterations / cryptoTime * 1000:F0} ops/sec)");

        // All tests should complete without errors
        Assert.Pass("Performance benchmarks completed successfully");
    }

    [Test]
    public static void BenchmarkUtilityPerformance()
    {
        Console.WriteLine("=== Utility Performance Benchmarks ===");

        // Run multiple iterations for more reliable results
        double powTotal = 0, hashTotal = 0, frexpTotal = 0;
        const int runs = 3;

        for (int run = 0; run < runs; run++)
        {
            // Benchmark myPow function
            powTotal += BenchmarkMyPow();

            // Benchmark hash functions
            hashTotal += BenchmarkHashFunctions();

            // Benchmark FrExp function
            frexpTotal += BenchmarkFrExp();
        }

        double powTime = powTotal / runs;
        double hashTime = hashTotal / runs;
        double frexpTime = frexpTotal / runs;

        Console.WriteLine($"myPow function: {powTime:F2} ms ({BenchmarkIterations / powTime * 1000:F0} ops/sec)");
        Console.WriteLine($"Hash functions: {hashTime:F2} ms ({BenchmarkIterations / hashTime * 1000:F0} ops/sec)");
        Console.WriteLine($"FrExp function: {frexpTime:F2} ms ({BenchmarkIterations / frexpTime * 1000:F0} ops/sec)");

        Assert.Pass("Utility performance benchmarks completed successfully");
    }

    [Test]
    public static void BenchmarkAdvancedMathPerformance()
    {
        Console.WriteLine("=== Advanced Math Performance Benchmarks ===");

        // Benchmark FastMath.FastPow
        var fastPowTime = BenchmarkFastPow();
        Console.WriteLine($"FastMath.FastPow: {fastPowTime:F2} ms ({BenchmarkIterations / fastPowTime * 1000:F0} ops/sec)");

        // Benchmark BoxMuller generation
        var boxMullerTime = BenchmarkBoxMuller();
        Console.WriteLine($"BoxMuller generation: {boxMullerTime:F2} ms ({BenchmarkIterations / boxMullerTime * 1000:F0} ops/sec)");

        Assert.Pass("Advanced math performance benchmarks completed successfully");
    }

    private static double BenchmarkFastPow()
    {
        var sw = Stopwatch.StartNew();
        var random = new Random(42);

        // Pre-generate test cases for FastMath.FastPow
        var testCases = new (float baseVal, int exp)[BenchmarkIterations];
        for (int i = 0; i < BenchmarkIterations; i++)
        {
            testCases[i] = ((float)(random.NextDouble() * 10), random.Next(-8, 16));
        }

        sw.Restart();
        for (int i = 0; i < BenchmarkIterations; i++)
        {
            var (baseVal, exp) = testCases[i];
            var result = FastMath.FastPow(baseVal, exp);
            // Use result to prevent optimization
            if (result > 1000000.0f) Console.Write("");
        }

        sw.Stop();
        return sw.Elapsed.TotalMilliseconds;
    }

    private static double BenchmarkBoxMuller()
    {
        var sw = Stopwatch.StartNew();
        var random = new Random(42);

        // Pre-generate uniform random values
        var u1Values = new double[BenchmarkIterations];
        var u2Values = new double[BenchmarkIterations];
        for (int i = 0; i < BenchmarkIterations; i++)
        {
            u1Values[i] = random.NextDouble();
            u2Values[i] = random.NextDouble() * 0.99 + 0.01; // Avoid values too close to 0
        }

        sw.Restart();
        for (int i = 0; i < BenchmarkIterations; i++)
        {
            var result = BoxMullerOptimized.GenerateGaussianPair(u1Values[i], u2Values[i]);
            // Use result to prevent optimization
            if (result.Item1 > 1000.0) Console.Write("");
        }

        sw.Stop();
        return sw.Elapsed.TotalMilliseconds;
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

        // Performance optimization: Pre-generate test cases that will benefit from optimizations
        var testCases = new (double baseVal, int exp)[BenchmarkIterations];
        for (int i = 0; i < BenchmarkIterations; i++)
        {
            // Mix of optimized cases and general cases
            switch (i % 10)
            {
                case 0: testCases[i] = (2.0, random.Next(0, 32)); break; // Powers of 2 (lookup table)
                case 1: testCases[i] = (10.0, random.Next(0, 16)); break; // Powers of 10 (lookup table)
                case 2: testCases[i] = (Math.E, random.Next(0, 16)); break; // Powers of e (lookup table)
                case 3: testCases[i] = (0.5, random.Next(0, 16)); break; // Powers of 0.5 (lookup table)
                case 4: testCases[i] = (random.NextDouble() * 10, 2); break; // Squares (fast path)
                case 5: testCases[i] = (random.NextDouble() * 10, 3); break; // Cubes (fast path)
                case 6: testCases[i] = (random.NextDouble() * 10, 4); break; // Fourth power (fast path)
                case 7: testCases[i] = (random.NextDouble() * 10, 0); break; // Power of 0 (fast path)
                case 8: testCases[i] = (random.NextDouble() * 10, 1); break; // Power of 1 (fast path)
                default: testCases[i] = (random.NextDouble() * 10, random.Next(-5, 16)); break; // General case
            }
        }

        sw.Restart();
        for (int i = 0; i < BenchmarkIterations; i++)
        {
            var (baseVal, exp) = testCases[i];
            var result = Utils.myPow(baseVal, exp);
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

        // Performance optimization: Pre-create test strings to avoid repeated concatenation
        var testStrings = new string[BenchmarkIterations / 3];
        for (int i = 0; i < testStrings.Length; i++)
        {
            testStrings[i] = testData + i.ToString();
        }

        for (int i = 0; i < testStrings.Length; i++)
        {
            var testString = testStrings[i];
            var md5 = Utils.GetMD5Hash(testString);
            var sha1 = Utils.GetSHA1Hash(testString);
            var sha256 = Utils.GetSHA256Hash(testString);
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