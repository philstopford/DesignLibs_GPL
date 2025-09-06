# MersenneTwister API Documentation

## Overview

The **MersenneTwister** library provides a high-quality pseudorandom number generator implementing the Mersenne Twister MT19937 algorithm. This generator is widely used in scientific computing and simulations due to its excellent statistical properties and very long period.

## Key Features

- **High-Quality Random Numbers** - Excellent statistical properties for scientific applications
- **Long Period** - Period of 2^19937 - 1, suitable for extensive simulations
- **Multiple Data Types** - Integer, float, and double precision random numbers
- **Flexible Seeding** - Single seed or array-based initialization
- **Thread-Safe Operations** - Parallel processing support for initialization

## Namespace
```csharp
using MersenneTwisterRNG;
```

## Important Notes

⚠️ **Security Warning**: MersenneTwister is designed for Monte Carlo simulations and statistical applications. It is **NOT secure for cryptographic purposes** and should never be used for security-sensitive applications like password generation or cryptographic keys.

## Core Class

### MersenneTwister
Main random number generator class implementing the MT19937 algorithm.

#### Constructors
```csharp
// Time-based seeding (uses current milliseconds)
public MersenneTwister()

// Single seed initialization  
public MersenneTwister(int seed)

// Array-based seeding (for more randomness)
public MersenneTwister(int[] init)
```

#### Properties
```csharp
// Maximum integer value that can be generated
public static int MaxRandomInt { get; }  // 0x7FFFFFFF (2,147,483,647)
```

#### Basic Usage Examples
```csharp
// Create generators with different seeding methods
var rng1 = new MersenneTwister();                    // Time-based seed
var rng2 = new MersenneTwister(12345);              // Fixed seed  
var rng3 = new MersenneTwister(new int[] { 1, 2, 3, 4, 5 }); // Array seed

// Generate basic random numbers
int randomInt = rng1.Next();                         // 0 to MaxRandomInt
float randomFloat = rng1.NextFloat();                // 0.0 to 1.0 (exclusive)
double randomDouble = rng1.NextDouble();             // 0.0 to 1.0 (exclusive)

Console.WriteLine($"Random integer: {randomInt}");
Console.WriteLine($"Random float: {randomFloat:F6}");
Console.WriteLine($"Random double: {randomDouble:F10}");
```

## Integer Random Numbers

### Basic Integer Generation
```csharp
// Generate integers in various ranges
public int Next()                           // [0, MaxRandomInt]
public int Next(int maxValue)              // [0, maxValue]  
public int Next(int minValue, int maxValue) // [minValue, maxValue]
```

#### Usage Examples
```csharp
var rng = new MersenneTwister(42);

// Generate integers in different ranges
int fullRange = rng.Next();                 // 0 to 2,147,483,647
int diceRoll = rng.Next(1, 7);             // 1 to 6 (inclusive)
int percentage = rng.Next(0, 101);          // 0 to 100 (inclusive)
int negativeRange = rng.Next(-100, 100);   // -100 to 100

// Generate arrays of random integers
var randomInts = new int[1000];
for (int i = 0; i < randomInts.Length; i++)
{
    randomInts[i] = rng.Next(0, 1000);
}

// Statistical analysis
double mean = randomInts.Average();
Console.WriteLine($"Mean of 1000 random integers [0,999]: {mean:F2}");
Console.WriteLine($"Expected mean: ~499.5");
```

### Specialized Integer Applications
```csharp
// Random sampling without replacement
public static List<int> RandomSample(MersenneTwister rng, int population, int sampleSize)
{
    if (sampleSize > population)
        throw new ArgumentException("Sample size cannot exceed population size");
    
    var indices = Enumerable.Range(0, population).ToList();
    var sample = new List<int>();
    
    for (int i = 0; i < sampleSize; i++)
    {
        int randomIndex = rng.Next(0, indices.Count);
        sample.Add(indices[randomIndex]);
        indices.RemoveAt(randomIndex);
    }
    
    return sample;
}

// Usage
var rng = new MersenneTwister();
var randomSample = RandomSample(rng, 1000, 50); // Sample 50 items from 1000
```

## Floating-Point Random Numbers

### Float Generation
```csharp
// Single-precision floating-point methods
public float NextFloat()                    // [0.0, 1.0) - excludes 1.0
public float NextFloat(bool includeOne)     // [0.0, 1.0] if true, [0.0, 1.0) if false
public float NextFloatPositive()           // (0.0, 1.0) - excludes both 0.0 and 1.0
```

### Double Generation
```csharp  
// Double-precision floating-point methods
public double NextDouble()                  // [0.0, 1.0) - excludes 1.0
public double NextDouble(bool includeOne)   // [0.0, 1.0] if true, [0.0, 1.0) if false
public double NextDoublePositive()         // (0.0, 1.0) - excludes both 0.0 and 1.0
```

#### Usage Examples
```csharp
var rng = new MersenneTwister(123);

// Different floating-point ranges
float standard = rng.NextFloat();           // [0.0, 1.0)
float inclusive = rng.NextFloat(true);      // [0.0, 1.0] 
float positive = rng.NextFloatPositive();   // (0.0, 1.0)

double preciseStd = rng.NextDouble();       // Higher precision
double preciseInc = rng.NextDouble(true);   
double precisePos = rng.NextDoublePositive();

Console.WriteLine($"Standard float: {standard:F10}");
Console.WriteLine($"Inclusive float: {inclusive:F10}");  
Console.WriteLine($"Positive float: {positive:F10}");
Console.WriteLine($"Precise double: {preciseStd:F15}");
```

### Custom Range Generation
```csharp
// Extension methods for custom ranges
public static class MersenneTwisterExtensions
{
    public static double NextDouble(this MersenneTwister rng, double min, double max)
    {
        return rng.NextDouble() * (max - min) + min;
    }
    
    public static float NextFloat(this MersenneTwister rng, float min, float max)
    {
        return rng.NextFloat() * (max - min) + min;
    }
    
    public static double NextGaussian(this MersenneTwister rng, double mean = 0.0, double stdDev = 1.0)
    {
        // Box-Muller transformation for normal distribution
        if (rng._hasSpareGaussian)
        {
            rng._hasSpareGaussian = false;
            return rng._spareGaussian * stdDev + mean;
        }
        
        rng._hasSpareGaussian = true;
        double u = rng.NextDoublePositive();
        double v = rng.NextDouble() * 2.0 * Math.PI;
        double mag = stdDev * Math.Sqrt(-2.0 * Math.Log(u));
        
        rng._spareGaussian = mag * Math.Sin(v);
        return mag * Math.Cos(v) + mean;
    }
    
    // Private fields for Gaussian generation (would need to be in actual extension)
    private static bool _hasSpareGaussian = false;
    private static double _spareGaussian = 0.0;
}

// Usage
var rng = new MersenneTwister();

double temperature = rng.NextDouble(-10.0, 40.0);     // Temperature range
float probability = rng.NextFloat(0.0f, 1.0f);       // Probability
double normalValue = rng.NextGaussian(100.0, 15.0);  // IQ distribution
```

## Statistical Applications

### Monte Carlo Simulations
```csharp
public class MonteCarloSimulation
{
    private readonly MersenneTwister _rng;
    
    public MonteCarloSimulation(int seed = 0)
    {
        _rng = seed == 0 ? new MersenneTwister() : new MersenneTwister(seed);
    }
    
    // Estimate π using random points in a circle
    public double EstimatePi(int iterations = 1000000)
    {
        int pointsInCircle = 0;
        
        for (int i = 0; i < iterations; i++)
        {
            double x = _rng.NextDouble(-1.0, 1.0);
            double y = _rng.NextDouble(-1.0, 1.0);
            
            if (x * x + y * y <= 1.0)
                pointsInCircle++;
        }
        
        return 4.0 * pointsInCircle / iterations;
    }
    
    // Simulate random walk
    public List<Point> RandomWalk2D(int steps, double stepSize = 1.0)
    {
        var path = new List<Point> { new Point(0, 0) };
        double x = 0, y = 0;
        
        for (int i = 0; i < steps; i++)
        {
            double angle = _rng.NextDouble() * 2.0 * Math.PI;
            x += stepSize * Math.Cos(angle);
            y += stepSize * Math.Sin(angle);
            path.Add(new Point(x, y));
        }
        
        return path;
    }
}

// Usage
var simulation = new MonteCarloSimulation(12345);
double piEstimate = simulation.EstimatePi(10000000);
Console.WriteLine($"π estimate: {piEstimate:F6} (actual: {Math.PI:F6})");

var walkPath = simulation.RandomWalk2D(1000);
Console.WriteLine($"Random walk ended at: ({walkPath.Last().X:F2}, {walkPath.Last().Y:F2})");
```

### Statistical Quality Testing
```csharp
public static class StatisticalTests
{
    public static void UniformityTest(MersenneTwister rng, int samples = 100000)
    {
        var values = new double[samples];
        for (int i = 0; i < samples; i++)
        {
            values[i] = rng.NextDouble();
        }
        
        // Basic statistics
        double mean = values.Average();
        double variance = values.Select(x => Math.Pow(x - mean, 2)).Average();
        double stdDev = Math.Sqrt(variance);
        
        Console.WriteLine("Uniformity Test Results:");
        Console.WriteLine($"  Mean: {mean:F6} (expected: 0.5)");
        Console.WriteLine($"  Std Dev: {stdDev:F6} (expected: ~0.289)");
        Console.WriteLine($"  Min: {values.Min():F6}");
        Console.WriteLine($"  Max: {values.Max():F6}");
        
        // Chi-square test for uniformity
        const int bins = 10;
        var histogram = new int[bins];
        
        foreach (double value in values)
        {
            int bin = (int)(value * bins);
            if (bin >= bins) bin = bins - 1; // Handle edge case
            histogram[bin]++;
        }
        
        double expectedCount = (double)samples / bins;
        double chiSquare = histogram.Sum(observed => Math.Pow(observed - expectedCount, 2) / expectedCount);
        
        Console.WriteLine($"  Chi-square statistic: {chiSquare:F2}");
        Console.WriteLine($"  Critical value (α=0.05, df=9): 16.919");
        Console.WriteLine($"  Test result: {(chiSquare < 16.919 ? "PASS" : "FAIL")} (uniform distribution)");
    }
    
    public static void CorrelationTest(MersenneTwister rng, int samples = 10000)
    {
        var x = new double[samples];
        var y = new double[samples];
        
        for (int i = 0; i < samples; i++)
        {
            x[i] = rng.NextDouble();
            y[i] = rng.NextDouble();
        }
        
        // Calculate correlation coefficient
        double meanX = x.Average();
        double meanY = y.Average();
        
        double numerator = x.Zip(y, (xi, yi) => (xi - meanX) * (yi - meanY)).Sum();
        double denomX = x.Sum(xi => Math.Pow(xi - meanX, 2));
        double denomY = y.Sum(yi => Math.Pow(yi - meanY, 2));
        
        double correlation = numerator / Math.Sqrt(denomX * denomY);
        
        Console.WriteLine("Correlation Test Results:");
        Console.WriteLine($"  Correlation coefficient: {correlation:F6} (expected: ~0.0)");
        Console.WriteLine($"  Test result: {(Math.Abs(correlation) < 0.1 ? "PASS" : "FAIL")} (independence)");
    }
}

// Usage
var rng = new MersenneTwister(42);
StatisticalTests.UniformityTest(rng);
StatisticalTests.CorrelationTest(rng);
```

## Performance Optimization

### Reproducible Results
```csharp
// Use fixed seeds for reproducible results in testing/debugging
public class ReproducibleSimulation
{
    public static void RunComparativeTest()
    {
        const int seed = 12345;
        const int samples = 1000;
        
        // Two identical generators
        var rng1 = new MersenneTwister(seed);
        var rng2 = new MersenneTwister(seed);
        
        // Generate identical sequences
        for (int i = 0; i < samples; i++)
        {
            double val1 = rng1.NextDouble();
            double val2 = rng2.NextDouble();
            
            if (Math.Abs(val1 - val2) > 1e-15)
            {
                Console.WriteLine($"Reproducibility failed at sample {i}");
                return;
            }
        }
        
        Console.WriteLine("Reproducibility test PASSED - identical sequences generated");
    }
}
```

### Batch Generation
```csharp
public static class BatchGeneration
{
    public static double[] GenerateDoubleArray(MersenneTwister rng, int count)
    {
        var result = new double[count];
        for (int i = 0; i < count; i++)
        {
            result[i] = rng.NextDouble();
        }
        return result;
    }
    
    public static int[] GenerateIntArray(MersenneTwister rng, int count, int minValue, int maxValue)
    {
        var result = new int[count];
        for (int i = 0; i < count; i++)
        {
            result[i] = rng.Next(minValue, maxValue);
        }
        return result;
    }
    
    // Parallel generation for large datasets (use with caution - creates multiple RNG instances)
    public static double[] GenerateDoubleArrayParallel(int seed, int count)
    {
        var result = new double[count];
        var partitionSize = Math.Max(1, count / Environment.ProcessorCount);
        
        Parallel.For(0, Environment.ProcessorCount, threadIndex =>
        {
            var threadRng = new MersenneTwister(seed + threadIndex); // Different seed per thread
            int start = threadIndex * partitionSize;
            int end = Math.Min(start + partitionSize, count);
            
            for (int i = start; i < end; i++)
            {
                result[i] = threadRng.NextDouble();
            }
        });
        
        return result;
    }
}

// Usage
var rng = new MersenneTwister();
double[] randomValues = BatchGeneration.GenerateDoubleArray(rng, 1000000);
int[] randomInts = BatchGeneration.GenerateIntArray(rng, 100000, 1, 100);
```

### Memory Management
```csharp
public static class EfficientGeneration
{
    // Reuse generator instance for better performance
    private static readonly MersenneTwister _sharedRng = new MersenneTwister();
    
    public static void ProcessLargeDataset(int[] dataset)
    {
        // Process in chunks to manage memory
        const int chunkSize = 10000;
        
        for (int i = 0; i < dataset.Length; i += chunkSize)
        {
            int currentChunkSize = Math.Min(chunkSize, dataset.Length - i);
            
            // Process chunk with random operations
            for (int j = 0; j < currentChunkSize; j++)
            {
                dataset[i + j] = _sharedRng.Next(0, 1000);
            }
            
            // Optional: Trigger GC for very large datasets
            if (i % (chunkSize * 100) == 0)
            {
                GC.Collect();
                GC.WaitForPendingFinalizers();
            }
        }
    }
}
```

## Integration Examples

### With Histogram Analysis
```csharp
using MersenneTwisterRNG;
using utility; // For Histo class

public static void AnalyzeRandomDistribution()
{
    var rng = new MersenneTwister(42);
    var histogram = new Histo(50, 0.0, 1.0);
    
    // Generate data
    const int samples = 100000;
    for (int i = 0; i < samples; i++)
    {
        double value = rng.NextDouble();
        histogram.Update(value);
    }
    
    Console.WriteLine($"Histogram Analysis ({samples} samples):");
    Console.WriteLine($"Mean: {histogram.Mean:F6}");
    Console.WriteLine($"Std Dev: {histogram.StandardDeviation:F6}");
    Console.WriteLine("\nDistribution:");
    Console.WriteLine(histogram.ToString());
}
```

### With geoWrangler Random Positioning
```csharp
using MersenneTwisterRNG;
using geoWrangler;
using Clipper2Lib;

public static PathsD GenerateRandomPolygons(int count, MersenneTwister rng)
{
    var polygons = new PathsD();
    
    for (int i = 0; i < count; i++)
    {
        // Random position and size
        double x = rng.NextDouble(-500, 500);
        double y = rng.NextDouble(-500, 500);  
        double width = rng.NextDouble(10, 100);
        double height = rng.NextDouble(10, 100);
        double rotation = rng.NextDouble(0, 360);
        
        // Create random polygon
        var rect = GeoWrangler.rectangle(width, height, x, y);
        var rotated = GeoWrangler.rotate(rect, rotation);
        
        polygons.Add(rotated);
    }
    
    return polygons;
}
```

## Best Practices

### Seeding Guidelines
```csharp
// Good practices for different use cases

// For reproducible scientific simulations
var scientificRng = new MersenneTwister(12345);

// For testing with known outcomes
var testRng = new MersenneTwister(0);

// For production simulations with good randomness
var entropy = Environment.TickCount ^ DateTime.Now.Millisecond ^ GetHashCode();
var productionRng = new MersenneTwister(entropy);

// For high-entropy initialization (best randomness)
var entropyArray = new int[8];
var systemRng = new Random();
for (int i = 0; i < entropyArray.Length; i++)
{
    entropyArray[i] = systemRng.Next();
}
var highEntropyRng = new MersenneTwister(entropyArray);
```

### Thread Safety
```csharp
// MersenneTwister is NOT thread-safe - use separate instances per thread
public class ThreadSafeRandomManager
{
    private static readonly ThreadLocal<MersenneTwister> _threadLocalRng = 
        new ThreadLocal<MersenneTwister>(() => new MersenneTwister(Thread.CurrentThread.ManagedThreadId));
    
    public static MersenneTwister ThreadRng => _threadLocalRng.Value!;
    
    // Thread-safe method for getting random values
    public static double GetRandomDouble() => ThreadRng.NextDouble();
    public static int GetRandomInt(int min, int max) => ThreadRng.Next(min, max);
}

// Usage in multi-threaded scenarios
Parallel.For(0, 10000, i =>
{
    double value = ThreadSafeRandomManager.GetRandomDouble();
    int index = ThreadSafeRandomManager.GetRandomInt(0, 1000);
    // Process with thread-safe random values
});
```

### Error Handling
```csharp
public static class SafeRandomOperations
{
    public static int SafeNext(MersenneTwister rng, int minValue, int maxValue)
    {
        if (minValue > maxValue)
        {
            throw new ArgumentException("minValue cannot be greater than maxValue");
        }
        
        try
        {
            return rng.Next(minValue, maxValue);
        }
        catch (Exception ex)
        {
            Console.WriteLine($"Random generation failed: {ex.Message}");
            return minValue; // Fallback value
        }
    }
    
    public static double SafeNextDouble(MersenneTwister rng, double min, double max)
    {
        if (double.IsNaN(min) || double.IsNaN(max) || min >= max)
        {
            throw new ArgumentException("Invalid range parameters");
        }
        
        double value = rng.NextDouble();
        return value * (max - min) + min;
    }
}
```

## Technical Specifications

### Algorithm Details
- **Algorithm**: MT19937 (Mersenne Twister)
- **Period**: 2^19937 - 1 (~10^6000)
- **State Size**: 624 × 32 bits = 19,968 bits
- **Equidistribution**: 623-dimensionally equidistributed
- **Statistical Quality**: Passes most statistical tests including DIEHARD

### Performance Characteristics
- **Generation Speed**: Approximately 2-4× slower than System.Random but much higher quality
- **Memory Usage**: ~2.5KB state (624 32-bit integers)
- **Initialization Cost**: Minimal for single seed, moderate for array seeding
- **Thread Safety**: Not thread-safe (requires separate instances per thread)

## Dependencies
- **.NET 8.0** - Target framework
- **System.Threading.Tasks** - Parallel processing for array initialization

## Related Libraries
- **[entropyRNG](entropyRNG-API.md)** - Thread-safe wrapper around MersenneTwister
- **[utility](utility-API.md)** - Uses MersenneTwister for histogram generation and statistical analysis
- **[Noise](Noise-API.md)** - May use MersenneTwister as underlying RNG for noise generation