# utility API Documentation

## Overview

The **utility** library provides general-purpose utility functions including mathematical operations, data compression, hashing, statistical analysis, and performance-optimized algorithms. It serves as a foundational library for common operations across the DesignLibs ecosystem.

## Key Features

- **Mathematical Operations** - Optimized power functions, fast trigonometry, precision operations
- **Hashing Functions** - MD5, SHA1, SHA256 with thread-safe implementations
- **Compression** - Data compression and decompression utilities
- **Statistical Analysis** - Histogram generation and analysis
- **Collection Extensions** - List shuffling and manipulation
- **Performance Optimization** - Fast math operations and lookup tables

## Namespace
```csharp
using utility;
```

## Core Classes and Functions

### Utils (Static Class)
Main utility class providing mathematical, hashing, and compression functions.

#### Number Formatting
```csharp
// Friendly number formatting
public static string friendlyNumber(int number)
public static string friendlyNumber(long number)
public static string friendlyNumber(double number)
```

#### Usage Example
```csharp
Console.WriteLine(Utils.friendlyNumber(1500));      // "1.5 thousand"
Console.WriteLine(Utils.friendlyNumber(2500000));   // "2.5 million"
Console.WriteLine(Utils.friendlyNumber(3.5e9));     // "3.5 billion"
```

#### Optimized Mathematical Operations
```csharp
// High-performance power function with lookup tables
public static double myPow(double num, int exp)
```

#### Usage Example
```csharp
// Optimized for common bases (2, 10, e, 0.5)
double result1 = Utils.myPow(2.0, 10);    // Fast lookup for powers of 2
double result2 = Utils.myPow(10.0, 6);    // Fast lookup for powers of 10
double result3 = Utils.myPow(Math.E, 3);  // Fast lookup for powers of e
double result4 = Utils.myPow(1.5, 4);     // General case calculation

// Performance comparison
var stopwatch = Stopwatch.StartNew();
for (int i = 0; i < 1000000; i++)
{
    double result = Utils.myPow(2.0, 10);  // Optimized
}
stopwatch.Stop();
Console.WriteLine($"Optimized: {stopwatch.ElapsedMilliseconds} ms");

stopwatch.Restart();
for (int i = 0; i < 1000000; i++)
{
    double result = Math.Pow(2.0, 10);     // Standard
}
stopwatch.Stop();
Console.WriteLine($"Standard: {stopwatch.ElapsedMilliseconds} ms");
```

#### Hash Functions
Thread-safe hashing with performance optimization.

```csharp
// Hash function methods
public static string getMD5Hash(string input)
public static string getSHA1Hash(string input)
public static string getSHA256Hash(string input)

// Byte array variants
public static byte[] getMD5HashBytes(byte[] input)
public static byte[] getSHA1HashBytes(byte[] input)
public static byte[] getSHA256HashBytes(byte[] input)
```

#### Usage Example
```csharp
string data = "Hello, World!";

// Generate different hash types
string md5Hash = Utils.getMD5Hash(data);
string sha1Hash = Utils.getSHA1Hash(data);
string sha256Hash = Utils.getSHA256Hash(data);

Console.WriteLine($"MD5: {md5Hash}");
Console.WriteLine($"SHA1: {sha1Hash}");
Console.WriteLine($"SHA256: {sha256Hash}");

// Hash byte arrays
byte[] binaryData = Encoding.UTF8.GetBytes(data);
byte[] hashBytes = Utils.getMD5HashBytes(binaryData);
string hashString = Convert.ToHexString(hashBytes);
```

#### Data Compression
```csharp
// Compression methods
public static byte[] compressData(byte[] data)
public static byte[] decompressData(byte[] compressedData)
public static string compressString(string text)
public static string decompressString(string compressedText)
```

#### Usage Example
```csharp
// Compress string data
string originalText = "This is a long text that will benefit from compression...";
string compressed = Utils.compressString(originalText);
string decompressed = Utils.decompressString(compressed);

Console.WriteLine($"Original length: {originalText.Length}");
Console.WriteLine($"Compressed length: {compressed.Length}");
Console.WriteLine($"Compression ratio: {(double)compressed.Length / originalText.Length:P2}");

// Compress binary data
byte[] originalData = File.ReadAllBytes("largefile.dat");
byte[] compressedData = Utils.compressData(originalData);
byte[] decompressedData = Utils.decompressData(compressedData);

Console.WriteLine($"Binary compression ratio: {(double)compressedData.Length / originalData.Length:P2}");
```

### FastMath (Static Class)
High-performance mathematical operations using lookup tables and optimized algorithms.

#### Fast Trigonometric Functions
```csharp
// Fast trigonometric operations using lookup tables
public static float FastSin(float x)    // Fast sine approximation
public static float FastCos(float x)    // Fast cosine approximation
public static float FastTan(float x)    // Fast tangent approximation
```

#### Fast Square Root Operations
```csharp
public static float FastSqrt(float x)           // Optimized square root
public static float FastInverseSqrt(float x)    // Fast inverse square root
```

#### Usage Example
```csharp
// Performance comparison for trigonometric functions
const int iterations = 1000000;
float[] angles = new float[iterations];
for (int i = 0; i < iterations; i++)
{
    angles[i] = (float)(i * Math.PI / 180.0); // Convert to radians
}

// Fast trigonometry (using lookup tables)
var stopwatch = Stopwatch.StartNew();
for (int i = 0; i < iterations; i++)
{
    float result = FastMath.FastSin(angles[i]);
}
stopwatch.Stop();
Console.WriteLine($"FastMath.FastSin: {stopwatch.ElapsedMilliseconds} ms");

// Standard trigonometry
stopwatch.Restart();
for (int i = 0; i < iterations; i++)
{
    float result = MathF.Sin(angles[i]);
}
stopwatch.Stop();
Console.WriteLine($"MathF.Sin: {stopwatch.ElapsedMilliseconds} ms");

// Fast square root
float value = 16.0f;
float fastSqrt = FastMath.FastSqrt(value);
float fastInvSqrt = FastMath.FastInverseSqrt(value);

Console.WriteLine($"FastSqrt({value}) = {fastSqrt}");
Console.WriteLine($"FastInverseSqrt({value}) = {fastInvSqrt}");
Console.WriteLine($"1/FastSqrt({value}) = {1.0f / fastSqrt}");
```

### Histo (Histogram Class)
Statistical histogram for data analysis and visualization.

#### Constructor
```csharp
// Create histogram with specified bins and range
Histo(int numBins, double minValue, double maxValue)
```

#### Methods
```csharp
// Data input
public void Update(double value)               // Add single data point
public void Update(double[] values)           // Add multiple data points
public void UpdateRange(IEnumerable<double> values) // Add collection of values

// Data analysis  
public int GetCount(int binIndex)             // Get count for specific bin
public double GetBinCenter(int binIndex)      // Get center value of bin
public double[] GetBinCenters()               // Get all bin centers
public int[] GetCounts()                      // Get all bin counts

// Statistics
public int TotalCount { get; }                // Total number of data points
public double Mean { get; }                   // Calculate mean value
public double StandardDeviation { get; }      // Calculate standard deviation

// Visualization
public void Reset()                           // Clear all data
public string ToString()                      // Text representation
```

#### Usage Example
```csharp
// Create histogram for analyzing random data
var histogram = new Histo(20, 0.0, 1.0);  // 20 bins from 0 to 1

// Generate and analyze random data
var random = new Random();
for (int i = 0; i < 10000; i++)
{
    double value = random.NextDouble();     // Uniform distribution [0,1)
    histogram.Update(value);
}

// Display results
Console.WriteLine("Histogram Analysis:");
Console.WriteLine($"Total samples: {histogram.TotalCount}");
Console.WriteLine($"Mean: {histogram.Mean:F4}");
Console.WriteLine($"Std Dev: {histogram.StandardDeviation:F4}");

// Display histogram
Console.WriteLine("\nHistogram:");
double[] centers = histogram.GetBinCenters();
int[] counts = histogram.GetCounts();

for (int i = 0; i < centers.Length; i++)
{
    string bar = new string('*', counts[i] / 50); // Scale for display
    Console.WriteLine($"[{centers[i]:F2}] {counts[i]:D4} {bar}");
}
```

#### Advanced Histogram Usage
```csharp
// Analyze performance measurement data
var perfHistogram = new Histo(50, 0.0, 1000.0);  // 50 bins, 0-1000ms range

// Collect timing data
var stopwatch = new Stopwatch();
for (int test = 0; test < 1000; test++)
{
    stopwatch.Restart();
    
    // Your operation to measure
    SomeOperation();
    
    stopwatch.Stop();
    perfHistogram.Update(stopwatch.Elapsed.TotalMilliseconds);
}

// Analyze performance distribution
Console.WriteLine("Performance Analysis:");
Console.WriteLine($"Average time: {perfHistogram.Mean:F2} ms");
Console.WriteLine($"Std deviation: {perfHistogram.StandardDeviation:F2} ms");

// Find performance percentiles
var sortedData = /* collect all data points */;
Array.Sort(sortedData);
Console.WriteLine($"50th percentile: {sortedData[sortedData.Length / 2]:F2} ms");
Console.WriteLine($"95th percentile: {sortedData[(int)(sortedData.Length * 0.95)]:F2} ms");
Console.WriteLine($"99th percentile: {sortedData[(int)(sortedData.Length * 0.99)]:F2} ms");
```

### IListExtensions
Extension methods for IList<T> collections.

#### Shuffle Method
```csharp
// Shuffle list elements using Fisher-Yates algorithm
public static void Shuffle<T>(this IList<T> ts)
```

#### Usage Example
```csharp
// Shuffle various collection types
var numbers = new List<int> { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
var names = new string[] { "Alice", "Bob", "Charlie", "David", "Eve" };

Console.WriteLine("Original numbers: " + string.Join(", ", numbers));
numbers.Shuffle();
Console.WriteLine("Shuffled numbers: " + string.Join(", ", numbers));

Console.WriteLine("Original names: " + string.Join(", ", names));
names.Shuffle();
Console.WriteLine("Shuffled names: " + string.Join(", ", names));

// Shuffle for random sampling
var dataSet = Enumerable.Range(1, 1000).ToList();
dataSet.Shuffle();
var randomSample = dataSet.Take(100).ToList(); // Get 100 random items

// Multiple shuffles for statistical analysis
var deck = Enumerable.Range(1, 52).ToList(); // Simulated deck of cards
for (int trial = 0; trial < 10000; trial++)
{
    deck.Shuffle();
    // Analyze card distribution...
}
```

### BoxMullerOptimized
Optimized Box-Muller transformation for normal distribution random numbers.

#### Usage Example
```csharp
// Generate normally distributed random numbers
var boxMuller = new BoxMullerOptimized();
var normalValues = new List<double>();

for (int i = 0; i < 10000; i++)
{
    double value = boxMuller.NextDouble(); // Standard normal (μ=0, σ=1)
    normalValues.Add(value);
}

// Convert to specific mean and standard deviation
double targetMean = 100.0;
double targetStdDev = 15.0;

var scaledValues = normalValues.Select(x => x * targetStdDev + targetMean).ToList();

// Analyze results with histogram
var normalHisto = new Histo(30, -4.0, 4.0);
normalValues.ForEach(normalHisto.Update);

Console.WriteLine("Normal Distribution Analysis:");
Console.WriteLine($"Generated mean: {normalValues.Average():F4} (expected: 0.0)");
Console.WriteLine($"Generated std dev: {Math.Sqrt(normalValues.Select(x => x * x).Average()):F4} (expected: 1.0)");
```

## Precision and Mathematical Utilities

### FrExp Function
Extracts mantissa and exponent from floating-point numbers.

```csharp
// Extract mantissa and exponent (like C frexp function)
public static double FrExp(double value, out int exponent)
```

#### Usage Example
```csharp
double[] testValues = { 0.0, 1.0, 2.0, 0.5, 1024.0, 0.125, Math.PI };

foreach (double value in testValues)
{
    double mantissa = Utils.FrExp(value, out int exponent);
    double reconstructed = mantissa * Math.Pow(2, exponent);
    
    Console.WriteLine($"{value:F6} = {mantissa:F6} * 2^{exponent} = {reconstructed:F6}");
}

// Use in scientific calculations
double scientificValue = 6.022e23; // Avogadro's number
double mantissa = Utils.FrExp(scientificValue, out int exp);
Console.WriteLine($"Avogadro's number: {mantissa:F6} * 2^{exp}");
```

### NextAfter Function
Find the next representable floating-point value.

```csharp
// Get next representable float value
public static double NextAfter(double from, double to)
```

#### Usage Example
```csharp
double value = 1.0;
double nextUp = Utils.NextAfter(value, double.PositiveInfinity);
double nextDown = Utils.NextAfter(value, double.NegativeInfinity);

Console.WriteLine($"Value: {value:R}");
Console.WriteLine($"Next up: {nextUp:R}");
Console.WriteLine($"Next down: {nextDown:R}");
Console.WriteLine($"ULP (Unit in Last Place): {nextUp - value:E}");

// Test floating-point precision limits
double x = 1.0;
for (int i = 0; i < 10; i++)
{
    x = Utils.NextAfter(x, double.PositiveInfinity);
    Console.WriteLine($"Step {i + 1}: {x:R}");
}
```

## Performance Optimization

### Thread-Safe Design
The utility library uses thread-safe patterns for optimal performance in multi-threaded scenarios.

```csharp
// Thread-local storage for hash algorithms
private static readonly ThreadLocal<MD5> _threadLocalMD5 = new(() => MD5.Create());
private static readonly ThreadLocal<SHA1> _threadLocalSHA1 = new(() => SHA1.Create());
private static readonly ThreadLocal<SHA256> _threadLocalSHA256 = new(() => SHA256.Create());

// Concurrent dictionary for serializer caching
private static readonly ConcurrentDictionary<Type, DataContractSerializer> _serializerCache = new();
```

### Lookup Table Optimizations
Fast mathematical operations using pre-calculated lookup tables.

```csharp
// Example: Powers lookup optimization
private static readonly double[] PowersOfTwo = new double[32];
private static readonly double[] PowersOfTen = new double[16];

// Usage in myPow function automatically benefits from lookup tables
double result = Utils.myPow(2.0, 15); // Uses PowersOfTwo lookup
double result2 = Utils.myPow(10.0, 8); // Uses PowersOfTen lookup
```

### Memory Management
```csharp
// Efficient memory usage patterns
public static void ProcessLargeDataSet(IEnumerable<double> data)
{
    const int batchSize = 10000;
    var histogram = new Histo(100, 0.0, 1000.0);
    
    int processed = 0;
    foreach (double value in data)
    {
        histogram.Update(value);
        processed++;
        
        // Periodic cleanup for very large datasets
        if (processed % (batchSize * 100) == 0)
        {
            GC.Collect();
            GC.WaitForPendingFinalizers();
        }
    }
}
```

## Integration Examples

### With Statistical Analysis
```csharp
using utility;

// Comprehensive statistical analysis pipeline
public class StatisticalAnalyzer
{
    public void AnalyzePerformanceData(IEnumerable<double> measurements)
    {
        // Convert to list for multiple passes
        var data = measurements.ToList();
        
        // Create histogram
        double min = data.Min();
        double max = data.Max();
        var histogram = new Histo(50, min, max);
        data.ForEach(histogram.Update);
        
        // Generate hash for data integrity
        string dataHash = Utils.getMD5Hash(string.Join(",", data));
        
        // Generate friendly statistics
        Console.WriteLine($"Data points: {Utils.friendlyNumber(data.Count)}");
        Console.WriteLine($"Range: {min:F3} to {max:F3}");
        Console.WriteLine($"Mean: {histogram.Mean:F3}");
        Console.WriteLine($"Std Dev: {histogram.StandardDeviation:F3}");
        Console.WriteLine($"Data hash: {dataHash}");
        
        // Visualize distribution
        Console.WriteLine("\nDistribution:");
        Console.WriteLine(histogram.ToString());
    }
}
```

### With Data Processing Pipeline
```csharp
using utility;

public class DataProcessor
{
    public byte[] ProcessAndCompress(string[] inputData)
    {
        // Shuffle data for random processing order
        var dataList = inputData.ToList();
        dataList.Shuffle();
        
        // Process data
        var processed = dataList.Select(ProcessSingleItem).ToArray();
        
        // Convert to binary and compress
        string combined = string.Join("\n", processed);
        byte[] compressed = Utils.compressData(Encoding.UTF8.GetBytes(combined));
        
        // Generate verification hash
        string hash = Utils.getSHA256Hash(combined);
        Console.WriteLine($"Processed {Utils.friendlyNumber(inputData.Length)} items");
        Console.WriteLine($"Compression ratio: {(double)compressed.Length / combined.Length:P2}");
        Console.WriteLine($"Verification hash: {hash}");
        
        return compressed;
    }
    
    private string ProcessSingleItem(string item)
    {
        // Use FastMath for performance-critical calculations
        float value = float.Parse(item);
        float result = FastMath.FastSin(value) * FastMath.FastCos(value * 2);
        return result.ToString("F6");
    }
}
```

## Best Practices

### Performance Guidelines
```csharp
// Use lookup table optimizations when possible
double power = Utils.myPow(2.0, 10);        // Fast - uses lookup table
double power2 = Utils.myPow(2.1, 10);       // Slower - uses Math.Pow

// Prefer FastMath for performance-critical loops  
for (int i = 0; i < 1000000; i++)
{
    float result = FastMath.FastSin(angles[i]); // Faster than MathF.Sin
}

// Reuse hash objects via thread-local storage (automatic)
string hash1 = Utils.getMD5Hash(data1);     // Gets thread-local MD5 instance
string hash2 = Utils.getMD5Hash(data2);     // Reuses same instance
```

### Memory Management
```csharp
// For large histograms, consider memory usage
var largeHistogram = new Histo(10000, 0.0, 1.0); // Uses significant memory

// Use appropriate bin counts
var efficientHistogram = new Histo(100, 0.0, 1.0); // Usually sufficient

// Reset histograms to free memory when done
histogram.Reset();
GC.Collect(); // Force cleanup if needed
```

### Error Handling
```csharp
// Validate inputs for mathematical operations
public static double SafePower(double baseValue, int exponent)
{
    if (double.IsInfinity(baseValue) || double.IsNaN(baseValue))
    {
        throw new ArgumentException("Invalid base value");
    }
    
    if (exponent > 1000) // Prevent overflow
    {
        throw new ArgumentException("Exponent too large");
    }
    
    return Utils.myPow(baseValue, exponent);
}

// Handle compression errors
public static byte[] SafeCompress(byte[] data)
{
    if (data == null || data.Length == 0)
    {
        return Array.Empty<byte>();
    }
    
    try
    {
        return Utils.compressData(data);
    }
    catch (Exception ex)
    {
        Console.WriteLine($"Compression failed: {ex.Message}");
        return data; // Return original data if compression fails
    }
}
```

## Dependencies
- **entropyRNG** - Used by IListExtensions.Shuffle for random number generation
- **System.Security.Cryptography** - Hash function implementations
- **System.IO.Compression** - Data compression functionality
- **.NET 8.0** - Target framework

## Related Libraries
- **[MersenneTwister](MersenneTwister-API.md)** - High-quality random number generation
- **[entropyRNG](entropyRNG-API.md)** - Thread-safe RNG system used by utility functions
- **[geoWrangler](geoWrangler-API.md)** - Uses utility functions for mathematical operations
- **[Noise](Noise-API.md)** - May use FastMath functions for performance optimization