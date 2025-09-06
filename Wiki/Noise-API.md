# Noise API Documentation

## Overview

The **Noise** library provides comprehensive implementations of procedural noise functions including Perlin, Simplex, and OpenSimplex noise algorithms. These are essential for generating natural-looking patterns, textures, terrain, and random variations in computer graphics, procedural content generation, and simulation applications.

## Key Features

- **Multiple Noise Types** - Perlin, Simplex, and OpenSimplex noise implementations
- **Multi-Dimensional Support** - 2D, 3D, and 4D noise generation
- **Octave-Based Generation** - Combine multiple frequencies for complex patterns
- **High Performance** - Optimized algorithms with lookup tables and vectorization
- **Seeded Generation** - Reproducible noise patterns with seed control
- **Natural Patterns** - Ideal for terrain, clouds, textures, and organic shapes

## Namespace
```csharp
using Noise;
```

## Core Classes

### PerlinNoise
Classic Perlin noise implementation providing smooth, natural-looking 3D noise.

#### Constructor
```csharp
public PerlinNoise(int seed)
```

#### Methods
```csharp
// Generate 3D noise value
public double Noise(double x, double y, double z)

// Octave-based noise (multiple frequencies combined)
public double OctaveNoise(double x, double y, double z, int octaves, double persistence)
```

#### Usage Examples
```csharp
// Create Perlin noise generator
var perlin = new PerlinNoise(12345);

// Generate single noise values
double value1 = perlin.Noise(10.5, 20.3, 5.0);    // 3D noise
double value2 = perlin.Noise(x, y, 0);             // 2D noise (z=0)

Console.WriteLine($"Perlin noise values: {value1:F6}, {value2:F6}");

// Generate 2D heightmap
int width = 256, height = 256;
double[,] heightmap = new double[width, height];

for (int x = 0; x < width; x++)
{
    for (int y = 0; y < height; y++)
    {
        // Scale coordinates to control noise frequency
        double nx = (double)x / 64.0;  // Lower values = larger features
        double ny = (double)y / 64.0;
        
        heightmap[x, y] = perlin.Noise(nx, ny, 0);
    }
}
```

#### Multi-Octave Terrain Generation
```csharp
public static double[,] GenerateTerrainHeightmap(int width, int height, int seed, 
    int octaves = 6, double persistence = 0.5, double scale = 64.0)
{
    var perlin = new PerlinNoise(seed);
    var heightmap = new double[width, height];
    
    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {
            double amplitude = 1.0;
            double frequency = 1.0;
            double noiseValue = 0.0;
            double maxValue = 0.0;
            
            // Combine multiple octaves
            for (int octave = 0; octave < octaves; octave++)
            {
                double nx = (double)x / scale * frequency;
                double ny = (double)y / scale * frequency;
                
                noiseValue += perlin.Noise(nx, ny, 0) * amplitude;
                maxValue += amplitude;
                
                amplitude *= persistence;  // Decrease amplitude
                frequency *= 2.0;         // Increase frequency
            }
            
            // Normalize to [0, 1]
            heightmap[x, y] = noiseValue / maxValue;
        }
    }
    
    return heightmap;
}

// Usage
var terrain = GenerateTerrainHeightmap(512, 512, 42, octaves: 8, persistence: 0.6);
```

### SimplexNoise
Advanced Simplex noise with better visual characteristics and performance than Perlin.

#### Constructor
```csharp
// Default settings (feature size ~30, persistence 0.7)
public SimplexNoise(int seed)

// Custom settings
public SimplexNoise(int largestFeature, double persistence, int seed)
```

#### Methods
```csharp
// Generate 2D noise values
public double GetNoise(int x, int y)       // Integer coordinates
public double GetNoise(double x, double y)  // Floating-point coordinates
```

#### Usage Examples
```csharp
// Create Simplex noise generators with different characteristics
var smoothNoise = new SimplexNoise(seed: 123);                    // Default smooth noise
var roughNoise = new SimplexNoise(largestFeature: 16, persistence: 0.8, seed: 123);  // Rougher terrain
var fineDetail = new SimplexNoise(largestFeature: 4, persistence: 0.3, seed: 123);   // Fine details

// Generate cloud texture
int size = 256;
double[,] cloudTexture = new double[size, size];

for (int x = 0; x < size; x++)
{
    for (int y = 0; y < size; y++)
    {
        // Multiple layers for realistic clouds
        double baseLayer = smoothNoise.GetNoise(x * 0.01, y * 0.01);
        double detailLayer = fineDetail.GetNoise(x * 0.05, y * 0.05) * 0.3;
        
        cloudTexture[x, y] = Math.Max(0, baseLayer + detailLayer);
    }
}

// Generate animated noise (time-based)
public static double GetAnimatedNoise(SimplexNoise noise, double x, double y, double time)
{
    // Use time as a third dimension for animation
    return (noise.GetNoise(x, y) + noise.GetNoise(x + time, y + time)) * 0.5;
}
```

#### Fractal Noise Patterns
```csharp
public class FractalNoise
{
    private readonly SimplexNoise[] _octaves;
    private readonly double[] _amplitudes;
    private readonly double[] _frequencies;
    
    public FractalNoise(int seed, int octaveCount = 6, double persistence = 0.5, double baseFrequency = 0.01)
    {
        _octaves = new SimplexNoise[octaveCount];
        _amplitudes = new double[octaveCount];
        _frequencies = new double[octaveCount];
        
        double amplitude = 1.0;
        double frequency = baseFrequency;
        
        for (int i = 0; i < octaveCount; i++)
        {
            _octaves[i] = new SimplexNoise(seed + i);
            _amplitudes[i] = amplitude;
            _frequencies[i] = frequency;
            
            amplitude *= persistence;
            frequency *= 2.0;
        }
    }
    
    public double GetNoise(double x, double y)
    {
        double result = 0.0;
        
        for (int i = 0; i < _octaves.Length; i++)
        {
            result += _octaves[i].GetNoise(x * _frequencies[i], y * _frequencies[i]) * _amplitudes[i];
        }
        
        return result;
    }
}

// Usage for complex terrain
var fractalNoise = new FractalNoise(seed: 42, octaveCount: 8, persistence: 0.6);
double terrainHeight = fractalNoise.GetNoise(worldX, worldZ) * 100.0; // Scale to world units
```

### OpenSimplexNoise  
High-quality noise with excellent visual properties and multi-dimensional support.

#### Constructor
```csharp
public OpenSimplexNoise(long seed)
```

#### Methods
```csharp
// 2D noise
public double Evaluate(double x, double y)

// 3D noise  
public double Evaluate(double x, double y, double z)

// 4D noise (includes time dimension)
public double Evaluate(double x, double y, double z, double w)
```

#### Usage Examples
```csharp
// Create OpenSimplex noise generator
var openSimplex = new OpenSimplexNoise(seed: 789);

// Generate various noise patterns
double value2D = openSimplex.Evaluate(10.0, 20.0);
double value3D = openSimplex.Evaluate(10.0, 20.0, 5.0);
double value4D = openSimplex.Evaluate(10.0, 20.0, 5.0, timeValue);

// Generate seamless texture
public static double[,] GenerateSeamlessTexture(OpenSimplexNoise noise, int size, double scale = 1.0)
{
    var texture = new double[size, size];
    
    for (int x = 0; x < size; x++)
    {
        for (int y = 0; y < size; y++)
        {
            // Map to circle in 4D space for seamless tiling
            double nx = (double)x / size;
            double ny = (double)y / size;
            
            double rdx = Math.Cos(nx * 2 * Math.PI) * scale;
            double rdy = Math.Sin(nx * 2 * Math.PI) * scale;
            double rdz = Math.Cos(ny * 2 * Math.PI) * scale;
            double rdw = Math.Sin(ny * 2 * Math.PI) * scale;
            
            texture[x, y] = noise.Evaluate(rdx, rdy, rdz, rdw);
        }
    }
    
    return texture;
}

// 3D volume texture for volumetric effects
public static double[,,] Generate3DVolumeTexture(OpenSimplexNoise noise, int width, int height, int depth)
{
    var volume = new double[width, height, depth];
    
    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {
            for (int z = 0; z < depth; z++)
            {
                double nx = (double)x / 32.0;  // Scale factor
                double ny = (double)y / 32.0;
                double nz = (double)z / 32.0;
                
                volume[x, y, z] = noise.Evaluate(nx, ny, nz);
            }
        }
    }
    
    return volume;
}
```

## Advanced Applications

### Procedural Terrain Generation
```csharp
public class TerrainGenerator
{
    private readonly PerlinNoise _heightNoise;
    private readonly SimplexNoise _moistureNoise;
    private readonly OpenSimplexNoise _temperatureNoise;
    
    public TerrainGenerator(int seed)
    {
        _heightNoise = new PerlinNoise(seed);
        _moistureNoise = new SimplexNoise(seed + 1);
        _temperatureNoise = new OpenSimplexNoise(seed + 2);
    }
    
    public TerrainData GenerateChunk(int chunkX, int chunkZ, int chunkSize)
    {
        var heights = new double[chunkSize, chunkSize];
        var moisture = new double[chunkSize, chunkSize];
        var temperature = new double[chunkSize, chunkSize];
        
        for (int x = 0; x < chunkSize; x++)
        {
            for (int z = 0; z < chunkSize; z++)
            {
                double worldX = chunkX * chunkSize + x;
                double worldZ = chunkZ * chunkSize + z;
                
                // Multi-octave height generation
                double height = 0;
                double amplitude = 100;
                double frequency = 0.01;
                
                for (int octave = 0; octave < 6; octave++)
                {
                    height += _heightNoise.Noise(worldX * frequency, worldZ * frequency, 0) * amplitude;
                    amplitude *= 0.5;
                    frequency *= 2.0;
                }
                
                heights[x, z] = Math.Max(0, height);
                
                // Climate generation
                moisture[x, z] = (_moistureNoise.GetNoise(worldX * 0.005, worldZ * 0.005) + 1.0) * 0.5;
                temperature[x, z] = (_temperatureNoise.Evaluate(worldX * 0.008, worldZ * 0.008) + 1.0) * 0.5;
            }
        }
        
        return new TerrainData 
        { 
            Heights = heights, 
            Moisture = moisture, 
            Temperature = temperature 
        };
    }
}

public class TerrainData
{
    public double[,] Heights { get; set; }
    public double[,] Moisture { get; set; }
    public double[,] Temperature { get; set; }
}
```

### Procedural Texture Generation
```csharp
public class TextureGenerator
{
    public static byte[,] GenerateMarbleTexture(int width, int height, int seed)
    {
        var noise1 = new PerlinNoise(seed);
        var noise2 = new SimplexNoise(seed + 1);
        var texture = new byte[width, height];
        
        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                double nx = (double)x / 64.0;
                double ny = (double)y / 64.0;
                
                // Create marble-like patterns
                double turbulence = noise2.GetNoise(x * 0.02, y * 0.02) * 10.0;
                double marble = Math.Sin((nx + turbulence) * Math.PI) * 0.5 + 0.5;
                
                // Add fine detail
                double detail = noise1.Noise(x * 0.05, y * 0.05, 0) * 0.1;
                marble = Math.Max(0, Math.Min(1, marble + detail));
                
                texture[x, y] = (byte)(marble * 255);
            }
        }
        
        return texture;
    }
    
    public static byte[,,] GenerateWoodTexture(int width, int height, int depth, int seed)
    {
        var noise = new OpenSimplexNoise(seed);
        var texture = new byte[width, height, depth];
        
        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                for (int z = 0; z < depth; z++)
                {
                    // Distance from center for ring pattern
                    double dx = (double)x - width / 2.0;
                    double dy = (double)y - height / 2.0;
                    double distance = Math.Sqrt(dx * dx + dy * dy);
                    
                    // Add noise to create natural variation
                    double disturbance = noise.Evaluate(x * 0.1, y * 0.1, z * 0.1) * 5.0;
                    
                    // Create wood rings
                    double rings = Math.Sin((distance + disturbance) * 0.3) * 0.5 + 0.5;
                    
                    texture[x, y, z] = (byte)(rings * 255);
                }
            }
        }
        
        return texture;
    }
}
```

### Animation and Time-Based Effects
```csharp
public class AnimatedNoise
{
    private readonly OpenSimplexNoise _noise4D;
    private readonly SimplexNoise _noise2D;
    
    public AnimatedNoise(int seed)
    {
        _noise4D = new OpenSimplexNoise(seed);
        _noise2D = new SimplexNoise(seed);
    }
    
    // Animated cloud movement
    public double GetAnimatedClouds(double x, double y, double time)
    {
        // Use 4D noise for smooth animation
        return _noise4D.Evaluate(x * 0.01, y * 0.01, time * 0.1, 0);
    }
    
    // Flowing water effect
    public double GetFlowingWater(double x, double y, double time)
    {
        double baseFlow = _noise2D.GetNoise(x * 0.02 + time * 0.5, y * 0.02);
        double ripples = _noise4D.Evaluate(x * 0.1, y * 0.1, time * 2.0, 0) * 0.3;
        
        return baseFlow + ripples;
    }
    
    // Fire/flame simulation
    public double GetFlameNoise(double x, double y, double time, double height)
    {
        // More turbulence at the base, smoother at the top
        double turbulence = Math.Max(0, 1.0 - height * 0.8);
        double flameNoise = _noise4D.Evaluate(
            x * 0.05, 
            y * 0.05 - time * 2.0,  // Upward movement
            time * 1.5, 
            0
        );
        
        return flameNoise * turbulence;
    }
}

// Usage in animation loop
var animatedNoise = new AnimatedNoise(123);
double currentTime = 0.0;

// In your game/animation loop
public void UpdateFrame(double deltaTime)
{
    currentTime += deltaTime;
    
    for (int x = 0; x < screenWidth; x++)
    {
        for (int y = 0; y < screenHeight; y++)
        {
            // Animated effects
            double clouds = animatedNoise.GetAnimatedClouds(x, y, currentTime);
            double water = animatedNoise.GetFlowingWater(x, y, currentTime);
            double flames = animatedNoise.GetFlameNoise(x, y, currentTime, (double)y / screenHeight);
            
            // Use for visual effects, particle systems, etc.
        }
    }
}
```

## Performance Optimization

### Lookup Table Optimization
```csharp
public class OptimizedNoiseGenerator
{
    private readonly OpenSimplexNoise _noise;
    private readonly double[,] _cachedNoise;
    private readonly int _cacheSize;
    private readonly double _scale;
    
    public OptimizedNoiseGenerator(int seed, int cacheSize = 256, double scale = 1.0)
    {
        _noise = new OpenSimplexNoise(seed);
        _cacheSize = cacheSize;
        _scale = scale;
        _cachedNoise = new double[cacheSize, cacheSize];
        
        // Pre-compute noise values
        GenerateCache();
    }
    
    private void GenerateCache()
    {
        for (int x = 0; x < _cacheSize; x++)
        {
            for (int y = 0; y < _cacheSize; y++)
            {
                _cachedNoise[x, y] = _noise.Evaluate(x * _scale, y * _scale);
            }
        }
    }
    
    public double GetCachedNoise(int x, int y)
    {
        // Use modulo for tiling
        int cx = Math.Abs(x) % _cacheSize;
        int cy = Math.Abs(y) % _cacheSize;
        return _cachedNoise[cx, cy];
    }
    
    // Bilinear interpolation for smooth scaling
    public double GetInterpolatedNoise(double x, double y)
    {
        double fx = x * _cacheSize;
        double fy = y * _cacheSize;
        
        int x1 = (int)Math.Floor(fx) % _cacheSize;
        int y1 = (int)Math.Floor(fy) % _cacheSize;
        int x2 = (x1 + 1) % _cacheSize;
        int y2 = (y1 + 1) % _cacheSize;
        
        double wx = fx - Math.Floor(fx);
        double wy = fy - Math.Floor(fy);
        
        double n1 = _cachedNoise[x1, y1];
        double n2 = _cachedNoise[x2, y1];
        double n3 = _cachedNoise[x1, y2];
        double n4 = _cachedNoise[x2, y2];
        
        double i1 = n1 * (1 - wx) + n2 * wx;
        double i2 = n3 * (1 - wx) + n4 * wx;
        
        return i1 * (1 - wy) + i2 * wy;
    }
}
```

### Parallel Processing
```csharp
public static class ParallelNoiseGeneration
{
    public static double[,] GenerateParallelTerrain(int width, int height, int seed, 
        int octaves = 6, double persistence = 0.5)
    {
        var terrain = new double[width, height];
        var noise = new PerlinNoise(seed);
        
        // Process in parallel chunks
        Parallel.For(0, height, y =>
        {
            for (int x = 0; x < width; x++)
            {
                double value = 0;
                double amplitude = 1.0;
                double frequency = 0.01;
                double maxValue = 0;
                
                for (int octave = 0; octave < octaves; octave++)
                {
                    value += noise.Noise(x * frequency, y * frequency, 0) * amplitude;
                    maxValue += amplitude;
                    amplitude *= persistence;
                    frequency *= 2.0;
                }
                
                terrain[x, y] = value / maxValue;
            }
        });
        
        return terrain;
    }
    
    // Batch processing for large datasets
    public static void ProcessLargeDataset(double[,] output, OpenSimplexNoise noise)
    {
        int width = output.GetLength(0);
        int height = output.GetLength(1);
        
        // Process in tiles for better cache locality
        const int tileSize = 64;
        
        Parallel.For(0, (width + tileSize - 1) / tileSize, tileX =>
        {
            Parallel.For(0, (height + tileSize - 1) / tileSize, tileY =>
            {
                int startX = tileX * tileSize;
                int startY = tileY * tileSize;
                int endX = Math.Min(startX + tileSize, width);
                int endY = Math.Min(startY + tileSize, height);
                
                for (int x = startX; x < endX; x++)
                {
                    for (int y = startY; y < endY; y++)
                    {
                        output[x, y] = noise.Evaluate(x * 0.01, y * 0.01);
                    }
                }
            });
        });
    }
}
```

## Integration Examples

### With geoWrangler for Shape Generation
```csharp
using Noise;
using geoWrangler;
using Clipper2Lib;

public static class NoiseBasedShapeGeneration
{
    public static PathD GenerateOrganicShape(int seed, int vertices = 32, double baseRadius = 50.0)
    {
        var noise = new PerlinNoise(seed);
        var shape = new PathD();
        
        for (int i = 0; i < vertices; i++)
        {
            double angle = (double)i / vertices * 2.0 * Math.PI;
            
            // Add noise to radius for organic variation
            double radiusVariation = noise.Noise(
                Math.Cos(angle) * 2.0, 
                Math.Sin(angle) * 2.0, 
                0
            ) * 15.0; // Max variation of 15 units
            
            double radius = baseRadius + radiusVariation;
            double x = Math.Cos(angle) * radius;
            double y = Math.Sin(angle) * radius;
            
            shape.Add(new PointD(x, y));
        }
        
        return shape;
    }
    
    public static PathsD GenerateNoiseBasedTerrain(int width, int segments, int seed)
    {
        var noise = new SimplexNoise(seed);
        var terrain = new PathsD();
        
        // Generate terrain profile
        var profile = new PathD();
        for (int i = 0; i <= segments; i++)
        {
            double x = (double)i / segments * width;
            double height = noise.GetNoise(x * 0.01, 0) * 100.0;
            profile.Add(new PointD(x, height));
        }
        
        terrain.Add(profile);
        return terrain;
    }
}
```

### With shapeEngine for Procedural Patterns
```csharp
using Noise;
using shapeEngine;

public static class ProceduralPatterns
{
    public static void ApplyNoiseToShape(ShapeLibrary shapeLib, int seed)
    {
        var noise = new OpenSimplexNoise(seed);
        
        // Generate base shape
        shapeLib.setShape((int)ShapeLibrary.shapeNames_all.rect);
        
        if (shapeLib.shapeValid)
        {
            // Apply noise-based modifications to vertices
            for (int i = 0; i < shapeLib.Vertex.Length; i++)
            {
                double noiseX = noise.Evaluate(shapeLib.Vertex[i].X * 0.1, shapeLib.Vertex[i].Y * 0.1);
                double noiseY = noise.Evaluate(shapeLib.Vertex[i].X * 0.1, shapeLib.Vertex[i].Y * 0.1, 1000);
                
                // Apply subtle distortion
                shapeLib.Vertex[i].X += noiseX * 5.0;
                shapeLib.Vertex[i].Y += noiseY * 5.0;
            }
        }
    }
}
```

## Best Practices

### Choosing the Right Noise Type
```csharp
// Perlin Noise - Good for:
// - Classic terrain generation
// - Traditional procedural textures  
// - When 3D noise is needed
var perlin = new PerlinNoise(seed);

// Simplex Noise - Good for:
// - Higher quality 2D patterns
// - Better visual characteristics than Perlin
// - Faster computation than Perlin in 2D
var simplex = new SimplexNoise(seed);

// OpenSimplex Noise - Good for:
// - Highest quality noise
// - 4D noise for seamless animation
// - Professional applications requiring the best quality
var openSimplex = new OpenSimplexNoise(seed);
```

### Parameter Guidelines
```csharp
public static class NoiseParameters
{
    // Terrain generation guidelines
    public const double TerrainScale = 0.01;      // Controls feature size
    public const int TerrainOctaves = 6;          // Detail levels
    public const double TerrainPersistence = 0.5; // Roughness control
    
    // Texture generation
    public const double TextureScale = 0.05;      // Higher frequency for details
    public const double TexturePersistence = 0.6; // More persistent detail
    
    // Animation
    public const double TimeScale = 0.1;          // Controls animation speed
    
    // Quality vs Performance trade-offs
    public static NoiseConfig GetConfig(QualityLevel quality)
    {
        return quality switch
        {
            QualityLevel.Low => new NoiseConfig { Octaves = 3, Scale = 0.02 },
            QualityLevel.Medium => new NoiseConfig { Octaves = 6, Scale = 0.01 },
            QualityLevel.High => new NoiseConfig { Octaves = 8, Scale = 0.005 },
            QualityLevel.Ultra => new NoiseConfig { Octaves = 12, Scale = 0.002 },
            _ => new NoiseConfig { Octaves = 6, Scale = 0.01 }
        };
    }
}

public enum QualityLevel { Low, Medium, High, Ultra }

public struct NoiseConfig
{
    public int Octaves;
    public double Scale;
}
```

### Error Handling and Validation
```csharp
public static class SafeNoiseOperations
{
    public static double SafeGetNoise(SimplexNoise noise, double x, double y)
    {
        try
        {
            if (double.IsNaN(x) || double.IsNaN(y) || 
                double.IsInfinity(x) || double.IsInfinity(y))
            {
                return 0.0; // Return neutral value for invalid input
            }
            
            return noise.GetNoise(x, y);
        }
        catch (Exception ex)
        {
            Console.WriteLine($"Noise generation failed: {ex.Message}");
            return 0.0;
        }
    }
    
    public static bool ValidateNoiseOutput(double value, double minExpected = -1.0, double maxExpected = 1.0)
    {
        return !double.IsNaN(value) && 
               !double.IsInfinity(value) && 
               value >= minExpected && 
               value <= maxExpected;
    }
}
```

## Dependencies
- **entropyRNG** - Used by SimplexNoise for random number generation
- **.NET 8.0** - Target framework
- **System.Threading.Tasks** - Parallel processing support

## Related Libraries
- **[MersenneTwister](MersenneTwister-API.md)** - High-quality random number generation used as basis for noise
- **[entropyRNG](entropyRNG-API.md)** - Thread-safe RNG system used by SimplexNoise
- **[geoWrangler](geoWrangler-API.md)** - Can use noise for organic shape generation
- **[shapeEngine](shapeEngine-API.md)** - Can incorporate noise for procedural shape variation