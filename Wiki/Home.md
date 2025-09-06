# DesignLibs_GPL API Documentation

Welcome to the comprehensive API documentation for DesignLibs_GPL, a collection of .NET 8.0 C# libraries for geometric design, mathematical operations, and graphics applications.

## Overview

DesignLibs_GPL comprises a multitude of libraries primarily built for use by the [Variance](https://github.com/philstopford/Variance_GPL) and [Quilt](https://github.com/philstopford/Quilt_GPL) tools. These libraries provide essential functionality for:

- **Geometric Operations**: Basic primitives, transformations, analysis, and layout file format support
- **Mathematical Functions**: Random number generation, noise functions, and mathematical parsing
- **Graphics and UI**: Cross-platform viewport libraries and error handling
- **Utility Functions**: Collections, data structures, and general-purpose utilities

## Library Categories

### Core Geometric Libraries
- **[geoLib](geoLib-API.md)** - Basic geometric primitives (points, arrays, matrices, rectangles)
- **[geoCore](geoCore-API.md)** - GDSII and Oasis file format parsing and writing
- **[geoWrangler](geoWrangler-API.md)** - Geometry transformations, raycasting, fragmentation, and decimation
- **[shapeEngine](shapeEngine-API.md)** - Shape generation and manipulation (rectangles, L-shapes, T-shapes, etc.)

### Polygon and Clipping Operations
- **[clipper](clipper-API.md)** - Advanced polygon clipping operations (Clipper2 C# port)
- **[LibTessDotNet](LibTessDotNet-API.md)** - Triangle tessellation for complex polygons

### Mathematical and Random Number Libraries
- **[MersenneTwister](MersenneTwister-API.md)** - High-quality pseudorandom number generation
- **[entropyRNG](entropyRNG-API.md)** - Thread-safe RNG system with entropy collection
- **[Noise](Noise-API.md)** - Perlin, Simplex, and OpenSimplex noise functions
- **[info.lundin.math](info.lundin.math-API.md)** - Mathematical expression parser and evaluator

### Graphics and UI Libraries
- **[Eto.VeldridSurface](Eto.VeldridSurface-API.md)** - Cross-platform 2D viewport for end-user applications
- **[errorReporter](errorReporter-API.md)** - Cross-platform error handling and reporting
- **[Color](Color-API.md)** - Color representations for viewports and graphics

### Utility and Data Structure Libraries
- **[utility](utility-API.md)** - General purpose utilities (compression, hashing, mathematical operations)
- **[MiscUtil](MiscUtil-API.md)** - Collections, ranges, and utility functions
- **[KDTree](KDTree-API.md)** - Spatial data structure for efficient nearest-neighbor searches
- **[SVGBuilder](SVGBuilder-API.md)** - SVG creation and manipulation utilities

### Additional Libraries
- **[Email](Email-API.md)** - Email functionality and utilities
- **[LHC](LHC-API.md)** - Specialized computational library  
- **[Tral.Randomness](Tral.Randomness-API.md)** - Additional randomness and statistical utilities

## Completed API Documentation

This comprehensive wiki provides detailed API documentation for the core DesignLibs_GPL libraries:

### ✅ Fully Documented Libraries
- **[geoLib](geoLib-API.md)** - Basic geometric primitives with cross-platform compatibility
- **[geoCore](geoCore-API.md)** - GDSII and Oasis file format parsing and writing  
- **[geoWrangler](geoWrangler-API.md)** - Advanced geometry transformations and operations
- **[shapeEngine](shapeEngine-API.md)** - Parametric shape generation (rectangles, L-shapes, T-shapes, etc.)
- **[utility](utility-API.md)** - Mathematical utilities, hashing, compression, and performance optimization
- **[clipper](clipper-API.md)** - Advanced polygon clipping operations (Clipper2 C# port)
- **[MersenneTwister](MersenneTwister-API.md)** - High-quality pseudorandom number generation
- **[Noise](Noise-API.md)** - Perlin, Simplex, and OpenSimplex noise functions
- **[Eto.VeldridSurface](Eto.VeldridSurface-API.md)** - Cross-platform 2D graphics viewport

Each documentation includes:
- 📋 Complete API reference with all classes, methods, and properties
- 💡 Comprehensive usage examples and code samples  
- 🔗 Integration examples with other DesignLibs libraries
- ⚡ Performance optimization guidelines
- 🛡️ Error handling and best practices
- 📚 Cross-references to related libraries

## Getting Started

### Prerequisites
- .NET 8.0 SDK or later
- Compatible operating system (Windows, Linux, macOS)

### Building the Libraries
```bash
# Clone the repository
git clone https://github.com/philstopford/DesignLibs_GPL.git
cd DesignLibs_GPL

# Restore dependencies
dotnet restore

# Build core libraries
dotnet build UnitTests/UnitTests.csproj

# Run tests to validate functionality
dotnet test UnitTests/UnitTests.csproj

# Run example application
dotnet run --project HistoTest/HistoTest.csproj
```

### Basic Usage Example

```csharp
using geoLib;
using geoWrangler;
using Clipper2Lib;

// Create basic geometric shapes
var rectangle = GeoWrangler.rectangle(100, 50, 10, 10); // width, height, x_offset, y_offset
var circle = GeoWrangler.ellipse(25, 25, 50, 50, 32); // rx, ry, cx, cy, resolution

// Perform boolean operations
var union = Clipper.Union(new PathsD { rectangle, circle }, FillRule.NonZero);

// Transform geometry
var rotated = GeoWrangler.rotate(union, 45.0); // rotate by 45 degrees
var translated = GeoWrangler.move(rotated, 100, 100); // translate by (100, 100)
```

## Key Concepts

### Coordinate Systems
- Most libraries use integer coordinates for precision (Point64, PathD with scaling)
- Floating-point operations are provided where necessary (PointD)
- GDSII/Oasis files use database units for precise manufacturing requirements

### Threading Support
- Many libraries are designed for multi-threaded environments
- Random number generators provide thread-safe implementations
- Geometric operations can be parallelized where appropriate

### Performance Considerations
- Integer-based geometry avoids floating-point precision issues
- Optimized algorithms for polygon operations and spatial queries
- Memory-efficient data structures for large datasets

## Contributing

This project follows semantic versioning and welcomes contributions. Please ensure:
- All changes include appropriate unit tests
- Documentation is updated for new features
- Code follows the established patterns and style

## License

DesignLibs_GPL is released under the GPLv3 license. See the LICENSE file for full details.