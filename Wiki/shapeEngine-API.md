# shapeEngine API Documentation

## Overview

The **shapeEngine** library provides comprehensive shape generation and manipulation capabilities for creating various geometric forms used in layout design and analysis. It supports both parametric shape generation and loading custom geometry from external sources.

## Key Features

- **Parametric Shape Generation** - Rectangles, L-shapes, T-shapes, X-shapes, U-shapes, S-shapes
- **Custom Shape Import** - Load shapes from GDSII/Oasis files or boolean operations
- **Shape Analysis** - Orthogonal geometry detection, bounding box calculation
- **Flexible Configuration** - Configurable shape parameters and positioning
- **Performance Optimization** - Caching and parallel processing support

## Namespace
```csharp
using shapeEngine;
```

## Core Classes

### ShapeLibrary
The main class for shape generation and management.

#### Supported Shape Types
```csharp
public enum shapeNames_all
{
    none,        // No shape selected
    rect,        // Rectangle or square
    Lshape,      // L-shaped polygon  
    Tshape,      // T-shaped polygon
    Xshape,      // X-shaped polygon (cross)
    Ushape,      // U-shaped polygon
    Sshape,      // S-shaped polygon
    GEOCORE,     // Shape loaded from GDS/Oasis file
    BOOLEAN,     // Shape from boolean operations
    text,        // Text-based shape
    bounding,    // Bounding box shape
    complex      // Complex layout shape
}
```

#### Shape Names List
```csharp
private static readonly List<string> availableShapes_all =
[
    "(None)", "Rectangle/Square", "L-shape", "T-shape", "X-shape", 
    "U-shape", "S-shape", "GDS/Oasis", "Boolean", "Text", "Bounding", "Layout"
];
```

#### Properties
```csharp
public int shapeIndex { get; set; }                    // Currently selected shape
public bool shapeValid { get; private set; }          // Shape validity flag
public bool geoCoreShapeOrthogonal { get; private set; } // Orthogonal geometry flag
public MyVertex[] Vertex { get; private set; }        // Shape vertices
public MyRound[] round1 { get; private set; }         // Corner rounding info
public bool[] tips { get; private set; }              // Edge tips configuration
public ShapeSettings LayerSettings { get; }           // Shape settings
public string last_error { get; set; }                // Last error message
```

#### Constructors
```csharp
// Basic constructor
ShapeLibrary(int[] shapes, ShapeSettings shapeSettings)

// Constructor with initial shape selection
ShapeLibrary(int[] shapes, int shapeIndex, ShapeSettings shapeSettings)
```

#### Usage Example
```csharp
// Define available shapes for client
int[] clientShapes = { 
    (int)ShapeLibrary.shapeNames_all.rect,
    (int)ShapeLibrary.shapeNames_all.Lshape,
    (int)ShapeLibrary.shapeNames_all.Tshape 
};

// Create shape settings
var settings = new ShapeSettings();
settings.setDouble(ShapeSettings.dbl.width, 100.0);
settings.setDouble(ShapeSettings.dbl.height, 50.0);

// Create shape library
var shapeLib = new ShapeLibrary(clientShapes, settings);

// Generate rectangle
shapeLib.setShape((int)ShapeLibrary.shapeNames_all.rect);
if (shapeLib.shapeValid)
{
    var vertices = shapeLib.Vertex;
    Console.WriteLine($"Generated {vertices.Length} vertices");
}
```

### ShapeSettings
Configuration class for shape parameters and positioning.

#### Parameter Types
```csharp
public enum dbl  // Double parameters
{
    width, height, cornerRadius, // Basic dimensions
    tipDepth, tipWidth,          // Tip parameters  
    offsetX, offsetY,            // Positioning
    rotation, scale              // Transformations
}

public enum intA // Integer parameters  
{
    tipLocations, subShapeCount, // Configuration
    layerIndex, datatype,        // Layer assignment
    arrayX, arrayY               // Array parameters
}

public enum bool_ // Boolean parameters
{
    enabled, orthogonal,         // Basic flags
    flipX, flipY,               // Mirroring
    roundCorners                // Rounding
}
```

#### Tip Locations
```csharp
public enum tipLocations 
{ 
    none, L, R, LR,              // Left, Right combinations
    T, B, TB,                    // Top, Bottom combinations  
    TL, TR, TLR,                // Top + side combinations
    BL, BR, BLR,                // Bottom + side combinations
    TBL, TBR, all               // Multi-side combinations
}
```

#### Sub-shape Positioning
```csharp
public enum subShapeLocations 
{ 
    TL, TR, BL, BR,             // Corner positions
    TS, RS, BS, LS,             // Side midpoints
    C                           // Center
}

public enum subShapeHorLocs { L, M, R }      // Left, Middle, Right
public enum subShapeVerLocs { B, M, T }      // Bottom, Middle, Top
```

#### Methods
```csharp
// Parameter setters
public void setDouble(dbl parameter, double value)
public void setInt(intA parameter, int value)  
public void setBool(bool_ parameter, bool value)

// Parameter getters
public double getDouble(dbl parameter)
public int getInt(intA parameter)
public bool getBool(bool_ parameter)

// Configuration helpers
public static List<string> getAvailableTipsLocations()
public static List<string> getAvailableSubShapePositions()
public static List<string> getPolyFillTypes()
```

#### Usage Example
```csharp
var settings = new ShapeSettings();

// Set basic dimensions
settings.setDouble(ShapeSettings.dbl.width, 200.0);
settings.setDouble(ShapeSettings.dbl.height, 100.0);
settings.setDouble(ShapeSettings.dbl.cornerRadius, 5.0);

// Configure positioning
settings.setDouble(ShapeSettings.dbl.offsetX, 50.0);
settings.setDouble(ShapeSettings.dbl.offsetY, 25.0);

// Set boolean options
settings.setBool(ShapeSettings.bool_.roundCorners, true);
settings.setBool(ShapeSettings.bool_.orthogonal, true);

// Configure tips
settings.setInt(ShapeSettings.intA.tipLocations, (int)ShapeSettings.tipLocations.LR);
```

## Shape Generation

### Basic Shapes

#### Rectangle
```csharp
// Generate rectangle
var settings = new ShapeSettings();
settings.setDouble(ShapeSettings.dbl.width, 100.0);
settings.setDouble(ShapeSettings.dbl.height, 50.0);

var shapeLib = new ShapeLibrary(clientShapes, settings);
shapeLib.setShape((int)ShapeLibrary.shapeNames_all.rect);

// Access vertices
foreach (var vertex in shapeLib.Vertex)
{
    Console.WriteLine($"Vertex: ({vertex.X:F2}, {vertex.Y:F2})");
}
```

#### L-Shape
```csharp
var settings = new ShapeSettings();
settings.setDouble(ShapeSettings.dbl.width, 100.0);      // Total width
settings.setDouble(ShapeSettings.dbl.height, 80.0);      // Total height
settings.setDouble(ShapeSettings.dbl.cornerWidth, 30.0); // Corner section width
settings.setDouble(ShapeSettings.dbl.cornerHeight, 25.0);// Corner section height

shapeLib.setShape((int)ShapeLibrary.shapeNames_all.Lshape);
```

#### T-Shape
```csharp
var settings = new ShapeSettings();
settings.setDouble(ShapeSettings.dbl.width, 120.0);      // Top bar width
settings.setDouble(ShapeSettings.dbl.height, 80.0);      // Total height
settings.setDouble(ShapeSettings.dbl.stemWidth, 40.0);   // Vertical stem width

shapeLib.setShape((int)ShapeLibrary.shapeNames_all.Tshape);
```

#### X-Shape (Cross)
```csharp
var settings = new ShapeSettings();
settings.setDouble(ShapeSettings.dbl.width, 100.0);      // Horizontal arm width
settings.setDouble(ShapeSettings.dbl.height, 100.0);     // Vertical arm height
settings.setDouble(ShapeSettings.dbl.armWidth, 20.0);    // Arm thickness

shapeLib.setShape((int)ShapeLibrary.shapeNames_all.Xshape);
```

#### U-Shape
```csharp
var settings = new ShapeSettings();
settings.setDouble(ShapeSettings.dbl.width, 100.0);      // Total width
settings.setDouble(ShapeSettings.dbl.height, 80.0);      // Total height
settings.setDouble(ShapeSettings.dbl.wallThickness, 15.0); // Wall thickness

shapeLib.setShape((int)ShapeLibrary.shapeNames_all.Ushape);
```

#### S-Shape
```csharp
var settings = new ShapeSettings();
settings.setDouble(ShapeSettings.dbl.width, 100.0);      // Total width
settings.setDouble(ShapeSettings.dbl.height, 120.0);     // Total height
settings.setDouble(ShapeSettings.dbl.segmentHeight, 30.0); // Segment height

shapeLib.setShape((int)ShapeLibrary.shapeNames_all.Sshape);
```

### Custom Shapes

#### From GeoCore (GDSII/Oasis)
```csharp
// Load custom shape from external geometry
PathD customGeometry = /* loaded from geoCore */;

shapeLib.setShape((int)ShapeLibrary.shapeNames_all.GEOCORE, customGeometry);

if (shapeLib.shapeValid)
{
    bool isOrthogonal = shapeLib.geoCoreShapeOrthogonal;
    Console.WriteLine($"Custom shape loaded, orthogonal: {isOrthogonal}");
}
```

#### From Boolean Operations
```csharp
// Create shape from boolean result
PathD booleanResult = /* result from boolean operation */;

shapeLib.setShape((int)ShapeLibrary.shapeNames_all.BOOLEAN, booleanResult);
```

## Advanced Features

### Corner Rounding
```csharp
// Configure corner rounding
var settings = new ShapeSettings();
settings.setBool(ShapeSettings.bool_.roundCorners, true);
settings.setDouble(ShapeSettings.dbl.cornerRadius, 10.0);

// Generate shape with rounded corners
shapeLib.setShape(shapeIndex);

// Access rounding information
foreach (var round in shapeLib.round1)
{
    Console.WriteLine($"Corner {round.index}: Max radius = {round.MaxRadius}");
}
```

### Edge Tips
```csharp
// Configure edge tips
var settings = new ShapeSettings();
settings.setInt(ShapeSettings.intA.tipLocations, (int)ShapeSettings.tipLocations.TLR);
settings.setDouble(ShapeSettings.dbl.tipDepth, 5.0);
settings.setDouble(ShapeSettings.dbl.tipWidth, 3.0);

shapeLib.setShape(shapeIndex);

// Check which edges have tips
for (int i = 0; i < shapeLib.tips.Length; i++)
{
    if (shapeLib.tips[i])
    {
        Console.WriteLine($"Edge {i} has tip");
    }
}
```

### Shape Transformations
```csharp
// Apply transformations
var settings = new ShapeSettings();
settings.setDouble(ShapeSettings.dbl.offsetX, 100.0);    // Translate X
settings.setDouble(ShapeSettings.dbl.offsetY, 50.0);     // Translate Y
settings.setDouble(ShapeSettings.dbl.rotation, 45.0);    // Rotate 45 degrees
settings.setDouble(ShapeSettings.dbl.scale, 1.5);        // Scale 150%
settings.setBool(ShapeSettings.bool_.flipX, true);       // Mirror X
```

## Shape Analysis

### Geometric Properties
```csharp
// Get shape pivot point (centroid)
PointD pivot = shapeLib.getPivotPoint();
Console.WriteLine($"Shape center: ({pivot.x:F2}, {pivot.y:F2})");

// Check orthogonal geometry
if (shapeLib.geoCoreShapeOrthogonal)
{
    Console.WriteLine("Shape contains only 90-degree angles");
}

// Validate shape
if (shapeLib.shapeValid)
{
    Console.WriteLine($"Shape has {shapeLib.Vertex.Length} vertices");
}
else
{
    Console.WriteLine($"Shape invalid: {shapeLib.last_error}");
}
```

### Vertex Analysis
```csharp
// Analyze vertices
foreach (var vertex in shapeLib.Vertex)
{
    Console.WriteLine($"Vertex ({vertex.X:F2}, {vertex.Y:F2}):");
    Console.WriteLine($"  Direction: {vertex.direction}");
    Console.WriteLine($"  Type: {vertex.type}");
    Console.WriteLine($"  Vertical: {vertex.vertical}");
    Console.WriteLine($"  Inner: {vertex.inner}");
}
```

### Bounding Analysis
```csharp
// Calculate bounding box
var vertices = shapeLib.Vertex;
var points = vertices.Select(v => new PointD(v.X, v.Y)).ToList();

double minX = points.Min(p => p.x);
double minY = points.Min(p => p.y);  
double maxX = points.Max(p => p.x);
double maxY = points.Max(p => p.y);

Console.WriteLine($"Bounds: ({minX:F2}, {minY:F2}) to ({maxX:F2}, {maxY:F2})");
Console.WriteLine($"Size: {maxX - minX:F2} x {maxY - minY:F2}");
```

## Integration with Other Libraries

### With geoWrangler
```csharp
using geoWrangler;

// Generate shape
shapeLib.setShape(shapeIndex);

// Convert to PathD for geoWrangler operations
var pathD = new PathD();
for (int i = 0; i < shapeLib.Vertex.Length - 1; i++) // Skip last vertex (closure)
{
    pathD.Add(new PointD(shapeLib.Vertex[i].X, shapeLib.Vertex[i].Y));
}

// Apply geoWrangler operations
var rotated = GeoWrangler.rotate(pathD, 30.0);
var scaled = GeoWrangler.scale(rotated, 2.0);
var translated = GeoWrangler.move(scaled, 100.0, 50.0);
```

### With geoCore
```csharp
using geoCoreLib;

// Generate shape
shapeLib.setShape(shapeIndex);

// Convert to Path64 for geoCore
var path64 = new Path64();
foreach (var vertex in shapeLib.Vertex)
{
    path64.Add(new Point64((long)(vertex.X * 1000), (long)(vertex.Y * 1000)));
}

// Add to geoCore cell
var cell = new GCCell();
cell.addPolygon(path64, layerNumber: 1, datatype: 0);
```

### With Clipper2Lib
```csharp
using Clipper2Lib;

// Generate multiple shapes
var shapes = new PathsD();
for (int i = 0; i < 3; i++)
{
    shapeLib.setShape(shapeIndex);
    var pathD = /* convert from shapeLib.Vertex */;
    shapes.Add(GeoWrangler.move(pathD, i * 150.0, 0.0));
}

// Apply boolean operations
var union = Clipper.Union(shapes, FillRule.NonZero);
var simplified = Clipper.SimplifyPaths(union, FillRule.NonZero);
```

## Performance Considerations

### Caching
```csharp
// Shape library caches computed geometry
// Repeated calls with same parameters are optimized

// Force recalculation if needed
shapeLib.setShape(shapeIndex); // Regenerates if parameters changed
```

### Threading
```csharp
// Shape generation supports parallel processing
// Vertex calculations can be parallelized internally

// Control threading via preprocessor directives
#define SHAPELIBSINGLETHREADED  // Disable threading for debugging
```

### Memory Management
```csharp
// Clean up after large operations
shapeLib.Vertex = null;
shapeLib.round1 = null;
shapeLib.tips = null;
GC.Collect(); // Force cleanup if needed
```

## Error Handling

### Validation
```csharp
// Always check shape validity
shapeLib.setShape(shapeIndex);
if (!shapeLib.shapeValid)
{
    Console.WriteLine($"Shape generation failed: {shapeLib.last_error}");
    return;
}

// Validate settings before shape generation
var settings = shapeLib.LayerSettings;
if (settings.getDouble(ShapeSettings.dbl.width) <= 0)
{
    Console.WriteLine("Invalid width specified");
    return;
}
```

### Common Errors
```csharp
// Handle common error scenarios
try
{
    shapeLib.setShape(shapeIndex, customGeometry);
}
catch (Exception ex)
{
    if (ex.Message.Contains("More shapes requested"))
    {
        Console.WriteLine("Shape index out of range");
    }
    else if (ex.Message.Contains("Invalid geometry"))
    {
        Console.WriteLine("Custom geometry is invalid");
    }
    else
    {
        Console.WriteLine($"Unexpected error: {ex.Message}");
    }
}
```

## Best Practices

### Shape Configuration
```csharp
// Use appropriate parameter types
settings.setDouble(ShapeSettings.dbl.width, 100.0);      // Use double for dimensions
settings.setInt(ShapeSettings.intA.layerIndex, 1);       // Use int for indices
settings.setBool(ShapeSettings.bool_.enabled, true);     // Use bool for flags

// Validate parameters before use
double width = settings.getDouble(ShapeSettings.dbl.width);
if (width <= 0)
{
    throw new ArgumentException("Width must be positive");
}
```

### Performance Optimization
```csharp
// Reuse ShapeLibrary instances
var shapeLib = new ShapeLibrary(clientShapes, settings);

// Generate multiple similar shapes efficiently
for (int i = 0; i < count; i++)
{
    settings.setDouble(ShapeSettings.dbl.offsetX, i * spacing);
    shapeLib.setShape(shapeIndex); // Reuses library instance
    
    // Process shape...
}
```

### Memory Efficiency
```csharp
// Process shapes in batches for large datasets
const int batchSize = 100;
for (int batch = 0; batch < totalShapes; batch += batchSize)
{
    var currentBatch = /* load batch */;
    
    foreach (var shapeConfig in currentBatch)
    {
        // Process shape
    }
    
    // Optional cleanup between batches
    if (batch % (batchSize * 10) == 0)
    {
        GC.Collect();
    }
}
```

## Dependencies
- **geoLib** - Basic geometric primitives and data types
- **geoWrangler** - Geometric operations and transformations
- **Clipper2Lib** - Polygon data types and operations
- **utility** - Mathematical and utility functions
- **.NET 8.0** - Target framework

## Related Libraries
- **[geoWrangler](geoWrangler-API.md)** - Geometric transformations applied to generated shapes
- **[geoLib](geoLib-API.md)** - Basic geometric primitives used in shape generation
- **[geoCore](geoCore-API.md)** - File I/O for loading custom shapes from GDSII/Oasis
- **[clipper](clipper-API.md)** - Boolean operations on generated shapes