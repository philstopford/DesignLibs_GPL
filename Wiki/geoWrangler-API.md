# geoWrangler API Documentation

## Overview

The **geoWrangler** library provides a comprehensive suite of geometric transformations, operations, and analysis tools. It serves as the primary geometric processing engine, offering everything from basic transformations to advanced operations like raycasting, fragmentation, and boolean operations.

## Key Features

- **Geometric Transformations** - Translation, rotation, scaling, mirroring
- **Boolean Operations** - Union, intersection, difference with advanced handling
- **Raycasting System** - Proximity analysis and collision detection  
- **Fragmentation & Decimation** - Polygon simplification and subdivision
- **Spatial Analysis** - Distance measurements, bounds calculation, containment tests
- **Array Generation** - Create regular arrays and patterns
- **Coordinate Conversion** - Seamless conversion between coordinate systems

## Namespace
```csharp
using geoWrangler;
```

## Core Classes

### GeoWrangler (Static Class)
The main static class providing all geometric operations.

## Coordinate System Operations

### Bounds Calculation
Calculate bounding boxes for geometry.

```csharp
// Get bounds of a polygon
PathD polygon = /* your polygon */;
PathD bounds = GeoWrangler.getBounds(polygon);
double minX = bounds[0].x, minY = bounds[0].y;
double maxX = bounds[1].x, maxY = bounds[1].y;

// Integer version
Path64 bounds64 = GeoWrangler.getBounds(polygonPath64);
```

### Distance Measurements
```csharp
// Distance between two points
PointD point1 = new(10, 20);
PointD point2 = new(30, 40);
double distance = GeoWrangler.distanceBetweenPoints(point1, point2);

// Shortest distance from point to polygon edge
double distToEdge = GeoWrangler.distanceBetweenPoints(point, polygonPath);
```

## Transformation Operations

### Translation (Move)
```csharp
// Move a single polygon
PathD moved = GeoWrangler.move(polygon, deltaX: 100, deltaY: 50);

// Move multiple polygons
PathsD movedPolygons = GeoWrangler.move(polygons, deltaX: 100, deltaY: 50);

// Move with PointD offset
PointD offset = new(100, 50);
PathD movedByPoint = GeoWrangler.move(polygon, offset);
```

### Rotation
```csharp
// Rotate around origin (degrees)
PathD rotated = GeoWrangler.rotate(polygon, angle: 45.0);

// Rotate around specific point
PointD center = new(50, 50);
PathD rotatedAroundCenter = GeoWrangler.rotate(polygon, angle: 45.0, center);

// Rotate multiple polygons
PathsD rotatedPolygons = GeoWrangler.rotate(polygons, angle: 90.0);
```

### Scaling  
```csharp
// Uniform scaling
PathD scaled = GeoWrangler.scale(polygon, factor: 2.0);

// Non-uniform scaling
PathD scaledXY = GeoWrangler.scale(polygon, xFactor: 2.0, yFactor: 1.5);

// Scale around specific center
PointD scaleCenter = new(100, 100);
PathD scaledFromCenter = GeoWrangler.scale(polygon, 2.0, scaleCenter);
```

### Mirroring
```csharp
// Mirror horizontally (flip X)
PathD mirroredX = GeoWrangler.flip(polygon, flipX: true);

// Mirror vertically (flip Y)  
PathD mirroredY = GeoWrangler.flip(polygon, flipY: true);

// Mirror both axes
PathD mirroredXY = GeoWrangler.flip(polygon, flipX: true, flipY: true);
```

## Shape Generation

### Basic Shapes
```csharp
// Create rectangle
PathD rectangle = GeoWrangler.rectangle(width: 100, height: 50, xOffset: 10, yOffset: 20);

// Create circle/ellipse
PathD circle = GeoWrangler.ellipse(radiusX: 25, radiusY: 25, centerX: 50, centerY: 50, resolution: 64);

// Create ellipse
PathD ellipse = GeoWrangler.ellipse(radiusX: 40, radiusY: 20, centerX: 0, centerY: 0, resolution: 32);
```

### Complex Shapes
```csharp
// L-shape
PathD lShape = GeoWrangler.lShape(width: 100, height: 80, cornerWidth: 30, cornerHeight: 25);

// T-shape  
PathD tShape = GeoWrangler.tShape(width: 120, height: 80, stemWidth: 40);

// Regular polygon
PathD hexagon = GeoWrangler.regularPolygon(sides: 6, radius: 50, centerX: 0, centerY: 0);
```

## Array Operations

### Create Arrays
```csharp
// Create 2D array from single polygon
PathsD array = GeoWrangler.makeArray(
    source: polygon,
    xCount: 5,        // 5 columns
    xPitch: 100.0,    // 100 units spacing in X
    yCount: 3,        // 3 rows
    yPitch: 80.0      // 80 units spacing in Y
);

// Create array from multiple polygons
PathsD complexArray = GeoWrangler.makeArray(
    source: polygons,
    xCount: 4,
    xPitch: 150.0,
    yCount: 2,
    yPitch: 200.0
);
```

## Boolean Operations

### Standard Boolean Operations
```csharp
using Clipper2Lib;

// Union (OR operation)
PathsD union = Clipper.Union(polygons, FillRule.NonZero);

// Intersection (AND operation)  
PathsD intersection = Clipper.Intersect(polygons1, polygons2, FillRule.NonZero);

// Difference (subtraction)
PathsD difference = Clipper.Difference(subject: polygons1, clip: polygons2, FillRule.NonZero);
```

### Advanced Boolean Operations
```csharp
// Custom boolean with layer flags
PathsD result = GeoWrangler.customBoolean(
    firstLayerOperator: (int)GeoWrangler.LayerFlag.none,
    firstLayer: polygons1,
    secondLayerOperator: (int)GeoWrangler.LayerFlag.NOT,  
    secondLayer: polygons2,
    booleanFlag: (int)GeoWrangler.booleanOperation.AND,
    resolution: 1000.0,
    extension: 0.0
);
```

### Polygon Inflation/Deflation
```csharp
// Expand polygons (positive sizing)
PathsD expanded = GeoWrangler.sizingOperation(polygons, sizing: 10.0, resolution: 1000);

// Shrink polygons (negative sizing)
PathsD shrunk = GeoWrangler.sizingOperation(polygons, sizing: -5.0, resolution: 1000);

// Precise integer sizing
Paths64 resized = GeoWrangler.sizingOperation(polygons64, sizing: 1000); // 1000 database units
```

## Raycasting System

### Basic Raycasting
```csharp
// Create raycast analysis
PathD emissionPolygon = /* source polygon */;
PathsD targetPolygons = /* collision targets */;

var rayCast = new RayCast(
    emissionPath: emissionPolygon,
    collisionPaths: targetPolygons,
    max: 10000,  // maximum ray length
    projectCorners: true
);

// Get results
PathsD rays = rayCast.getRays();           // All cast rays
PathsD clippedRays = rayCast.getClippedRays(); // Rays clipped to collisions
double rayLength = rayCast.getRayLength(0);     // Length of specific ray
```

### Advanced Raycasting Options
```csharp
var advancedRayCast = new RayCast(
    emissionPath: emissionPolygon,
    collisionPaths: targetPolygons,
    max: 50000,
    projectCorners: false,                              // Use averaged normals
    invert: RayCast.inversionMode.x,                   // Invert X direction
    multisampleRayCount: 3,                            // Multiple rays per vertex
    runOuterLoopThreaded: true,                        // Enable threading
    startOffset: new PointD(5.0, 0.0),               // Start offset
    endOffset: new PointD(-2.0, 0.0),                // End offset  
    sideRayFallOff: RayCast.Falloff.gaussian,        // Side ray falloff
    sideRayFallOffMultiplier: 0.8,                    // Falloff strength
    dirOverride: RayCast.forceSingleDirection.vertical // Force vertical rays
);
```

### Raycast Analysis
```csharp
// Analyze raycast results
PathsD clippedRays = rayCast.getClippedRays();
for (int i = 0; i < clippedRays.Count; i++)
{
    double length = rayCast.getRayLength(i);
    Console.WriteLine($"Ray {i}: Length = {length}");
    
    if (clippedRays[i].Count >= 2)
    {
        PointD start = clippedRays[i][0];
        PointD end = clippedRays[i][^1];
        Console.WriteLine($"  From ({start.x:F2}, {start.y:F2}) to ({end.x:F2}, {end.y:F2})");
    }
}
```

## Geometric Analysis

### Area Calculations
```csharp
// Calculate polygon area
double area = GeoWrangler.area(polygon);

// Area of multiple polygons
double totalArea = 0;
foreach (var poly in polygons)
{
    totalArea += GeoWrangler.area(poly);
}
```

### Containment and Intersection Tests
```csharp
// Point-in-polygon test
bool isInside = GeoWrangler.pointInPolygon(testPoint, polygon);

// Polygon intersection test
bool intersects = GeoWrangler.polygonsIntersect(polygon1, polygon2);

// Enclosure test
bool enclosed = GeoWrangler.enclosed(innerPolygons, outerPolygons);
```

### Angle and Direction Analysis
```csharp
// Calculate angle between points
double angle = GeoWrangler.angle(center, point1, point2);

// Check if polygon is clockwise
bool isClockwise = GeoWrangler.clockwise(polygon);

// Get polygon orientation
bool clockwiseOriented = GeoWrangler.clockwise(polygon);
if (!clockwiseOriented)
{
    polygon.Reverse(); // Make clockwise
}
```

## Polygon Processing

### Simplification and Decimation
```csharp
// Simplify polygon (reduce vertex count)
PathD simplified = GeoWrangler.simplify(polygon, tolerance: 1.0);

// Decimate with specific algorithm
PathD decimated = GeoWrangler.decimate(polygon, factor: 0.5);

// Remove collinear points
PathD cleaned = GeoWrangler.removeCollinear(polygon, tolerance: 0.1);
```

### Fragmentation
```csharp
// Fragment polygons for processing
var fragmenter = new Fragmenter(resolution: 1000.0);
PathsD fragments = fragmenter.fragmentPolygons(polygons);

// Process fragments individually
foreach (var fragment in fragments)
{
    // Perform operations on each fragment
    var processed = GeoWrangler.someOperation(fragment);
}
```

### Reordering and Sanitization
```csharp
// Ensure consistent vertex ordering
PathD reordered = GeoWrangler.reorderPolygon(polygon);

// Sanitize polygon (remove invalid points, fix winding)
PathD sanitized = GeoWrangler.sanitizePolygon(polygon);

// Strip duplicate points
PathD stripped = GeoWrangler.stripDuplicates(polygon, tolerance: 0.001);
```

## Coordinate Conversion

### Type Conversions
```csharp
// Convert PathD to Path64 (with scaling)
Path64 scaledPath = GeoWrangler.path64FromPathD(pathD, scaling: 1000.0);

// Convert Path64 to PathD (with scaling)  
PathD scaledPathD = GeoWrangler.PathDFromPath64(path64, scaling: 0.001);

// Batch conversions
PathsD convertedPaths = GeoWrangler.pathsDFromPaths64(paths64, scaling: 0.001);
```

### Precision Handling
```csharp
// High-precision operations using integers
Paths64 precisePaths = /* your integer paths */;
Paths64 preciseResult = GeoWrangler.booleanOperation64(precisePaths, operation);

// Convert back to floating point for display
PathsD displayPaths = GeoWrangler.pathsDFromPaths64(preciseResult, scaling: 1e-3);
```

## Advanced Features

### Keyhole Processing  
```csharp
// Handle keyhole (concave) polygons
PathsD keyholed = GeoWrangler.processKeyholes(complexPolygons);

// Keyhole-aware boolean operations
PathsD result = GeoWrangler.keyholeBooleans(polygons1, polygons2, operation);
```

### Proximity Analysis
```csharp
// Find nearest neighbors
var proximityResults = GeoWrangler.proximityAnalysis(
    sourcePoints: points,
    targetPolygons: polygons,
    maxDistance: 100.0
);

// Distance fields
PathsD distanceField = GeoWrangler.generateDistanceField(
    sources: polygons,
    resolution: 1.0,
    maxDistance: 50.0
);
```

### Edge Extension
```csharp
// Extend polygon edges
PathsD extendedEdges = GeoWrangler.extendEdges(edges, extension: 10.0);

// Extend single edge
PathD extendedEdge = GeoWrangler.extendEdge(edge, extension: 5.0);
```

## Performance Optimization

### Threading Control
```csharp
// Control threading for operations
#define GWSINGLETHREADED  // Disable threading for debugging

// Many operations automatically parallelize:
// - Array generation
// - Coordinate conversion  
// - Geometric transformations
// - Boolean operations (when beneficial)
```

### Memory Management
```csharp
// Efficient large-scale processing
const int batchSize = 1000;
var results = new List<PathsD>();

for (int i = 0; i < largeDataSet.Count; i += batchSize)
{
    var batch = largeDataSet.Skip(i).Take(batchSize);
    var processed = GeoWrangler.processBatch(batch);
    results.Add(processed);
    
    // Optional: Force garbage collection between batches
    if (i % (batchSize * 10) == 0)
    {
        GC.Collect();
    }
}
```

## Common Usage Patterns

### Design Rule Checking
```csharp
// Minimum width check
double minWidth = 10.0;
foreach (var polygon in polygons)
{
    // Shrink and expand to find violations
    var shrunk = GeoWrangler.sizingOperation([polygon], -minWidth/2, 1000);
    var expanded = GeoWrangler.sizingOperation(shrunk, minWidth/2, 1000);
    
    // Differences indicate width violations
    var violations = Clipper.Difference([polygon], expanded, FillRule.NonZero);
    if (violations.Count > 0)
    {
        Console.WriteLine($"Width violation detected in polygon");
    }
}
```

### Proximity Analysis
```csharp
// Find polygons within specific distance
double checkDistance = 25.0;
var sourcePolygon = polygons[0];

foreach (var targetPolygon in polygons.Skip(1))
{
    double distance = GeoWrangler.minimumDistance(sourcePolygon, targetPolygon);
    if (distance < checkDistance)
    {
        Console.WriteLine($"Polygons too close: {distance:F2} < {checkDistance}");
    }
}
```

### Pattern Generation
```csharp
// Create complex repeated patterns
PathD baseUnit = GeoWrangler.rectangle(10, 10, 0, 0);

// Create staggered array
PathsD evenRows = GeoWrangler.makeArray(baseUnit, 10, 20.0, 5, 20.0);
PathsD oddRows = GeoWrangler.makeArray(baseUnit, 10, 20.0, 5, 20.0);
oddRows = GeoWrangler.move(oddRows, 10.0, 10.0); // Offset odd rows

// Combine into final pattern
var pattern = new PathsD();
pattern.AddRange(evenRows);
pattern.AddRange(oddRows);
```

## Error Handling and Debugging

### Validation
```csharp
// Validate polygon before processing
bool isValid = GeoWrangler.validatePolygon(polygon);
if (!isValid)
{
    Console.WriteLine("Invalid polygon detected");
    polygon = GeoWrangler.sanitizePolygon(polygon);
}

// Check for self-intersections
bool hasSelfIntersections = GeoWrangler.hasSelfIntersections(polygon);
if (hasSelfIntersections)
{
    // Fix or reject polygon
    polygon = GeoWrangler.removeSelfIntersections(polygon);
}
```

### Tolerance Management
```csharp
// Set appropriate tolerances for operations
const double geometryTolerance = 0.001;    // For coordinate comparisons
const double areaTolerancePercent = 1.0;   // 1% area difference tolerance

// Use consistent resolution for operations  
const double operationResolution = 1000.0; // 1000 units per user unit
```

## Integration Examples

### With geoCore (File I/O)
```csharp
using geoCoreLib;

// Load geometry from GDSII
var geoCore = /* load from file */;
var drawing = geoCore.getDrawing();
var cell = drawing.cellList[0];

// Extract polygons for processing
var polygons = new PathsD();
foreach (var element in cell.elementList.OfType<GCPolygon>())
{
    polygons.Add(element.pointarray);
}

// Process with geoWrangler
var processed = GeoWrangler.sizingOperation(polygons, 5.0, 1000);

// Add back to geoCore
foreach (var polygon in processed)
{
    cell.addPolygon(GeoWrangler.path64FromPathD(polygon, 1000), 2, 0);
}
```

### With shapeEngine
```csharp
using shapeEngine;

// Generate shapes with shapeEngine
var shapeLib = new ShapeLibrary();
var shape = shapeLib.generateShape(/* parameters */);

// Process with geoWrangler
var processed = GeoWrangler.move(shape, 100, 100);
processed = GeoWrangler.rotate(processed, 45.0);

// Apply boolean operations
var final = Clipper.Union(processed, FillRule.NonZero);
```

## Constants and Configuration

### Important Constants
```csharp
public static class Constants
{
    public const double tolerance = 1e-10;     // Geometric tolerance
    public const double pi = Math.PI;          // Pi constant
    public const double deg2rad = Math.PI / 180.0;  // Degree to radian conversion
    public const double rad2deg = 180.0 / Math.PI;  // Radian to degree conversion
}
```

## Dependencies
- **Clipper2Lib** - Core polygon operations and data types
- **geoLib** - Basic geometric primitives  
- **utility** - Mathematical and utility functions
- **.NET 8.0** - Target framework

## Related Libraries
- **[geoLib](geoLib-API.md)** - Basic geometric primitives used by geoWrangler
- **[clipper](clipper-API.md)** - Advanced polygon clipping operations
- **[shapeEngine](shapeEngine-API.md)** - Shape generation using geoWrangler operations
- **[geoCore](geoCore-API.md)** - File I/O that works with geoWrangler-processed geometry