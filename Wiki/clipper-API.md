# clipper API Documentation

## Overview

The **clipper** library is a comprehensive C# port of the Clipper2 polygon clipping library by Angus Johnson. It provides high-performance, robust polygon boolean operations, path offsetting, and geometric analysis. This library is essential for advanced polygon processing in CAD, graphics, and computational geometry applications.

## Key Features

- **Boolean Operations** - Union, intersection, difference, XOR with high precision
- **Path Offsetting** - Inflate/deflate polygons with various join and end types  
- **Robust Algorithms** - Handles complex polygons, self-intersections, and edge cases
- **High Performance** - Optimized for speed with integer and floating-point variants
- **Point-in-Polygon** - Fast containment testing and polygon analysis
- **Multiple Precision** - Integer (Point64/Path64) and floating-point (PointD/PathD) support

## Namespace
```csharp
using Clipper2Lib;
```

## Core Data Types

### Point Types
```csharp
// Integer precision point
public struct Point64
{
    public long X { get; set; }
    public long Y { get; set; }
    
    // Constructors
    public Point64(long x, long y)
    public Point64(double x, double y)  // With rounding
    public Point64(Point64 pt)          // Copy constructor
    
    // Operators
    public static Point64 operator +(Point64 lhs, Point64 rhs)
    public static Point64 operator -(Point64 lhs, Point64 rhs)
    public static bool operator ==(Point64 lhs, Point64 rhs)
    public static bool operator !=(Point64 lhs, Point64 rhs)
}

// Floating-point precision point  
public struct PointD
{
    public double x { get; set; }
    public double y { get; set; }
    
    // Constructors
    public PointD(double x, double y)
    public PointD(PointD pt)           // Copy constructor  
    public PointD(Point64 pt)          // From integer point
    
    // Operators (similar to Point64)
}
```

### Path Types
```csharp
// Collections of points forming polygons or polylines
public class Path64 : List<Point64> { }      // Integer precision path
public class Paths64 : List<Path64> { }      // Collection of integer paths

public class PathD : List<PointD> { }        // Floating-point precision path  
public class PathsD : List<PathD> { }        // Collection of floating-point paths
```

### Enumerations

#### ClipType - Boolean Operation Types
```csharp
public enum ClipType
{
    None,         // No operation
    Intersection, // A ∩ B - Common area
    Union,        // A ∪ B - Combined area  
    Difference,   // A - B - Subject minus clip
    Xor           // A ⊕ B - Exclusive or (non-overlapping areas)
}
```

#### FillRule - Polygon Fill Rules
```csharp
public enum FillRule
{
    EvenOdd,    // Alternating fill (standard graphics rule)
    NonZero,    // Non-zero winding (CAD standard)
    Positive,   // Positive winding only
    Negative    // Negative winding only  
}
```

#### JoinType - Path Offsetting Corner Styles
```csharp
public enum JoinType
{
    Miter,  // Sharp corners with miter limit
    Square, // Square corners  
    Bevel,  // Beveled (cut-off) corners
    Round   // Rounded corners
}
```

#### EndType - Path End Styles
```csharp
public enum EndType
{
    Polygon, // Closed polygon (default)
    Joined,  // Open path with joined ends
    Butt,    // Flat end caps
    Square,  // Square end caps (extend beyond path)
    Round    // Rounded end caps
}
```

## Main Classes

### Clipper (Static Class)
Simplified interface providing the most commonly used operations.

#### Boolean Operations
```csharp
// Intersection - Find common area
public static Paths64 Intersect(Paths64 subject, Paths64 clip, FillRule fillRule)
public static PathsD Intersect(PathsD subject, PathsD clip, FillRule fillRule, int precision = 2)

// Union - Combine areas  
public static Paths64 Union(Paths64 subject, FillRule fillRule)                    // Self-union
public static Paths64 Union(Paths64 subject, Paths64 clip, FillRule fillRule)     // Two-operand
public static PathsD Union(PathsD subject, FillRule fillRule)
public static PathsD Union(PathsD subject, PathsD clip, FillRule fillRule, int precision = 2)

// Difference - Subtract clip from subject
public static Paths64 Difference(Paths64 subject, Paths64 clip, FillRule fillRule)
public static PathsD Difference(PathsD subject, PathsD clip, FillRule fillRule, int precision = 2)

// XOR - Exclusive or (non-overlapping areas)
public static Paths64 Xor(Paths64 subject, Paths64 clip, FillRule fillRule)
public static PathsD Xor(PathsD subject, PathsD clip, FillRule fillRule, int precision = 2)
```

#### Usage Examples
```csharp
// Create two overlapping rectangles
var rect1 = new Path64 
{
    new Point64(0, 0), new Point64(100, 0), 
    new Point64(100, 100), new Point64(0, 100)
};

var rect2 = new Path64
{
    new Point64(50, 50), new Point64(150, 50),
    new Point64(150, 150), new Point64(50, 150)
};

var subjects = new Paths64 { rect1 };
var clips = new Paths64 { rect2 };

// Boolean operations
var intersection = Clipper.Intersect(subjects, clips, FillRule.NonZero);
var union = Clipper.Union(subjects, clips, FillRule.NonZero);
var difference = Clipper.Difference(subjects, clips, FillRule.NonZero);
var xor = Clipper.Xor(subjects, clips, FillRule.NonZero);

Console.WriteLine($"Intersection: {intersection.Count} polygons");
Console.WriteLine($"Union: {union.Count} polygons");
Console.WriteLine($"Difference: {difference.Count} polygons");
Console.WriteLine($"XOR: {xor.Count} polygons");
```

#### Path Offsetting (Inflation/Deflation)
```csharp
// Offset paths (positive = inflate, negative = deflate)
public static Paths64 InflatePaths(Paths64 paths, double delta, JoinType joinType,
    EndType endType, double miterLimit = 2.0, double arcTolerance = 0.0)
    
public static PathsD InflatePaths(PathsD paths, double delta, JoinType joinType, 
    EndType endType, double miterLimit = 2.0, double arcTolerance = 0.0, int precision = 2)
```

#### Usage Examples
```csharp
// Create a polygon to offset
var polygon = new Paths64 
{
    new Path64 
    {
        new Point64(50, 50), new Point64(150, 50),
        new Point64(150, 150), new Point64(50, 150)
    }
};

// Inflate (expand) with rounded corners
var inflated = Clipper.InflatePaths(polygon, 20.0, JoinType.Round, EndType.Polygon);

// Deflate (shrink) with sharp corners
var deflated = Clipper.InflatePaths(polygon, -10.0, JoinType.Miter, EndType.Polygon);

// Offset open paths with different end styles  
var openPath = new Paths64 
{
    new Path64 { new Point64(0, 50), new Point64(100, 50), new Point64(200, 100) }
};

var offsetOpen = Clipper.InflatePaths(openPath, 15.0, JoinType.Round, EndType.Round);
```

### ClipperOffset
Advanced path offsetting with detailed control over parameters.

#### Constructor and Configuration
```csharp
public ClipperOffset(double miterLimit = 2.0, double arcTolerance = 0.0)

// Configuration methods
public void AddPath(Path64 path, JoinType joinType, EndType endType)
public void AddPaths(Paths64 paths, JoinType joinType, EndType endType)
public void Execute(double delta, Paths64 solution)
public void Execute(double delta, PolyTree64 polytree)
public void Clear()
```

#### Advanced Usage Example
```csharp
var offsetter = new ClipperOffset(miterLimit: 3.0, arcTolerance: 0.25);

// Add different path types with specific join/end styles
offsetter.AddPath(outerPath, JoinType.Round, EndType.Polygon);
offsetter.AddPath(innerPath, JoinType.Miter, EndType.Polygon);  
offsetter.AddPath(openPath, JoinType.Square, EndType.Square);

// Execute with multiple offsets
var results1 = new Paths64();
offsetter.Execute(10.0, results1);

var results2 = new Paths64();
offsetter.Execute(20.0, results2);

// Clear and reuse
offsetter.Clear();
```

### Clipper64 and ClipperD 
Advanced clipper engines for complex operations.

#### Clipper64 - Integer Precision Engine
```csharp
public class Clipper64 : ClipperBase
{
    // Path management
    public void AddPath(Path64 path, PathType pathType, bool isOpen = false)
    public void AddPaths(Paths64 paths, PathType pathType, bool isOpen = false)  
    public void AddSubject(Paths64 paths)
    public void AddClip(Paths64 paths)
    public void AddOpenSubject(Paths64 paths)
    
    // Execution
    public bool Execute(ClipType clipType, FillRule fillRule, Paths64 closedPaths)
    public bool Execute(ClipType clipType, FillRule fillRule, Paths64 closedPaths, Paths64 openPaths)
    public bool Execute(ClipType clipType, FillRule fillRule, PolyTree64 polytree)
    public bool Execute(ClipType clipType, FillRule fillRule, PolyTree64 polytree, Paths64 openPaths)
    
    // Configuration
    public bool PreserveCollinear { get; set; }
    public bool ReverseSolution { get; set; }
    
    // Cleanup
    public void Clear()
}
```

#### Usage Example
```csharp
var clipper = new Clipper64();

// Add multiple subject and clip paths
clipper.AddSubject(subjectPolygons);
clipper.AddClip(clipPolygons);
clipper.AddOpenSubject(openPaths);

// Configure options
clipper.PreserveCollinear = true;
clipper.ReverseSolution = false;

// Execute with separate open and closed results
var closedResults = new Paths64();
var openResults = new Paths64();

bool success = clipper.Execute(ClipType.Union, FillRule.NonZero, closedResults, openResults);

if (success)
{
    Console.WriteLine($"Closed paths: {closedResults.Count}");
    Console.WriteLine($"Open paths: {openResults.Count}");
}
else
{
    Console.WriteLine("Clipping operation failed");
}

// Execute with hierarchical result (polytree)
var polytree = new PolyTree64();
success = clipper.Execute(ClipType.Union, FillRule.NonZero, polytree);

if (success)
{
    Console.WriteLine($"Polytree children: {polytree.Count}");
    TraversePolytree(polytree);
}

clipper.Clear();
```

#### ClipperD - Floating-Point Engine
```csharp
public class ClipperD : ClipperBase
{
    public ClipperD(int precision = 2) // Decimal places for rounding
    
    // Similar methods to Clipper64 but using PathsD, PathD, PolyTreeD
    public void AddPath(PathD path, PathType pathType, bool isOpen = false)
    public void AddPaths(PathsD paths, PathType pathType, bool isOpen = false)
    public bool Execute(ClipType clipType, FillRule fillRule, PathsD closedPaths)
    // ... etc
}
```

## Advanced Features

### PolyTree - Hierarchical Polygon Results
Represents polygon hierarchies (holes within polygons, islands within holes, etc.)

#### PolyTree64 and PolyTreeD
```csharp
public class PolyTree64 : PolyPath64
{
    public PolyPath64? GetFirst()        // Get first polygon
    public void Clear()                  // Clear all data
    public int Count { get; }            // Number of top-level polygons
}

public class PolyPath64
{
    public PolyPath64? Parent { get; }   // Parent polygon (null for top-level)
    public Path64 Polygon { get; }       // The actual polygon path
    public List<PolyPath64> Childs { get; } // Child polygons (holes or islands)
    public bool IsHole { get; }          // True if this is a hole
    public PolyPath64? Next { get; }     // Next sibling polygon
    public double Area { get; }          // Signed area (negative for holes)
}
```

#### Usage Example
```csharp
var clipper = new Clipper64();
clipper.AddSubject(complexPolygons);

var polytree = new PolyTree64();
clipper.Execute(ClipType.Union, FillRule.NonZero, polytree);

// Traverse hierarchy
void TraversePolytree(PolyPath64 polypath, int level = 0)
{
    string indent = new string(' ', level * 2);
    
    foreach (var child in polypath.Childs)
    {
        Console.WriteLine($"{indent}Polygon: {child.Polygon.Count} vertices, " +
                         $"Area: {child.Area:F2}, IsHole: {child.IsHole}");
        
        if (child.Childs.Count > 0)
        {
            TraversePolytree(child, level + 1);
        }
    }
}

TraversePolytree(polytree);
```

### Point-in-Polygon Testing
```csharp
public enum PointInPolygonResult
{
    IsOn,      // Point lies exactly on polygon boundary
    IsInside,  // Point is inside polygon
    IsOutside  // Point is outside polygon  
}

// Test point location
public static PointInPolygonResult PointInPolygon(Point64 pt, Path64 polygon)
```

#### Usage Example
```csharp
var polygon = new Path64 
{
    new Point64(0, 0), new Point64(100, 0),
    new Point64(100, 100), new Point64(0, 100)
};

var testPoints = new Point64[]
{
    new Point64(50, 50),   // Inside
    new Point64(150, 50),  // Outside  
    new Point64(0, 50),    // On boundary
    new Point64(-10, 50)   // Outside
};

foreach (var point in testPoints)
{
    var result = InternalClipper.PointInPolygon(point, polygon);
    Console.WriteLine($"Point ({point.X}, {point.Y}): {result}");
}
```

### Area Calculation
```csharp
// Calculate signed area (positive = counterclockwise, negative = clockwise)
public static double Area(Path64 path)
public static double Area(PathD path)

// Usage
var polygon = new Path64 { /* vertices */ };
double area = Clipper.Area(polygon);

if (area > 0)
    Console.WriteLine($"Counterclockwise polygon, area: {area:F2}");
else if (area < 0)  
    Console.WriteLine($"Clockwise polygon, area: {Math.Abs(area):F2}");
else
    Console.WriteLine("Degenerate polygon (zero area)");
```

### Path Simplification and Cleanup
```csharp
// Remove duplicate points
public static Path64 StripDuplicates(Path64 path, bool isClosedPath = true)
public static PathD StripDuplicates(PathD path, bool isClosedPath = true)

// Simplify paths (remove unnecessary vertices)
public static Paths64 SimplifyPaths(Paths64 paths, FillRule fillRule, bool strictlySimple = false)
public static PathsD SimplifyPaths(PathsD paths, FillRule fillRule, bool strictlySimple = false)

// Usage
var noisyPath = new Path64 { /* path with many close points */ };
var cleanPath = Clipper.StripDuplicates(noisyPath);

var complexPaths = new Paths64 { /* complex polygons */ };
var simplified = Clipper.SimplifyPaths(complexPaths, FillRule.NonZero, strictlySimple: true);
```

## Performance Optimization

### Precision Considerations
```csharp
// Integer precision (fastest, most precise)
var intPaths = new Paths64();
var intResult = Clipper.Union(intPaths, FillRule.NonZero);

// Floating-point precision (convenient but slower)
var floatPaths = new PathsD(); 
var floatResult = Clipper.Union(floatPaths, FillRule.NonZero, precision: 2);

// Convert between precisions
Paths64 toInt = Clipper.ScalePaths64(floatPaths, 100.0);    // Scale up to preserve precision
PathsD toFloat = Clipper.ScalePathsD(intPaths, 0.01);       // Scale down to floating-point
```

### Batch Operations
```csharp
// Reuse clipper instances for better performance
var clipper = new Clipper64();

foreach (var batch in largeBatches)
{
    clipper.Clear();
    clipper.AddSubject(batch.subjects);
    clipper.AddClip(batch.clips);
    
    var results = new Paths64();
    clipper.Execute(ClipType.Union, FillRule.NonZero, results);
    
    ProcessResults(results);
}
```

### Memory Management
```csharp
// Efficient large-scale processing
public static void ProcessLargeDataset(IEnumerable<Paths64> datasets)
{
    const int batchSize = 100;
    var clipper = new Clipper64();
    
    var batch = new List<Paths64>(batchSize);
    
    foreach (var dataset in datasets)
    {
        batch.Add(dataset);
        
        if (batch.Count >= batchSize)
        {
            ProcessBatch(clipper, batch);
            batch.Clear();
            
            // Optional: Force garbage collection for very large datasets
            if (batch.Count % (batchSize * 10) == 0)
            {
                GC.Collect();
                GC.WaitForPendingFinalizers();
            }
        }
    }
    
    if (batch.Count > 0)
    {
        ProcessBatch(clipper, batch);
    }
}

private static void ProcessBatch(Clipper64 clipper, List<Paths64> batch)
{
    clipper.Clear();
    
    foreach (var paths in batch)
    {
        clipper.AddSubject(paths);
    }
    
    var results = new Paths64();
    clipper.Execute(ClipType.Union, FillRule.NonZero, results);
    
    // Process results...
}
```

## Common Usage Patterns

### CAD Boolean Operations
```csharp
public static Paths64 PerformCADBoolean(Paths64 layerA, Paths64 layerB, string operation)
{
    return operation.ToUpper() switch
    {
        "AND" => Clipper.Intersect(layerA, layerB, FillRule.NonZero),
        "OR" => Clipper.Union(layerA, layerB, FillRule.NonZero),
        "NOT" => Clipper.Difference(layerA, layerB, FillRule.NonZero),
        "XOR" => Clipper.Xor(layerA, layerB, FillRule.NonZero),
        _ => throw new ArgumentException($"Unknown operation: {operation}")
    };
}
```

### Polygon Offsetting for Tolerances
```csharp
public static Paths64 ApplyManufacturingTolerances(Paths64 originalPaths, double tolerance)
{
    // Apply both positive and negative tolerance to find min/max boundaries
    var expanded = Clipper.InflatePaths(originalPaths, tolerance, JoinType.Round, EndType.Polygon);
    var contracted = Clipper.InflatePaths(originalPaths, -tolerance, JoinType.Round, EndType.Polygon);
    
    // Return envelope that contains both
    return Clipper.Union(expanded, contracted, FillRule.NonZero);
}
```

### Complex Shape Decomposition
```csharp
public static List<Path64> DecomposeComplexPolygon(Path64 complexPolygon)
{
    var paths = new Paths64 { complexPolygon };
    
    // Simplify to remove self-intersections
    var simplified = Clipper.SimplifyPaths(paths, FillRule.NonZero, strictlySimple: true);
    
    // Use polytree to handle holes properly
    var clipper = new Clipper64();
    clipper.AddSubject(simplified);
    
    var polytree = new PolyTree64();
    clipper.Execute(ClipType.Union, FillRule.NonZero, polytree);
    
    var result = new List<Path64>();
    ExtractPolygonsFromTree(polytree, result);
    
    return result;
}

private static void ExtractPolygonsFromTree(PolyPath64 polypath, List<Path64> result)
{
    foreach (var child in polypath.Childs)
    {
        if (!child.IsHole && child.Polygon.Count > 0)
        {
            result.Add(child.Polygon);
        }
        
        ExtractPolygonsFromTree(child, result);
    }
}
```

## Error Handling and Validation

### Robust Operations
```csharp
public static Paths64 SafeBooleanOperation(Paths64 subject, Paths64 clip, ClipType operation)
{
    try
    {
        // Validate inputs
        if (subject == null || subject.Count == 0) 
            return new Paths64();
            
        // Clean inputs
        var cleanSubject = subject.Select(p => Clipper.StripDuplicates(p)).Where(p => p.Count >= 3).ToList();
        var cleanClip = clip?.Select(p => Clipper.StripDuplicates(p)).Where(p => p.Count >= 3).ToList();
        
        if (cleanSubject.Count == 0) 
            return new Paths64();
            
        // Perform operation
        var result = operation switch
        {
            ClipType.Union => cleanClip != null ? 
                Clipper.Union(new Paths64(cleanSubject), new Paths64(cleanClip), FillRule.NonZero) :
                Clipper.Union(new Paths64(cleanSubject), FillRule.NonZero),
            ClipType.Intersection => cleanClip != null ?
                Clipper.Intersect(new Paths64(cleanSubject), new Paths64(cleanClip), FillRule.NonZero) :
                new Paths64(),
            ClipType.Difference => cleanClip != null ?
                Clipper.Difference(new Paths64(cleanSubject), new Paths64(cleanClip), FillRule.NonZero) :
                new Paths64(cleanSubject),
            ClipType.Xor => cleanClip != null ?
                Clipper.Xor(new Paths64(cleanSubject), new Paths64(cleanClip), FillRule.NonZero) :
                new Paths64(cleanSubject),
            _ => new Paths64()
        };
        
        // Validate result
        return result.Where(p => p.Count >= 3 && Math.Abs(Clipper.Area(p)) > 1e-10).ToList();
    }
    catch (Exception ex)
    {
        Console.WriteLine($"Boolean operation failed: {ex.Message}");
        return new Paths64();
    }
}
```

## Integration Examples

### With geoWrangler
```csharp
using Clipper2Lib;
using geoWrangler;

// Combine Clipper operations with geoWrangler transformations
var polygon1 = GeoWrangler.rectangle(100, 50, 0, 0);
var polygon2 = GeoWrangler.ellipse(30, 30, 50, 25, 32);

// Convert to Clipper-compatible format
var paths1 = new PathsD { polygon1 };
var paths2 = new PathsD { polygon2 };

// Perform boolean operations
var union = Clipper.Union(paths1, paths2, FillRule.NonZero);
var intersection = Clipper.Intersect(paths1, paths2, FillRule.NonZero);

// Transform results back with geoWrangler
var transformedUnion = union.Select(p => GeoWrangler.rotate(p, 45.0)).ToList();
```

### With geoCore
```csharp
using Clipper2Lib;
using geoCoreLib;

// Extract polygons from geoCore, process with Clipper, add back
var geoCore = new GeoCore(/* load file */);
var drawing = geoCore.getDrawing();
var cell = drawing.cellList[0];

// Extract polygons
var polygons = new PathsD();
foreach (var element in cell.elementList.OfType<GCPolygon>())
{
    var pathD = element.pointarray.Select(pt => new PointD(pt.X * 1e-6, pt.Y * 1e-6)).ToList();
    polygons.Add(new PathD(pathD));
}

// Process with Clipper
var processed = Clipper.Union(polygons, FillRule.NonZero);

// Add back to geoCore  
foreach (var polygon in processed)
{
    var path64 = polygon.Select(pt => new Point64((long)(pt.x * 1e6), (long)(pt.y * 1e6))).ToList();
    cell.addPolygon(new Path64(path64), 2, 0);
}
```

## Dependencies
- **.NET 8.0** - Target framework
- **System.Collections.Generic** - Collection types
- **System.Runtime.CompilerServices** - Performance optimizations

## Related Libraries
- **[geoWrangler](geoWrangler-API.md)** - Uses Clipper2Lib data types and operations
- **[geoCore](geoCore-API.md)** - File I/O that works with Clipper-compatible polygon data  
- **[shapeEngine](shapeEngine-API.md)** - Shape generation that outputs Clipper-compatible geometry
- **[geoLib](geoLib-API.md)** - Basic geometric primitives that integrate with Clipper data types