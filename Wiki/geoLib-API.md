# geoLib API Documentation

## Overview

The **geoLib** library provides basic geometric primitives designed for cross-platform compatibility. It was created to avoid dependencies on platform-specific libraries like System.Drawing, making it suitable for headless environments and multi-platform operations.

## Key Features

- **Cross-platform compatibility** - Works on Windows, Linux, and macOS
- **Integer and floating-point support** - Provides both precision modes
- **Lightweight design** - Minimal dependencies
- **Clipper2Lib integration** - Uses modern polygon library data types

## Namespace
```csharp
using geoLib;
```

## Core Classes and Structures

### MyVertex
Represents a vertex in 2D space with additional geometric properties for shape definition.

#### Properties
```csharp
public double X { get; set; }                    // X coordinate
public double Y { get; set; }                    // Y coordinate
public typeDirection direction { get; set; }     // Direction type
public bool vertical { get; set; }              // True if attached to vertical edge
public bool inner { get; set; }                 // True if inner vertex
public typeVertex type { get; set; }            // Vertex type (corner or center)
public bool xBiasApplied { get; set; }          // X bias tracking
public bool yBiasApplied { get; set; }          // Y bias tracking
```

#### Constructors
```csharp
// Create vertex with full properties
MyVertex(double X, double Y, typeDirection direction, bool vertical, bool inner, typeVertex type)

// Copy constructor
MyVertex(MyVertex source)
```

#### Usage Example
```csharp
// Create a corner vertex
var vertex = new MyVertex(100.0, 50.0, typeDirection.up1, false, false, typeVertex.corner);

// Copy existing vertex
var vertexCopy = new MyVertex(vertex);
```

### MyRound
Represents rounding information for geometric elements.

#### Properties
```csharp
public int index { get; set; }           // Index identifier
public int verFace { get; set; }         // Vertical face identifier  
public int horFace { get; set; }         // Horizontal face identifier
public double MaxRadius { get; set; }    // Maximum radius for rounding
public typeRound direction { get; set; } // Rounding direction (inner/exterior)
```

### GeoLibArray
Array configuration with integer coordinates for geometric operations.

#### Properties
```csharp
public Point64 point { get; set; }      // Base point for the array
public Point64 pitch { get; init; }     // Spacing between elements (read-only after init)
public Point64 count { get; init; }     // Number of elements in each direction
```

#### Usage Example
```csharp
var array = new GeoLibArray 
{
    point = new Point64(0, 0),
    pitch = new Point64(100, 100),  // 100 unit spacing in both directions
    count = new Point64(5, 3)       // 5x3 array
};
```

### GeoLibArrayF
Array configuration with floating-point coordinates for precise positioning.

#### Properties
```csharp
public PointD point { get; set; }       // Base point (floating-point)
public PointD pitch { get; set; }       // Spacing between elements
public Point64 count { get; set; }      // Number of elements in each direction
```

#### Usage Example
```csharp
var preciseArray = new GeoLibArrayF
{
    point = new PointD(0.5, 0.5),
    pitch = new PointD(10.25, 15.75),
    count = new Point64(4, 6)
};
```

### GeoLibVector3
3D vector implementation for geometric calculations.

#### Properties
```csharp
public double x { get; set; }  // X component
public double y { get; set; }  // Y component  
public double z { get; set; }  // Z component
```

#### Constructors
```csharp
// Copy constructor
GeoLibVector3(GeoLibVector3 source)

// Integer coordinates
GeoLibVector3(int X, int Y, int Z)

// Floating-point coordinates
GeoLibVector3(double X, double Y, double Z)
```

#### Usage Example
```csharp
// Create vectors
var origin = new GeoLibVector3(0, 0, 0);
var point = new GeoLibVector3(10.5, 20.3, 5.0);
var copy = new GeoLibVector3(point);
```

### GeoLibMatrix
Matrix operations for 2D transformations.

#### Properties
```csharp
public double[] m { get; set; }  // 6-element transformation matrix [m11, m12, m21, m22, dx, dy]
```

#### Constructor
```csharp
GeoLibMatrix(float m11, float m12, float m21, float m22, float dx, float dy)
```

#### Methods
```csharp
public void Rotate(double ang)              // Rotate by angle in radians
public PointD transform(PointD inputPt)     // Transform a point
public Point64 transform(Point64 inputPt)   // Transform integer point
```

#### Usage Example
```csharp
// Create identity matrix with translation
var matrix = new GeoLibMatrix(1.0f, 0.0f, 0.0f, 1.0f, 10.0f, 20.0f);

// Rotate by 45 degrees
matrix.Rotate(Math.PI / 4);

// Transform points
var transformed = matrix.transform(new PointD(100, 100));
```

### GeoLibRectangle
Cross-platform rectangle implementation avoiding System.Drawing dependencies.

#### Properties
```csharp
public Point64 Location { get; set; }  // Top-left corner location
public int Width { get; set; }         // Rectangle width
public int Height { get; set; }        // Rectangle height
```

#### Constructors
```csharp
GeoLibRectangle()                                          // Default (0,0,0,0)
GeoLibRectangle(int x, int y, int width, int height)      // With dimensions
```

#### Usage Example
```csharp
// Create rectangle
var rect = new GeoLibRectangle(10, 20, 100, 50);

// Access properties
Console.WriteLine($"Area: {rect.Width * rect.Height}");
Console.WriteLine($"Top-left: ({rect.Location.X}, {rect.Location.Y})");
```

## Enumerations

### typeDirection
Direction types for geometry elements.
```csharp
public enum typeDirection { left1, right1, up1, down1, tilt1 }
```

### typeVertex  
Vertex type classification.
```csharp
public enum typeVertex { corner, center }
```

### typeRound
Rounding types for geometric elements.
```csharp
public enum typeRound { inner, exter }
```

## Integration with Clipper2Lib

geoLib is designed to work seamlessly with Clipper2Lib data types:

- Uses `Point64` for integer coordinates
- Uses `PointD` for floating-point coordinates
- Compatible with `PathD` and `PathsD` for polygon operations

```csharp
// Convert geoLib structures to Clipper2Lib paths
var vertices = new List<MyVertex> { /* vertices */ };
var path = new PathD();
foreach (var vertex in vertices)
{
    path.Add(new PointD(vertex.X, vertex.Y));
}
```

## Best Practices

### Coordinate System Choice
- Use **integer coordinates** (`Point64`) for precise manufacturing/layout applications
- Use **floating-point coordinates** (`PointD`) for graphics and approximate calculations
- Consider scaling factors when converting between coordinate systems

### Memory Management
- geoLib structures are lightweight value types where possible
- Copy constructors create deep copies to avoid reference issues
- Arrays use read-only properties after initialization to prevent accidental modification

### Thread Safety
- Most geoLib structures are immutable or provide value semantics
- Matrix operations should be performed on separate instances in multi-threaded scenarios
- No global state is maintained by the library

## Common Usage Patterns

### Creating Shape Vertices
```csharp
var vertices = new List<MyVertex>
{
    new MyVertex(0, 0, typeDirection.right1, false, false, typeVertex.corner),
    new MyVertex(100, 0, typeDirection.up1, true, false, typeVertex.corner),
    new MyVertex(100, 50, typeDirection.left1, false, false, typeVertex.corner),
    new MyVertex(0, 50, typeDirection.down1, true, false, typeVertex.corner)
};
```

### Array Generation
```csharp
// Create regular grid
var grid = new GeoLibArray
{
    point = new Point64(0, 0),
    pitch = new Point64(100, 100),
    count = new Point64(10, 10)
};

// Generate array positions
var positions = new List<Point64>();
for (int x = 0; x < grid.count.X; x++)
{
    for (int y = 0; y < grid.count.Y; y++)
    {
        positions.Add(new Point64(
            grid.point.X + x * grid.pitch.X,
            grid.point.Y + y * grid.pitch.Y
        ));
    }
}
```

### Transformation Matrices
```csharp
// Create transformation pipeline
var translate = new GeoLibMatrix(1, 0, 0, 1, 50, 100);  // Translate by (50, 100)
var rotate = new GeoLibMatrix(1, 0, 0, 1, 0, 0);        // Identity matrix
rotate.Rotate(Math.PI / 6);                              // Rotate by 30 degrees

// Apply transformations
var point = new PointD(10, 10);
point = translate.transform(point);
point = rotate.transform(point);
```

## Dependencies
- **Clipper2Lib** - For Point64 and PointD data types
- **.NET 8.0** - Target framework

## Related Libraries
- **[geoWrangler](geoWrangler-API.md)** - Advanced geometry operations using geoLib primitives
- **[shapeEngine](shapeEngine-API.md)** - Shape generation using geoLib structures
- **[clipper](clipper-API.md)** - Polygon operations on geoLib-compatible data types