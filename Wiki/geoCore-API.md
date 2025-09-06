# geoCore API Documentation

## Overview

The **geoCore** library is a comprehensive layout parsing and generation library for GDSII and OASIS file formats. It provides integer-based geometry representation to avoid precision issues and is inspired by various industry-standard libraries including layouteditor and KLayout.

## Key Features

- **Dual Format Support** - Handles both GDSII Stream and OASIS formats
- **Integer-Based Geometry** - Avoids floating-point precision issues 
- **Hierarchical Design** - Full support for cell references and arrays
- **Manufacturing Precision** - Uses database units for exact manufacturing specifications
- **Error Handling** - Comprehensive error reporting during file operations

## Namespace
```csharp
using geoCoreLib;
```

## Core Classes

### GeoCore
Main entry point for file operations and geometry management.

#### Properties
```csharp
public static double tolerance { get; set; }    // Default tolerance for operations (0.001)
public List<string> error_msgs { get; set; }    // Error messages from operations
```

#### File Format Support
```csharp
public enum fileType 
{ 
    gds,     // GDSII Stream format
    oasis    // OASIS format
}
```

#### Key Methods
```csharp
// File operations
public void setValid(bool val)              // Set validity state
public bool isValid()                       // Check if data is valid
public List<string> getStructureList()     // Get all structure names
public List<string> getActiveStructureLDList()  // Get layer/datatype list

// Data access
public GCDrawingfield getDrawing()         // Get the main drawing data
public void updateGeometry(int structure, int layer)  // Update geometry for specific layer
```

#### Usage Example
```csharp
var geoCore = new GeoCore();

// Load a GDSII file
if (geoCore.isValid())
{
    var structures = geoCore.getStructureList();
    Console.WriteLine($"Found {structures.Count} structures:");
    foreach (var structure in structures)
    {
        Console.WriteLine($"  - {structure}");
    }
}
```

### GCDrawingfield
Represents the main drawing container with all geometric data.

#### Properties
```csharp
// Core data
public List<GCCell> cellList { get; set; }      // All cells in the drawing
public int active_cell { get; set; }            // Currently active cell index
public string libname { get; set; }             // Library name

// Units and scaling
public double databaseunits { get; set; }       // Database units (default: 1E-9)
public double userunits { get; set; }           // User units (default: 1E-3)

// Timestamps
public short modyear, modmonth, modday { get; set; }     // Modification date
public short modhour, modmin, modsec { get; set; }       // Modification time
public short accyear, accmonth, accday { get; set; }     // Access date  
public short acchour, accmin, accsec { get; set; }       // Access time
```

#### Constants
```csharp
public const double default_databaseunits = 1E-9;  // 1 nanometer
public const double default_userunits = 1E-3;      // 1 millimeter
```

#### Methods
```csharp
public GCDrawingfield copy()                    // Create a copy
public void reset()                             // Reset to default state
public double getDrawingScale()                 // Get current scale factor
```

#### Usage Example
```csharp
var drawingField = new GCDrawingfield("MyLibrary");
Console.WriteLine($"Database units: {drawingField.databaseunits}");
Console.WriteLine($"User units: {drawingField.userunits}");
Console.WriteLine($"Scale factor: {drawingField.getDrawingScale()}");
```

### GCCell
Represents a single cell (structure) containing geometric elements.

#### Properties
```csharp
public string cellName { get; set; }                // Cell name
public bool saved { get; set; }                     // Saved state flag
public List<GCElement> elementList { get; set; }    // Geometric elements

// Timestamps (same as GCDrawingfield)
public short modyear, modmonth, modday { get; set; }
public short modhour, modmin, modsec { get; set; }
public short accyear, accmonth, accday { get; set; }
public short acchour, accmin, accsec { get; set; }
```

#### Methods
```csharp
// Basic operations
public void setName(string name)               // Set cell name
public string getName()                        // Get cell name

// Add geometric elements
public void addBox(int x, int y, int width, int height, int layer, int datatype)
public void addBox(Path64 points, int layer, int datatype)
public void addPolygon(Path64 points, int layer, int datatype) 
public void addPath(Path64 points, int pathType, int layer, int datatype)
public void addText(string text, int x, int y, int layer, int datatype)

// Cell references
public void addCellref(GCCellref cellRef)      // Add cell reference
public void addCellrefArray(GCCellrefArray cellRefArray)  // Add cell reference array
```

#### Usage Example
```csharp
var cell = new GCCell();
cell.setName("INVERTER");

// Add a rectangle
cell.addBox(0, 0, 1000, 500, 1, 0);  // 1000x500 rectangle on layer 1

// Add a polygon
var polygonPoints = new Path64
{
    new Point64(0, 0),
    new Point64(500, 0), 
    new Point64(500, 250),
    new Point64(0, 250)
};
cell.addPolygon(polygonPoints, 2, 0);
```

### GCElement
Base class for all geometric elements in a cell.

#### Derived Classes
- **GCBox** - Rectangular elements
- **GCPolygon** - Arbitrary polygon shapes  
- **GCPath** - Path elements with width
- **GCText** - Text labels
- **GCCellref** - References to other cells
- **GCCellrefArray** - Arrays of cell references

#### Common Properties
```csharp
public int layer_nr { get; set; }      // Layer number
public int datatype_nr { get; set; }   // Datatype number
public bool select { get; set; }       // Selection state
public bool pos { get; set; }          // Positive geometry
```

### GCBox
Rectangular geometry element.

#### Properties
```csharp
public Point64 x { get; set; }         // Bottom-left corner
public Point64 y { get; set; }         // Top-right corner  
public int width { get; set; }         // Width
public int height { get; set; }        // Height
```

#### Constructor
```csharp
GCBox(int x, int y, int width, int height, int layer, int datatype)
```

### GCPolygon
Arbitrary polygon geometry element.

#### Properties
```csharp
public Path64 pointarray { get; set; } // Polygon vertices
```

#### Constructor
```csharp
GCPolygon(Path64 points, int layer, int datatype)
```

### GCPath
Path element with specified width.

#### Properties
```csharp
public Path64 pointarray { get; set; } // Path points
public int width { get; set; }         // Path width
public int cap { get; set; }           // End cap style
```

#### Cap Styles
```csharp
public const int cap_butt = 0;         // Flush end caps
public const int cap_round = 1;        // Round end caps  
public const int cap_square = 2;       // Square end caps
```

### GCCellref
Reference to another cell with transformation.

#### Properties
```csharp
public string cell_ref { get; set; }   // Referenced cell name
public Point64 x { get; set; }         // Reference position
public double angle { get; set; }      // Rotation angle
public double mag { get; set; }        // Magnification factor
public bool x_flipped { get; set; }    // X-axis mirror
public bool y_flipped { get; set; }    // Y-axis mirror
```

### GCCellrefArray
Array of cell references.

#### Properties  
```csharp
public GCCellref cellRef { get; set; } // Base cell reference
public Point64 count { get; set; }     // Array dimensions (nx, ny)
public Point64 x_pitch { get; set; }   // X spacing
public Point64 y_pitch { get; set; }   // Y spacing
```

## File Operations

### Loading Files

#### GDSII Files
```csharp
var gdsReader = new gdsReader();
var geoCore = gdsReader.load("design.gds");

if (geoCore.isValid())
{
    var drawing = geoCore.getDrawing();
    Console.WriteLine($"Loaded library: {drawing.libname}");
}
else
{
    foreach (var error in geoCore.error_msgs)
    {
        Console.WriteLine($"Error: {error}");
    }
}
```

#### OASIS Files
```csharp
var oasReader = new oasReader();
var geoCore = oasReader.load("design.oas");
```

### Writing Files

#### GDSII Export
```csharp
var gdsWriter = new gdsWriter(geoCore, "output.gds");
gdsWriter.write();
```

#### OASIS Export
```csharp
var oasWriter = new oasWriter(geoCore, "output.oas");
oasWriter.write();
```

## Coordinate System and Units

### Database Units
- Default: **1 nanometer** (1E-9 meters)
- Used for precise manufacturing coordinates
- All internal coordinates are integers in database units

### User Units  
- Default: **1 millimeter** (1E-3 meters)
- Used for display and user interaction
- Conversion factor: `userunits / databaseunits`

### Coordinate Conversion
```csharp
// Convert user coordinates to database coordinates
double userCoord = 10.5;  // 10.5 mm
long dbCoord = (long)(userCoord * drawingField.userunits / drawingField.databaseunits);

// Convert database coordinates to user coordinates  
long dbCoord = 10500000;  // 10.5 million nanometers
double userCoord = dbCoord * drawingField.databaseunits / drawingField.userunits;
```

## Hierarchical Design

### Cell References
```csharp
// Create a basic cell
var transistor = new GCCell();
transistor.setName("NMOS");
transistor.addBox(0, 0, 1000, 500, 1, 0);  // Gate
transistor.addBox(-200, -100, 1400, 100, 2, 0);  // Source/Drain

// Create a higher-level cell that references the transistor
var inverter = new GCCell();
inverter.setName("INV_X1");

// Add transistor instance
var nmos_ref = new GCCellref();
nmos_ref.cell_ref = "NMOS";
nmos_ref.x = new Point64(0, 0);
nmos_ref.angle = 0;
nmos_ref.mag = 1.0;
inverter.addCellref(nmos_ref);

// Add PMOS (flipped version)
var pmos_ref = new GCCellref();  
pmos_ref.cell_ref = "NMOS";
pmos_ref.x = new Point64(0, 1000);
pmos_ref.y_flipped = true;
inverter.addCellref(pmos_ref);
```

### Cell Arrays
```csharp
// Create array of inverters
var array = new GCCellrefArray();
array.cellRef = new GCCellref { cell_ref = "INV_X1" };
array.count = new Point64(8, 4);        // 8x4 array
array.x_pitch = new Point64(2000, 0);   // 2µm spacing in X
array.y_pitch = new Point64(0, 3000);   // 3µm spacing in Y

var arrayCell = new GCCell();
arrayCell.setName("INV_ARRAY");
arrayCell.addCellrefArray(array);
```

## Layer Management

### Layer/Datatype Combinations
```csharp
// Get all layer/datatype combinations in active structure
var ldList = geoCore.getActiveStructureLDList();
foreach (var ld in ldList)
{
    Console.WriteLine($"Layer/Datatype: {ld}");
}

// Typical layer assignments
const int ACTIVE = 1;      // Active silicon
const int POLY = 2;        // Polysilicon  
const int METAL1 = 3;      // First metal layer
const int VIA1 = 4;        // Via between metal layers
const int METAL2 = 5;      // Second metal layer
```

### Layer Names
```csharp
// Set up layer name mapping for better readability
var layerNames = new Dictionary<string, string>
{
    {"1:0", "ACTIVE"},
    {"2:0", "POLY"},
    {"3:0", "METAL1"},
    {"4:0", "VIA1"},
    {"5:0", "METAL2"}
};
```

## Error Handling

### Validation
```csharp
if (!geoCore.isValid())
{
    Console.WriteLine("GeoCore data is invalid!");
    foreach (var error in geoCore.error_msgs)
    {
        Console.WriteLine($"  - {error}");
    }
    return;
}
```

### Common Error Scenarios
- **File format errors** - Invalid GDSII/OASIS structure
- **Missing cell references** - References to undefined cells
- **Invalid coordinates** - Coordinates outside valid range
- **Layer conflicts** - Invalid layer/datatype combinations

## Integration with Other Libraries

### With geoWrangler
```csharp
using geoWrangler;

// Extract geometry from geoCore for processing
var cell = drawing.cellList[0];
foreach (var element in cell.elementList)
{
    if (element is GCPolygon poly)
    {
        // Use geoWrangler for geometric operations
        var bounds = GeoWrangler.getBounds(poly.pointarray);
        var area = GeoWrangler.area(poly.pointarray);
    }
}
```

### With Clipper2Lib
```csharp
using Clipper2Lib;

// geoCore uses Clipper2Lib data types natively
Path64 polygon1 = /* from GCPolygon */;
Path64 polygon2 = /* from another element */;

// Perform boolean operations
var union = Clipper.Union(new PathsD { polygon1, polygon2 }, FillRule.NonZero);
```

## Performance Considerations

### Large File Handling
- Use streaming for very large files
- Process cells individually to manage memory
- Consider parallel processing for independent operations

### Memory Management
```csharp
// Reset drawing field to free memory
drawingField.reset();

// Clear error messages after processing
geoCore.error_msgs.Clear();
```

### Coordinate Precision
- Use integer coordinates throughout for exact precision
- Avoid floating-point calculations in critical paths
- Scale coordinates appropriately for the target process

## Best Practices

### File Organization
- Use meaningful cell names
- Organize hierarchy logically (bottom-up)
- Keep cell complexity manageable
- Use consistent layer assignments

### Error Handling
- Always check `isValid()` after operations
- Process error messages for debugging
- Validate coordinates before creating elements

### Performance  
- Batch similar operations together
- Use cell references for repeated patterns
- Minimize coordinate transformations
- Cache frequently accessed data

## Dependencies
- **geoLib** - Basic geometric primitives
- **geoWrangler** - Geometric operations and utilities
- **Clipper2Lib** - Polygon data types and operations
- **.NET 8.0** - Target framework

## Related Libraries
- **[geoWrangler](geoWrangler-API.md)** - Geometric transformations and analysis
- **[geoLib](geoLib-API.md)** - Basic geometric primitives
- **[clipper](clipper-API.md)** - Advanced polygon operations