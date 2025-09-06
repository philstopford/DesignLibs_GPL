# Eto.VeldridSurface API Documentation

## Overview

The **Eto.VeldridSurface** library provides a cross-platform 2D graphics viewport for end-user applications using the Eto.Forms framework and Veldrid rendering engine. This library enables high-performance, hardware-accelerated graphics rendering across Windows, macOS, and Linux platforms.

## Key Features

- **Cross-Platform Graphics** - Works on Windows, macOS, and Linux with native performance
- **Hardware Acceleration** - Uses Veldrid for GPU-accelerated rendering
- **Eto.Forms Integration** - Seamlessly integrates with Eto.Forms UI applications
- **2D Viewport Operations** - Pan, zoom, selection, and grid display
- **Polygon Rendering** - High-performance rendering of complex geometry
- **Customizable Appearance** - Configurable colors, grids, and visual settings

## Namespace
```csharp
using VeldridEto;
```

## Core Classes

### VeldridDriver
Main graphics driver that controls rendering to a VeldridSurface control.

#### Constructor
```csharp
public VeldridDriver(ref OVPSettings settings, ref VeldridSurface surface)
```

#### Key Methods
```csharp
public void SetUpVeldrid()              // Initialize Veldrid resources
public void updateViewport()           // Refresh the viewport display  
public void reset()                    // Reset viewport to default state
public void zoomExtents()              // Zoom to fit all geometry
public void setSelection(PointF location)  // Set selection point
```

#### Usage Example
```csharp
// Create viewport settings
var settings = new OVPSettings();

// Create Veldrid surface control
var surface = new VeldridSurface();

// Initialize driver
var driver = new VeldridDriver(ref settings, ref surface);
driver.SetUpVeldrid();

// Add to form
var form = new Form();
form.Content = surface;
```

### OVPSettings
Configuration class for viewport appearance and behavior.

#### Viewport Bounds
```csharp
public float minX { get; set; }     // Minimum X coordinate
public float maxX { get; set; }     // Maximum X coordinate  
public float minY { get; set; }     // Minimum Y coordinate
public float maxY { get; set; }     // Maximum Y coordinate
```

#### Colors and Appearance
```csharp
public Color minorGridColor { get; set; }        // Minor grid line color
public Color majorGridColor { get; set; }        // Major grid line color
public Color axisColor { get; set; }             // Axis line color
public Color backColor { get; set; }             // Background color
public Color selectionColor { get; set; }        // Selection highlight color
public Color inverSelectionColor { get; set; }   // Inverse selection color
```

#### Geometry Data
```csharp
public List<ovp_Poly>? polyList { get; set; }           // Main polygon list
public List<int>? polyListPtCount { get; set; }         // Point counts per polygon
public List<int>? polySourceIndex { get; set; }         // Source tracking
public List<bool>? polyMask { get; set; }               // Visibility mask

public List<ovp_Poly>? bgPolyList { get; set; }         // Background polygons
public List<ovp_Poly>? lineList { get; set; }           // Line geometry
public List<ovp_Poly>? tessPolyList { get; set; }       // Tessellated triangles
```

#### Configuration Methods
```csharp
// Grid and axes
public bool drawGrid()                  // Get grid visibility
public void drawGrid(bool val)          // Set grid visibility
public bool drawAxes()                  // Get axes visibility  
public void drawAxes(bool val)          // Set axes visibility

// Zoom and pan
public bool isZoomAndPanAllowed()       // Get zoom/pan state
public void allowZoomAndPan(bool val)   // Enable/disable zoom/pan
public void setCameraPos(PointF pos)    // Set camera position
public PointF getCameraPos()            // Get camera position

// Visual settings
public bool aA()                        // Get anti-aliasing state
public void aA(bool val)                // Set anti-aliasing
public bool filled()                    // Get filled polygon state
public void filled(bool val)            // Set filled polygon rendering
```

#### Usage Example
```csharp
var settings = new OVPSettings();

// Configure appearance
settings.backColor = Colors.White;
settings.minorGridColor = Colors.LightGray;
settings.majorGridColor = Colors.Gray;
settings.axisColor = Colors.Black;
settings.selectionColor = Colors.Blue;

// Enable features
settings.drawGrid(true);
settings.drawAxes(true);
settings.allowZoomAndPan(true);
settings.aA(true);
settings.filled(true);

// Set viewport bounds
settings.minX = -1000;
settings.maxX = 1000;
settings.minY = -1000;
settings.maxY = 1000;
```

### ovp_Poly
Polygon representation class for viewport rendering.

#### Properties
```csharp
public PathD poly { get; set; }         // Polygon geometry
public Color color { get; set; }        // Polygon color
public float alpha { get; set; }        // Transparency (0.0-1.0)
```

#### Constructor
```csharp
public ovp_Poly(PathD inputPoly, Color inputColor, float inputAlpha = 1.0f)
```

#### Usage Example
```csharp
// Create polygon from geometric data
var polygonPath = new PathD
{
    new PointD(0, 0),
    new PointD(100, 0), 
    new PointD(100, 100),
    new PointD(0, 100)
};

var polygon = new ovp_Poly(polygonPath, Colors.Red, 0.8f);

// Add to viewport settings
settings.polyList = new List<ovp_Poly> { polygon };
settings.polyListPtCount = new List<int> { polygonPath.Count };
```

## Integration with Geometry Libraries

### With geoWrangler
```csharp
using VeldridEto;
using geoWrangler;
using Clipper2Lib;

public static class GeometryViewport
{
    public static void DisplayGeometry(VeldridDriver driver, OVPSettings settings, PathsD geometry)
    {
        var polyList = new List<ovp_Poly>();
        var ptCounts = new List<int>();
        
        foreach (var path in geometry)
        {
            var poly = new ovp_Poly(path, Colors.Blue, 0.7f);
            polyList.Add(poly);
            ptCounts.Add(path.Count);
        }
        
        settings.polyList = polyList;
        settings.polyListPtCount = ptCounts;
        settings.changed = true;
        
        driver.updateViewport();
    }
    
    public static void DisplayBooleanResult(VeldridDriver driver, OVPSettings settings, 
        PathsD subject, PathsD clip, PathsD result)
    {
        var polyList = new List<ovp_Poly>();
        var ptCounts = new List<int>();
        
        // Subject polygons in blue
        foreach (var path in subject)
        {
            polyList.Add(new ovp_Poly(path, Colors.Blue, 0.5f));
            ptCounts.Add(path.Count);
        }
        
        // Clip polygons in green  
        foreach (var path in clip)
        {
            polyList.Add(new ovp_Poly(path, Colors.Green, 0.5f));
            ptCounts.Add(path.Count);
        }
        
        // Result polygons in red (filled)
        foreach (var path in result)
        {
            polyList.Add(new ovp_Poly(path, Colors.Red, 0.8f));
            ptCounts.Add(path.Count);
        }
        
        settings.polyList = polyList;
        settings.polyListPtCount = ptCounts;
        settings.changed = true;
        
        driver.updateViewport();
        driver.zoomExtents(); // Fit all geometry
    }
}
```

### With geoCore File Display
```csharp
using VeldridEto;
using geoCoreLib;

public static class FileViewer
{
    public static void DisplayGDSIIFile(VeldridDriver driver, OVPSettings settings, 
        GeoCore geoCore, int cellIndex, int layerFilter = -1)
    {
        var drawing = geoCore.getDrawing();
        if (cellIndex >= drawing.cellList.Count) return;
        
        var cell = drawing.cellList[cellIndex];
        var polyList = new List<ovp_Poly>();
        var ptCounts = new List<int>();
        
        // Color map for different layers
        var layerColors = new Dictionary<int, Color>
        {
            {1, Colors.Red}, {2, Colors.Blue}, {3, Colors.Green},
            {4, Colors.Yellow}, {5, Colors.Purple}, {6, Colors.Orange}
        };
        
        foreach (var element in cell.elementList)
        {
            if (layerFilter != -1 && element.layer_nr != layerFilter) continue;
            
            Color layerColor = layerColors.ContainsKey(element.layer_nr) 
                ? layerColors[element.layer_nr] 
                : Colors.Gray;
            
            switch (element)
            {
                case GCPolygon polygon:
                    var pathD = polygon.pointarray.Select(pt => 
                        new PointD(pt.X * drawing.getDrawingScale(), 
                                  pt.Y * drawing.getDrawingScale())).ToList();
                    
                    polyList.Add(new ovp_Poly(new PathD(pathD), layerColor));
                    ptCounts.Add(pathD.Count);
                    break;
                    
                case GCBox box:
                    var boxPath = new PathD
                    {
                        new PointD(box.x.X * drawing.getDrawingScale(), box.x.Y * drawing.getDrawingScale()),
                        new PointD(box.y.X * drawing.getDrawingScale(), box.x.Y * drawing.getDrawingScale()),
                        new PointD(box.y.X * drawing.getDrawingScale(), box.y.Y * drawing.getDrawingScale()),
                        new PointD(box.x.X * drawing.getDrawingScale(), box.y.Y * drawing.getDrawingScale())
                    };
                    
                    polyList.Add(new ovp_Poly(boxPath, layerColor));
                    ptCounts.Add(4);
                    break;
            }
        }
        
        settings.polyList = polyList;
        settings.polyListPtCount = ptCounts;
        settings.changed = true;
        
        driver.updateViewport();
        driver.zoomExtents();
    }
}
```

## Interactive Features

### Mouse and Keyboard Handling
```csharp
public class ViewportController
{
    private VeldridDriver driver;
    private OVPSettings settings;
    private bool isDragging = false;
    private PointF lastMousePos;
    
    public ViewportController(VeldridDriver driver, OVPSettings settings)
    {
        this.driver = driver;
        this.settings = settings;
        
        // Set up event handlers (if surface provides mouse events)
        // surface.MouseDown += OnMouseDown;
        // surface.MouseMove += OnMouseMove;
        // surface.MouseUp += OnMouseUp;
    }
    
    private void OnMouseDown(object sender, MouseEventArgs e)
    {
        if (e.Buttons == MouseButtons.Primary)
        {
            if (settings.isZoomAndPanAllowed())
            {
                isDragging = true;
                lastMousePos = e.Location;
            }
            else
            {
                // Selection mode
                var worldPos = ScreenToWorld(e.Location);
                driver.setSelection(worldPos);
            }
        }
    }
    
    private void OnMouseMove(object sender, MouseEventArgs e)
    {
        if (isDragging && settings.isZoomAndPanAllowed())
        {
            var deltaX = e.Location.X - lastMousePos.X;
            var deltaY = e.Location.Y - lastMousePos.Y;
            
            var currentPos = settings.getCameraPos();
            var newPos = new PointF(
                currentPos.X - deltaX * settings.getZoomFactor(),
                currentPos.Y + deltaY * settings.getZoomFactor() // Invert Y
            );
            
            settings.setCameraPos(newPos);
            lastMousePos = e.Location;
            driver.updateViewport();
        }
    }
    
    private void OnMouseUp(object sender, MouseEventArgs e)
    {
        isDragging = false;
    }
    
    private PointF ScreenToWorld(PointF screenPos)
    {
        // Convert screen coordinates to world coordinates
        // Implementation depends on viewport transformation matrix
        var cameraPos = settings.getCameraPos();
        var zoomFactor = settings.getZoomFactor();
        
        return new PointF(
            screenPos.X * zoomFactor + cameraPos.X,
            screenPos.Y * zoomFactor + cameraPos.Y
        );
    }
}
```

### Animation and Real-Time Updates
```csharp
public class AnimatedViewport
{
    private VeldridDriver driver;
    private OVPSettings settings;
    private Timer animationTimer;
    private double animationTime = 0.0;
    
    public AnimatedViewport(VeldridDriver driver, OVPSettings settings)
    {
        this.driver = driver;
        this.settings = settings;
        
        // Set up animation timer
        animationTimer = new Timer(UpdateAnimation, null, 0, 16); // ~60 FPS
    }
    
    private void UpdateAnimation(object state)
    {
        animationTime += 0.016; // Add frame time
        
        // Update animated geometry
        if (settings.polyList != null)
        {
            for (int i = 0; i < settings.polyList.Count; i++)
            {
                var poly = settings.polyList[i];
                
                // Animate color or alpha
                float alpha = (float)(0.5 + 0.5 * Math.Sin(animationTime + i));
                poly.alpha = alpha;
                
                // Could also animate geometry positions
            }
            
            settings.changed = true;
            
            // Update viewport on UI thread
            Application.Instance.AsyncInvoke(() => driver.updateViewport());
        }
    }
    
    public void StopAnimation()
    {
        animationTimer?.Dispose();
    }
}
```

## Performance Optimization

### Efficient Geometry Management
```csharp
public static class ViewportOptimization
{
    public static void OptimizeGeometry(OVPSettings settings)
    {
        if (settings.polyList == null) return;
        
        var optimizedPolys = new List<ovp_Poly>();
        var optimizedPtCounts = new List<int>();
        
        for (int i = 0; i < settings.polyList.Count; i++)
        {
            var poly = settings.polyList[i];
            
            // Skip very small polygons that won't be visible
            var bounds = GetPolygonBounds(poly.poly);
            if (bounds.Width < 1.0 && bounds.Height < 1.0) continue;
            
            // Simplify complex polygons
            var simplified = SimplifyPolygon(poly.poly, tolerance: 0.1);
            if (simplified.Count >= 3)
            {
                optimizedPolys.Add(new ovp_Poly(simplified, poly.color, poly.alpha));
                optimizedPtCounts.Add(simplified.Count);
            }
        }
        
        settings.polyList = optimizedPolys;
        settings.polyListPtCount = optimizedPtCounts;
    }
    
    private static RectangleF GetPolygonBounds(PathD polygon)
    {
        if (polygon.Count == 0) return RectangleF.Empty;
        
        double minX = polygon.Min(p => p.x);
        double maxX = polygon.Max(p => p.x);
        double minY = polygon.Min(p => p.y);
        double maxY = polygon.Max(p => p.y);
        
        return new RectangleF((float)minX, (float)minY, 
            (float)(maxX - minX), (float)(maxY - minY));
    }
    
    private static PathD SimplifyPolygon(PathD polygon, double tolerance)
    {
        // Implement Douglas-Peucker or similar algorithm
        // For now, return original
        return polygon;
    }
}
```

### Level-of-Detail Rendering
```csharp
public static class LODRenderer
{
    public static void SetupLevelOfDetail(OVPSettings settings, float currentZoom)
    {
        if (settings.polyList == null) return;
        
        // Adjust rendering based on zoom level
        if (currentZoom < 0.1f)
        {
            // Very zoomed out - show simplified geometry
            settings.filled(false);  // Show outlines only
            settings.aA(false);      // Disable anti-aliasing for performance
        }
        else if (currentZoom < 1.0f)
        {
            // Medium zoom - moderate detail
            settings.filled(true);
            settings.aA(true);
        }
        else
        {
            // Zoomed in - full detail
            settings.filled(true);
            settings.aA(true);
            settings.drawGrid(true); // Show grid when zoomed in
        }
    }
}
```

## Error Handling and Diagnostics

### Viewport Validation
```csharp
public static class ViewportValidation
{
    public static bool ValidateSettings(OVPSettings settings, out string errorMessage)
    {
        errorMessage = string.Empty;
        
        // Check bounds
        if (settings.minX >= settings.maxX || settings.minY >= settings.maxY)
        {
            errorMessage = "Invalid viewport bounds";
            return false;
        }
        
        // Check polygon data consistency
        if (settings.polyList != null && settings.polyListPtCount != null)
        {
            if (settings.polyList.Count != settings.polyListPtCount.Count)
            {
                errorMessage = "Polygon count mismatch";
                return false;
            }
            
            for (int i = 0; i < settings.polyList.Count; i++)
            {
                if (settings.polyList[i].poly.Count != settings.polyListPtCount[i])
                {
                    errorMessage = $"Point count mismatch for polygon {i}";
                    return false;
                }
            }
        }
        
        return true;
    }
    
    public static void LogRenderingInfo(VeldridDriver driver, OVPSettings settings)
    {
        int totalPolygons = settings.polyList?.Count ?? 0;
        int totalPoints = settings.polyListPtCount?.Sum() ?? 0;
        
        Console.WriteLine("Viewport Rendering Info:");
        Console.WriteLine($"  Total polygons: {totalPolygons}");
        Console.WriteLine($"  Total points: {totalPoints}");
        Console.WriteLine($"  Viewport bounds: ({settings.minX}, {settings.minY}) to ({settings.maxX}, {settings.maxY})");
        Console.WriteLine($"  Grid enabled: {settings.drawGrid()}");
        Console.WriteLine($"  Anti-aliasing: {settings.aA()}");
    }
}
```

## Best Practices

### Application Integration
```csharp
public class ViewportForm : Form
{
    private VeldridDriver driver;
    private OVPSettings settings;
    private VeldridSurface surface;
    
    public ViewportForm()
    {
        InitializeComponent();
        SetupViewport();
    }
    
    private void SetupViewport()
    {
        // Initialize viewport components
        settings = new OVPSettings();
        surface = new VeldridSurface();
        driver = new VeldridDriver(ref settings, ref surface);
        
        // Configure initial settings
        ConfigureDefaultSettings();
        
        // Set up the viewport
        driver.SetUpVeldrid();
        
        // Add to form
        Content = surface;
        
        // Handle form events
        Closing += OnFormClosing;
    }
    
    private void ConfigureDefaultSettings()
    {
        settings.backColor = Colors.White;
        settings.minorGridColor = Color.FromArgb(240, 240, 240);
        settings.majorGridColor = Color.FromArgb(200, 200, 200);
        settings.axisColor = Colors.Black;
        settings.selectionColor = Colors.Blue;
        
        settings.drawGrid(true);
        settings.drawAxes(true);
        settings.allowZoomAndPan(true);
        settings.aA(true);
        settings.filled(true);
        
        // Set reasonable default bounds
        settings.minX = -1000;
        settings.maxX = 1000;
        settings.minY = -1000;
        settings.maxY = 1000;
    }
    
    public void DisplayGeometry(PathsD geometry)
    {
        GeometryViewport.DisplayGeometry(driver, settings, geometry);
    }
    
    private void OnFormClosing(object sender, CancelEventArgs e)
    {
        // Clean up resources
        driver?.Dispose();
    }
}

// Usage
var form = new ViewportForm();
form.Show();

// Display some geometry
var geometry = new PathsD { /* your polygon data */ };
form.DisplayGeometry(geometry);
```

## Dependencies
- **Eto.Forms** - Cross-platform UI framework
- **Eto.Veldrid** - Veldrid integration for Eto.Forms
- **Veldrid** - Cross-platform graphics API
- **Clipper2Lib** - Polygon data types
- **LibTessDotNet** - Polygon tessellation for filled rendering
- **.NET 8.0** - Target framework

## Related Libraries
- **[geoWrangler](geoWrangler-API.md)** - Provides geometry that can be displayed in the viewport
- **[geoCore](geoCore-API.md)** - File I/O for geometry that can be visualized
- **[clipper](clipper-API.md)** - Polygon operations whose results can be displayed
- **[errorReporter](errorReporter-API.md)** - Error handling integration for viewport errors