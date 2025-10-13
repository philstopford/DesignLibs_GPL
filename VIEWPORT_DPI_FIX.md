# Viewport High DPI Scaling Fix

## Problem Statement
On scaled high DPI displays, the viewport geometry was not fitted correctly and the view was not centered. Objects would appear smaller than expected, and zoom-to-extents would not properly fit the geometry to the viewport.

## Root Cause
The viewport code was inconsistently using physical pixel dimensions (`RenderWidth`/`RenderHeight`) and logical dimensions (`Width`/`Height`) throughout the rendering pipeline.

On high DPI displays:
- **Logical dimensions** (`Width`/`Height`): The control size in logical pixels (e.g., 800x600)
- **Physical dimensions** (`RenderWidth`/`RenderHeight`): The actual framebuffer size in physical pixels, scaled by DPI (e.g., 1600x1200 on a 2x DPI display)

The issue was that some parts of the code used physical dimensions while others used logical dimensions, causing:
1. The orthographic projection to be calculated based on physical pixels, making the view frustum too large in world space
2. The zoom calculations to divide world space extents by physical pixels, making objects appear smaller
3. Mouse coordinate conversions to scale by DPI when they should use logical coordinates

## Solution
Changed the viewport to consistently use **logical dimensions** throughout:

### Files Modified

#### 1. `Eto/Eto.VeldridSurface/VeldridDriver_Draw.cs`
- **Draw() method**: Changed orthographic projection to use `Surface.Width`/`Surface.Height` instead of `Surface.RenderWidth`/`Surface.RenderHeight`
- **drawGrid() method**: Changed grid extent calculations to use logical dimensions
- **drawAxes() method**: Changed axis extent calculations to use logical dimensions

#### 2. `Eto/Eto.VeldridSurface/VeldridDriver_Public.cs`
- **zoomExtents() method**: Changed zoom level calculations to use `Surface.Width`/`Surface.Height` instead of `Surface.RenderWidth`/`Surface.RenderHeight`

#### 3. `Eto/Eto.VeldridSurface/VeldridDriver_Conversions.cs`
- **ScreenToWorld() method**: Changed coordinate conversions to use logical dimensions
- **WorldToScreen() method**: Changed coordinate conversions to use logical dimensions

#### 4. `Eto/Eto.VeldridSurface/VeldridDriver_Handlers.cs`
- **dragHandler() method**: Removed scaling of mouse coordinates by `LogicalPixelSize`
- **selectByClick() method**: Removed scaling of mouse coordinates by `LogicalPixelSize`

## Technical Details

### Before (Incorrect)
```csharp
// Orthographic projection using physical pixels
float left = ovpSettings.getCameraX() - (float)Surface!.RenderWidth / 2 * zoom;
float right = ovpSettings.getCameraX() + (float)Surface!.RenderWidth / 2 * zoom;

// Zoom calculation using physical pixels
float zoomLevel_x = dX / Surface!.RenderWidth;

// Mouse coordinates scaled by DPI
PointF scaledLocation = e.Location * Surface!.ParentWindow.LogicalPixelSize;
```

### After (Correct)
```csharp
// Orthographic projection using logical dimensions
float left = ovpSettings.getCameraX() - (float)Surface!.Width / 2 * zoom;
float right = ovpSettings.getCameraX() + (float)Surface!.Width / 2 * zoom;

// Zoom calculation using logical dimensions
float zoomLevel_x = dX / Surface!.Width;

// Mouse coordinates in logical pixels (no scaling needed)
PointF scaledLocation = e.Location;
```

## Why This Works

1. **Orthographic Projection**: By using logical dimensions, the view frustum size in world space remains consistent regardless of DPI scaling. The GPU automatically handles the mapping to the physical framebuffer.

2. **Zoom Calculations**: Using logical dimensions ensures that the same world space extent results in the same zoom level, regardless of DPI.

3. **Coordinate Conversions**: Mouse events report coordinates in logical pixels, and by using logical dimensions throughout, we maintain consistency between input coordinates and rendered coordinates.

4. **Framebuffer**: The framebuffer size (`RenderWidth`/`RenderHeight`) is still correctly set by the platform handlers (GTK/WPF) to match the physical pixel count, ensuring sharp rendering on high DPI displays.

## Testing
- All 437 existing unit tests pass
- Build completes without errors
- Code is ready for manual testing on different DPI settings (1x, 1.5x, 2x, etc.)

## Expected Behavior After Fix
- Geometry is properly fitted to the viewport on zoom-to-extents
- View is correctly centered on the geometry
- Zoom levels are consistent across different DPI settings
- Mouse panning and selection work correctly on high DPI displays
- Grid and axes render at appropriate scales
