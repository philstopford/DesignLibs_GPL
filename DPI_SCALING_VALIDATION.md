# DPI Scaling Logic Validation

## Scenario 1: Normal DPI (1x)
- Window logical size: 800x600
- RenderWidth: 800, RenderHeight: 600 (same as logical)
- LogicalPixelSize: 1.0

### Before Fix (Incorrect):
- Orthographic projection used RenderWidth/RenderHeight (800x600)
- Mouse coords scaled by LogicalPixelSize (1.0) - no effect
- Result: Correct behavior by accident

### After Fix (Correct):
- Orthographic projection uses Width/Height (800x600)
- Mouse coords NOT scaled
- Result: Correct behavior by design

## Scenario 2: High DPI (2x)
- Window logical size: 800x600
- RenderWidth: 1600, RenderHeight: 1200 (2x scaled)
- LogicalPixelSize: 2.0

### Before Fix (Incorrect):
- Orthographic projection used RenderWidth/RenderHeight (1600x1200)
  - View frustum: 1600 * zoom in world units wide
  - Objects appear 2x smaller than intended
- Mouse coords at (400, 300) logical pixels scaled to (800, 600)
  - ScreenToWorld used RenderWidth/RenderHeight (1600x1200)
  - Incorrect world coordinate calculated
- ZoomExtents divides world extent by 1600, making objects appear smaller
- Result: Geometry too small, not centered

### After Fix (Correct):
- Orthographic projection uses Width/Height (800x600)
  - View frustum: 800 * zoom in world units wide
  - Objects appear correct size
- Mouse coords at (400, 300) logical pixels NOT scaled
  - ScreenToWorld uses Width/Height (800x600)
  - Correct world coordinate calculated
- ZoomExtents divides world extent by 800, correct zoom level
- Result: Geometry properly sized and centered

## Mathematical Example

World space geometry bounds: 0 to 1000 units in X

### Normal DPI (1x):
- Before: zoomLevel = 1000 / 800 = 1.25
- After: zoomLevel = 1000 / 800 = 1.25
- Same result ✓

### High DPI (2x):
- Before: zoomLevel = 1000 / 1600 = 0.625 (too zoomed out!)
- After: zoomLevel = 1000 / 800 = 1.25 (correct!)
- Fixed ✓

## Key Insight

The framebuffer size (RenderWidth/RenderHeight) should only be used when:
1. Creating the swapchain (already done correctly by platform handlers)
2. Allocating GPU resources

For all viewport calculations (projection, zoom, coordinate conversion), we should use the logical dimensions (Width/Height) to maintain consistent behavior across different DPI settings.
