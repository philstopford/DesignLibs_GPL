# Viewport Performance Diagnostics Guide

This guide explains how to use the comprehensive diagnostic logging added to help troubleshoot viewport update delays on Linux/GTK.

## Overview

The viewport rendering system now includes detailed timing and diagnostic information at multiple levels:
- Overall viewport update timing
- Individual drawing task performance (axes, grid, lines, polygons)
- GTK-specific rendering pipeline events
- Platform detection and state information

## Enabling Diagnostics

There are two ways to enable diagnostic output:

### Method 1: Environment Variable (Recommended for client tools)

Set the `VIEWPORT_DIAGNOSTICS` environment variable before running your application:

```bash
# On Linux/GTK
export VIEWPORT_DIAGNOSTICS=1
./YourApplication

# Or run directly
VIEWPORT_DIAGNOSTICS=1 ./YourApplication
```

```powershell
# On Windows
$env:VIEWPORT_DIAGNOSTICS="1"
.\YourApplication.exe
```

This method automatically enables diagnostics for both the VeldridDriver and GTK-specific handlers.

### Method 2: Programmatic Control

For applications that embed the viewport, you can enable diagnostics programmatically:

```csharp
// After creating the VeldridDriver
driver.EnableDiagnostics = true;
```

This method only affects the VeldridDriver diagnostics, not GTK-specific output. Use environment variable for complete diagnostics.

## Understanding the Output

### Viewport Update Flow

When diagnostics are enabled, you'll see output like this:

```
[VIEWPORT DIAG] updateViewport() called on GTK platform
[VIEWPORT DIAG] ovpSettings.changed=True, drawing=False, done_drawing=False
[VIEWPORT DIAG] Starting geometry processing...
[VIEWPORT DIAG]   Polygons: fg=100, bg=50, tess=200
[VIEWPORT DIAG]   Lines: 75
[VIEWPORT DIAG]   drawAxes() took 2ms
[VIEWPORT DIAG]   drawGrid() took 5ms
[VIEWPORT DIAG]   drawLines() took 12ms
[VIEWPORT DIAG]   drawPolygons() took 45ms
[VIEWPORT DIAG] All drawing tasks completed in 47ms
[VIEWPORT DIAG] pUpdateViewportAsync completed in 48ms
[VIEWPORT DIAG] After async: done_drawing=True
[VIEWPORT DIAG] updateHostFunc callback took 1ms
[VIEWPORT DIAG] Surface.Invalidate() took 0ms
[VIEWPORT DIAG] Total updateViewport() time: 50ms
[VIEWPORT DIAG] ----------------------------------------
```

### GTK-Specific Events

On GTK, you'll also see platform-specific rendering events:

```
[GTK DIAG] Invalidate() called, skipDraw=False
[GTK DIAG] glArea_Render called, skipDraw=False
[GTK DIAG] Calling OnDraw callback
[GTK DIAG] OnDraw callback completed in 25ms
[GTK DIAG] glArea_Render completed
```

## Key Metrics to Monitor

### 1. Total Update Time
- **Normal**: < 100ms for typical geometry
- **Concerning**: > 500ms 
- **Critical**: > 1000ms (1 second)

Look for: `Total updateViewport() time: XXms`

### 2. Drawing Task Distribution
Individual task times help identify bottlenecks:
- `drawPolygons()` - Usually the longest (polygon tessellation and vertex buffer creation)
- `drawLines()` - Medium duration
- `drawGrid()` and `drawAxes()` - Usually fastest

### 3. Platform-Specific Delays

#### On GTK:
- Check time between `Invalidate()` and `glArea_Render`
- Large gaps (> 100ms) indicate GTK event loop delays
- Monitor `OnDraw callback completed` time

#### On Windows/WPF:
- Similar but without GTK-specific messages
- Generally faster due to better thread affinity

### 4. State Issues

Watch for warning signs:
- `ovpSettings.changed=False` when it should be True
- `drawing=True` when calling updateViewport (indicates concurrent calls)
- `done_drawing=False` after async completion (shouldn't happen)
- Early returns from pUpdateViewportAsync (check the reason)

## Common Performance Issues

### Issue 1: Long Polygon Processing Time
**Symptom**: `drawPolygons() took XXXms` where XXX > 100

**Possible Causes**:
- Too many polygons
- Complex tessellation
- Inefficient parallel processing

**Debug**: Check polygon counts in output:
```
[VIEWPORT DIAG]   Polygons: fg=10000, bg=5000, tess=20000
```

### Issue 2: GTK Event Loop Delays
**Symptom**: Long gap between `Invalidate()` and `glArea_Render`

**Possible Causes**:
- GTK main thread busy with other tasks
- Heavy UI updates blocking event processing
- Thread affinity issues

**Debug**: Look for patterns in timing between these events

### Issue 3: Callback Overhead
**Symptom**: `updateHostFunc callback took XXXms` where XXX > 50

**Possible Causes**:
- Client application doing heavy work in callback
- Blocking operations in update handler
- Recursive viewport updates

**Debug**: Review what the client's updateHostFunc does

## Performance Comparison: Windows vs Linux/GTK

When comparing performance between platforms, look for:

1. **Total time differences**:
   - Windows: Typically 20-50ms
   - Linux/GTK: Should be similar (30-60ms)
   - Differences > 100ms indicate platform-specific issues

2. **Event loop responsiveness**:
   - Windows: Immediate invalidation → render
   - GTK: May have slight delay (10-30ms normal)
   - Delays > 100ms on GTK indicate event loop issues

3. **Drawing task performance**:
   - Should be similar across platforms
   - Significant differences indicate thread/parallelization issues

## Reporting Performance Issues

When reporting performance problems, include:

1. **Full diagnostic output** from both platforms if comparing
2. **Geometry complexity**: Number of polygons, lines, points
3. **Platform information**: 
   - OS version
   - GTK version (on Linux)
   - .NET version
4. **Timing patterns**: 
   - Average times
   - Worst-case times
   - Any patterns (e.g., first update slow, subsequent fast)

## Example: Debugging Slow Updates on GTK

1. Enable diagnostics:
```bash
export VIEWPORT_DIAGNOSTICS=1
./quilt
```

2. Perform the slow operation and capture output

3. Look for:
   - Which phase is slow? (geometry processing vs rendering vs invalidation)
   - Are there GTK event loop delays?
   - How does it compare to Windows?

4. Common fixes:
   - If polygon processing is slow: Reduce geometry complexity or optimize tessellation
   - If GTK event delays: Check for UI thread blocking
   - If invalidation is slow: May be driver/compositor issue

## Disabling Diagnostics

To disable diagnostics:

```bash
# Remove or unset the environment variable
unset VIEWPORT_DIAGNOSTICS

# Or set to 0
export VIEWPORT_DIAGNOSTICS=0
```

Or programmatically:
```csharp
driver.EnableDiagnostics = false;
```

## Performance Tips

1. **Reduce geometry updates**: Only call `updateViewport()` when geometry actually changes
2. **Batch updates**: Group multiple geometry changes before updating
3. **Monitor callback performance**: Keep `updateHostFunc` lightweight
4. **Check platform-specific code**: GTK may need different optimization strategies than WPF

## Technical Notes

- Diagnostics use `System.Diagnostics.Stopwatch` for high-precision timing
- Output goes to `Console.WriteLine` (stderr by default)
- Minimal performance impact when disabled
- Thread-safe (diagnostics from parallel tasks are serialized)
- Environment variable checked on each access (can enable/disable without restart)
