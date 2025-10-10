# Viewport Performance Fix for Linux/GTK

## Problem Statement
The viewport on Linux/GTK was 10-20 seconds slower to update with geometry compared to Windows for the same data.

## Root Cause Analysis

### Issue 1: Async/Await Race Condition (CRITICAL)
The `pUpdateViewport()` method was declared as `async void`, which caused a critical race condition:

```csharp
// BEFORE (BROKEN):
private async void pUpdateViewport() {
    // ... async work ...
    done_drawing = true;
}

public void updateViewport() {
    pUpdateViewport();  // Returns immediately, doesn't wait!
    if (done_drawing) {  // Always false here - async work not done yet
        Surface!.Invalidate();  // NEVER CALLED
    }
}
```

**Impact**: The viewport would never invalidate after geometry updates because `done_drawing` was checked before the async work completed. Only timer-based updates worked, causing the perceived 10-20 second delays.

### Issue 2: Inefficient Nested Parallelization
Tessellated polygons used a nested `Parallel.For` loop:

```csharp
// BEFORE (INEFFICIENT):
Parallel.For(0, tessPolyListCount, poly => {
    Parallel.For(0, 3, pt => {  // Parallelizing 3 iterations!
        // Process triangle point
    });
});
```

**Impact**: The overhead of spawning threads for just 3 iterations far exceeded any benefit, especially on GTK/Linux where thread scheduling has more overhead.

## Solution

### Fix 1: Synchronous Completion Pattern
The initial async void approach still had issues because client tools couldn't wait for completion. Changed to use `.GetAwaiter().GetResult()` to force synchronous completion:

```csharp
// AFTER (FIXED):
private async Task pUpdateViewportAsync() {
    // ... async work ...
    done_drawing = true;
}

public void updateViewport() {
    pUpdateViewportAsync().GetAwaiter().GetResult();  // Blocks until complete!
    if (done_drawing) {
        updateHostFunc?.Invoke();
        Surface!.Invalidate();
        ovpSettings.changed = false;
        drawing = false;
        done_drawing = false;
    }
}
```

This ensures that when client tools call `updateViewport()`, the method doesn't return until all geometry processing is complete and the viewport is invalidated.

### Fix 2: Sequential Inner Loop
```csharp
// AFTER (OPTIMIZED):
Parallel.For(0, tessPolyListCount, poly => {
    for (int pt = 0; pt < 3; pt++) {  // Sequential for small iterations
        // Process triangle point
    }
});
```

## Files Modified
- `Eto/Eto.VeldridSurface/VeldridDriver_Draw.cs` - Fixed async pattern and parallelization
- `Eto/Eto.VeldridSurface/VeldridDriver_Public.cs` - Updated `updateViewport()` to await properly
- `Eto/Eto.VeldridSurface/VeldridDriver_Handlers.cs` - Updated `Clock_Elapsed()` to await properly

## Testing
- All 437 unit tests pass
- No regressions introduced
- GTK test application builds successfully

## Expected Performance Improvement
The viewport should now update immediately after client tools call `updateViewport()` on Linux/GTK, matching Windows performance. The method now blocks until the update is complete, ensuring proper sequencing and eliminating the 10-20 second delays caused by the async void pattern.
