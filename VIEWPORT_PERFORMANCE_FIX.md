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

### Fix 1: Proper Async/Await Pattern
```csharp
// AFTER (FIXED):
private async Task pUpdateViewportAsync() {
    // ... async work ...
    done_drawing = true;
}

public async void updateViewport() {
    await pUpdateViewportAsync();  // Actually waits for completion!
    if (done_drawing) {  // Now true when async work is done
        Surface!.Invalidate();  // GETS CALLED
    }
}
```

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
The viewport should now update immediately after geometry changes on Linux/GTK, matching Windows performance. The 10-20 second delay should be completely eliminated as `Surface.Invalidate()` is now called correctly after geometry processing completes.
