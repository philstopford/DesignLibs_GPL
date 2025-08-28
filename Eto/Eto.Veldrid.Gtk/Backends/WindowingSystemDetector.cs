using System;
using System.Runtime.InteropServices;

namespace Eto.Veldrid.Gtk.Backends;

/// <summary>
/// Utility for detecting the current windowing system (X11 or Wayland)
/// </summary>
internal static class WindowingSystemDetector
{
    public enum WindowingSystem
    {
        Unknown,
        X11,
        Wayland
    }

    private static WindowingSystem? _cachedSystem;

    /// <summary>
    /// Detects the current windowing system
    /// </summary>
    public static WindowingSystem DetectWindowingSystem()
    {
        if (_cachedSystem.HasValue)
            return _cachedSystem.Value;

        _cachedSystem = DetectWindowingSystemImpl();
        return _cachedSystem.Value;
    }

    private static WindowingSystem DetectWindowingSystemImpl()
    {
        Console.WriteLine("WindowingSystemDetector: Starting windowing system detection");
        
        // First check environment variables
        string? waylandDisplay = Environment.GetEnvironmentVariable("WAYLAND_DISPLAY");
        string? x11Display = Environment.GetEnvironmentVariable("DISPLAY");
        
        Console.WriteLine($"WindowingSystemDetector: WAYLAND_DISPLAY={waylandDisplay}, DISPLAY={x11Display}");

        // If WAYLAND_DISPLAY is set and DISPLAY is not, we're likely on Wayland
        if (!string.IsNullOrEmpty(waylandDisplay) && string.IsNullOrEmpty(x11Display))
        {
            Console.WriteLine("WindowingSystemDetector: Detected pure Wayland (no X11)");
            return WindowingSystem.Wayland;
        }

        // If DISPLAY is set, we might be on X11 or XWayland
        if (!string.IsNullOrEmpty(x11Display))
        {
            // Check if we're running under XWayland
            if (IsRunningUnderXWayland())
            {
                Console.WriteLine("WindowingSystemDetector: Detected XWayland");
                return WindowingSystem.Wayland;
            }
            Console.WriteLine("WindowingSystemDetector: Detected X11");
            return WindowingSystem.X11;
        }

        // Try to detect based on GDK display type
        try
        {
            var display = Gdk.Display.Default;
            if (display != null)
            {
                string displayTypeName = display.GetType().Name;
                Console.WriteLine($"WindowingSystemDetector: GDK display type: {displayTypeName}");
                
                if (displayTypeName.Contains("Wayland", StringComparison.OrdinalIgnoreCase))
                {
                    Console.WriteLine("WindowingSystemDetector: Detected Wayland via GDK");
                    return WindowingSystem.Wayland;
                }
                else if (displayTypeName.Contains("X11", StringComparison.OrdinalIgnoreCase))
                {
                    Console.WriteLine("WindowingSystemDetector: Detected X11 via GDK");
                    return WindowingSystem.X11;
                }
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine($"WindowingSystemDetector: Error in GDK display detection: {ex.Message}");
        }

        // Fall back to checking if Wayland compositor is available
        if (!string.IsNullOrEmpty(waylandDisplay))
        {
            Console.WriteLine("WindowingSystemDetector: Fallback to Wayland (WAYLAND_DISPLAY present)");
            return WindowingSystem.Wayland;
        }

        // Default to X11 if we can't determine
        Console.WriteLine("WindowingSystemDetector: Fallback to X11 (default)");
        return WindowingSystem.X11;
    }

    private static bool IsRunningUnderXWayland()
    {
        // Check for XWayland-specific environment variables or process indicators
        try
        {
            // XWayland typically sets this environment variable
            string? waylandDisplay = Environment.GetEnvironmentVariable("WAYLAND_DISPLAY");
            
            // If both DISPLAY and WAYLAND_DISPLAY are set, we're likely in XWayland
            if (!string.IsNullOrEmpty(waylandDisplay) && 
                !string.IsNullOrEmpty(Environment.GetEnvironmentVariable("DISPLAY")))
            {
                return true;
            }

            // Check for XWayland process or compositor hints
            string? xdgSessionType = Environment.GetEnvironmentVariable("XDG_SESSION_TYPE");
            if (string.Equals(xdgSessionType, "wayland", StringComparison.OrdinalIgnoreCase))
            {
                return true;
            }
        }
        catch
        {
            // Ignore errors
        }

        return false;
    }
}