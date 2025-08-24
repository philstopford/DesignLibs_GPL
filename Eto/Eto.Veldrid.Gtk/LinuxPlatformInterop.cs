using System;
using System.Runtime.InteropServices;

namespace Eto.Veldrid.Gtk;

/// <summary>
/// Provides interop functions for Linux display servers (X11, Wayland, XWayland)
/// </summary>
public static class LinuxPlatformInterop
{
	private const string
		linux_libgdk_x11_name = "libgdk-3.so.0",
		linux_libgdk_wayland_name = "libgdk-3.so.0",
		linux_libGL_name = "libGL.so.1";

	public enum DisplayServerType
	{
		Unknown,
		X11,
		Wayland,
		XWayland
	}

	// X11 interop functions
	[DllImport(linux_libgdk_x11_name)]
	public static extern IntPtr gdk_x11_display_get_xdisplay(IntPtr gdkDisplay);

	[DllImport(linux_libgdk_x11_name)]
	public static extern int gdk_x11_screen_get_screen_number(IntPtr gdkScreen);

	[DllImport(linux_libgdk_x11_name)]
	public static extern IntPtr gdk_x11_window_get_xid(IntPtr gdkWindow);

	// Wayland interop functions
	[DllImport(linux_libgdk_wayland_name)]
	public static extern IntPtr gdk_wayland_display_get_wl_display(IntPtr gdkDisplay);

	[DllImport(linux_libgdk_wayland_name)]
	public static extern IntPtr gdk_wayland_window_get_wl_surface(IntPtr gdkWindow);

	// OpenGL proc address
	[DllImport(linux_libGL_name)]
	public static extern IntPtr glXGetProcAddress(string name);

	// Platform detection functions
	[DllImport(linux_libgdk_x11_name)]
	public static extern bool gdk_display_is_x11(IntPtr display);

	[DllImport(linux_libgdk_wayland_name)]
	public static extern bool gdk_display_is_wayland(IntPtr display);

	// Note: gdk_window_is_x11 and gdk_window_is_wayland were added in GDK 4
	// For GDK 3, we use display-level detection and try-catch approaches

	/// <summary>
	/// Determines the type of display server currently being used
	/// </summary>
	/// <param name="gdkDisplay">The GDK display handle</param>
	/// <param name="gdkWindow">The GDK window handle</param>
	/// <returns>The type of display server</returns>
	public static DisplayServerType GetDisplayServerType(IntPtr gdkDisplay, IntPtr gdkWindow)
	{
		if (gdkDisplay == IntPtr.Zero)
		{
			Console.WriteLine("[LinuxPlatformInterop] GDK display handle is null");
			return DisplayServerType.Unknown;
		}

		if (gdkWindow == IntPtr.Zero)
		{
			Console.WriteLine("[LinuxPlatformInterop] GDK window handle is null");
			return DisplayServerType.Unknown;
		}

		try
		{
			// Use display-level detection since window-level functions don't exist in GDK 3
			bool displayIsWayland = gdk_display_is_wayland(gdkDisplay);
			Console.WriteLine($"[LinuxPlatformInterop] gdk_display_is_wayland(): {displayIsWayland}");
			
			if (displayIsWayland)
			{
				Console.WriteLine("[LinuxPlatformInterop] Detected Wayland display");
				return DisplayServerType.Wayland;
			}

			bool displayIsX11 = gdk_display_is_x11(gdkDisplay);
			Console.WriteLine($"[LinuxPlatformInterop] gdk_display_is_x11(): {displayIsX11}");
			
			if (displayIsX11)
			{
				// Check if we're running under XWayland by looking at the session type
				string? xdgSessionType = Environment.GetEnvironmentVariable("XDG_SESSION_TYPE");
				Console.WriteLine($"[LinuxPlatformInterop] XDG_SESSION_TYPE: {xdgSessionType}");
				
				if (xdgSessionType == "wayland")
				{
					Console.WriteLine("[LinuxPlatformInterop] Detected XWayland (X11 display in Wayland session)");
					return DisplayServerType.XWayland;
				}
				Console.WriteLine("[LinuxPlatformInterop] Detected native X11 display");
				return DisplayServerType.X11;
			}
		}
		catch (DllNotFoundException ex)
		{
			Console.WriteLine($"[LinuxPlatformInterop] DLL not found, falling back to environment detection: {ex.Message}");
			// If we can't load the libraries, fall back to environment detection
			return GetDisplayServerTypeFromEnvironment();
		}
		catch (EntryPointNotFoundException ex)
		{
			Console.WriteLine($"[LinuxPlatformInterop] Entry point not found, falling back to environment detection: {ex.Message}");
			// If the functions don't exist, fall back to environment detection
			return GetDisplayServerTypeFromEnvironment();
		}
		catch (Exception ex)
		{
			Console.WriteLine($"[LinuxPlatformInterop] Error during detection, falling back to environment detection: {ex.Message}");
			return GetDisplayServerTypeFromEnvironment();
		}

		Console.WriteLine("[LinuxPlatformInterop] Unknown display server type");
		return DisplayServerType.Unknown;
	}

	/// <summary>
	/// Determines the type of display server currently being used (legacy method)
	/// </summary>
	/// <param name="gdkDisplay">The GDK display handle</param>
	/// <returns>The type of display server</returns>
	public static DisplayServerType GetDisplayServerType(IntPtr gdkDisplay)
	{
		// Legacy method - falls back to environment detection for safety
		Console.WriteLine("[LinuxPlatformInterop] Using legacy detection method, falling back to environment detection");
		return GetDisplayServerTypeFromEnvironment();
	}

	/// <summary>
	/// Fallback method to detect display server type from environment variables
	/// This is used when GDK functions are not available
	/// </summary>
	/// <returns>The detected display server type</returns>
	public static DisplayServerType GetDisplayServerTypeFromEnvironment()
	{
		// Check environment variables to determine the session type
		string? waylandDisplay = Environment.GetEnvironmentVariable("WAYLAND_DISPLAY");
		string? xdgSessionType = Environment.GetEnvironmentVariable("XDG_SESSION_TYPE");
		string? display = Environment.GetEnvironmentVariable("DISPLAY");

		// If the session type is explicitly set, use that
		if (xdgSessionType == "wayland")
		{
			// In a Wayland session, if we have DISPLAY, it could be XWayland
			// But since we can't use GDK functions here, assume native Wayland
			return DisplayServerType.Wayland;
		}

		if (xdgSessionType == "x11")
		{
			return DisplayServerType.X11;
		}

		// Fallback to checking display variables
		if (!string.IsNullOrEmpty(waylandDisplay))
		{
			return DisplayServerType.Wayland;
		}

		if (!string.IsNullOrEmpty(display))
		{
			return DisplayServerType.X11;
		}

		return DisplayServerType.Unknown;
	}

	/// <summary>
	/// Safely checks if a window supports Wayland operations by attempting to use Wayland functions
	/// </summary>
	/// <param name="display">The GDK display handle</param>
	/// <param name="window">The GDK window handle</param>
	/// <returns>True if the window supports Wayland operations, false otherwise</returns>
	public static bool SafelyCheckWaylandWindow(IntPtr display, IntPtr window)
	{
		try
		{
			// First check if the display is Wayland
			if (!gdk_display_is_wayland(display))
			{
				Console.WriteLine("[LinuxPlatformInterop] Display is not Wayland");
				return false;
			}

			// Try to get Wayland surface - if this works, the window supports Wayland
			IntPtr wlSurface = gdk_wayland_window_get_wl_surface(window);
			bool result = wlSurface != IntPtr.Zero;
			Console.WriteLine($"[LinuxPlatformInterop] Wayland surface check result: {result}");
			return result;
		}
		catch (Exception ex)
		{
			Console.WriteLine($"[LinuxPlatformInterop] Failed to check Wayland window support: {ex.Message}");
			return false;
		}
	}

	/// <summary>
	/// Safely checks if a window supports X11 operations by attempting to use X11 functions
	/// </summary>
	/// <param name="display">The GDK display handle</param>
	/// <param name="window">The GDK window handle</param>
	/// <returns>True if the window supports X11 operations, false otherwise</returns>
	public static bool SafelyCheckX11Window(IntPtr display, IntPtr window)
	{
		try
		{
			// First check if the display is X11
			if (!gdk_display_is_x11(display))
			{
				Console.WriteLine("[LinuxPlatformInterop] Display is not X11");
				return false;
			}

			// Try to get X11 window ID - if this works, the window supports X11
			IntPtr xid = gdk_x11_window_get_xid(window);
			bool result = xid != IntPtr.Zero;
			Console.WriteLine($"[LinuxPlatformInterop] X11 window ID check result: {result}");
			return result;
		}
		catch (Exception ex)
		{
			Console.WriteLine($"[LinuxPlatformInterop] Failed to check X11 window support: {ex.Message}");
			return false;
		}
	}
}