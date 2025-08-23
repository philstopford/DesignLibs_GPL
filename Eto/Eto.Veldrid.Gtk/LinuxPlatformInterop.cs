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

	/// <summary>
	/// Determines the type of display server currently being used
	/// </summary>
	/// <param name="gdkDisplay">The GDK display handle</param>
	/// <returns>The type of display server</returns>
	public static DisplayServerType GetDisplayServerType(IntPtr gdkDisplay)
	{
		if (gdkDisplay == IntPtr.Zero)
			return DisplayServerType.Unknown;

		try
		{
			// Check if it's a Wayland display
			if (gdk_display_is_wayland(gdkDisplay))
			{
				// Check if we're actually running under XWayland
				// XWayland sets the DISPLAY environment variable
				string? display = Environment.GetEnvironmentVariable("DISPLAY");
				if (!string.IsNullOrEmpty(display))
				{
					return DisplayServerType.XWayland;
				}
				return DisplayServerType.Wayland;
			}

			// Check if it's an X11 display
			if (gdk_display_is_x11(gdkDisplay))
			{
				return DisplayServerType.X11;
			}
		}
		catch (DllNotFoundException)
		{
			// If we can't load the libraries, fall back to environment detection
			return GetDisplayServerTypeFromEnvironment();
		}
		catch (EntryPointNotFoundException)
		{
			// If the functions don't exist, fall back to environment detection
			return GetDisplayServerTypeFromEnvironment();
		}

		return DisplayServerType.Unknown;
	}

	/// <summary>
	/// Fallback method to detect display server type from environment variables
	/// </summary>
	/// <returns>The detected display server type</returns>
	public static DisplayServerType GetDisplayServerTypeFromEnvironment()
	{
		// Check for Wayland session
		string? waylandDisplay = Environment.GetEnvironmentVariable("WAYLAND_DISPLAY");
		string? xdgSessionType = Environment.GetEnvironmentVariable("XDG_SESSION_TYPE");
		string? display = Environment.GetEnvironmentVariable("DISPLAY");

		if (!string.IsNullOrEmpty(waylandDisplay) || xdgSessionType == "wayland")
		{
			// If we have both WAYLAND_DISPLAY and DISPLAY, we're likely in XWayland
			if (!string.IsNullOrEmpty(display))
			{
				return DisplayServerType.XWayland;
			}
			return DisplayServerType.Wayland;
		}

		if (!string.IsNullOrEmpty(display) || xdgSessionType == "x11")
		{
			return DisplayServerType.X11;
		}

		return DisplayServerType.Unknown;
	}
}