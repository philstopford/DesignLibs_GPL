using System;
using System.Runtime.InteropServices;

namespace Eto.Veldrid.Gtk;

internal static class X11Interop
{
	private const string
		linux_libgdk_name = "libgdk-3.so.0",
		linux_libGL_name = "libGL.so.1",
		linux_libglib_name = "libglib-2.0.so.0";

	// GDK backend detection functions
	[DllImport(linux_libgdk_name)]
	public static extern IntPtr gdk_display_get_default();

	[DllImport(linux_libgdk_name)]
	public static extern IntPtr gdk_display_get_name(IntPtr display);

	[DllImport(linux_libgdk_name)]
	public static extern IntPtr gdk_display_get_default_display();
	
	// Helper to convert IntPtr to string
	private static string? PtrToString(IntPtr ptr)
	{
		return ptr == IntPtr.Zero ? null : Marshal.PtrToStringAnsi(ptr);
	}

	// X11-specific functions
	[DllImport(linux_libgdk_name)]
	public static extern IntPtr gdk_x11_display_get_xdisplay(IntPtr gdkDisplay);

	[DllImport(linux_libgdk_name)]
	public static extern int gdk_x11_screen_get_screen_number(IntPtr gdkScreen);

	[DllImport(linux_libgdk_name)]
	public static extern IntPtr gdk_x11_window_get_xid(IntPtr gdkWindow);

	// Wayland-specific functions
	[DllImport(linux_libgdk_name)]
	public static extern IntPtr gdk_wayland_display_get_wl_display(IntPtr gdkDisplay);

	[DllImport(linux_libgdk_name)]
	public static extern IntPtr gdk_wayland_window_get_wl_surface(IntPtr gdkWindow);

	// OpenGL function
	[DllImport(linux_libGL_name)]
	public static extern IntPtr glXGetProcAddress(string name);

	/// <summary>
	/// Checks if the current GDK display is running on X11
	/// </summary>
	public static bool IsX11Display(IntPtr gdkDisplay)
	{
		try
		{
			// First try: Check display name for X11 indicators
			var displayNamePtr = gdk_display_get_name(gdkDisplay);
			var displayName = PtrToString(displayNamePtr);
			
			if (!string.IsNullOrEmpty(displayName))
			{
				// X11 displays typically have format like ":0", ":1", "hostname:0", etc.
				if (displayName.Contains(':') && !displayName.StartsWith("wayland"))
				{
					return true;
				}
				// Explicit wayland indicator
				if (displayName.StartsWith("wayland"))
				{
					return false;
				}
			}

			// Second try: Attempt to call X11-specific function as fallback
			// This is risky but we wrap it in try-catch
			var xDisplay = gdk_x11_display_get_xdisplay(gdkDisplay);
			return xDisplay != IntPtr.Zero;
		}
		catch
		{
			return false;
		}
	}

	/// <summary>
	/// Checks if the current GDK display is running on Wayland
	/// </summary>
	public static bool IsWaylandDisplay(IntPtr gdkDisplay)
	{
		try
		{
			// First try: Check display name for Wayland indicators
			var displayNamePtr = gdk_display_get_name(gdkDisplay);
			var displayName = PtrToString(displayNamePtr);
			
			if (!string.IsNullOrEmpty(displayName))
			{
				// Explicit wayland indicator
				if (displayName.StartsWith("wayland"))
				{
					return true;
				}
				// X11 displays typically have format like ":0", ":1", "hostname:0", etc.
				if (displayName.Contains(':') && !displayName.StartsWith("wayland"))
				{
					return false;
				}
			}

			// Second try: Attempt to call Wayland-specific function as fallback
			// This is risky but we wrap it in try-catch
			var waylandDisplay = gdk_wayland_display_get_wl_display(gdkDisplay);
			return waylandDisplay != IntPtr.Zero;
		}
		catch
		{
			return false;
		}
	}
}