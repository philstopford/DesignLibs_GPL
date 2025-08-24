using System;
using System.Runtime.InteropServices;

namespace Eto.Veldrid.Gtk;

internal static class X11Interop
{
	private const string
		linux_libgdk_name = "libgdk-3.so.0",
		linux_libGL_name = "libGL.so.1";

	// Common GDK functions
	[DllImport(linux_libgdk_name)]
	public static extern IntPtr gdk_display_get_default();

	[DllImport(linux_libgdk_name)]
	public static extern IntPtr gdk_display_get_name(IntPtr display);

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
			// Try to get X11 display handle - if this succeeds, we're on X11
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
			// Try to get Wayland display handle - if this succeeds, we're on Wayland
			var waylandDisplay = gdk_wayland_display_get_wl_display(gdkDisplay);
			return waylandDisplay != IntPtr.Zero;
		}
		catch
		{
			return false;
		}
	}
}