using System;
using System.Runtime.InteropServices;

namespace Eto.Veldrid.Gtk;

internal static class X11Interop
{
	private const string
		linux_libgdk_name = "libgdk-3.so.0",
		linux_libGL_name = "libGL.so.1";

	// GDK backend detection functions
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

	// Additional GDK functions for native window management
	[DllImport(linux_libgdk_name)]
	public static extern void gdk_window_ensure_native(IntPtr gdkWindow);

	// OpenGL function
	[DllImport(linux_libGL_name)]
	public static extern IntPtr glXGetProcAddress(string name);

	/// <summary>
	/// Determines if the current display is X11-based by checking the display name.
	/// </summary>
	public static bool IsX11Display(IntPtr gdkDisplay)
	{
		try
		{
			var displayNamePtr = gdk_display_get_name(gdkDisplay);
			if (displayNamePtr == IntPtr.Zero)
				return false;
				
			var displayName = System.Runtime.InteropServices.Marshal.PtrToStringAnsi(displayNamePtr);
			if (string.IsNullOrEmpty(displayName))
				return false;

			// X11 displays typically have names like ":0", ":0.0", "localhost:10.0", etc.
			// They don't contain "wayland" in the name
			return displayName.Contains(":") && !displayName.ToLower().Contains("wayland");
		}
		catch
		{
			return false;
		}
	}

	/// <summary>
	/// Determines if the current display is Wayland-based by checking the display name.
	/// </summary>
	public static bool IsWaylandDisplay(IntPtr gdkDisplay)
	{
		try
		{
			var displayNamePtr = gdk_display_get_name(gdkDisplay);
			if (displayNamePtr == IntPtr.Zero)
				return false;
				
			var displayName = System.Runtime.InteropServices.Marshal.PtrToStringAnsi(displayNamePtr);
			if (string.IsNullOrEmpty(displayName))
				return false;

			// Wayland displays typically have names like "wayland-0", "wayland-1", etc.
			return displayName.ToLower().Contains("wayland");
		}
		catch
		{
			return false;
		}
	}
}