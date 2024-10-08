﻿using System;
using System.Runtime.InteropServices;

namespace Eto.Veldrid.Gtk;

internal static class X11Interop
{
	private const string
		linux_libgdk_x11_name = "libgdk-3.so.0",
		linux_libGL_name = "libGL.so.1";

	[DllImport(linux_libgdk_x11_name)]
	public static extern IntPtr gdk_x11_display_get_xdisplay(IntPtr gdkDisplay);

	[DllImport(linux_libgdk_x11_name)]
	public static extern int gdk_x11_screen_get_screen_number(IntPtr gdkScreen);

	[DllImport(linux_libgdk_x11_name)]
	public static extern IntPtr gdk_x11_window_get_xid(IntPtr gdkDisplay);

	[DllImport(linux_libGL_name)]
	public static extern IntPtr glXGetProcAddress(string name);
}