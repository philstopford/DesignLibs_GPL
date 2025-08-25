using Eto.Drawing;
using Eto.GtkSharp.Forms;
using Eto.Veldrid;
using Eto.Veldrid.Gtk;
using Gtk;
using System;
using Veldrid;

[assembly: Eto.ExportHandler(typeof(VeldridSurface), typeof(GtkVeldridSurfaceHandler))]

namespace Eto.Veldrid.Gtk;

public class GtkVeldridSurfaceHandler : GtkControl<global::Gtk.Widget, VeldridSurface, VeldridSurface.ICallback>, VeldridSurface.IHandler, VeldridSurface.IOpenGL
{
	private GLArea? glArea;
	private System.Action _makeCurrent;
	private System.Action _clearCurrent;
	public Size RenderSize => Size.Round((SizeF)Widget.Size * Scale);

	private float Scale => Widget.ParentWindow?.Screen?.LogicalPixelSize ?? 1;

	public override global::Gtk.Widget ContainerContentControl => glArea ?? base.ContainerContentControl;

	public GtkVeldridSurfaceHandler()
	{
		_makeCurrent = MakeCurrent;
		_clearCurrent = ClearCurrent;
	}

	public Swapchain? CreateSwapchain()
	{
		if (Widget.Backend == GraphicsBackend.OpenGL)
		{
			return Widget.GraphicsDevice?.MainSwapchain;
		}

		// For Vulkan backend, create swapchain immediately using simplified approach
		return CreateSwapchainNow();
	}

	private Swapchain? CreateSwapchainNow()
	{
		// To embed Veldrid in an Eto control, these platform-specific
		// versions of CreateSwapchain use the technique outlined here:
		//
		//   https://github.com/mellinoe/veldrid/issues/155
		//
		SwapchainSource source;

		// Detect whether we're running on X11 or Wayland and create appropriate SwapchainSource
		var gdkDisplay = Control.Display.Handle;

		bool isX11 = X11Interop.IsX11Display(gdkDisplay);
		bool isWayland = X11Interop.IsWaylandDisplay(gdkDisplay);

		// Add debugging information for diagnostic purposes
		var displayNamePtr = X11Interop.gdk_display_get_name(gdkDisplay);
		var displayName = displayNamePtr == IntPtr.Zero ? "unknown" : System.Runtime.InteropServices.Marshal.PtrToStringAnsi(displayNamePtr) ?? "unknown";
		Console.WriteLine($"[DEBUG] Display name: {displayName}, IsX11: {isX11}, IsWayland: {isWayland}");

		if (isX11 && !isWayland)
		{
			// X11 path - use Xlib SwapchainSource
			Console.WriteLine("[DEBUG] Using X11/Xlib SwapchainSource");
			source = SwapchainSource.CreateXlib(
				X11Interop.gdk_x11_display_get_xdisplay(gdkDisplay),
				X11Interop.gdk_x11_window_get_xid(Control.Window.Handle));
		}
		else if (isWayland && !isX11)
		{
			// Wayland path - use Wayland SwapchainSource with simplified approach
			Console.WriteLine("[DEBUG] Using Wayland SwapchainSource");
			source = CreateWaylandSwapchainSource(gdkDisplay);
		}
		else
		{
			throw new NotSupportedException($"Unsupported or ambiguous windowing system detected. Display: {displayName}, IsX11: {isX11}, IsWayland: {isWayland}. Only X11 and Wayland are supported for Vulkan backend.");
		}

		Size renderSize = RenderSize;
		var swapchain = Widget.GraphicsDevice?.ResourceFactory.CreateSwapchain(
			new SwapchainDescription(
				source,
				(uint)renderSize.Width,
				(uint)renderSize.Height,
				Widget.GraphicsDeviceOptions.SwapchainDepthFormat,
				Widget.GraphicsDeviceOptions.SyncToVerticalBlank,
				Widget.GraphicsDeviceOptions.SwapchainSrgbFormat));

		return swapchain;
	}

	private SwapchainSource CreateWaylandSwapchainSource(IntPtr gdkDisplay)
	{
		Console.WriteLine("[DEBUG] Creating Wayland SwapchainSource");
		
		// Ensure widget is properly realized before accessing Wayland handles
		if (!Control.IsRealized)
		{
			throw new InvalidOperationException("Widget must be realized before creating Wayland surface");
		}
		
		// Basic validation that we have a valid window
		if (Control.Window == null)
		{
			throw new InvalidOperationException("Widget window is null");
		}
		
		Console.WriteLine($"[DEBUG] Window state - Visible: {Control.Window.IsVisible}, Realized: {Control.IsRealized}, Handle: {Control.Window.Handle}");
		
		// Get Wayland handles
		var waylandDisplay = X11Interop.gdk_wayland_display_get_wl_display(gdkDisplay);
		var waylandSurface = X11Interop.gdk_wayland_window_get_wl_surface(Control.Window.Handle);
		
		if (waylandDisplay == IntPtr.Zero || waylandSurface == IntPtr.Zero)
		{
			throw new InvalidOperationException($"Failed to get Wayland handles - Display: {waylandDisplay}, Surface: {waylandSurface}");
		}
		
		Console.WriteLine($"[DEBUG] Successfully obtained Wayland handles - Display: {waylandDisplay}, Surface: {waylandSurface}");
		
		return SwapchainSource.CreateWayland(waylandDisplay, waylandSurface);
	}

	private void glArea_InitializeGraphicsBackend(object? sender, EventArgs e)
	{
		if (glArea == null)
			return;
		// Make context current to manually initialize a Veldrid GraphicsDevice.
		glArea.Context.MakeCurrent();
		Callback.OnInitializeBackend(Widget, new InitializeEventArgs(RenderSize));
		// Veldrid clears the context at the end of initialization and sets it current in the worker thread.

		// Clear context in the worker thread for now to make Mesa happy.
		if (Widget.GraphicsDevice?.GetOpenGLInfo(out BackendInfoOpenGL glInfo) == true)
		{
			// This action has to wait so GTK can manage the context after this method.
			glInfo.ExecuteOnGLThread(_clearCurrent);//, wait: true);
		}

		glArea.Render += glArea_Render;
		glArea.Resize += glArea_Resize;
	}

	private void Control_InitializeGraphicsBackend(object? sender, EventArgs e)
	{
		// Simple immediate backend initialization for both X11 and Wayland
		Console.WriteLine("[DEBUG] Initializing graphics backend");
		Callback.OnInitializeBackend(Widget, new InitializeEventArgs(RenderSize));
	}

	private bool skipDraw;

	private void glArea_Resize(object o, ResizeArgs args)
	{
		skipDraw = false;
		Callback.OnResize(Widget, new ResizeEventArgs(RenderSize));
	}

	private void glArea_Render(object o, RenderArgs args)
	{
		if (!skipDraw)
		{
			skipDraw = true;

			// GTK makes the context current for us, so we need to clear it to hand it over to the Veldrid worker.
			Gdk.GLContext.ClearCurrent();
			if (Widget.GraphicsDevice == null)
				return;

			// Make context current on the Veldrid worker.
			if (Widget.GraphicsDevice.GetOpenGLInfo(out BackendInfoOpenGL glInfo))
			{
				// No need for this action to wait, we just need it done some time before issuing commands.
				glInfo.ExecuteOnGLThread(_makeCurrent);//, wait: false);
			}

			// It's important to only issue Veldrid commands in OnDraw,
			// since we only have a GL context current in the worker here.
			Callback.OnDraw(Widget, EventArgs.Empty);

			// Clear the context from the worker so GTK can use it again.
			// This action needs to wait so the context is cleared before we continue
			// out of this method (either for the coming MakeCurrent or for GTK).
			glInfo?.ExecuteOnGLThread(_clearCurrent);//, wait: true);

			// GTK does not seem to need the context current after the Render event,
			// but setting it back is safer if this assumption changes in the future.
			glArea?.MakeCurrent();
		}
		skipDraw = false;
	}

	private void MakeCurrent()
	{
		glArea?.MakeCurrent();
	}

	private void ClearCurrent()
	{
		Gdk.GLContext.ClearCurrent();
	}

	// TODO: Figure this one out! The docstring for this property in Veldrid's OpenGLPlatformInfo is ambiguous.
	IntPtr VeldridSurface.IOpenGL.OpenGLContextHandle => glArea?.Context.Handle ?? IntPtr.Zero;

	IntPtr VeldridSurface.IOpenGL.GetProcAddress(string name) => X11Interop.glXGetProcAddress(name);

	void VeldridSurface.IOpenGL.MakeCurrent(IntPtr context) => MakeCurrent();

	IntPtr VeldridSurface.IOpenGL.GetCurrentContext() => Gdk.GLContext.Current.Handle;

	void VeldridSurface.IOpenGL.ClearCurrentContext() => ClearCurrent();

	void VeldridSurface.IOpenGL.DeleteContext(IntPtr context)
	{
	}

	void VeldridSurface.IOpenGL.SwapBuffers()
	{
		// GLArea doesn't support drawing directly, so we queue a render but don't actually call OnDraw
		if (skipDraw)
			return;

		skipDraw = true;
		glArea?.QueueRender();
	}

	void VeldridSurface.IOpenGL.SetSyncToVerticalBlank(bool on)
	{
	}

	void VeldridSurface.IOpenGL.SetSwapchainFramebuffer()
	{
	}

	void VeldridSurface.IOpenGL.ResizeSwapchain(uint width, uint height)
	{
	}

	void Forms.Control.IHandler.Invalidate(Rectangle rect, bool invalidateChildren)
	{
		skipDraw = false;
		glArea?.QueueRender();
	}

	void Forms.Control.IHandler.Invalidate(bool invalidateChildren)
	{
		skipDraw = false;
		glArea?.QueueRender();
	}

	protected override global::Gtk.Widget CreateControl()
	{
		if (Widget.Backend == GraphicsBackend.OpenGL)
		{
			glArea = new GLArea();
			glArea.CanFocus = true;
			glArea.CanDefault = true;

			// Veldrid technically supports as low as OpenGL 3.0, but the full
			// complement of features is only available with 3.3 and higher.
			glArea.SetRequiredVersion(3, 3);

			glArea.HasDepthBuffer = true;
			glArea.HasStencilBuffer = true;
			// Control.Child = glArea;
			glArea.Realized += glArea_InitializeGraphicsBackend;
			return glArea;
		}
		else
		{
			EtoEventBox box = new();
			box.CanFocus = true;
			box.CanDefault = true;
			
			// Use Realized event for both X11 and Wayland - keep it simple
			box.Realized += Control_InitializeGraphicsBackend;
			
			return box;
		}
	}
}