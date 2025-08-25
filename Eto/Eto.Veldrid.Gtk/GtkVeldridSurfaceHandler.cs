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
		Swapchain? swapchain;

		if (Widget.Backend == GraphicsBackend.OpenGL)
		{
			swapchain = Widget.GraphicsDevice?.MainSwapchain;
		}
		else
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
				// Wayland path - use Wayland SwapchainSource
				Console.WriteLine("[DEBUG] Using Wayland SwapchainSource");
				source = CreateWaylandSwapchainSource(gdkDisplay);
			}
			else
			{
				throw new NotSupportedException($"Unsupported or ambiguous windowing system detected. Display: {displayName}, IsX11: {isX11}, IsWayland: {isWayland}. Only X11 and Wayland are supported for Vulkan backend.");
			}

			Size renderSize = RenderSize;
			swapchain = Widget.GraphicsDevice?.ResourceFactory.CreateSwapchain(
				new SwapchainDescription(
					source,
					(uint)renderSize.Width,
					(uint)renderSize.Height,
					Widget.GraphicsDeviceOptions.SwapchainDepthFormat,
					Widget.GraphicsDeviceOptions.SyncToVerticalBlank,
					Widget.GraphicsDeviceOptions.SwapchainSrgbFormat));
		}

		return swapchain;
	}

	private SwapchainSource CreateWaylandSwapchainSource(IntPtr gdkDisplay)
	{
		Console.WriteLine("[DEBUG] Creating Wayland SwapchainSource");
		
		// Ensure the widget is realized and properly configured
		if (!Control.IsRealized)
		{
			Console.WriteLine("[DEBUG] Widget not realized, calling Realize()...");
			Control.Realize();
		}
		
		// Ensure we have a window
		if (Control.Window == null)
		{
			throw new InvalidOperationException("Control window is null - widget not properly initialized");
		}
		
		// Ensure the window is visible and mapped
		if (!Control.Window.IsVisible || !Control.IsMapped)
		{
			Console.WriteLine("[DEBUG] Window not visible/mapped, showing widget...");
			Control.ShowAll();
			
			// Process events to ensure the window is properly mapped
			while (global::Gtk.Application.EventsPending())
			{
				global::Gtk.Application.RunIteration(false);
			}
			
			// Wait for the widget to be properly mapped
			int attempts = 0;
			while (!Control.IsMapped && attempts < 50)
			{
				System.Threading.Thread.Sleep(10);
				while (global::Gtk.Application.EventsPending())
				{
					global::Gtk.Application.RunIteration(false);
				}
				attempts++;
			}
			
			if (!Control.IsMapped)
			{
				throw new InvalidOperationException("Widget could not be properly mapped for Wayland surface creation");
			}
		}
		
		Console.WriteLine($"[DEBUG] Window state - Visible: {Control.Window.IsVisible}, Mapped: {Control.IsMapped}, Handle: {Control.Window.Handle}");
		
		// Additional check: ensure the widget has been allocated proper size
		var allocation = Control.Allocation;
		if (allocation.Width <= 1 || allocation.Height <= 1)
		{
			Console.WriteLine($"[DEBUG] Warning: Widget has small allocation ({allocation.Width}x{allocation.Height})");
		}
		
		// Get Wayland handles with enhanced error checking
		var waylandDisplay = X11Interop.gdk_wayland_display_get_wl_display(gdkDisplay);
		if (waylandDisplay == IntPtr.Zero)
		{
			throw new InvalidOperationException("Failed to get Wayland display handle - display may not be ready");
		}
		
		var waylandSurface = X11Interop.gdk_wayland_window_get_wl_surface(Control.Window.Handle);
		if (waylandSurface == IntPtr.Zero)
		{
			throw new InvalidOperationException("Failed to get Wayland window surface handle - surface may not be ready");
		}
		
		Console.WriteLine($"[DEBUG] Successfully obtained Wayland handles - Display: {waylandDisplay}, Surface: {waylandSurface}");
		
		// Ensure all pending surface operations are committed before creating swapchain
		Console.WriteLine("[DEBUG] Committing surface state before swapchain creation");
		while (global::Gtk.Application.EventsPending())
		{
			global::Gtk.Application.RunIteration(false);
		}
		
		// Force a surface commit to ensure the compositor has fully processed the surface
		if (Control.Window.IsVisible)
		{
			Control.Window.ProcessUpdates(false);
		}
		
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
		Callback.OnInitializeBackend(Widget, new InitializeEventArgs(RenderSize));
	}

	private bool _waylandInitialized = false;
	
	private void Control_InitializeGraphicsBackendWayland(object? sender, EventArgs e)
	{
		// For Wayland, we only want to initialize once and only after the widget has been
		// properly mapped and visible to avoid protocol errors
		if (_waylandInitialized)
			return;
		
		if (sender is EtoEventBox box && box.IsVisible && box.IsMapped && 
		    box.Allocation.Width > 1 && box.Allocation.Height > 1)
		{
			Console.WriteLine("[DEBUG] Wayland widget fully mapped and visible, proceeding with initialization");
			_waylandInitialized = true;
			
			// Ensure the widget's window is also properly configured
			if (box.Window != null && box.Window.IsVisible)
			{
				// Process all pending events to ensure surface is fully committed
				while (global::Gtk.Application.EventsPending())
				{
					global::Gtk.Application.RunIteration(false);
				}
				
				// Additional delay to ensure compositor has fully processed the surface
				global::GLib.Timeout.Add(100, () => {
					if (box.Window != null && box.Window.IsVisible && box.IsMapped)
					{
						Callback.OnInitializeBackend(Widget, new InitializeEventArgs(RenderSize));
					}
					return false; // Don't repeat
				});
			}
		}
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
			
			// For Wayland with Vulkan, we need to wait until the window is properly mapped
			// and visible before initializing the graphics backend
			var gdkDisplay = box.Display?.Handle ?? X11Interop.gdk_display_get_default();
			bool isWayland = X11Interop.IsWaylandDisplay(gdkDisplay);
			
			if (isWayland)
			{
				Console.WriteLine("[DEBUG] Detected Wayland, using maximum deferred initialization");
				// For Wayland, wait until the widget is fully mapped and visible
				// This ensures the surface is completely ready for graphics operations
				box.Mapped += Control_InitializeGraphicsBackendWayland;
			}
			else
			{
				// For X11, use the regular initialization timing
				box.Realized += Control_InitializeGraphicsBackend;
			}
			
			return box;
		}
	}
}