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
	private DrawingArea? drawingArea;
	private System.Action _makeCurrent;
	private System.Action _clearCurrent;
	public Size RenderSize => Size.Round((SizeF)Widget.Size * Scale);

	private float Scale => Widget.ParentWindow?.Screen?.LogicalPixelSize ?? 1;

	public override global::Gtk.Widget ContainerContentControl => glArea ?? drawingArea ?? base.ContainerContentControl;

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

		// For Vulkan backend, attempt direct creation without conservative validation
		return CreateVulkanSwapchain();
	}

	private Swapchain? CreateVulkanSwapchain()
	{
		// Ensure we have a proper native window for Vulkan operations
		if (drawingArea?.Window == null)
		{
			Console.WriteLine("[DEBUG] DrawingArea window not available for Vulkan swapchain creation");
			return null;
		}

		// Ensure the widget is properly realized and has a native window
		if (!drawingArea.IsRealized)
		{
			Console.WriteLine("[DEBUG] DrawingArea not realized, cannot create Vulkan swapchain");
			return null;
		}

		// Get display info for debugging
		var gdkDisplay = drawingArea.Display.Handle;
		bool isWayland = X11Interop.IsWaylandDisplay(gdkDisplay);
		
		Console.WriteLine($"[DEBUG] Creating Vulkan swapchain - Wayland: {isWayland}, Realized: {drawingArea.IsRealized}");

		try
		{
			SwapchainSource source = CreateSwapchainSource();
			
			Size renderSize = RenderSize;
			var swapchain = Widget.GraphicsDevice?.ResourceFactory.CreateSwapchain(
				new SwapchainDescription(
					source,
					(uint)renderSize.Width,
					(uint)renderSize.Height,
					Widget.GraphicsDeviceOptions.SwapchainDepthFormat,
					Widget.GraphicsDeviceOptions.SyncToVerticalBlank,
					Widget.GraphicsDeviceOptions.SwapchainSrgbFormat));

			if (swapchain != null)
			{
				Console.WriteLine("[DEBUG] Vulkan swapchain created successfully");
			}
			else
			{
				Console.WriteLine("[DEBUG] Failed to create Vulkan swapchain");
			}

			return swapchain;
		}
		catch (Exception ex)
		{
			Console.WriteLine($"[DEBUG] Exception creating Vulkan swapchain: {ex.Message}");
			return null;
		}
	}

	private SwapchainSource CreateSwapchainSource()
	{
		// Use drawingArea for Vulkan operations instead of Control
		var gdkDisplay = drawingArea!.Display.Handle;
		
		bool isX11 = X11Interop.IsX11Display(gdkDisplay);
		bool isWayland = X11Interop.IsWaylandDisplay(gdkDisplay);

		// Add debugging information
		var displayNamePtr = X11Interop.gdk_display_get_name(gdkDisplay);
		var displayName = displayNamePtr == IntPtr.Zero ? "unknown" : System.Runtime.InteropServices.Marshal.PtrToStringAnsi(displayNamePtr) ?? "unknown";
		Console.WriteLine($"[DEBUG] Display name: {displayName}, IsX11: {isX11}, IsWayland: {isWayland}");

		if (isX11 && !isWayland)
		{
			Console.WriteLine("[DEBUG] Creating X11/Xlib SwapchainSource");
			return SwapchainSource.CreateXlib(
				X11Interop.gdk_x11_display_get_xdisplay(gdkDisplay),
				X11Interop.gdk_x11_window_get_xid(drawingArea.Window.Handle));
		}
		else if (isWayland && !isX11)
		{
			Console.WriteLine("[DEBUG] Creating Wayland SwapchainSource");
			return CreateWaylandSwapchainSource(gdkDisplay);
		}
		else
		{
			throw new NotSupportedException($"Unsupported or ambiguous windowing system. Display: {displayName}, IsX11: {isX11}, IsWayland: {isWayland}");
		}
	}

	private SwapchainSource CreateWaylandSwapchainSource(IntPtr gdkDisplay)
	{
		Console.WriteLine("[DEBUG] Creating Wayland SwapchainSource from DrawingArea");
		
		// Use drawingArea for Wayland operations
		if (drawingArea?.Window == null)
		{
			throw new InvalidOperationException("DrawingArea window is null");
		}
		
		Console.WriteLine($"[DEBUG] DrawingArea state - Visible: {drawingArea.Visible}, Realized: {drawingArea.IsRealized}, Mapped: {drawingArea.IsMapped}, Window.IsVisible: {drawingArea.Window.IsVisible}");
		
		// Get Wayland handles from DrawingArea
		var waylandDisplay = X11Interop.gdk_wayland_display_get_wl_display(gdkDisplay);
		var waylandSurface = X11Interop.gdk_wayland_window_get_wl_surface(drawingArea.Window.Handle);
		
		if (waylandDisplay == IntPtr.Zero || waylandSurface == IntPtr.Zero)
		{
			throw new InvalidOperationException($"Failed to get Wayland handles - Display: {waylandDisplay}, Surface: {waylandSurface}");
		}
		
		Console.WriteLine($"[DEBUG] Successfully obtained Wayland handles - Display: {waylandDisplay}, Surface: {waylandSurface}");
		
		try
		{
			var swapchainSource = SwapchainSource.CreateWayland(waylandDisplay, waylandSurface);
			Console.WriteLine("[DEBUG] Wayland SwapchainSource created successfully");
			return swapchainSource;
		}
		catch (Exception ex)
		{
			Console.WriteLine($"[DEBUG] Failed to create Wayland SwapchainSource: {ex.Message}");
			throw;
		}
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

	private void DrawingArea_InitializeGraphicsBackend(object? sender, EventArgs e)
	{
		Console.WriteLine("[DEBUG] DrawingArea mapped and ready for graphics initialization");
		
		// Ensure native window for proper Vulkan surface access
		if (drawingArea?.Window != null)
		{
			X11Interop.gdk_window_ensure_native(drawingArea.Window.Handle);
			Console.WriteLine("[DEBUG] Ensured native window for DrawingArea");
		}
		
		// Ensure the widget is shown - this is critical for Wayland
		if (drawingArea != null && !drawingArea.Visible)
		{
			Console.WriteLine("[DEBUG] Making DrawingArea visible");
			drawingArea.ShowAll();
		}
		
		// Get display info for debugging
		var gdkDisplay = drawingArea!.Display.Handle;
		bool isWayland = X11Interop.IsWaylandDisplay(gdkDisplay);
		
		if (isWayland)
		{
			Console.WriteLine("[DEBUG] Wayland detected, processing events before initialization");
			// Process pending events to ensure surface is ready
			while (global::Gtk.Application.EventsPending())
			{
				global::Gtk.Application.RunIteration();
			}
			
			// Small delay to ensure compositor has committed the surface
			System.Threading.Thread.Sleep(50);
		}
		
		Console.WriteLine($"[DEBUG] DrawingArea state - Visible: {drawingArea.Visible}, Realized: {drawingArea.IsRealized}, Mapped: {drawingArea.IsMapped}");
		
		// Initialize graphics backend
		Console.WriteLine("[DEBUG] Initializing graphics backend");
		
		try
		{
			Callback.OnInitializeBackend(Widget, new InitializeEventArgs(RenderSize));
			
			// Add a draw event handler for continuous rendering
			if (drawingArea != null)
			{
				drawingArea.Drawn += DrawingArea_Draw;
				Console.WriteLine("[DEBUG] Added DrawingArea draw event handler");
			}
		}
		catch (Exception ex)
		{
			Console.WriteLine($"[DEBUG] Graphics backend initialization failed: {ex.Message}");
			throw;
		}
	}

	private void DrawingArea_Draw(object o, DrawnArgs args)
	{
		if (!skipDraw)
		{
			skipDraw = true;
			
			try
			{
				// Trigger Veldrid drawing for Vulkan backend
				Console.WriteLine("[DEBUG] DrawingArea draw event triggered");
				Callback.OnDraw(Widget, EventArgs.Empty);
			}
			catch (Exception ex)
			{
				Console.WriteLine($"[DEBUG] Exception in DrawingArea_Draw: {ex.Message}");
			}
		}
		skipDraw = false;
	}

	private void Control_InitializeGraphicsBackend(object? sender, EventArgs e)
	{
		// Legacy method for EtoEventBox - kept for compatibility
		Console.WriteLine("[DEBUG] Initializing graphics backend (legacy)");
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
		if (glArea != null)
		{
			glArea.QueueRender();
		}
		else if (drawingArea != null)
		{
			// For Vulkan backend, trigger drawing directly
			TriggerVulkanDraw();
		}
	}

	void Forms.Control.IHandler.Invalidate(bool invalidateChildren)
	{
		skipDraw = false;
		if (glArea != null)
		{
			glArea.QueueRender();
		}
		else if (drawingArea != null)
		{
			// For Vulkan backend, trigger drawing directly
			TriggerVulkanDraw();
		}
	}

	private void TriggerVulkanDraw()
	{
		try
		{
			// Queue a draw on the DrawingArea widget to trigger the Drawn event
			drawingArea?.QueueDraw();
		}
		catch (Exception ex)
		{
			Console.WriteLine($"[DEBUG] Exception in TriggerVulkanDraw: {ex.Message}");
		}
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
			glArea.Realized += glArea_InitializeGraphicsBackend;
			return glArea;
		}
		else
		{
			// Use DrawingArea for Vulkan backend instead of EtoEventBox for better native window support
			drawingArea = new DrawingArea();
			drawingArea.CanFocus = true;
			drawingArea.CanDefault = true;
			
			// Set minimum size to ensure proper widget allocation
			drawingArea.SetSizeRequest(100, 100);
			
			// Ensure native window creation for proper Vulkan surface access
			// This is critical for Wayland support
			drawingArea.AppPaintable = true;
			
			// Make sure the widget is visible
			drawingArea.Visible = true;
			
			// Use Map event for initialization - ensures widget is visible and ready for graphics operations
			drawingArea.Mapped += DrawingArea_InitializeGraphicsBackend;
			
			Console.WriteLine("[DEBUG] Created DrawingArea for Vulkan backend");
			return drawingArea;
		}
	}
}