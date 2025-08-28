using Eto.Drawing;
using Eto.GtkSharp.Forms;
using Eto.Veldrid;
using Eto.Veldrid.Gtk.Backends;
using Gtk;
using System;
using Veldrid;

namespace Eto.Veldrid.Gtk;

public class GtkVeldridSurfaceHandler : GtkControl<global::Gtk.Widget, VeldridSurface, VeldridSurface.ICallback>, VeldridSurface.IHandler, VeldridSurface.IOpenGL
{
	private IVeldridBackendHandler? _backendHandler;
	private global::Gtk.Widget? _widget;
	
	public Size RenderSize => Size.Round((SizeF)Widget.Size * Scale);

	private float Scale => Widget.ParentWindow?.Screen?.LogicalPixelSize ?? 1;

	public override global::Gtk.Widget ContainerContentControl => _widget ?? base.ContainerContentControl;

	public GtkVeldridSurfaceHandler()
	{
	}

	public Swapchain? CreateSwapchain()
	{
		return _backendHandler?.CreateSwapchain(Widget, RenderSize);
	}

	public void InitializeGraphicsDevice(VeldridSurface surface, InitializeEventArgs e)
	{
		_backendHandler?.InitializeGraphicsDevice(surface, new Size(e.Width, e.Height));
	}

	protected override global::Gtk.Widget CreateControl()
	{
		// Create appropriate backend handler based on the graphics backend
		_backendHandler = VeldridBackendFactory.CreateBackendHandler(Widget.Backend);
		_backendHandler.SetupEventHandlers(Callback, Widget);
		
		_widget = _backendHandler.CreateWidget();
		return _widget;
	}

	// IOpenGL interface implementation - delegate to OpenGL backend handler if available
	IntPtr VeldridSurface.IOpenGL.OpenGLContextHandle => 
		_backendHandler is OpenGLBackendHandler openGLHandler ? openGLHandler.OpenGLContextHandle : IntPtr.Zero;

	IntPtr VeldridSurface.IOpenGL.GetProcAddress(string name) => 
		_backendHandler is OpenGLBackendHandler openGLHandler ? openGLHandler.GetProcAddress(name) : IntPtr.Zero;

	void VeldridSurface.IOpenGL.MakeCurrent(IntPtr context) => 
		(_backendHandler as OpenGLBackendHandler)?.MakeCurrent(context);

	IntPtr VeldridSurface.IOpenGL.GetCurrentContext() => 
		_backendHandler is OpenGLBackendHandler openGLHandler ? openGLHandler.GetCurrentContext() : IntPtr.Zero;

	void VeldridSurface.IOpenGL.ClearCurrentContext() => 
		(_backendHandler as OpenGLBackendHandler)?.ClearCurrentContext();

	void VeldridSurface.IOpenGL.DeleteContext(IntPtr context) => 
		(_backendHandler as OpenGLBackendHandler)?.DeleteContext(context);

	void VeldridSurface.IOpenGL.SwapBuffers() => 
		(_backendHandler as OpenGLBackendHandler)?.SwapBuffers();

	void VeldridSurface.IOpenGL.SetSyncToVerticalBlank(bool enable) => 
		(_backendHandler as OpenGLBackendHandler)?.SetSyncToVerticalBlank(enable);

	void VeldridSurface.IOpenGL.SetSwapchainFramebuffer() => 
		(_backendHandler as OpenGLBackendHandler)?.SetSwapchainFramebuffer();

	void VeldridSurface.IOpenGL.ResizeSwapchain(uint width, uint height) => 
		(_backendHandler as OpenGLBackendHandler)?.ResizeSwapchain(width, height);

	void Forms.Control.IHandler.Invalidate(Rectangle rect, bool invalidateChildren)
	{
		_backendHandler?.Invalidate();
	}

	void Forms.Control.IHandler.Invalidate(bool invalidateChildren)
	{
		_backendHandler?.Invalidate();
	}

	protected override void Dispose(bool disposing)
	{
		if (disposing)
		{
			_backendHandler?.Dispose();
			_backendHandler = null;
		}
		base.Dispose(disposing);
	}
}