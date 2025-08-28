using Eto.Drawing;
using Eto.Veldrid;
using Gtk;
using System;
using Veldrid;
using Veldrid.OpenGL;

namespace Eto.Veldrid.Gtk.Backends;

/// <summary>
/// Backend handler for OpenGL using GTK's GLArea
/// </summary>
internal class OpenGLBackendHandler : IVeldridBackendHandler, VeldridSurface.IOpenGL
{
    private GLArea? _glArea;
    private VeldridSurface.ICallback? _callback;
    private VeldridSurface? _surface;
    private bool _skipDraw;
    private readonly System.Action _makeCurrent;
    private readonly System.Action _clearCurrent;

    public GraphicsBackend Backend => GraphicsBackend.OpenGL;

    public OpenGLBackendHandler()
    {
        _makeCurrent = MakeCurrent;
        _clearCurrent = ClearCurrent;
    }

    public global::Gtk.Widget CreateWidget()
    {
        _glArea = new GLArea
        {
            CanFocus = true,
            CanDefault = true,
            HasDepthBuffer = true,
            HasStencilBuffer = true
        };

        // Veldrid technically supports as low as OpenGL 3.0, but the full
        // complement of features is only available with 3.3 and higher.
        _glArea.SetRequiredVersion(3, 3);

        return _glArea;
    }

    public Swapchain? CreateSwapchain(VeldridSurface surface, Size renderSize)
    {
        // For OpenGL, use the main swapchain from the graphics device
        return surface.GraphicsDevice?.MainSwapchain;
    }

    public void InitializeGraphicsDevice(VeldridSurface surface, Size renderSize)
    {
        if (_glArea == null)
            return;

        // Make context current to manually initialize a Veldrid GraphicsDevice
        _glArea.Context.MakeCurrent();

        // Create the OpenGL graphics device
        surface.GraphicsDevice = GraphicsDevice.CreateOpenGL(
            surface.GraphicsDeviceOptions,
            new OpenGLPlatformInfo(
                OpenGLContextHandle,
                GetProcAddress,
                MakeCurrent,
                GetCurrentContext,
                ClearCurrentContext,
                DeleteContext,
                SwapBuffers,
                SetSyncToVerticalBlank,
                SetSwapchainFramebuffer,
                ResizeSwapchain),
            (uint)renderSize.Width,
            (uint)renderSize.Height);

        // Clear context in the worker thread for now to make Mesa happy
        if (surface.GraphicsDevice?.GetOpenGLInfo(out BackendInfoOpenGL glInfo) == true)
        {
            // This action has to wait so GTK can manage the context after this method
            glInfo.ExecuteOnGLThread(_clearCurrent);
        }
    }

    public void HandleResize(Size newSize)
    {
        _skipDraw = false;
        _callback?.OnResize(_surface!, new ResizeEventArgs(newSize));
    }

    public void Invalidate()
    {
        _skipDraw = false;
        _glArea?.QueueRender();
    }

    public void SetupEventHandlers(VeldridSurface.ICallback callback, VeldridSurface surface)
    {
        _callback = callback;
        _surface = surface;

        if (_glArea != null)
        {
            _glArea.Realized += OnGLAreaRealized;
            _glArea.Render += OnGLAreaRender;
            _glArea.Resize += OnGLAreaResize;
        }
    }

    private void OnGLAreaRealized(object? sender, EventArgs e)
    {
        if (_glArea == null || _callback == null || _surface == null)
            return;

        var size = new Size((int)_glArea.AllocatedWidth, (int)_glArea.AllocatedHeight);
        _callback.OnInitializeBackend(_surface, new InitializeEventArgs(size));
    }

    private void OnGLAreaRender(object o, RenderArgs args)
    {
        if (_skipDraw || _surface?.GraphicsDevice == null)
        {
            _skipDraw = false;
            return;
        }

        _skipDraw = true;

        // GTK makes the context current for us, so we need to clear it to hand it over to the Veldrid worker
        Gdk.GLContext.ClearCurrent();

        // Make context current on the Veldrid worker
        if (_surface.GraphicsDevice.GetOpenGLInfo(out BackendInfoOpenGL glInfo))
        {
            // No need for this action to wait, we just need it done some time before issuing commands
            glInfo.ExecuteOnGLThread(_makeCurrent);
        }

        // It's important to only issue Veldrid commands in OnDraw,
        // since we only have a GL context current in the worker here
        _callback?.OnDraw(_surface, EventArgs.Empty);

        // Clear the context from the worker so GTK can use it again
        glInfo?.ExecuteOnGLThread(_clearCurrent);

        // GTK does not seem to need the context current after the Render event,
        // but setting it back is safer if this assumption changes in the future
        _glArea?.MakeCurrent();

        _skipDraw = false;
    }

    private void OnGLAreaResize(object o, ResizeArgs args)
    {
        var size = new Size((int)args.Width, (int)args.Height);
        HandleResize(size);
    }

    private void MakeCurrent()
    {
        _glArea?.MakeCurrent();
    }

    private void ClearCurrent()
    {
        Gdk.GLContext.ClearCurrent();
    }

    // IOpenGL implementation
    public IntPtr OpenGLContextHandle => _glArea?.Context.Handle ?? IntPtr.Zero;

    public IntPtr GetProcAddress(string name) => X11Interop.glXGetProcAddress(name);

    public void MakeCurrent(IntPtr context) => MakeCurrent();

    public IntPtr GetCurrentContext() => Gdk.GLContext.Current.Handle;

    public void ClearCurrentContext() => ClearCurrent();

    public void DeleteContext(IntPtr context)
    {
        // GTK manages the context lifecycle
    }

    public void SwapBuffers()
    {
        // GLArea doesn't support drawing directly, so we queue a render but don't actually call OnDraw
        if (_skipDraw)
            return;

        _skipDraw = true;
        _glArea?.QueueRender();
    }

    public void SetSyncToVerticalBlank(bool enable)
    {
        // Handled by GTK/driver
    }

    public void SetSwapchainFramebuffer()
    {
        // Not needed for GLArea
    }

    public void ResizeSwapchain(uint width, uint height)
    {
        // Handled automatically by GLArea
    }

    public void Dispose()
    {
        if (_glArea != null)
        {
            _glArea.Realized -= OnGLAreaRealized;
            _glArea.Render -= OnGLAreaRender;
            _glArea.Resize -= OnGLAreaResize;
            _glArea = null;
        }
        _callback = null;
        _surface = null;
    }
}