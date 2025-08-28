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
        try
        {
            Console.WriteLine("OpenGLBackendHandler: Creating GLArea widget");
            
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
            
            // Try to set auto render for better compatibility in headless environments
            try
            {
                _glArea.AutoRender = true;
            }
            catch (Exception ex)
            {
                Console.WriteLine($"OpenGLBackendHandler: Warning - could not set AutoRender: {ex.Message}");
            }

            Console.WriteLine("OpenGLBackendHandler: GLArea widget created successfully");
            return _glArea;
        }
        catch (Exception ex)
        {
            Console.WriteLine($"OpenGLBackendHandler: Error creating GLArea widget: {ex.Message}");
            throw;
        }
    }

    public Swapchain? CreateSwapchain(VeldridSurface surface, Size renderSize)
    {
        // For OpenGL, use the main swapchain from the graphics device
        return surface.GraphicsDevice?.MainSwapchain;
    }

    public void InitializeGraphicsDevice(VeldridSurface surface, Size renderSize)
    {
        try
        {
            Console.WriteLine($"OpenGLBackendHandler: InitializeGraphicsDevice called with size {renderSize}");
            
            if (_glArea == null)
            {
                Console.WriteLine("OpenGLBackendHandler: GLArea is null, cannot initialize");
                return;
            }

            Console.WriteLine($"OpenGLBackendHandler: GLArea context available: {_glArea.Context != null}");
            
            // Make context current to manually initialize a Veldrid GraphicsDevice
            try
            {
                _glArea.Context.MakeCurrent();
                Console.WriteLine("OpenGLBackendHandler: Made OpenGL context current");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"OpenGLBackendHandler: Failed to make context current: {ex.Message}");
                throw new VeldridException($"Failed to make OpenGL context current: {ex.Message}", ex);
            }

            // Create the OpenGL graphics device
            Console.WriteLine("OpenGLBackendHandler: Creating OpenGL graphics device...");
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

            Console.WriteLine($"OpenGLBackendHandler: Successfully created graphics device: {surface.GraphicsDevice}");

            // Clear context in the worker thread for now to make Mesa happy
            if (surface.GraphicsDevice?.GetOpenGLInfo(out BackendInfoOpenGL glInfo) == true)
            {
                Console.WriteLine("OpenGLBackendHandler: Clearing context on GL thread");
                // This action has to wait so GTK can manage the context after this method
                glInfo.ExecuteOnGLThread(_clearCurrent);
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine($"OpenGLBackendHandler: Error in InitializeGraphicsDevice: {ex.Message}");
            Console.WriteLine($"OpenGLBackendHandler: Stack trace: {ex.StackTrace}");
            throw new VeldridException($"Failed to initialize OpenGL graphics device: {ex.Message}", ex);
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
            
            // Schedule a manual initialization check in case Realized event doesn't fire
            // This is common in headless/CI environments
            GLib.Timeout.Add(50, () => {
                Console.WriteLine($"OpenGLBackendHandler: Manual initialization check - IsRealized: {_glArea.IsRealized}");
                if (!_glArea.IsRealized && _surface?.GraphicsDevice == null)
                {
                    Console.WriteLine("OpenGLBackendHandler: GLArea not realized after timeout, attempting manual initialization");
                    // Try to trigger manual initialization
                    TryManualInitialization();
                }
                return false; // Don't repeat
            });
        }
    }

    private void OnGLAreaRealized(object? sender, EventArgs e)
    {
        Console.WriteLine("OpenGLBackendHandler: OnGLAreaRealized called");
        TryInitialization();
    }
    
    private void TryManualInitialization()
    {
        if (_glArea == null || _callback == null || _surface == null)
        {
            Console.WriteLine("OpenGLBackendHandler: Manual initialization failed - missing components");
            return;
        }
        
        // For headless environments, we may need to force the size
        var width = _glArea.AllocatedWidth > 0 ? _glArea.AllocatedWidth : 800;
        var height = _glArea.AllocatedHeight > 0 ? _glArea.AllocatedHeight : 600;
        
        Console.WriteLine($"OpenGLBackendHandler: Manual initialization with size {width}x{height}");
        var size = new Size(width, height);
        _callback.OnInitializeBackend(_surface, new InitializeEventArgs(size));
    }
    
    private void TryInitialization()
    {
        if (_glArea == null || _callback == null || _surface == null)
        {
            Console.WriteLine("OpenGLBackendHandler: Initialization failed - missing components");
            return;
        }

        var size = new Size((int)_glArea.AllocatedWidth, (int)_glArea.AllocatedHeight);
        Console.WriteLine($"OpenGLBackendHandler: Initializing with size {size}");
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