using Eto.Drawing;
using Eto.GtkSharp.Forms;
using Eto.Veldrid;
using Gtk;
using System;
using System.Runtime.InteropServices;
using Veldrid;

namespace Eto.Veldrid.Gtk.Backends;

/// <summary>
/// Backend handler for Vulkan with support for both X11 and Wayland
/// </summary>
internal class VulkanBackendHandler : IVeldridBackendHandler
{
    private EtoEventBox? _eventBox;
    private VeldridSurface.ICallback? _callback;
    private VeldridSurface? _surface;
    private WindowingSystemDetector.WindowingSystem _windowingSystem;

    public GraphicsBackend Backend => GraphicsBackend.Vulkan;

    public VulkanBackendHandler()
    {
        _windowingSystem = WindowingSystemDetector.DetectWindowingSystem();
        Console.WriteLine($"VulkanBackendHandler: Detected windowing system: {_windowingSystem}");
    }

    public global::Gtk.Widget CreateWidget()
    {
        _eventBox = new EtoEventBox
        {
            CanFocus = true,
            CanDefault = true
        };

        return _eventBox;
    }

    public Swapchain? CreateSwapchain(VeldridSurface surface, Size renderSize)
    {
        try
        {
            Console.WriteLine($"VulkanBackendHandler: Attempting to create swapchain with size {renderSize}");
            
            if (_eventBox == null || surface.GraphicsDevice == null)
            {
                Console.WriteLine($"VulkanBackendHandler: Cannot create swapchain - eventBox: {_eventBox != null}, graphicsDevice: {surface.GraphicsDevice != null}");
                return null;
            }

            SwapchainSource source = CreateSwapchainSource();
            Console.WriteLine($"VulkanBackendHandler: Created swapchain source for {_windowingSystem}");

            var swapchain = surface.GraphicsDevice.ResourceFactory.CreateSwapchain(
                new SwapchainDescription(
                    source,
                    (uint)renderSize.Width,
                    (uint)renderSize.Height,
                    surface.GraphicsDeviceOptions.SwapchainDepthFormat,
                    surface.GraphicsDeviceOptions.SyncToVerticalBlank,
                    surface.GraphicsDeviceOptions.SwapchainSrgbFormat));
            
            Console.WriteLine($"VulkanBackendHandler: Successfully created swapchain: {swapchain}");
            return swapchain;
        }
        catch (Exception ex)
        {
            Console.WriteLine($"VulkanBackendHandler: Failed to create swapchain: {ex.Message}");
            throw;
        }
    }

    private SwapchainSource CreateSwapchainSource()
    {
        if (_eventBox == null)
            throw new InvalidOperationException("EventBox widget is not available");

        switch (_windowingSystem)
        {
            case WindowingSystemDetector.WindowingSystem.Wayland:
                return CreateWaylandSwapchainSource();
            
            case WindowingSystemDetector.WindowingSystem.X11:
            default:
                return CreateX11SwapchainSource();
        }
    }

    private SwapchainSource CreateWaylandSwapchainSource()
    {
        // Get Wayland display and surface
        IntPtr waylandDisplay = WaylandInterop.gdk_wayland_display_get_wl_display(_eventBox!.Display.Handle);
        IntPtr waylandSurface = WaylandInterop.gdk_wayland_window_get_wl_surface(_eventBox.Window.Handle);

        return SwapchainSource.CreateWayland(waylandDisplay, waylandSurface);
    }

    private SwapchainSource CreateX11SwapchainSource()
    {
        // Use existing X11 interop
        IntPtr x11Display = X11Interop.gdk_x11_display_get_xdisplay(_eventBox!.Display.Handle);
        IntPtr x11Window = X11Interop.gdk_x11_window_get_xid(_eventBox.Window.Handle);

        return SwapchainSource.CreateXlib(x11Display, x11Window);
    }

    public void InitializeGraphicsDevice(VeldridSurface surface, Size renderSize)
    {
        try
        {
            Console.WriteLine($"VulkanBackendHandler: Attempting to create Vulkan graphics device with size {renderSize}");
            
            // Create the Vulkan graphics device
            surface.GraphicsDevice = GraphicsDevice.CreateVulkan(surface.GraphicsDeviceOptions);
            
            Console.WriteLine($"VulkanBackendHandler: Successfully created Vulkan graphics device: {surface.GraphicsDevice}");
        }
        catch (Exception ex)
        {
            Console.WriteLine($"VulkanBackendHandler: Failed to create Vulkan graphics device: {ex.Message}");
            
            // Don't throw immediately - let the calling code handle the fallback
            // This allows the application to continue and potentially switch to OpenGL
            throw new VeldridException($"Failed to initialize Vulkan graphics device: {ex.Message}", ex);
        }
    }

    public void HandleResize(Size newSize)
    {
        _callback?.OnResize(_surface!, new ResizeEventArgs(newSize));
    }

    public void Invalidate()
    {
        // For Vulkan, we don't need to queue renders like with GLArea
        // The application will handle rendering through the Veldrid surface
    }

    public void SetupEventHandlers(VeldridSurface.ICallback callback, VeldridSurface surface)
    {
        _callback = callback;
        _surface = surface;

        if (_eventBox != null)
        {
            _eventBox.Realized += OnEventBoxRealized;
            _eventBox.SizeAllocated += OnEventBoxSizeAllocated;
            
            Console.WriteLine("VulkanBackendHandler: Setting up immediate initialization");
            
            // Try immediate initialization - don't wait for GTK events
            TryManualInitialization();
        }
    }

    private void OnEventBoxRealized(object? sender, EventArgs e)
    {
        Console.WriteLine("VulkanBackendHandler: OnEventBoxRealized called");
        TryInitialization();
    }
    
    private void TryManualInitialization()
    {
        if (_eventBox == null || _callback == null || _surface == null)
        {
            Console.WriteLine("VulkanBackendHandler: Manual initialization failed - missing components");
            return;
        }
        
        if (_surface.GraphicsDevice != null)
        {
            Console.WriteLine("VulkanBackendHandler: Graphics device already exists, skipping manual initialization");
            return;
        }
        
        try
        {
            Console.WriteLine("VulkanBackendHandler: Starting manual initialization");
            
            // Force EventBox to show and map if needed
            if (!_eventBox.IsMapped)
            {
                Console.WriteLine("VulkanBackendHandler: EventBox not mapped, attempting to map");
                _eventBox.Show();
            }
            
            // For headless environments, we may need to force the size
            var width = _eventBox.AllocatedWidth > 0 ? _eventBox.AllocatedWidth : 800;
            var height = _eventBox.AllocatedHeight > 0 ? _eventBox.AllocatedHeight : 600;
            
            Console.WriteLine($"VulkanBackendHandler: Manual initialization with size {width}x{height}");
            Console.WriteLine($"VulkanBackendHandler: EventBox state - IsRealized: {_eventBox.IsRealized}, IsMapped: {_eventBox.IsMapped}, Visible: {_eventBox.Visible}");
            
            var size = new Size(width, height);
            _callback.OnInitializeBackend(_surface, new InitializeEventArgs(size));
            Console.WriteLine("VulkanBackendHandler: Manual initialization completed successfully");
        }
        catch (Exception ex)
        {
            Console.WriteLine($"VulkanBackendHandler: Error in manual initialization: {ex.Message}");
            Console.WriteLine($"VulkanBackendHandler: Stack trace: {ex.StackTrace}");
        }
    }
    
    private void TryInitialization()
    {
        try
        {
            Console.WriteLine("VulkanBackendHandler: EventBox realized");
            
            if (_eventBox == null || _callback == null || _surface == null)
            {
                Console.WriteLine($"VulkanBackendHandler: Missing components - eventBox: {_eventBox != null}, callback: {_callback != null}, surface: {_surface != null}");
                return;
            }

            var size = new Size(_eventBox.AllocatedWidth, _eventBox.AllocatedHeight);
            Console.WriteLine($"VulkanBackendHandler: EventBox size: {size}");
            
            _callback.OnInitializeBackend(_surface, new InitializeEventArgs(size));
            Console.WriteLine("VulkanBackendHandler: Called OnInitializeBackend");
        }
        catch (Exception ex)
        {
            Console.WriteLine($"VulkanBackendHandler: Error in TryInitialization: {ex.Message}");
            Console.WriteLine($"VulkanBackendHandler: Stack trace: {ex.StackTrace}");
        }
    }

    private void OnEventBoxSizeAllocated(object o, SizeAllocatedArgs args)
    {
        var size = new Size(args.Allocation.Width, args.Allocation.Height);
        HandleResize(size);
    }

    public void Dispose()
    {
        if (_eventBox != null)
        {
            _eventBox.Realized -= OnEventBoxRealized;
            _eventBox.SizeAllocated -= OnEventBoxSizeAllocated;
            _eventBox = null;
        }
        _callback = null;
        _surface = null;
    }
}

/// <summary>
/// Wayland interop functions for creating Vulkan surfaces
/// </summary>
internal static class WaylandInterop
{
    private const string linux_libgdk_wayland_name = "libgdk-3.so.0";

    [DllImport(linux_libgdk_wayland_name)]
    public static extern IntPtr gdk_wayland_display_get_wl_display(IntPtr gdkDisplay);

    [DllImport(linux_libgdk_wayland_name)]
    public static extern IntPtr gdk_wayland_window_get_wl_surface(IntPtr gdkWindow);
}