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
        if (_eventBox == null || surface.GraphicsDevice == null)
            return null;

        SwapchainSource source = CreateSwapchainSource();

        return surface.GraphicsDevice.ResourceFactory.CreateSwapchain(
            new SwapchainDescription(
                source,
                (uint)renderSize.Width,
                (uint)renderSize.Height,
                surface.GraphicsDeviceOptions.SwapchainDepthFormat,
                surface.GraphicsDeviceOptions.SyncToVerticalBlank,
                surface.GraphicsDeviceOptions.SwapchainSrgbFormat));
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
        // Create the Vulkan graphics device
        surface.GraphicsDevice = GraphicsDevice.CreateVulkan(surface.GraphicsDeviceOptions);
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
        }
    }

    private void OnEventBoxRealized(object? sender, EventArgs e)
    {
        if (_eventBox == null || _callback == null || _surface == null)
            return;

        var size = new Size(_eventBox.AllocatedWidth, _eventBox.AllocatedHeight);
        _callback.OnInitializeBackend(_surface, new InitializeEventArgs(size));
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