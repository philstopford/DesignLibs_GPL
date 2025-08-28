using Eto.Drawing;
using Eto.Veldrid;
using Gtk;
using System;
using Veldrid;

namespace Eto.Veldrid.Gtk.Backends;

/// <summary>
/// Interface for handling different Veldrid graphics backends in GTK
/// </summary>
internal interface IVeldridBackendHandler : IDisposable
{
    /// <summary>
    /// The graphics backend this handler supports
    /// </summary>
    GraphicsBackend Backend { get; }

    /// <summary>
    /// Creates the appropriate GTK widget for this backend
    /// </summary>
    global::Gtk.Widget CreateWidget();

    /// <summary>
    /// Creates a swapchain for this backend
    /// </summary>
    Swapchain? CreateSwapchain(VeldridSurface surface, Size renderSize);

    /// <summary>
    /// Initializes the graphics device for this backend
    /// </summary>
    void InitializeGraphicsDevice(VeldridSurface surface, Size renderSize);

    /// <summary>
    /// Handles resize events
    /// </summary>
    void HandleResize(Size newSize);

    /// <summary>
    /// Forces a render/invalidation
    /// </summary>
    void Invalidate();

    /// <summary>
    /// Sets up event handlers
    /// </summary>
    void SetupEventHandlers(VeldridSurface.ICallback callback, VeldridSurface surface);
}