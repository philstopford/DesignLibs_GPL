using System;
using Veldrid;

namespace Eto.Veldrid.Gtk.Backends;

/// <summary>
/// Factory for creating appropriate backend handlers
/// </summary>
internal static class VeldridBackendFactory
{
    /// <summary>
    /// Creates a backend handler for the specified graphics backend
    /// </summary>
    public static IVeldridBackendHandler CreateBackendHandler(GraphicsBackend backend)
    {
        return backend switch
        {
            GraphicsBackend.OpenGL => new OpenGLBackendHandler(),
            GraphicsBackend.Vulkan => new VulkanBackendHandler(),
            _ => throw new NotSupportedException($"Graphics backend {backend} is not supported on GTK platform")
        };
    }

    /// <summary>
    /// Checks if the specified backend is supported on the current system
    /// </summary>
    public static bool IsBackendSupported(GraphicsBackend backend)
    {
        switch (backend)
        {
            case GraphicsBackend.OpenGL:
                return GraphicsDevice.IsBackendSupported(GraphicsBackend.OpenGL);
            
            case GraphicsBackend.Vulkan:
                return GraphicsDevice.IsBackendSupported(GraphicsBackend.Vulkan);
            
            default:
                return false;
        }
    }

    /// <summary>
    /// Gets the preferred backend for the current system
    /// </summary>
    public static GraphicsBackend GetPreferredBackend()
    {
        // Prefer Vulkan if available, otherwise fall back to OpenGL
        if (IsBackendSupported(GraphicsBackend.Vulkan))
        {
            return GraphicsBackend.Vulkan;
        }
        
        if (IsBackendSupported(GraphicsBackend.OpenGL))
        {
            return GraphicsBackend.OpenGL;
        }

        throw new NotSupportedException("No supported Veldrid backend found for GTK platform");
    }
}