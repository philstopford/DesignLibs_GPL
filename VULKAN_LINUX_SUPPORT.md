# Vulkan Support for Linux

This repository now includes comprehensive Vulkan support for Linux with automatic detection and compatibility for all major display server configurations.

## Supported Display Servers

✅ **X11**: Traditional X Window System  
✅ **Wayland**: Modern compositor-based display server  
✅ **XWayland**: X11 applications running under Wayland  

## Automatic Backend Selection

The VeldridSurface now automatically selects the best available graphics backend on Linux:

1. **Vulkan** (preferred) - For modern graphics performance and features
2. **OpenGL** (fallback) - For compatibility when Vulkan is unavailable

## Platform Detection

The system automatically detects your display server configuration using:

- **GDK Display Detection**: Primary method using GTK/GDK APIs
- **Environment Variables**: Fallback detection using standard Linux environment variables
  - `DISPLAY` - X11 display identifier
  - `WAYLAND_DISPLAY` - Wayland display identifier  
  - `XDG_SESSION_TYPE` - Session type hint

## Usage

No code changes required! Simply create a VeldridSurface as usual:

```csharp
// Automatically selects Vulkan if available, OpenGL otherwise
var surface = new VeldridSurface();

// Or explicitly specify Vulkan
var surface = new VeldridSurface(GraphicsBackend.Vulkan);
```

## Environment Examples

### Pure X11 Session
```bash
DISPLAY=:0
XDG_SESSION_TYPE=x11
# Uses X11 interop with Vulkan/OpenGL
```

### Pure Wayland Session  
```bash
WAYLAND_DISPLAY=wayland-0
XDG_SESSION_TYPE=wayland
# Uses Wayland interop with Vulkan/OpenGL
```

### XWayland (X11 apps under Wayland)
```bash
DISPLAY=:0
WAYLAND_DISPLAY=wayland-0
XDG_SESSION_TYPE=wayland
# Automatically detected as XWayland, uses X11 interop
```

## Requirements

- **Vulkan**: Install Vulkan drivers and loader
  ```bash
  # Ubuntu/Debian
  sudo apt install vulkan-tools libvulkan1 mesa-vulkan-drivers
  
  # Fedora
  sudo dnf install vulkan-tools vulkan-loader mesa-vulkan-drivers
  ```

- **Development**: Vulkan development headers (if building from source)
  ```bash
  # Ubuntu/Debian  
  sudo apt install libvulkan-dev
  
  # Fedora
  sudo dnf install vulkan-headers vulkan-loader-devel
  ```

## Troubleshooting

### Check Vulkan Support
```bash
# Verify Vulkan installation
vulkaninfo --summary

# List available devices
vkcube # Should show a spinning cube
```

### Check Display Server
```bash
# Check current session type
echo $XDG_SESSION_TYPE

# Check display variables
echo "X11: $DISPLAY"
echo "Wayland: $WAYLAND_DISPLAY"
```

### Force OpenGL Fallback
If you encounter Vulkan issues, you can force OpenGL:

```csharp
var surface = new VeldridSurface(GraphicsBackend.OpenGL);
```

## Technical Details

The implementation includes:

- **LinuxPlatformInterop**: Platform detection and interop functions
- **Dynamic SwapchainSource**: Automatic selection between Xlib and Wayland
- **Error Handling**: Graceful fallbacks with informative error messages
- **Environment Detection**: Robust detection across different Linux distributions

For more details, see the source code in `Eto/Eto.Veldrid.Gtk/`.