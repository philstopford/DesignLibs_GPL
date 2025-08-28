using Eto.Forms;
using Eto.Veldrid;
using System;
using Veldrid;

namespace Eto.Veldrid.Gtk.Examples
{
    /// <summary>
    /// Example demonstrating backend switching functionality
    /// </summary>
    public class BackendSwitchingExample
    {
        public static void Main(string[] args)
        {
            var platform = new Eto.GtkSharp.Platform();
            
            // Parse backend from command line
            GraphicsBackend backend = VeldridSurface.PreferredBackend;
            if (args.Length > 0 && Enum.TryParse<GraphicsBackend>(args[0], true, out var requestedBackend))
            {
                if (GraphicsDevice.IsBackendSupported(requestedBackend))
                {
                    backend = requestedBackend;
                    Console.WriteLine($"Using requested backend: {backend}");
                }
                else
                {
                    Console.WriteLine($"Requested backend {requestedBackend} not supported, using: {backend}");
                }
            }
            
            new Application(platform).Run(new ExampleForm(backend));
        }
    }
    
    public class ExampleForm : Form
    {
        private VeldridSurface surface;
        
        public ExampleForm(GraphicsBackend backend)
        {
            Title = $"Veldrid Example - {backend} Backend";
            Size = new Eto.Drawing.Size(800, 600);
            
            // Create surface with specified backend
            var options = new GraphicsDeviceOptions(
                debug: false,
                swapchainDepthFormat: PixelFormat.R32Float,
                syncToVerticalBlank: true,
                resourceBindingModel: ResourceBindingModel.Improved);
                
            surface = new VeldridSurface(backend, options);
            
            surface.VeldridInitialized += (sender, e) =>
            {
                Console.WriteLine($"Veldrid initialized with {backend} backend");
                Console.WriteLine($"Graphics Device: {surface.GraphicsDevice?.GetType().Name}");
                Console.WriteLine($"Backend Features: {surface.GraphicsDevice?.Features}");
            };
            
            surface.Draw += OnDraw;
            surface.Resize += OnResize;
            
            Content = surface;
            
            // Add menu for backend info
            var aboutCommand = new Command { MenuText = "About Backend", ToolBarText = "About" };
            aboutCommand.Executed += (sender, e) => ShowBackendInfo();
            
            Menu = new MenuBar
            {
                Items = { new ButtonMenuItem { Text = "&Help", Items = { aboutCommand } } }
            };
        }
        
        private void OnDraw(object sender, EventArgs e)
        {
            // Basic rendering example - clear to blue
            var device = surface.GraphicsDevice;
            var commandList = device.ResourceFactory.CreateCommandList();
            
            commandList.Begin();
            commandList.SetFramebuffer(device.SwapchainFramebuffer);
            commandList.ClearColorTarget(0, new RgbaFloat(0.2f, 0.4f, 0.8f, 1.0f));
            commandList.End();
            
            device.SubmitCommands(commandList);
            device.SwapBuffers();
            
            commandList.Dispose();
        }
        
        private void OnResize(object sender, ResizeEventArgs e)
        {
            Console.WriteLine($"Surface resized: {e.Width}x{e.Height}");
        }
        
        private void ShowBackendInfo()
        {
            var device = surface.GraphicsDevice;
            var info = $"Backend: {surface.Backend}\n" +
                      $"Device: {device?.DeviceName}\n" +
                      $"Vendor: {device?.VendorName}\n" +
                      $"API Version: {device?.ApiVersion}\n" +
                      $"Driver: {device?.DriverName} {device?.DriverInfo}\n" +
                      $"Uniform Buffer Alignment: {device?.UniformBufferMinOffsetAlignment}\n" +
                      $"Max Texture Size: {device?.GetPixelFormatSupport(PixelFormat.R8G8B8A8UNorm, TextureType.Texture2D, TextureUsage.Sampled)}";
                      
            MessageBox.Show(this, info, "Backend Information");
        }
    }
}