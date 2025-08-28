using Eto.Forms;
using Eto.Veldrid;
using System;
using Veldrid;

namespace TestEtoVeldrid.Gtk
{
	public static class MainClass
	{
		[STAThread]
		public static void Main(string[] args)
		{
			// Check for headless environment - but be less strict about it
			var display = Environment.GetEnvironmentVariable("DISPLAY");
			var waylandDisplay = Environment.GetEnvironmentVariable("WAYLAND_DISPLAY");
			
			if (string.IsNullOrEmpty(display) && string.IsNullOrEmpty(waylandDisplay))
			{
				Console.WriteLine("ERROR: No display available. This application requires a graphical environment.");
				Console.WriteLine("For headless testing, consider using Xvfb:");
				Console.WriteLine("  xvfb-run -a dotnet run");
				Console.WriteLine("Or set up a virtual display.");
				Environment.Exit(1);
			}

			var platform = new Eto.GtkSharp.Platform();
			var app = new Application(platform);

			// Parse command line arguments for backend selection
			GraphicsBackend? requestedBackend = null;
			
			if (args.Length > 0)
			{
				if (Enum.TryParse<GraphicsBackend>(args[0], true, out var parsedBackend))
				{
					if (GraphicsDevice.IsBackendSupported(parsedBackend))
					{
						requestedBackend = parsedBackend;
						Console.WriteLine($"Using requested backend: {requestedBackend}");
					}
					else
					{
						Console.WriteLine($"Requested backend {parsedBackend} is not supported, will use default with fallback");
					}
				}
				else
				{
					Console.WriteLine($"Invalid backend specified: {args[0]}, will use default with fallback");
				}
			}

			// Determine backend to use with fallback logic
			GraphicsBackend backend = requestedBackend ?? VeldridSurface.PreferredBackend;
			
			if (requestedBackend == null)
			{
				Console.WriteLine($"Using default backend: {backend}");
			}

			// Try to create main form with the selected backend
			MainForm? mainForm = null;
			
			try
			{
				mainForm = new MainForm(backend);
			}
			catch (VeldridException vex) when (backend == GraphicsBackend.Vulkan && requestedBackend == null)
			{
				Console.WriteLine($"Vulkan backend failed: {vex.Message}");
				Console.WriteLine("Attempting fallback to OpenGL...");
				
				// Try with OpenGL as fallback
				try
				{
					mainForm = new MainForm(GraphicsBackend.OpenGL);
				}
				catch (Exception fallbackEx)
				{
					Console.WriteLine($"OpenGL fallback also failed: {fallbackEx.Message}");
					Console.WriteLine("No working graphics backend found. Please check your graphics drivers.");
					Environment.Exit(1);
				}
			}
			catch (Exception ex)
			{
				Console.WriteLine($"Application failed to start: {ex.Message}");
				
				// Check if this is likely a display/graphics issue
				if (ex.Message.Contains("display") || ex.Message.Contains("GL") || ex.Message.Contains("graphics"))
				{
					Console.WriteLine();
					Console.WriteLine("This appears to be a graphics/display issue. Common solutions:");
					Console.WriteLine("1. Ensure you have working graphics drivers installed");
					Console.WriteLine("2. For headless environments, use: xvfb-run -a dotnet run");
					Console.WriteLine("3. Try running with software rendering: LIBGL_ALWAYS_SOFTWARE=1 dotnet run");
					Console.WriteLine("4. Check that your X11/Wayland display is working with: glxinfo");
				}
				
				Environment.Exit(1);
			}

			if (mainForm != null)
			{
				app.Run(mainForm);
			}
		}
	}
}
