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
				Environment.Exit(1);
			}

			if (mainForm != null)
			{
				app.Run(mainForm);
			}
		}
	}
}
