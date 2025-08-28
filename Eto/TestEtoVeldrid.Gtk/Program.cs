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

			// Parse command line arguments for backend selection
			GraphicsBackend backend = VeldridSurface.PreferredBackend;
			
			if (args.Length > 0)
			{
				if (Enum.TryParse<GraphicsBackend>(args[0], true, out var requestedBackend))
				{
					if (GraphicsDevice.IsBackendSupported(requestedBackend))
					{
						backend = requestedBackend;
						Console.WriteLine($"Using requested backend: {backend}");
					}
					else
					{
						Console.WriteLine($"Requested backend {requestedBackend} is not supported, using default: {backend}");
					}
				}
				else
				{
					Console.WriteLine($"Invalid backend specified: {args[0]}, using default: {backend}");
				}
			}
			else
			{
				Console.WriteLine($"Using default backend: {backend}");
			}

			new Application(platform).Run(new MainForm(backend));
		}
	}
}
