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

			// Print available backends and platform information
			Console.WriteLine("=== Vulkan Viewport Linux Support Test ===");
			Console.WriteLine($"Operating System: {Environment.OSVersion}");
			Console.WriteLine($"Environment Variables:");
			Console.WriteLine($"  DISPLAY = {Environment.GetEnvironmentVariable("DISPLAY")}");
			Console.WriteLine($"  WAYLAND_DISPLAY = {Environment.GetEnvironmentVariable("WAYLAND_DISPLAY")}");
			Console.WriteLine($"  XDG_SESSION_TYPE = {Environment.GetEnvironmentVariable("XDG_SESSION_TYPE")}");
			
			Console.WriteLine("\n=== Available Graphics Backends ===");
			foreach (GraphicsBackend backend in Enum.GetValues<GraphicsBackend>())
			{
				bool supported = GraphicsDevice.IsBackendSupported(backend);
				Console.WriteLine($"  {backend}: {(supported ? "✓ Supported" : "✗ Not Supported")}");
			}

			Console.WriteLine($"\n=== Preferred Backend ===");
			try
			{
				var preferredBackend = VeldridSurface.PreferredBackend;
				Console.WriteLine($"  Selected: {preferredBackend}");
			}
			catch (Exception e)
			{
				Console.WriteLine($"  Error: {e.Message}");
			}

			Console.WriteLine("\n=== Starting Application ===");
			new Application(platform).Run(new MainForm());
		}
	}
}
