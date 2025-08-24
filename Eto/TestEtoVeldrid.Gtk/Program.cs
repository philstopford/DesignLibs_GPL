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

			// Check if user wants to force Vulkan backend for testing
			bool forceVulkan = args.Length > 0 && (args[0] == "--vulkan" || args[0] == "-v");
			
			if (forceVulkan)
			{
				Console.WriteLine("Forcing Vulkan backend for testing...");
				new Application(platform).Run(new MainForm(GraphicsBackend.Vulkan));
			}
			else
			{
				new Application(platform).Run(new MainForm());
			}
		}
	}
}
