using Eto.Forms;
using Eto.Veldrid;
using System;
using Eto.GtkSharp;
using Veldrid;

namespace TestEtoVeldrid.Gtk;

public static class MainClass
{
	[STAThread]
	public static void Main(string[] args)
	{
		Platform platform = new();

		int backend_arg = Array.IndexOf(args, "--graphicsMode");

		GraphicsBackend backend = VeldridSurface.PreferredBackend;
			
		if (backend_arg != -1)
		{
			string backend_string = args[backend_arg + 1];
			switch (backend_string.ToLower())
			{
				case "opengl":
					backend = GraphicsBackend.OpenGL;
					break;
				case "d3d":
					backend = GraphicsBackend.Direct3D11;
					break;
				case "metal":
					backend = GraphicsBackend.Metal;
					break;
				case "vulkan":
				default:
					backend = GraphicsBackend.Vulkan;
					break;
			}
		}

		new Application(platform).Run(new MainForm(backend));
	}
}