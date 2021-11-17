using Eto.Forms;
using System;
using Eto.GtkSharp;
using TestEtoGl;

namespace TestEtoVeldrid.Gtk;

public static class MainClass
{
	[STAThread]
	public static void Main(string[] args)
	{
		Platform platform = new();

		new Application(platform).Run(new MainForm());
	}
}