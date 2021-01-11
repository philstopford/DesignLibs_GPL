using Eto.Forms;
using System;
using TestEtoGl;

namespace TestEtoVeldrid.Gtk
{
	public static class MainClass
	{
		[STAThread]
		public static void Main(string[] args)
		{
			var platform = new Eto.GtkSharp.Platform();

			new Application(platform).Run(new MainForm());
		}
	}
}
