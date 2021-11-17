using Eto.Forms;
using Eto.Veldrid;
using System;
using Eto.Wpf;

namespace TestEtoVeldrid.Wpf;

public static class MainClass
{
	[STAThread]
	public static void Main(string[] args)
	{
		Platform platform = new();

		new Application(platform).Run(new MainForm());
	}
}