using Eto.Forms;
using System;
using Eto.WinForms;

namespace TestEtoVeldrid.WinForms;

public static class Program
{
	[STAThread]
	public static void Main(string[] args)
	{
		Platform platform = new();

		new Application(platform).Run(new MainForm());
	}
}