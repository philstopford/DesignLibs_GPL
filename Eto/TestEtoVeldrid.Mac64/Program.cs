using Eto.Forms;
using System;
using TestEtoVeldrid;

namespace TestEtoVeldrid.Mac64
{
	public static class MainClass
	{
		[STAThread]
		public static void Main(string[] args)
		{
			var platform = new Eto.Mac.Platform();

			new Application(platform).Run(new MainForm());
		}
	}
}
