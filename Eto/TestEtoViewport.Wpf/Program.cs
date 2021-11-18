using Eto.OpenTK;
using Eto.Forms;
using System;
using Eto.Wpf;

namespace TestEtoGl.WPF_WinformsHost;

internal static class Program
{
    /// <summary>
    /// The main entry point for the application.
    /// </summary>
    [STAThread]
    private static void Main()
    {
        Platform platform = new();
        platform.Add<GLSurface.IHandler>(() => new Eto.OpenTK.Wpf.WpfWinGLSurfaceHandler());

        new Application(platform).Run(new MainForm());
    }
}