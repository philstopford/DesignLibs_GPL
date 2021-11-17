using Eto.OpenTK;
using Eto.Forms;
using System;
using Eto.WinForms;

namespace TestEtoGl.WPF_Framebuffer;

internal static class Program
{
    /// <summary>
    /// The main entry point for the application.
    /// </summary>
    [STAThread]
    private static void Main()
    {
        Platform platform = new();
        platform.Add<GLSurface.IHandler>(() => new Eto.OpenTK.WinForms.WinGLSurfaceHandler());

        new Application(platform).Run(new MainForm());
    }
}