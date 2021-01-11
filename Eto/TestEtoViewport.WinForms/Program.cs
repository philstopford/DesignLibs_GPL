using Eto.OpenTK;
using Eto.Forms;
using System;

namespace TestEtoGl.WPF_Framebuffer
{
    static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            var platform = new Eto.WinForms.Platform();
            platform.Add<GLSurface.IHandler>(() => new Eto.OpenTK.WinForms.WinGLSurfaceHandler());

            new Application(platform).Run(new MainForm());
        }
    }
}