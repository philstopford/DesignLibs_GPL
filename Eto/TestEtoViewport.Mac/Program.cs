using System;
using Eto.Forms;
using Eto.Mac;
using Eto.OpenTK;
using Eto.OpenTK.Mac;

namespace TestEtoGl.Mac;

public static class Program
{
    [STAThread]
    public static void Main(string[] args)
    {
        Platform gen = new();

        // shouldn't be needed, Eto needs fixing
        gen.Add<GLSurface.IHandler>(() => new MacGLSurfaceHandler());

        // run application with our main form
        new Application(gen).Run(new MainForm());
    }
}