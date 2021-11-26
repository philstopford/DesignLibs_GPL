using System;
using System.Globalization;
using Burkardt.NavierStokesNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace NavierStokes2DExactTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    navier_stokes_2d_exact_test tests navier_stokes_2d_exact().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("navier_stokes_2d_exact_test");
        Console.WriteLine("  Test navier_stokes_2d_exact().");
        //
        //  GMS Flow.
        //
        uvp_gms_test();
        uvp_gms_test2();
        rhs_gms_test();
        resid_gms_test();
        gnuplot_gms_test();
        //
        //  Lukas Bystricky Flow.
        //
        uvp_lukas_test();
        uvp_lukas_test2();
        rhs_lukas_test();
        resid_lukas_test();
        gnuplot_lukas_test();
        //
        //  Poiseuille Flow.
        //
        uvp_poiseuille_test();
        uvp_poiseuille_test2();
        rhs_poiseuille_test();
        resid_poiseuille_test();
        gnuplot_poiseuille_test();
        parameter_poiseuille_test();
        //
        //  Spiral Flow.
        //
        uvp_spiral_test();
        uvp_spiral_test2();
        rhs_spiral_test();
        resid_spiral_test();
        gnuplot_spiral_test();
        parameter_spiral_test();
        //
        //  Taylor Flow.
        //
        uvp_taylor_test();
        uvp_taylor_test2();
        rhs_taylor_test();
        resid_taylor_test();
        gnuplot_taylor_test();
        parameter_taylor_test();
        //
        //  Vortex Flow.
        //
        uvp_vortex_test();
        uvp_vortex_test2();
        rhs_vortex_test();
        resid_vortex_test();
        gnuplot_vortex_test();
        parameter_vortex_test();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("navier_stokes_2d_exact_test");
        Console.WriteLine("  Normal end of execution.");
    }

    private static void uvp_gms_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    uvp_gms_test samples the GMS flow at a specific time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2020
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;
        const double t = 1.0;

        Console.WriteLine("");
        Console.WriteLine("uvp_gms_test");
        Console.WriteLine("  GMS flow:");
        Console.WriteLine("  Estimate the range of velocity and pressure");
        Console.WriteLine("  at time T = 1,");
        Console.WriteLine("  over the [-1,+1]x[-1,+1] square.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];

        const double r8_lo = -1.0;
        const double r8_hi = 1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);

        GMS.uvp_gms(nu, rho, n, x, y, t, ref u, ref v, ref p);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  U:"
                          + "  " + typeMethods.r8vec_min(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  V:"
                          + "  " + typeMethods.r8vec_min(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  P:"
                          + "  " + typeMethods.r8vec_min(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void uvp_gms_test2()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    uvp_gms_test2 samples the GMS flow on the boundary.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2020
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        const double nu = 1.0;
        const double rho = 1.0;
        const double t = 1.0;

        Console.WriteLine("");
        Console.WriteLine("uvp_gms_test2");
        Console.WriteLine("  GMS flow:");
        Console.WriteLine("  Estimate the range of velocity and pressure");
        Console.WriteLine("  at time T = 1,");
        Console.WriteLine("  over the boundary of the [-1,+1]x[-1,+1] square.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 400;

        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];
        double[] x = new double[n];
        double[] y = new double[n];

        const double r8_lo = -1.0;
        const double r8_hi = 1.0;

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref x);
        for (i = 0; i < 100; i++)
        {
            y[i] = r8_lo;
        }

        for (i = 100; i < 200; i++)
        {
            x[i] = r8_hi;
        }

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref y, index: +100);

        typeMethods.r8vec_linspace(100, r8_hi, r8_lo, ref x, index: +200);
        for (i = 200; i < 300; i++)
        {
            y[i] = r8_hi;
        }

        for (i = 300; i < 400; i++)
        {
            x[i] = r8_lo;
        }

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref y, index: +300);

        GMS.uvp_gms(nu, rho, n, x, y, t, ref u, ref v, ref p);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  U:"
                          + "  " + typeMethods.r8vec_min(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  V:"
                          + "  " + typeMethods.r8vec_min(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  P:"
                          + "  " + typeMethods.r8vec_min(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void rhs_gms_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    rhs_gms_test samples the GMS right hand side at a specific time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2020
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;
        const double t = 1.0;

        Console.WriteLine("");
        Console.WriteLine("rhs_gms_test");
        Console.WriteLine("  GMS flow:");
        Console.WriteLine("  Sample the Navier-Stokes right hand sides");
        Console.WriteLine("  at time T = 1,");
        Console.WriteLine("  over the [-1,+1]x[-1,+1] square.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] f = new double[n];
        double[] g = new double[n];
        double[] h = new double[n];

        const double r8_lo = -1.0;
        const double r8_hi = 1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);

        GMS.rhs_gms(nu, rho, n, x, y, t, ref f, ref g, ref h);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  F:"
                          + "  " + typeMethods.r8vec_min(n, f).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, f).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  G:"
                          + "  " + typeMethods.r8vec_min(n, g).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, g).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  H:"
                          + "  " + typeMethods.r8vec_min(n, h).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, h).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void resid_gms_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    resid_gms_test samples the GMS residual at a specific time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2020
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;
        const double t = 1.0;

        Console.WriteLine("");
        Console.WriteLine("resid_gms_test");
        Console.WriteLine("  GMS flow:");
        Console.WriteLine("  Sample the Navier-Stokes residuals");
        Console.WriteLine("  at time T = 1,");
        Console.WriteLine("  over the [-1,+1]x[-1,+1] square.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        int n = 1000;

        double[] pr = new double[n];
        double[] ur = new double[n];
        double[] vr = new double[n];

        const double r8_lo = -1.0;
        const double r8_hi = 1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);

        GMS.resid_gms(nu, rho, n, x, y, t, ref ur, ref vr, ref pr);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  Ur:"
                          + "  " + typeMethods.r8vec_amin(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Vr:"
                          + "  " + typeMethods.r8vec_amin(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Pr:"
                          + "  " + typeMethods.r8vec_amin(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void gnuplot_gms_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    gnuplot_gms_test plots the GMS flow.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2020
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int x_num = 21;
        const int y_num = 21;

        const double nu = 1.0;
        const double rho = 1.0;
        const double t = 1.0;

        Console.WriteLine("");
        Console.WriteLine("gnuplot_gms_test:");
        Console.WriteLine("  GMS flow:");
        Console.WriteLine("  Generate a velocity field on a regular grid.");
        Console.WriteLine("  Store in GNUPLOT data and command files.");

        const double x_lo = -1.0;
        const double x_hi = 1.0;

        const double y_lo = -1.0;
        const double y_hi = 1.0;

        double[] x = new double[x_num * y_num];
        double[] y = new double[x_num * y_num];

        NavierStokes2DExact.grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi, ref x, ref y);

        int n = x_num * y_num;

        double[] u = new double[x_num * y_num];
        double[] v = new double[x_num * y_num];
        double[] p = new double[x_num * y_num];

        GMS.uvp_gms(nu, rho, n, x, y, t, ref u, ref v, ref p);

        string header = "gms";
        const double s = 0.25;
        NavierStokes2DExact.ns2de_gnuplot(header, n, x, y, u, v, p, s);

    }

    private static void uvp_lukas_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    uvp_lukas_test samples Lukas Bystricky's flow at the initial time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 March 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("uvp_lukas_test");
        Console.WriteLine("  Lukas Bystricky Flow:");
        Console.WriteLine("  Estimate the range of velocity and pressure");
        Console.WriteLine("  at the initial time T = 0, over the unit square.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];

        const double r8_lo = 0.0;
        const double r8_hi = 1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        const double t = 0.0;

        LukasBystricky.uvp_lukas(nu, rho, n, x, y, t, ref u, ref v, ref p);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  U:"
                          + "  " + typeMethods.r8vec_min(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  V:"
                          + "  " + typeMethods.r8vec_min(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  P:"
                          + "  " + typeMethods.r8vec_min(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void uvp_lukas_test2()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    uvp_lukas_test2 samples Lukas Bystricky's flow on the boundary.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 March 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        const double nu = 1.0;
        const double rho = 1.0;
        const double t = 0.0;

        Console.WriteLine("");
        Console.WriteLine("uvp_lukas_test2");
        Console.WriteLine("  Lukas Bystricky Flow:");
        Console.WriteLine("  Estimate the range of velocity and pressure");
        Console.WriteLine("  on the boundary");
        Console.WriteLine("  at the initial time T = 0, over the unit square.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 400;

        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];
        double[] x = new double[n];
        double[] y = new double[n];

        const double r8_lo = 0.0;
        const double r8_hi = 1.0;

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref x);
        for (i = 0; i < 100; i++)
        {
            y[i] = r8_lo;
        }

        for (i = 100; i < 200; i++)
        {
            x[i] = r8_hi;
        }

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref y, index: +100);

        typeMethods.r8vec_linspace(100, r8_hi, r8_lo, ref x, index: +200);
        for (i = 200; i < 300; i++)
        {
            y[i] = r8_hi;
        }

        for (i = 300; i < 400; i++)
        {
            x[i] = r8_lo;
        }

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref y, index: +300);

        LukasBystricky.uvp_lukas(nu, rho, n, x, y, t, ref u, ref v, ref p);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  U:"
                          + "  " + typeMethods.r8vec_min(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  V:"
                          + "  " + typeMethods.r8vec_min(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  P:"
                          + "  " + typeMethods.r8vec_min(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void rhs_lukas_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    rhs_lukas_test samples Lukas Bystricky's right hand side at the initial time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 March 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("rhs_lukas_test");
        Console.WriteLine("  Lukas Bystricky Flow:");
        Console.WriteLine("  Sample the Navier-Stokes right hand sides");
        Console.WriteLine("  at the initial time T = 0, using the unit square.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] f = new double[n];
        double[] g = new double[n];
        double[] h = new double[n];

        const double r8_lo = 0.0;
        const double r8_hi = 1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        const double t = 0.0;

        LukasBystricky.rhs_lukas(nu, rho, n, x, y, t, ref f, ref g, ref h);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  F:"
                          + "  " + typeMethods.r8vec_min(n, f).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, f).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  G:"
                          + "  " + typeMethods.r8vec_min(n, g).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, g).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  H:"
                          + "  " + typeMethods.r8vec_min(n, h).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, h).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void resid_lukas_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    resid_lukas_test samples Lukas Bystricky's residual at the initial time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 March 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("resid_lukas_test");
        Console.WriteLine("  Lukas Bystricky Flow:");
        Console.WriteLine("  Sample the Navier-Stokes residuals");
        Console.WriteLine("  at the initial time T = 0, on the unit square.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] pr = new double[n];
        double[] ur = new double[n];
        double[] vr = new double[n];

        const double r8_lo = 0.0;
        const double r8_hi = 1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        const double t = 0.0;

        LukasBystricky.resid_lukas(nu, rho, n, x, y, t, ref ur, ref vr, ref pr);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  Ur:"
                          + "  " + typeMethods.r8vec_amin(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Vr:"
                          + "  " + typeMethods.r8vec_amin(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Pr:"
                          + "  " + typeMethods.r8vec_amin(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void gnuplot_lukas_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    gnuplot_lukas_test plots Lukas Bystricky's flow.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 March 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int x_num = 21;
        const int y_num = 21;

        Console.WriteLine("");
        Console.WriteLine("gnuplot_lukas_test:");
        Console.WriteLine("  Lukas Bystricky Flow:");
        Console.WriteLine("  Generate a velocity field on a regular grid.");
        Console.WriteLine("  Store in GNUPLOT data and command files.");

        const double x_lo = 0.0;
        const double x_hi = 1.0;

        const double y_lo = 0.0;
        const double y_hi = 1.0;

        double[] x = new double[x_num * y_num];
        double[] y = new double[x_num * y_num];

        NavierStokes2DExact.grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi, ref x, ref y);

        const double nu = 1.0;
        const double rho = 1.0;
        const int n = x_num * y_num;
        const double t = 0.0;

        double[] u = new double[x_num * y_num];
        double[] v = new double[x_num * y_num];
        double[] p = new double[x_num * y_num];

        LukasBystricky.uvp_lukas(nu, rho, n, x, y, t, ref u, ref v, ref p);

        const string header = "lukas";
        const double s = 0.25;
        NavierStokes2DExact.ns2de_gnuplot(header, n, x, y, u, v, p, s);

    }

    private static void uvp_poiseuille_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    uvp_poiseuille_test samples Poiseuille flow at the initial time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("uvp_poiseuille_test");
        Console.WriteLine("  Poiseuille Flow:");
        Console.WriteLine("  Estimate the range of velocity and pressure");
        Console.WriteLine("  at the initial time T = 0, over a channel region.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];

        const double x_lo = 0.0;
        const double x_hi = 6.0;
        const double y_lo = -1.0;
        const double y_hi = +1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, x_lo, x_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, y_lo, y_hi, ref seed);
        const double t = 0.0;

        Poiseuille.uvp_poiseuille(nu, rho, n, x, y, t, ref u, ref v, ref p);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  U:"
                          + "  " + typeMethods.r8vec_min(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  V:"
                          + "  " + typeMethods.r8vec_min(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  P:"
                          + "  " + typeMethods.r8vec_min(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void uvp_poiseuille_test2()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    uvp_poiseuille_test2 samples Poiseuille flow on the boundary.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        const double nu = 1.0;
        const double rho = 1.0;
        const double t = 0.0;

        Console.WriteLine("");
        Console.WriteLine("uvp_poiseuille_test2");
        Console.WriteLine("  Poiseuille Flow:");
        Console.WriteLine("  Estimate the range of velocity and pressure");
        Console.WriteLine("  on the boundary");
        Console.WriteLine("  at the initial time T = 0, over a channel region.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 400;

        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];
        double[] x = new double[n];
        double[] y = new double[n];

        const double x_lo = 0.0;
        const double x_hi = 6.0;
        const double y_lo = -1.0;
        const double y_hi = +1.0;

        typeMethods.r8vec_linspace(100, x_lo, x_hi, ref x);
        for (i = 0; i < 100; i++)
        {
            y[i] = y_lo;
        }

        for (i = 100; i < 200; i++)
        {
            x[i] = x_hi;
        }

        typeMethods.r8vec_linspace(100, y_lo, y_hi, ref y, index: +100);

        typeMethods.r8vec_linspace(100, x_hi, x_lo, ref x, index: +200);
        for (i = 200; i < 300; i++)
        {
            y[i] = y_hi;
        }

        for (i = 300; i < 400; i++)
        {
            x[i] = x_lo;
        }

        typeMethods.r8vec_linspace(100, y_hi, y_lo, ref y, index: +300);

        Poiseuille.uvp_poiseuille(nu, rho, n, x, y, t, ref u, ref v, ref p);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  U:"
                          + "  " + typeMethods.r8vec_min(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  V:"
                          + "  " + typeMethods.r8vec_min(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  P:"
                          + "  " + typeMethods.r8vec_min(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void rhs_poiseuille_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    rhs_poiseuille_test samples Poiseuille right hand sides.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("rhs_poiseuille_test");
        Console.WriteLine("  Poiseuille Flow:");
        Console.WriteLine("  Sample the Navier-Stokes right hand sides");
        Console.WriteLine("  at the initial time T = 0, over a channel region.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] f = new double[n];
        double[] g = new double[n];
        double[] h = new double[n];

        const double x_lo = 0.0;
        const double x_hi = 6.0;
        const double y_lo = -1.0;
        const double y_hi = +1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, x_lo, x_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, y_lo, y_hi, ref seed);
        const double t = 0.0;

        Poiseuille.rhs_poiseuille(nu, rho, n, x, y, t, ref f, ref g, ref h);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  F:"
                          + "  " + typeMethods.r8vec_min(n, f).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, f).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  G:"
                          + "  " + typeMethods.r8vec_min(n, g).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, g).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  H:"
                          + "  " + typeMethods.r8vec_min(n, h).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, h).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void resid_poiseuille_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    resid_poiseuille_test samples Poiseuille residuals.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("resid_poiseuille_test");
        Console.WriteLine("  Poiseuille Flow:");
        Console.WriteLine("  Sample the Navier-Stokes residuals");
        Console.WriteLine("  at the initial time T = 0, over a channel region.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] pr = new double[n];
        double[] ur = new double[n];
        double[] vr = new double[n];

        const double x_lo = 0.0;
        const double x_hi = 6.0;
        const double y_lo = -1.0;
        const double y_hi = +1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, x_lo, x_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, y_lo, y_hi, ref seed);
        const double t = 0.0;

        Poiseuille.resid_poiseuille(nu, rho, n, x, y, t, ref ur, ref vr, ref pr);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  Ur:"
                          + "  " + typeMethods.r8vec_amin(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Vr:"
                          + "  " + typeMethods.r8vec_amin(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Pr:"
                          + "  " + typeMethods.r8vec_amin(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void gnuplot_poiseuille_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    gnuplot_poiseuille_test plots the Poiseuille flow.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int x_num = 21;
        const int y_num = 21;

        Console.WriteLine("");
        Console.WriteLine("gnuplot_poiseuille_test:");
        Console.WriteLine("  Poiseuille Flow:");
        Console.WriteLine("  Generate a velocity field on a regular grid.");
        Console.WriteLine("  Store in GNUPLOT data and command files.");

        const double x_lo = 0.0;
        const double x_hi = 6.0;

        const double y_lo = -1.0;
        const double y_hi = +1.0;

        double[] x = new double[x_num * y_num];
        double[] y = new double[x_num * y_num];

        NavierStokes2DExact.grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi, ref x, ref y);

        const double nu = 1.0;
        const double rho = 1.0;
        const int n = x_num * y_num;
        const double t = 0.0;

        double[] u = new double[x_num * y_num];
        double[] v = new double[x_num * y_num];
        double[] p = new double[x_num * y_num];

        Poiseuille.uvp_poiseuille(nu, rho, n, x, y, t, ref u, ref v, ref p);

        const string header = "poiseuille";
        const double s = 5.00;
        NavierStokes2DExact.ns2de_gnuplot(header, n, x, y, u, v, p, s);

    }

    private static void parameter_poiseuille_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARAMETER_poiseuille_test monitors Poiseuille solution norms for NU, RHO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int k;
        double p_norm;
        double t;
        double u_norm;
        double v_norm;

        Console.WriteLine("");
        Console.WriteLine("PARAMETER_poiseuille_test");
        Console.WriteLine("  Poiseuille Flow:");
        Console.WriteLine("  Monitor solution norms over time for various");
        Console.WriteLine("  values of NU, RHO.");

        const int n = 1000;

        double[] u = new double[n];
        double[] v = new double[n];
        double[] p = new double[n];

        const double x_lo = 0.0;
        const double x_hi = 6.0;
        const double y_lo = -1.0;
        const double y_hi = +1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, x_lo, x_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, y_lo, y_hi, ref seed);
        //
        //  Vary RHO.
        //
        Console.WriteLine("");
        Console.WriteLine("  RHO affects the pressure scaling.");
        Console.WriteLine("");
        Console.WriteLine("     RHO         NU           T     ||U||       ||V||       ||P||");
        Console.WriteLine("");

        double nu = 1.0;
        double rho = 1.0;

        for (j = 1; j <= 3; j++)
        {
            for (k = 0; k <= 5; k++)
            {
                t = k / 5.0;

                Poiseuille.uvp_poiseuille(nu, rho, n, x, y, t, ref u, ref v, ref p);

                u_norm = typeMethods.r8vec_norm_l2(n, u) / n;
                v_norm = typeMethods.r8vec_norm_l2(n, v) / n;
                p_norm = typeMethods.r8vec_norm_l2(n, p) / n;

                Console.WriteLine("  " + rho.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + nu.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + v_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + p_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }

            Console.WriteLine("");
            rho /= 100.0;
        }

        /*
        Vary NU.
        */
        Console.WriteLine("");
        Console.WriteLine("  NU affects the time scaling.");
        Console.WriteLine("");
        Console.WriteLine("     RHO         NU           T     ||U||       ||V||       ||P||");
        Console.WriteLine("");

        nu = 1.0;
        rho = 1.0;

        for (i = 1; i <= 4; i++)
        {
            for (k = 0; k <= 5; k++)
            {
                t = k / 5.0;

                Poiseuille.uvp_poiseuille(nu, rho, n, x, y, t, ref u, ref v, ref p);

                u_norm = typeMethods.r8vec_norm_l2(n, u) / n;
                v_norm = typeMethods.r8vec_norm_l2(n, v) / n;
                p_norm = typeMethods.r8vec_norm_l2(n, p) / n;

                Console.WriteLine("  " + rho.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + nu.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + v_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + p_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }

            Console.WriteLine("");

            nu /= 10.0;
        }

    }

    private static void uvp_spiral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    uvp_spiral_test samples the Spiral flow at the initial time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("uvp_spiral_test");
        Console.WriteLine("  Spiral Flow:");
        Console.WriteLine("  Estimate the range of velocity and pressure");
        Console.WriteLine("  at the initial time T = 0, over the unit square.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];

        const double r8_lo = 0.0;
        const double r8_hi = 1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        const double t = 0.0;

        Spiral.uvp_spiral(nu, rho, n, x, y, t, ref u, ref v, ref p);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  U:"
                          + "  " + typeMethods.r8vec_min(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  V:"
                          + "  " + typeMethods.r8vec_min(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  P:"
                          + "  " + typeMethods.r8vec_min(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void uvp_spiral_test2()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    uvp_spiral_test2 samples the Spiral flow on the boundary.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 March 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        const double nu = 1.0;
        const double rho = 1.0;
        const double t = 0.0;

        Console.WriteLine("");
        Console.WriteLine("uvp_spiral_test2");
        Console.WriteLine("  Spiral Flow:");
        Console.WriteLine("  Estimate the range of velocity and pressure");
        Console.WriteLine("  on the boundary");
        Console.WriteLine("  at the initial time T = 0, over the unit square.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 400;

        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];
        double[] x = new double[n];
        double[] y = new double[n];

        double r8_lo = 0.0;
        double r8_hi = 1.0;

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref x);
        for (i = 0; i < 100; i++)
        {
            y[i] = r8_lo;
        }

        for (i = 100; i < 200; i++)
        {
            x[i] = r8_hi;
        }

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref y, index: +100);

        typeMethods.r8vec_linspace(100, r8_hi, r8_lo, ref x, index: +200);
        for (i = 200; i < 300; i++)
        {
            y[i] = r8_hi;
        }

        for (i = 300; i < 400; i++)
        {
            x[i] = r8_lo;
        }

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref y, index: +300);

        Spiral.uvp_spiral(nu, rho, n, x, y, t, ref u, ref v, ref p);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  U:"
                          + "  " + typeMethods.r8vec_min(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  V:"
                          + "  " + typeMethods.r8vec_min(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  P:"
                          + "  " + typeMethods.r8vec_min(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void rhs_spiral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    rhs_spiral_test samples the Spiral right hand side.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("rhs_spiral_test");
        Console.WriteLine("  Spiral Flow:");
        Console.WriteLine("  Sample the Navier-Stokes right hand sides");
        Console.WriteLine("  at the initial time T = 0, using the unit square.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] f = new double[n];
        double[] g = new double[n];
        double[] h = new double[n];

        const double r8_lo = 0.0;
        const double r8_hi = 1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        const double t = 0.0;

        Spiral.rhs_spiral(nu, rho, n, x, y, t, ref f, ref g, ref h);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  F:"
                          + "  " + typeMethods.r8vec_min(n, f).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, f).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  G:"
                          + "  " + typeMethods.r8vec_min(n, g).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, g).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  H:"
                          + "  " + typeMethods.r8vec_min(n, h).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, h).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void resid_spiral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    resid_spiral_test samples the Spiral residual.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("resid_spiral_test");
        Console.WriteLine("  Spiral Flow:");
        Console.WriteLine("  Sample the Navier-Stokes residuals");
        Console.WriteLine("  at the initial time T = 0, on the unit square.");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] pr = new double[n];
        double[] ur = new double[n];
        double[] vr = new double[n];

        const double r8_lo = 0.0;
        const double r8_hi = 1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        const double t = 0.0;

        Spiral.resid_spiral(nu, rho, n, x, y, t, ref ur, ref vr, ref pr);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  Ur:"
                          + "  " + typeMethods.r8vec_amin(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Vr:"
                          + "  " + typeMethods.r8vec_amin(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Pr:"
                          + "  " + typeMethods.r8vec_amin(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void gnuplot_spiral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    gnuplot_spiral_test plots the Spiral flow.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int x_num = 21;
        const int y_num = 21;

        Console.WriteLine("");
        Console.WriteLine("gnuplot_spiral_test:");
        Console.WriteLine("  Spiral Flow:");
        Console.WriteLine("  Generate a velocity field on a regular grid.");
        Console.WriteLine("  Store in GNUPLOT data and command files.");

        const double x_lo = 0.0;
        const double x_hi = 1.0;

        const double y_lo = 0.0;
        const double y_hi = 1.0;

        double[] x = new double[x_num * y_num];
        double[] y = new double[x_num * y_num];

        NavierStokes2DExact.grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi, ref x, ref y);

        const double nu = 1.0;
        const double rho = 1.0;
        const int n = x_num * y_num;
        const double t = 0.0;

        double[] u = new double[x_num * y_num];
        double[] v = new double[x_num * y_num];
        double[] p = new double[x_num * y_num];

        Spiral.uvp_spiral(nu, rho, n, x, y, t, ref u, ref v, ref p);

        const string header = "spiral";
        const double s = 5.00;
        NavierStokes2DExact.ns2de_gnuplot(header, n, x, y, u, v, p, s);

    }

    private static void parameter_spiral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARAMETER_spiral_test monitors Spiral solution norms for NU, RHO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int k;
        double p_norm;
        double t;
        double u_norm;
        double v_norm;

        Console.WriteLine("");
        Console.WriteLine("PARAMETER_spiral_test");
        Console.WriteLine("  Spiral Flow:");
        Console.WriteLine("  Monitor solution norms over time for various");
        Console.WriteLine("  values of NU, RHO.");

        const int n = 1000;

        double[] u = new double[n];
        double[] v = new double[n];
        double[] p = new double[n];

        const double r8_lo = 0.0;
        const double r8_hi = 1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        //
        //  Vary RHO.
        //
        Console.WriteLine("");
        Console.WriteLine("  RHO affects the pressure scaling.");
        Console.WriteLine("");
        Console.WriteLine("     RHO         NU           T     ||U||       ||V||       ||P||");
        Console.WriteLine("");

        double nu = 1.0;
        double rho = 1.0;

        for (j = 1; j <= 3; j++)
        {
            for (k = 0; k <= 5; k++)
            {
                t = k / 5.0;

                Spiral.uvp_spiral(nu, rho, n, x, y, t, ref u, ref v, ref p);

                u_norm = typeMethods.r8vec_norm_l2(n, u) / n;
                v_norm = typeMethods.r8vec_norm_l2(n, v) / n;
                p_norm = typeMethods.r8vec_norm_l2(n, p) / n;

                Console.WriteLine("  " + rho.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + nu.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + v_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + p_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }

            Console.WriteLine("");
            rho /= 100.0;
        }

        /*
        Vary NU.
        */
        Console.WriteLine("");
        Console.WriteLine("  NU affects the time scaling.");
        Console.WriteLine("");
        Console.WriteLine("     RHO         NU           T     ||U||       ||V||       ||P||");
        Console.WriteLine("");

        nu = 1.0;
        rho = 1.0;

        for (i = 1; i <= 4; i++)
        {
            for (k = 0; k <= 5; k++)
            {
                t = k / 5.0;

                Spiral.uvp_spiral(nu, rho, n, x, y, t, ref u, ref v, ref p);

                u_norm = typeMethods.r8vec_norm_l2(n, u) / n;
                v_norm = typeMethods.r8vec_norm_l2(n, v) / n;
                p_norm = typeMethods.r8vec_norm_l2(n, p) / n;

                Console.WriteLine("  " + rho.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + nu.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + v_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + p_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }

            Console.WriteLine("");

            nu /= 10.0;
        }

    }

    private static void uvp_taylor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    uvp_taylor_test samples the Taylor solution at the initial time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("uvp_taylor_test");
        Console.WriteLine("  Taylor Flow:");
        Console.WriteLine("  Estimate the range of velocity and pressure");
        Console.WriteLine("  at the initial time T = 0, using a region that is");
        Console.WriteLine("  the square centered at (1.5,1.5) with 'radius' 1.0,");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];

        const double r8_lo = 0.5;
        const double r8_hi = +2.5;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        const double t = 0.0;

        Taylor.uvp_taylor(nu, rho, n, x, y, t, ref u, ref v, ref p);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  U:"
                          + "  " + typeMethods.r8vec_min(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  V:"
                          + "  " + typeMethods.r8vec_min(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  P:"
                          + "  " + typeMethods.r8vec_min(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void uvp_taylor_test2()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    uvp_taylor_test2 samples the Taylor solution on the boundary at the initial time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 March 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        const double nu = 1.0;
        const double rho = 1.0;
        const double t = 0.0;

        Console.WriteLine("");
        Console.WriteLine("uvp_taylor_test2");
        Console.WriteLine("  Taylor Flow:");
        Console.WriteLine("  Estimate the range of velocity and pressure");
        Console.WriteLine("  on the boundary");
        Console.WriteLine("  at the initial time T = 0, using a region that is");
        Console.WriteLine("  the square centered at (1.5,1.5) with 'radius' 1.0,");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 400;

        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];
        double[] x = new double[n];
        double[] y = new double[n];

        const double r8_lo = 0.5;
        const double r8_hi = +2.5;

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref x);
        for (i = 0; i < 100; i++)
        {
            y[i] = r8_lo;
        }

        for (i = 100; i < 200; i++)
        {
            x[i] = r8_hi;
        }

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref y, index: +100);

        typeMethods.r8vec_linspace(100, r8_hi, r8_lo, ref x, index: +200);
        for (i = 200; i < 300; i++)
        {
            y[i] = r8_hi;
        }

        for (i = 300; i < 400; i++)
        {
            x[i] = r8_lo;
        }

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref y, index: +300);

        Taylor.uvp_taylor(nu, rho, n, x, y, t, ref u, ref v, ref p);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  U:"
                          + "  " + typeMethods.r8vec_min(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  V:"
                          + "  " + typeMethods.r8vec_min(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  P:"
                          + "  " + typeMethods.r8vec_min(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void rhs_taylor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    rhs_taylor_test samples the Taylor right hand side at the initial time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("rhs_taylor_test");
        Console.WriteLine("  Taylor Flow:");
        Console.WriteLine("  Sample the Navier-Stokes right hand sides");
        Console.WriteLine("  at the initial time T = 0, using a region that is");
        Console.WriteLine("  the square centered at (1.5,1.5) with 'radius' 1.0,");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] f = new double[n];
        double[] g = new double[n];
        double[] h = new double[n];

        const double r8_lo = 0.5;
        const double r8_hi = +2.5;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        const double t = 0.0;

        Taylor.rhs_taylor(nu, rho, n, x, y, t, ref f, ref g, ref h);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  F:"
                          + "  " + typeMethods.r8vec_min(n, f).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, f).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  G:"
                          + "  " + typeMethods.r8vec_min(n, g).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, g).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  H:"
                          + "  " + typeMethods.r8vec_min(n, h).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, h).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void resid_taylor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    resid_taylor_test samples the Taylor residual at the initial time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("resid_taylor_test");
        Console.WriteLine("  Taylor Flow:");
        Console.WriteLine("  Sample the Navier-Stokes residuals");
        Console.WriteLine("  at the initial time T = 0, using a region that is");
        Console.WriteLine("  the square centered at (1.5,1.5) with 'radius' 1.0,");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] pr = new double[n];
        double[] ur = new double[n];
        if (ur == null)
        {
            throw new ArgumentNullException(nameof(ur));
        }

        double[] vr = new double[n];

        const double r8_lo = 0.5;
        const double r8_hi = +2.5;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        const double t = 0.0;

        Taylor.resid_taylor(nu, rho, n, x, y, t, ref ur, ref vr, ref pr);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  Ur:"
                          + "  " + typeMethods.r8vec_amin(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Vr:"
                          + "  " + typeMethods.r8vec_amin(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Pr:"
                          + "  " + typeMethods.r8vec_amin(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void gnuplot_taylor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    gnuplot_taylor_test plots the Taylor flow.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int x_num = 21;
        const int y_num = 21;

        Console.WriteLine("");
        Console.WriteLine("gnuplot_taylor_test:");
        Console.WriteLine("  Taylor Flow:");
        Console.WriteLine("  Generate a Taylor vortex velocity field on a regular grid.");
        Console.WriteLine("  Store in GNUPLOT data and command files.");

        const double x_lo = 0.5;
        const double x_hi = 2.5;

        const double y_lo = 0.5;
        const double y_hi = 2.5;

        double[] x = new double[x_num * y_num];
        double[] y = new double[x_num * y_num];

        NavierStokes2DExact.grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi, ref x, ref y);

        const double nu = 1.0;
        const double rho = 1.0;
        const int n = x_num * y_num;
        const double t = 0.0;

        double[] u = new double[x_num * y_num];
        double[] v = new double[x_num * y_num];
        double[] p = new double[x_num * y_num];

        Taylor.uvp_taylor(nu, rho, n, x, y, t, ref u, ref v, ref p);

        const string header = "taylor";
        const double s = 0.10;
        NavierStokes2DExact.ns2de_gnuplot(header, n, x, y, u, v, p, s);

    }

    private static void parameter_taylor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARAMETER_taylor_test monitors Taylor solution norms for various values of NU, RHO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int k;
        double p_norm;
        double t;
        double u_norm;
        double v_norm;

        Console.WriteLine("");
        Console.WriteLine("PARAMETER_taylor_test");
        Console.WriteLine("  Taylor Flow:");
        Console.WriteLine("  Monitor solution norms over time for various");
        Console.WriteLine("  values of NU, RHO.");

        const int n = 1000;

        double[] u = new double[n];
        double[] v = new double[n];
        double[] p = new double[n];

        const double r8_lo = 0.5;
        const double r8_hi = +2.5;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        //
        //  Vary RHO.
        //
        Console.WriteLine("");
        Console.WriteLine("  RHO affects the pressure scaling.");
        Console.WriteLine("");
        Console.WriteLine("     RHO         NU           T     ||U||       ||V||       ||P||");
        Console.WriteLine("");

        double nu = 1.0;
        double rho = 1.0;

        for (j = 1; j <= 3; j++)
        {
            for (k = 0; k <= 5; k++)
            {
                t = k / 5.0;

                Taylor.uvp_taylor(nu, rho, n, x, y, t, ref u, ref v, ref p);

                u_norm = typeMethods.r8vec_norm_l2(n, u) / n;
                v_norm = typeMethods.r8vec_norm_l2(n, v) / n;
                p_norm = typeMethods.r8vec_norm_l2(n, p) / n;

                Console.WriteLine("  " + rho.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + nu.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + v_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + p_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }

            Console.WriteLine("");
            rho /= 100.0;
        }

        /*
        Vary NU.
        */
        Console.WriteLine("");
        Console.WriteLine("  NU affects the time scaling.");
        Console.WriteLine("");
        Console.WriteLine("     RHO         NU           T     ||U||       ||V||       ||P||");
        Console.WriteLine("");

        nu = 1.0;
        rho = 1.0;

        for (i = 1; i <= 4; i++)
        {
            for (k = 0; k <= 5; k++)
            {
                t = k / 5.0;

                Taylor.uvp_taylor(nu, rho, n, x, y, t, ref u, ref v, ref p);

                u_norm = typeMethods.r8vec_norm_l2(n, u) / n;
                v_norm = typeMethods.r8vec_norm_l2(n, v) / n;
                p_norm = typeMethods.r8vec_norm_l2(n, p) / n;

                Console.WriteLine("  " + rho.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + nu.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + v_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + p_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }

            Console.WriteLine("");

            nu /= 10.0;
        }

    }

    private static void uvp_vortex_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    uvp_vortex_test samples the Vortex solution at the initial time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("uvp_vortex_test");
        Console.WriteLine("  Vortex Flow:");
        Console.WriteLine("  Estimate the range of velocity and pressure");
        Console.WriteLine("  at the initial time T = 0, using a region that is");
        Console.WriteLine("  the square centered at (1.0,1.0) with 'radius' 0.5,");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];

        const double r8_lo = 0.5;
        const double r8_hi = +1.5;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        const double t = 0.0;

        Vortex.uvp_vortex(nu, rho, n, x, y, t, ref u, ref v, ref p);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  U:"
                          + "  " + typeMethods.r8vec_min(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  V:"
                          + "  " + typeMethods.r8vec_min(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  P:"
                          + "  " + typeMethods.r8vec_min(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void uvp_vortex_test2()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    uvp_vortex_test2 samples the Vortex solution on the boundary at the initial time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        const double nu = 1.0;
        const double rho = 1.0;
        const double t = 0.0;

        Console.WriteLine("");
        Console.WriteLine("uvp_vortex_test2");
        Console.WriteLine("  Vortex Flow:");
        Console.WriteLine("  Estimate the range of velocity and pressure");
        Console.WriteLine("  on the boundary");
        Console.WriteLine("  at the initial time T = 0, using a region that is");
        Console.WriteLine("  the square centered at (1.0,1.0) with 'radius' 0.5,");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 400;

        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];
        double[] x = new double[n];
        double[] y = new double[n];

        const double r8_lo = 0.5;
        const double r8_hi = +1.5;

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref x);
        for (i = 0; i < 100; i++)
        {
            y[i] = r8_lo;
        }

        for (i = 100; i < 200; i++)
        {
            x[i] = r8_hi;
        }

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref y, index: +100);

        typeMethods.r8vec_linspace(100, r8_hi, r8_lo, ref x, index: +200);
        for (i = 200; i < 300; i++)
        {
            y[i] = r8_hi;
        }

        for (i = 300; i < 400; i++)
        {
            x[i] = r8_lo;
        }

        typeMethods.r8vec_linspace(100, r8_lo, r8_hi, ref y, index: +300);

        Vortex.uvp_vortex(nu, rho, n, x, y, t, ref u, ref v, ref p);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  U:"
                          + "  " + typeMethods.r8vec_min(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  V:"
                          + "  " + typeMethods.r8vec_min(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  P:"
                          + "  " + typeMethods.r8vec_min(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void rhs_vortex_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    rhs_vortex_test samples the Vortex right hand side at the initial time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("rhs_vortex_test");
        Console.WriteLine("  Vortex Flow:");
        Console.WriteLine("  Sample the Navier-Stokes right hand sides");
        Console.WriteLine("  at the initial time T = 0, using a region that is");
        Console.WriteLine("  the square centered at (1.0,1.0) with 'radius' 0.5,");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] f = new double[n];
        double[] g = new double[n];
        double[] h = new double[n];

        const double r8_lo = 0.5;
        const double r8_hi = +1.5;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        const double t = 0.0;

        Vortex.rhs_vortex(nu, rho, n, x, y, t, ref f, ref g, ref h);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  F:"
                          + "  " + typeMethods.r8vec_min(n, f).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, f).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  G:"
                          + "  " + typeMethods.r8vec_min(n, g).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, g).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  H:"
                          + "  " + typeMethods.r8vec_min(n, h).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_max(n, h).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void resid_vortex_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    resid_vortex_test samples the Vortex residual at the initial time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 1.0;
        const double rho = 1.0;

        Console.WriteLine("");
        Console.WriteLine("resid_vortex_test");
        Console.WriteLine("  Vortex Flow:");
        Console.WriteLine("  Sample the Navier-Stokes residuals");
        Console.WriteLine("  at the initial time T = 0, using a region that is");
        Console.WriteLine("  the square centered at (1.0,1.0) with 'radius' 0.5,");
        Console.WriteLine("  Kinematic viscosity NU = " + nu + "");
        Console.WriteLine("  Fluid density RHO = " + rho + "");

        const int n = 1000;

        double[] pr = new double[n];
        double[] ur = new double[n];
        double[] vr = new double[n];

        const double r8_lo = 0.5;
        const double r8_hi = +1.5;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        const double t = 0.0;

        Vortex.resid_vortex(nu, rho, n, x, y, t, ref ur, ref vr, ref pr);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  Ur:"
                          + "  " + typeMethods.r8vec_amin(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14).ToString(CultureInfo.InvariantCulture).PadLeft(14) +
                          "");
        Console.WriteLine("  Vr:"
                          + "  " + typeMethods.r8vec_amin(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14).ToString(CultureInfo.InvariantCulture).PadLeft(14) +
                          "");
        Console.WriteLine("  Pr:"
                          + "  " + typeMethods.r8vec_amin(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14).ToString(CultureInfo.InvariantCulture).PadLeft(14) +
                          "");

    }

    private static void gnuplot_vortex_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    gnuplot_vortex_test plots the Vortex flow.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int x_num = 21;
        const int y_num = 21;

        Console.WriteLine("");
        Console.WriteLine("gnuplot_vortex_test:");
        Console.WriteLine("  Vortex Flow:");
        Console.WriteLine("  Generate a Vortex vortex velocity field on a regular grid.");
        Console.WriteLine("  Store in GNUPLOT data and command files.");

        const double x_lo = 0.5;
        const double x_hi = 1.5;

        const double y_lo = 0.5;
        const double y_hi = 1.5;

        double[] x = new double[x_num * y_num];
        double[] y = new double[x_num * y_num];

        NavierStokes2DExact.grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi, ref x, ref y);

        const double nu = 1.0;
        const double rho = 1.0;
        const int n = x_num * y_num;
        const double t = 0.0;

        double[] u = new double[x_num * y_num];
        double[] v = new double[x_num * y_num];
        double[] p = new double[x_num * y_num];

        Vortex.uvp_vortex(nu, rho, n, x, y, t, ref u, ref v, ref p);

        const string header = "vortex";
        const double s = 0.10;
        NavierStokes2DExact.ns2de_gnuplot(header, n, x, y, u, v, p, s);

    }

    private static void parameter_vortex_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARAMETER_vortex_test monitors Vortex solution norms for various values of NU, RHO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int k;
        double p_norm;
        double t;
        double u_norm;
        double v_norm;

        Console.WriteLine("");
        Console.WriteLine("PARAMETER_vortex_test");
        Console.WriteLine("  Vortex Flow:");
        Console.WriteLine("  Monitor solution norms over time for various");
        Console.WriteLine("  values of NU, RHO.");

        const int n = 1000;

        double[] u = new double[n];
        double[] v = new double[n];
        double[] p = new double[n];

        const double r8_lo = 0.5;
        const double r8_hi = +1.5;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, r8_lo, r8_hi, ref seed);
        //
        //  Vary RHO.
        //
        Console.WriteLine("");
        Console.WriteLine("  RHO affects the pressure scaling.");
        Console.WriteLine("");
        Console.WriteLine("     RHO         NU           T     ||U||       ||V||       ||P||");
        Console.WriteLine("");

        double nu = 1.0;
        double rho = 1.0;

        for (j = 1; j <= 3; j++)
        {
            for (k = 0; k <= 5; k++)
            {
                t = k / 5.0;

                Vortex.uvp_vortex(nu, rho, n, x, y, t, ref u, ref v, ref p);

                u_norm = typeMethods.r8vec_norm_l2(n, u) / n;
                v_norm = typeMethods.r8vec_norm_l2(n, v) / n;
                p_norm = typeMethods.r8vec_norm_l2(n, p) / n;

                Console.WriteLine("  " + rho.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + nu.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + v_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + p_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }

            Console.WriteLine("");
            rho /= 100.0;
        }

        /*
        Vary NU.
        */
        Console.WriteLine("");
        Console.WriteLine("  NU affects the time scaling.");
        Console.WriteLine("");
        Console.WriteLine("     RHO         NU           T     ||U||       ||V||       ||P||");
        Console.WriteLine("");

        nu = 1.0;
        rho = 1.0;

        for (i = 1; i <= 4; i++)
        {
            for (k = 0; k <= 5; k++)
            {
                t = k / 5.0;

                Vortex.uvp_vortex(nu, rho, n, x, y, t, ref u, ref v, ref p);

                u_norm = typeMethods.r8vec_norm_l2(n, u) / n;
                v_norm = typeMethods.r8vec_norm_l2(n, v) / n;
                p_norm = typeMethods.r8vec_norm_l2(n, p) / n;

                Console.WriteLine("  " + rho.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + nu.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + v_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + p_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }

            Console.WriteLine("");

            nu /= 10.0;
        }
    }
}