using System;
using System.Globalization;
using Burkardt.NavierStokesNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace NavierStokes3DExactTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NS3DE_TEST tests the NS3DE library.
        //
        //  Location:
        //
        //    http://people.sc.fsu.edu/~jburkardt/cpp_src/navier_stokes_3d_exact/ns3de_test.cpp
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("NS3DE_TEST");
        Console.WriteLine("  Test the NS3DE library.");

        uvwp_burgers_test();
        resid_burgers_test();

        uvwp_ethier_test();
        resid_ethier_test();

        Console.WriteLine("");
        Console.WriteLine("NS3DE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void uvwp_burgers_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UVWP_BURGERS_TEST tests UVWP_BURGERS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 0.25;

        Console.WriteLine("");
        Console.WriteLine("UVWP_BURGERS_TEST");
        Console.WriteLine("  UVWP_BURGERS evaluates the Burgers solution.");
        Console.WriteLine("  Estimate the range of velocity and pressure");
        Console.WriteLine("  at the initial time T = 0, in a region that is the");
        Console.WriteLine("  cube centered at (0,0,0) with 'radius' 1.0.");
        Console.WriteLine("  Viscosity = " + nu + "");

        const int n = 1000;

        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];
        double[] w = new double[n];

        const double xyz_lo = -1.0;
        const double xyz_hi = +1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, xyz_lo, xyz_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, xyz_lo, xyz_hi, ref seed);
        double[] z = UniformRNG.r8vec_uniform_ab_new(n, xyz_lo, xyz_hi, ref seed);
        double[] t = typeMethods.r8vec_zeros_new(n);

        Burgers.uvwp_burgers(nu, n, x, y, z, t, ref u, ref v, ref w, ref p);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  U:  "
                          + "  " + typeMethods.r8vec_amin(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  V:  "
                          + "  " + typeMethods.r8vec_amin(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  W:  "
                          + "  " + typeMethods.r8vec_amin(n, w).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, w).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  P:  "
                          + "  " + typeMethods.r8vec_amin(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
    }

    private static void resid_burgers_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RESID_BURGERS_TEST tests RESID_BURGERS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double nu = 0.25;

        Console.WriteLine("");
        Console.WriteLine("RESID_BURGERS_TEST");
        Console.WriteLine("  RESID_BURGERS evaluates the Burgers residual.");
        Console.WriteLine("  Sample the Navier-Stokes residuals");
        Console.WriteLine("  at the initial time T = 0, using a region that is");
        Console.WriteLine("  the cube centered at (0,0,0) with 'radius' 1.0,");
        Console.WriteLine("  Viscosity = " + nu + "");

        const int n = 1000;

        double[] pr = new double[n];
        double[] ur = new double[n];
        double[] vr = new double[n];
        double[] wr = new double[n];

        const double xyz_lo = -1.0;
        const double xyz_hi = +1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, xyz_lo, xyz_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, xyz_lo, xyz_hi, ref seed);
        double[] z = UniformRNG.r8vec_uniform_ab_new(n, xyz_lo, xyz_hi, ref seed);
        double[] t = typeMethods.r8vec_zeros_new(n);

        Burgers.resid_burgers(nu, n, x, y, z, t, ref ur, ref vr, ref wr, ref pr);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  Ur:  "
                          + "  " + typeMethods.r8vec_amin(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Vr:  "
                          + "  " + typeMethods.r8vec_amin(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Wr:  "
                          + "  " + typeMethods.r8vec_amin(n, wr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, wr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Pr:  "
                          + "  " + typeMethods.r8vec_amin(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
    }

    private static void uvwp_ethier_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UVWP_ETHIER_TEST tests UVWP_ETHIER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = Math.PI / 4.0;
        double d = Math.PI / 2.0;

        Console.WriteLine("");
        Console.WriteLine("UVWP_ETHIER_TEST");
        Console.WriteLine("  UVWP_ETHIER evaluates the Ethier solution.");
        Console.WriteLine("  Estimate the range of velocity and pressure");
        Console.WriteLine("  at the initial time T = 0, in a region that is the");
        Console.WriteLine("  cube centered at (0,0,0) with 'radius' 1.0.");
        Console.WriteLine("  Parameter A = " + a + "");
        Console.WriteLine("  Parameter D = " + d + "");

        const int n = 1000;

        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];
        double[] w = new double[n];

        const double xyz_lo = -1.0;
        const double xyz_hi = +1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, xyz_lo, xyz_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, xyz_lo, xyz_hi, ref seed);
        double[] z = UniformRNG.r8vec_uniform_ab_new(n, xyz_lo, xyz_hi, ref seed);
        double[] t = typeMethods.r8vec_zeros_new(n);

        Ethier.uvwp_ethier(a, d, n, x, y, z, t, ref u, ref v, ref w, ref p);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  U:  "
                          + "  " + typeMethods.r8vec_amin(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, u).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  V:  "
                          + "  " + typeMethods.r8vec_amin(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, v).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  W:  "
                          + "  " + typeMethods.r8vec_amin(n, w).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, w).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  P:  "
                          + "  " + typeMethods.r8vec_amin(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, p).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
    }

    private static void resid_ethier_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RESID_ETHIER_TEST tests RESID_ETHIER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = Math.PI / 4.0;
        double d = Math.PI / 2.0;

        Console.WriteLine("");
        Console.WriteLine("RESID_ETHIER_TEST");
        Console.WriteLine("  RESID_ETHIER evaluates the Ethier residual.");
        Console.WriteLine("  Sample the Navier-Stokes residuals");
        Console.WriteLine("  at the initial time T = 0, using a region that is");
        Console.WriteLine("  the cube centered at (0,0,0) with 'radius' 1.0,");
        Console.WriteLine("  Parameter A = " + a + "");
        Console.WriteLine("  Parameter D = " + d + "");

        const int n = 1000;

        double[] pr = new double[n];
        double[] ur = new double[n];
        double[] vr = new double[n];
        double[] wr = new double[n];

        const double xyz_lo = -1.0;
        const double xyz_hi = +1.0;
        int seed = 123456789;

        double[] x = UniformRNG.r8vec_uniform_ab_new(n, xyz_lo, xyz_hi, ref seed);
        double[] y = UniformRNG.r8vec_uniform_ab_new(n, xyz_lo, xyz_hi, ref seed);
        double[] z = UniformRNG.r8vec_uniform_ab_new(n, xyz_lo, xyz_hi, ref seed);
        double[] t = typeMethods.r8vec_zeros_new(n);

        Ethier.resid_ethier(a, d, n, x, y, z, t, ref ur, ref vr, ref wr, ref pr);

        Console.WriteLine("");
        Console.WriteLine("           Minimum       Maximum");
        Console.WriteLine("");
        Console.WriteLine("  Ur:  "
                          + "  " + typeMethods.r8vec_amin(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, ur).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Vr:  "
                          + "  " + typeMethods.r8vec_amin(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, vr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Wr:  "
                          + "  " + typeMethods.r8vec_amin(n, wr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, wr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Pr:  "
                          + "  " + typeMethods.r8vec_amin(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + typeMethods.r8vec_amax(n, pr).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }
}