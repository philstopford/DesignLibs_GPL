using System;
using Burkardt.HyperGeometry.HypersphereNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace HypersphereTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for HYPERSPHERE_PROPERTIES_TEST.
        //
        //  Discussion:
        //
        //    HYPERSPHERE_PROPERTIES_TEST tests the HYPERSPHERE_PROPERTIES library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("HYPERSPHERE_PROPERTIES_TEST:");
        Console.WriteLine("  Test the HYPERSPHERE_PROPERTIES library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();

        Console.WriteLine("");
        Console.WriteLine("HYPERSPHERE_PROPERTIES_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests the coordinate conversion routines.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c;
        double err;
        int m;
        int n;
        double[] r;
        int seed;
        int test;
        double[] theta;
        double[] x;
        double[] x2;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Test the coordinate conversion routines:");
        Console.WriteLine("  CARTESIAN_TO_HYPERSPHERE: X       -> R,Theta");
        Console.WriteLine("  HYPERSPHERE_TO_CARTESIAN: R,Theta -> X.");
        Console.WriteLine("");
        Console.WriteLine("  Pick a random X, and compute X2 by converting X");
        Console.WriteLine("  to hypersphere and back.  Consider norm of difference.");
        Console.WriteLine("");
        Console.WriteLine("  M    || X - X2 ||");

        seed = 123456789;

        n = 1;
        r = new double[n];

        for (m = 1; m <= 5; m++)
        {
            Console.WriteLine("");

            theta = new double[(m - 1) * n];

            for (test = 1; test <= 5; test++)
            {
                x = UniformRNG.r8mat_uniform_01_new(m, n, ref seed);
                c = UniformRNG.r8vec_uniform_01_new(m, ref seed);
                Hypersphere.cartesian_to_hypersphere(m, n, c, x, ref r, ref theta);
                x2 = Hypersphere.hypersphere_to_cartesian(m, n, c, r, theta);
                err = typeMethods.r8mat_norm_fro_affine(m, n, x, x2);
                Console.WriteLine("  " + m.ToString().PadLeft(2)
                                       + "  " + err.ToString().PadLeft(14) + "");
            }
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests HYPERSPHERE_01_SURFACE_UNIFORM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m;
        int n;
        int seed;
        int test;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  HYPERSPHERE_01_SURFACE_UNIFORM samples uniformly from the");
        Console.WriteLine("  surface of the unit hypersphere");

        seed = 123456789;
        typeMethods.r8vecNormalData data = new();

        n = 1;
        for (m = 1; m <= 5; m++)
        {
            for (test = 1; test <= 3; test++)
            {
                x = Hypersphere.hypersphere_01_surface_uniform(m, n, ref data, ref seed);
                typeMethods.r8vec_transpose_print(m, x, "  Random hypersphere point:");
            }
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests HYPERSPHERE_01_AREA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double area = 0;
        double area2;
        int m = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("TEST03:");
        Console.WriteLine("  HYPERSPHERE_01_AREA evaluates the area of the unit");
        Console.WriteLine("  hypersphere in M dimensions.");
        Console.WriteLine("");
        Console.WriteLine("       M      Exact       Computed");
        Console.WriteLine("              Area        Area");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Hypersphere.hypersphere_01_area_values(ref n_data, ref m, ref area);

            if (n_data == 0)
            {
                break;
            }

            area2 = Hypersphere.hypersphere_01_area(m);

            Console.WriteLine("  " + m.ToString().PadLeft(6)
                                   + "  " + area.ToString().PadLeft(10)
                                   + "  " + area2.ToString().PadLeft(10) + "");
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests HYPERSPHERE_01_VOLUME.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m = 0;
        int n_data;
        double volume = 0;
        double volume2;

        Console.WriteLine("");
        Console.WriteLine("TEST04:");
        Console.WriteLine("  HYPERSPHERE_01_VOLUME evaluates the area of the unit");
        Console.WriteLine("  hypersphere in M dimensions.");
        Console.WriteLine("  HYPERSPHERE_01_VOLUME_VALUES returns some test values.");
        Console.WriteLine("");
        Console.WriteLine("       M      Exact       Computed");
        Console.WriteLine("              Volume      Volume");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Hypersphere.hypersphere_01_volume_values(ref n_data, ref m, ref volume);

            if (n_data == 0)
            {
                break;
            }

            volume2 = Hypersphere.hypersphere_01_volume(m);

            Console.WriteLine("  " + m.ToString().PadLeft(6)
                                   + "  " + volume.ToString().PadLeft(10)
                                   + "  " + volume2.ToString().PadLeft(10) + "");
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests HYPERSPHERE_AREA, HYPERSPHERE_VOLUME.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double area;
        int m;
        double r;
        double volume;

        r = 1.5;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  For a hypersphere in M dimensions:");
        Console.WriteLine("  HYPERSPHERE_AREA computes the area");
        Console.WriteLine("  HYPERSPHERE_VOLUME computes the volume.");
        Console.WriteLine("");
        Console.WriteLine("  Notice that both quantities eventually decrease");
        Console.WriteLine("");
        Console.WriteLine("  We use a radius of R = " + r + "");
        Console.WriteLine("");
        Console.WriteLine("    M        Area          Volume    Area / Volume ");
        Console.WriteLine("");

        for (m = 1; m <= 20; m++)
        {
            area = Hypersphere.hypersphere_area(m, r);
            volume = Hypersphere.hypersphere_volume(m, r);
            Console.WriteLine("  " + m.ToString().PadLeft(3)
                                   + "  " + area.ToString().PadLeft(14)
                                   + "  " + volume.ToString().PadLeft(14)
                                   + "  " + (area / volume).ToString().PadLeft(14) + "");
        }
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests the stereographic mapping.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double err;
        int m;
        int n;
        int seed;
        int test;
        double[] x1;
        double[] x2;
        double[] x3;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  Test the stereographic mapping:");
        Console.WriteLine("  HYPERSPHERE_STEREOGRAPH maps hypersphere points to the plane.");
        Console.WriteLine("  HYPERSPHERE_STEREOGRAPH_INVERSE inverts the mapping.");
        Console.WriteLine("");
        Console.WriteLine("  Pick a random X1 on the hypersphere.");
        Console.WriteLine("  Map it to a point X2 on the plane.");
        Console.WriteLine("  Map it back to a point X3 on the hypersphere.");
        Console.WriteLine("  Consider norm of difference.");
        Console.WriteLine("");
        Console.WriteLine("  M    || X1 - X3 ||");

        seed = 123456789;
        typeMethods.r8vecNormalData data = new();

        n = 1;
        for (m = 2; m <= 5; m++)
        {
            Console.WriteLine("");
            for (test = 1; test <= 5; test++)
            {
                x1 = Hypersphere.hypersphere_01_surface_uniform(m, n, ref data, ref seed);
                x2 = Hypersphere.hypersphere_stereograph(m, n, x1);
                x3 = Hypersphere.hypersphere_stereograph_inverse(m, n, x2);
                err = typeMethods.r8mat_norm_fro_affine(m, n, x1, x3);
                Console.WriteLine("  " + m.ToString().PadLeft(2)
                                       + "  " + err.ToString().PadLeft(14) + "");
            }
        }
    }
}