﻿using System;
using System.Globalization;
using Burkardt.Lagrange;
using Burkardt.Types;
using Burkardt.Uniform;

namespace LagrangeInterpnDTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LAGRANGE_INTERP_ND_TEST.
        //
        //  Discussion:
        //
        //    LAGRANGE_INTERP_ND_TEST tests the LAGRANGE_INTERP_ND library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LAGRANGE_INTERP_ND_TEST:");
            
        Console.WriteLine("  Test the LAGRANGE_INTERP_ND library.");
        //
        //  Use the interface that passes in the orders directly.
        //
        test01();
        test02();
        test03();
        test04();
        //
        //  Repeat tests 1, 2, 3, and 4,
        //  using the interface that passes in the orders indirectly,
        //  based on the "level".
        //
        test05();
        test06();
        test07();
        test08();
        //
        //  Experiment with anisotropic orders.
        //
        test09();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("LAGRANGE_INTERP_ND_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 interpolates in 1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Interpolate in 1D, using orders.");
        Console.WriteLine("  LAGRANGE_INTERP_ND_GRID sets the interpolant.");
        Console.WriteLine("  LAGRANGE_INTERP_ND_VALUE evaluates it.");

        int m = 1;

        int[] n_1d = new int[m];
        for (i = 0; i < m; i++)
        {
            n_1d[i] = 5;
        }

        double[] a = new double[m];
        double[] b = new double[m];
        for (i = 0; i < m; i++)
        {
            a[i] = 0.0;
            b[i] = 1.0;
        }

        int nd = LagrangenD.lagrange_interp_nd_size(m, n_1d);
        double[] xd = LagrangenD.lagrange_interp_nd_grid(m, n_1d, a, b, nd);
        double[] zd = f_sinr(m, nd, xd);
        //
        //  Evaluate.
        //
        Console.WriteLine("");
        Console.WriteLine("         Zinterp          Zexact      Error");
        Console.WriteLine("");

        int ni = 5;
        int seed = 123456789;
        double[] xi = UniformRNG.r8mat_uniform_01_new(m, ni, ref seed);
        double[] ze = f_sinr(m, ni, xi);
        double[] zi = LagrangenD.lagrange_interp_nd_value(m, n_1d, a, b, nd, zd, ni, xi);

        for (j = 0; j < ni; j++)
        {
            Console.WriteLine("  " + zi[j].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + ze[j].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(zi[j] - ze[j]).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 interpolates in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  Interpolate in 2D, using orders.");
        Console.WriteLine("  LAGRANGE_INTERP_ND_GRID sets the interpolant.");
        Console.WriteLine("  LAGRANGE_INTERP_ND_VALUE evaluates it.");

        const int m = 2;

        int[] n_1d = new int[m];
        for (i = 0; i < m; i++)
        {
            n_1d[i] = 5;
        }

        double[] a = new double[m];
        double[] b = new double[m];
        for (i = 0; i < m; i++)
        {
            a[i] = 0.0;
            b[i] = 1.0;
        }

        int nd = LagrangenD.lagrange_interp_nd_size(m, n_1d);
        double[] xd = LagrangenD.lagrange_interp_nd_grid(m, n_1d, a, b, nd);
        double[] zd = f_sinr(m, nd, xd);
        //
        //  Evaluate.
        //
        Console.WriteLine("");
        Console.WriteLine("         Zinterp          Zexact      Error");
        Console.WriteLine("");

        int ni = 5;
        int seed = 123456789;
        double[] xi = UniformRNG.r8mat_uniform_01_new(m, ni, ref seed);
        double[] ze = f_sinr(m, ni, xi);
        double[] zi = LagrangenD.lagrange_interp_nd_value(m, n_1d, a, b, nd, zd, ni, xi);

        for (j = 0; j < ni; j++)
        {
            Console.WriteLine("  " + zi[j].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + ze[j].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(zi[j] - ze[j]).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 interpolates in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST03:");
        Console.WriteLine("  Interpolate in 3D, using orders.");
        Console.WriteLine("  LAGRANGE_INTERP_ND_GRID sets the interpolant.");
        Console.WriteLine("  LAGRANGE_INTERP_ND_VALUE evaluates it.");

        const int m = 3;

        int[] n_1d = new int[m];
        for (i = 0; i < m; i++)
        {
            n_1d[i] = 5;
        }

        double[] a = new double[m];
        double[] b = new double[m];
        for (i = 0; i < m; i++)
        {
            a[i] = 0.0;
            b[i] = 1.0;
        }

        int nd = LagrangenD.lagrange_interp_nd_size(m, n_1d);
        double[] xd = LagrangenD.lagrange_interp_nd_grid(m, n_1d, a, b, nd);
        double[] zd = f_sinr(m, nd, xd);
        //
        //  Evaluate.
        //
        Console.WriteLine("");
        Console.WriteLine("         Zinterp          Zexact      Error");
        Console.WriteLine("");

        int ni = 5;
        int seed = 123456789;
        double[] xi = UniformRNG.r8mat_uniform_01_new(m, ni, ref seed);
        double[] ze = f_sinr(m, ni, xi);
        double[] zi = LagrangenD.lagrange_interp_nd_value(m, n_1d, a, b, nd, zd, ni, xi);

        for (j = 0; j < ni; j++)
        {
            Console.WriteLine("  " + zi[j].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + ze[j].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(zi[j] - ze[j]).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 interpolates in 3D, using increasing resolution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int l;

        Console.WriteLine("");
        Console.WriteLine("TEST04:");
        Console.WriteLine("  Interpolate in 3D, using orders.");
        Console.WriteLine("  Use a sequence of increasing orders.");

        const int m = 3;

        int[] n_1d = new int[m];

        double[] a = new double[m];
        double[] b = new double[m];
        for (i = 0; i < m; i++)
        {
            a[i] = 0.0;
            b[i] = 1.0;
        }

        int ni = 20;
        int seed = 123456789;
        double[] xi = UniformRNG.r8mat_uniform_01_new(m, ni, ref seed);
        double[] ze = f_sinr(m, ni, xi);

        Console.WriteLine("");
        Console.WriteLine("  Level     Order   Average Error");
        Console.WriteLine("");

        for (l = 1; l <= 5; l++)
        {
            int order = l switch
            {
                0 => 1,
                _ => (int) Math.Pow(2, l) + 1
            };

            for (i = 0; i < m; i++)
            {
                n_1d[i] = order;
            }

            int nd = LagrangenD.lagrange_interp_nd_size(m, n_1d);
            double[] xd = LagrangenD.lagrange_interp_nd_grid(m, n_1d, a, b, nd);
            double[] zd = f_sinr(m, nd, xd);
            //
            //  Evaluate.
            //
            double[] zi = LagrangenD.lagrange_interp_nd_value(m, n_1d, a, b, nd, zd, ni, xi);

            double e = typeMethods.r8vec_norm_affine(ni, zi, ze) / ni;

            Console.WriteLine("  " + l.ToString(CultureInfo.InvariantCulture).PadLeft(5)
                                   + "  " + nd.ToString(CultureInfo.InvariantCulture).PadLeft(5)
                                   + "  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 repeats test 1 using levels.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST05:");
        Console.WriteLine("  Repeat test #1, using levels.");
        Console.WriteLine("  LAGRANGE_INTERP_ND_GRID2 sets the interpolant.");
        Console.WriteLine("  LAGRANGE_INTERP_ND_VALUE2 evaluates it.");

        const int m = 1;

        int[] ind = new int[m];
        for (i = 0; i < m; i++)
        {
            ind[i] = 2;
        }

        double[] a = new double[m];
        double[] b = new double[m];
        for (i = 0; i < m; i++)
        {
            a[i] = 0.0;
            b[i] = 1.0;
        }

        int nd = LagrangenD.lagrange_interp_nd_size2(m, ind);
        double[] xd = LagrangenD.lagrange_interp_nd_grid2(m, ind, a, b, nd);
        double[] zd = f_sinr(m, nd, xd);
        //
        //  Evaluate.
        //
        Console.WriteLine("");
        Console.WriteLine("         Zinterp          Zexact      Error");
        Console.WriteLine("");

        int ni = 5;
        int seed = 123456789;
        double[] xi = UniformRNG.r8mat_uniform_01_new(m, ni, ref seed);
        double[] ze = f_sinr(m, ni, xi);
        double[] zi = LagrangenD.lagrange_interp_nd_value2(m, ind, a, b, nd, zd, ni, xi);

        for (j = 0; j < ni; j++)
        {
            Console.WriteLine("  " + zi[j].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + ze[j].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(zi[j] - ze[j]).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 repeats test 2 using levels.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST06:");
        Console.WriteLine("  Repeat test #2, using levels.");
        Console.WriteLine("  LAGRANGE_INTERP_ND_GRID2 sets the interpolant.");
        Console.WriteLine("  LAGRANGE_INTERP_ND_VALUE2 evaluates it.");

        const int m = 2;

        int[] ind = new int[m];
        for (i = 0; i < m; i++)
        {
            ind[i] = 2;
        }

        double[] a = new double[m];
        double[] b = new double[m];
        for (i = 0; i < m; i++)
        {
            a[i] = 0.0;
            b[i] = 1.0;
        }

        int nd = LagrangenD.lagrange_interp_nd_size2(m, ind);
        double[] xd = LagrangenD.lagrange_interp_nd_grid2(m, ind, a, b, nd);
        double[] zd = f_sinr(m, nd, xd);
        //
        //  Evaluate.
        //
        Console.WriteLine("");
        Console.WriteLine("         Zinterp          Zexact      Error");
        Console.WriteLine("");

        int ni = 5;
        int seed = 123456789;
        double[] xi = UniformRNG.r8mat_uniform_01_new(m, ni, ref seed);
        double[] ze = f_sinr(m, ni, xi);
        double[] zi = LagrangenD.lagrange_interp_nd_value2(m, ind, a, b, nd, zd, ni, xi);

        for (j = 0; j < ni; j++)
        {
            Console.WriteLine("  " + zi[j].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + ze[j].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(zi[j] - ze[j]).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 repeats test 3 using levels.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST07:");
        Console.WriteLine("  Repeat test #3,  using levels.");
        Console.WriteLine("  LAGRANGE_INTERP_ND_GRID2 sets the interpolant.");
        Console.WriteLine("  LAGRANGE_INTERP_ND_VALUE2 evaluates it.");

        int m = 3;

        int[] ind = new int[m];
        for (i = 0; i < m; i++)
        {
            ind[i] = 2;
        }

        double[] a = new double[m];
        double[] b = new double[m];
        for (i = 0; i < m; i++)
        {
            a[i] = 0.0;
            b[i] = 1.0;
        }

        int nd = LagrangenD.lagrange_interp_nd_size2(m, ind);
        double[] xd = LagrangenD.lagrange_interp_nd_grid2(m, ind, a, b, nd);
        double[] zd = f_sinr(m, nd, xd);
        //
        //  Evaluate.
        //
        Console.WriteLine("");
        Console.WriteLine("         Zinterp          Zexact      Error");
        Console.WriteLine("");

        int ni = 5;
        int seed = 123456789;
        double[] xi = UniformRNG.r8mat_uniform_01_new(m, ni, ref seed);
        double[] ze = f_sinr(m, ni, xi);
        double[] zi = LagrangenD.lagrange_interp_nd_value2(m, ind, a, b, nd, zd, ni, xi);

        for (j = 0; j < ni; j++)
        {
            Console.WriteLine("  " + zi[j].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + ze[j].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(zi[j] - ze[j]).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 repeats test 4 using levels.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int l;

        Console.WriteLine("");
        Console.WriteLine("TEST08:");
        Console.WriteLine("  Interpolate in 3D, using levels.");
        Console.WriteLine("  Use a sequence of increasing levels.");

        const int m = 3;

        int[] ind = new int[m];

        double[] a = new double[m];
        double[] b = new double[m];
        for (i = 0; i < m; i++)
        {
            a[i] = 0.0;
            b[i] = 1.0;
        }

        int ni = 20;
        int seed = 123456789;
        double[] xi = UniformRNG.r8mat_uniform_01_new(m, ni, ref seed);
        double[] ze = f_sinr(m, ni, xi);

        Console.WriteLine("");
        Console.WriteLine("  Level     Order   Average Error");
        Console.WriteLine("");

        for (l = 0; l <= 5; l++)
        {
            for (i = 0; i < m; i++)
            {
                ind[i] = l;
            }

            int nd = LagrangenD.lagrange_interp_nd_size2(m, ind);
            double[] xd = LagrangenD.lagrange_interp_nd_grid2(m, ind, a, b, nd);
            double[] zd = f_sinr(m, nd, xd);
            //
            //  Evaluate.
            //
            double[] zi = LagrangenD.lagrange_interp_nd_value2(m, ind, a, b, nd, zd, ni, xi);

            double e = typeMethods.r8vec_norm_affine(ni, zi, ze) / ni;

            Console.WriteLine("  " + l.ToString(CultureInfo.InvariantCulture).PadLeft(5)
                                   + "  " + nd.ToString(CultureInfo.InvariantCulture).PadLeft(5)
                                   + "  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        }
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 interpolates in 3D, using anisotropic resolution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int l;

        Console.WriteLine("");
        Console.WriteLine("TEST09:");
        Console.WriteLine("  Interpolate in 3D, using orders.");
        Console.WriteLine("  Use a sequence of increasing orders.");
        Console.WriteLine("  Use anisotropic resolution.");
        Console.WriteLine("  The interpoland is a polynomial of degrees 3, 5, 2");
        Console.WriteLine("  so our orders need to be at least 4, 6, 3 to match it.");

        const int m = 3;

        int[] n_1d = new int[m];

        double[] a = new double[m];
        double[] b = new double[m];
        for (i = 0; i < m; i++)
        {
            a[i] = 0.0;
            b[i] = 1.0;
        }

        int ni = 20;
        int seed = 123456789;
        double[] xi = UniformRNG.r8mat_uniform_01_new(m, ni, ref seed);
        double[] ze = f_poly352(m, ni, xi);

        Console.WriteLine("");
        Console.WriteLine("  Level     Orders   Average Error");
        Console.WriteLine("");

        for (l = 0; l <= 10; l++)
        {
            switch (l)
            {
                case 0:
                    n_1d[0] = 1;
                    n_1d[1] = 1;
                    n_1d[2] = 1;
                    break;
                case 1:
                    n_1d[0] = 2;
                    n_1d[1] = 1;
                    n_1d[2] = 1;
                    break;
                case 2:
                    n_1d[0] = 1;
                    n_1d[1] = 2;
                    n_1d[2] = 1;
                    break;
                case 3:
                    n_1d[0] = 1;
                    n_1d[1] = 1;
                    n_1d[2] = 2;
                    break;
                case 4:
                    n_1d[0] = 4;
                    n_1d[1] = 2;
                    n_1d[2] = 2;
                    break;
                case 5:
                    n_1d[0] = 2;
                    n_1d[1] = 4;
                    n_1d[2] = 2;
                    break;
                case 6:
                    n_1d[0] = 2;
                    n_1d[1] = 2;
                    n_1d[2] = 4;
                    break;
                case 8:
                    n_1d[0] = 6;
                    n_1d[1] = 4;
                    n_1d[2] = 4;
                    break;
                case 9:
                    n_1d[0] = 4;
                    n_1d[1] = 6;
                    n_1d[2] = 4;
                    break;
                case 10:
                    n_1d[0] = 4;
                    n_1d[1] = 4;
                    n_1d[2] = 6;
                    break;
            }

            int nd = LagrangenD.lagrange_interp_nd_size(m, n_1d);
            double[] xd = LagrangenD.lagrange_interp_nd_grid(m, n_1d, a, b, nd);
            double[] zd = f_poly352(m, nd, xd);
            //
            //  Evaluate.
            //
            double[] zi = LagrangenD.lagrange_interp_nd_value(m, n_1d, a, b, nd, zd, ni, xi);

            double e = typeMethods.r8vec_norm_affine(ni, zi, ze) / ni;

            Console.WriteLine("  " + l.ToString(CultureInfo.InvariantCulture).PadLeft(5)
                                   + "  " + n_1d[0].ToString(CultureInfo.InvariantCulture).PadLeft(5)
                                   + "  " + n_1d[1].ToString(CultureInfo.InvariantCulture).PadLeft(5)
                                   + "  " + n_1d[2].ToString(CultureInfo.InvariantCulture).PadLeft(5)
                                   + "  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        }
    }

    private static double[] f_sinr(int m, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_SINR is a scalar function of an M-dimensional argument, to be interpolated.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[M*N], the points.
        //
        //    Output, double F_SINR[N], the value of the function at each point.
        //
    {
        int j;

        double[] z = new double[n];

        for (j = 0; j < n; j++)
        {
            double r = 0.0;
            int i;
            for (i = 0; i < m; i++)
            {
                r += Math.Pow(x[i + j * m], 2);
            }

            r = Math.Sqrt(r);
            z[j] = Math.Sin(r);
        }

        return z;
    }

    private static double[] f_poly352(int m, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_POLY253 is a scalar function of a 3-dimensional argument, to be interpolated.
        //
        //  Discussion:
        //
        //    The polynomial is of maximum degrees 3, 5, and 2, in the three variables.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[M*N], the points.
        //
        //    Output, double F_POLY253[N], the value of the function at each point.
        //
    {
        int j;

        double[] z = new double[n];

        for (j = 0; j < n; j++)
        {
            z[j] =
                1.0
                + Math.Pow(x[0 + j * m], 2) * Math.Pow(x[1 + j * m], 5) * Math.Pow(x[2 + j * m], 2)
                + x[0 + j * m] * Math.Pow(x[1 + j * m], 2) * x[2 + j * m]
                + Math.Pow(x[0 + j * m], 3);
        }

        return z;
    }
}