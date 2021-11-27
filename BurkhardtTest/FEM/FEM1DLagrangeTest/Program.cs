﻿using System;
using System.Globalization;
using Burkardt.FEM;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace FEM1DLagrangeTest;

internal static class Program
{
    private static void Main()
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM1D_LAGRANGE_TEST.
//
//  Discussion:
//
//    FEM1D_LAGRANGE_TEST tests FEM1D_LAGRANGE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2014
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("FEM1D_LAGRANGE_TEST");
        Console.WriteLine("  Test the FEM1D_LAGRANGE library.");

        legendre_set_test();
        lagrange_value_test();
        lagrange_derivative_test();

        int x_num = 11;
        int q_num = 5;
        fem1d_lagrange_stiffness_test(x_num, q_num);

        x_num = 11;
        q_num = 10;
        fem1d_lagrange_stiffness_test(x_num, q_num);

        Console.WriteLine("");
        Console.WriteLine("FEM1D_LAGRANGE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void legendre_set_test()
//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_TEST tests LEGENDRE_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2014
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_SET_TEST");
        Console.WriteLine("  LEGENDRE_SET returns points and weights of");
        Console.WriteLine("  Gauss-Legendre quadrature rules.");
        Console.WriteLine("");
        Console.WriteLine("   N               1             X^4           Runge");
        Console.WriteLine("");

        for (int n = 1; n <= 10; n++)
        {
            double[] x = new double[n];
            double[] w = new double[n];

            LegendreQuadrature.legendre_set(n, ref x, ref w);
            double e1 = 0.0;
            double e2 = 0.0;
            double e3 = 0.0;
            for (int i = 0; i < n; i++)
            {
                e1 += w[i];
                e2 += w[i] * Math.Pow(x[i], 4);
                e3 += w[i] / (1.0 + 25.0 * x[i] * x[i]);
            }

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + e1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + e2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + e3.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void lagrange_value_test()
//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_VALUE_TEST tests LAGRANGE_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2014
//
//  Author:
//
//    John Burkardt
//
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("LAGRANGE_VALUE_TEST");
        Console.WriteLine("  LAGRANGE_VALUE evaluates the Lagrange basis polynomials.");

        const int nd = 5;
        const double xlo = 0.0;
        const double xhi = nd - 1;
        double[] xd = typeMethods.r8vec_linspace_new(nd, xlo, xhi);

        typeMethods.r8vec_print(nd, xd, "  Lagrange basis points:");
//
//  Evaluate the polynomials.
//
        Console.WriteLine("");
        Console.WriteLine("   I      X          L1(X)       L2(X)       L3(X)" +
                          "       L4(X)       L5(X)");
        Console.WriteLine("");

        int ni = 2 * nd - 1;
        double[] xi = typeMethods.r8vec_linspace_new(ni, xlo, xhi);

        double[] li = FEM_1D_Lagrange.lagrange_value(nd, xd, ni, xi);

        for (i = 0; i < ni; i++)
        {
            string cout = "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                               + "  " + xi[i].ToString(CultureInfo.InvariantCulture).PadLeft(10);
            int j;
            for (j = 0; j < nd; j++)
            {
                cout += "  " + li[i + j * ni].ToString(CultureInfo.InvariantCulture).PadLeft(10);
            }

            Console.WriteLine(cout);
        }
    }

    private static void lagrange_derivative_test()
//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_DERIVATIVE_TEST tests LAGRANGE_DERIVATIVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2014
//
//  Author:
//
//    John Burkardt
//
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("LAGRANGE_DERIVATIVE_TEST");
        Console.WriteLine("  LAGRANGE_DERIVATIVE evaluates the Lagrange basis derivative.");

        const int nd = 5;
        const double xlo = 0.0;
        const double xhi = nd - 1;
        double[] xd = typeMethods.r8vec_linspace_new(nd, xlo, xhi);

        typeMethods.r8vec_print(nd, xd, "  Lagrange basis points:");
//
//  Evaluate the polynomials.
//
        Console.WriteLine("");
        Console.WriteLine("   I      X         L1'(X)      L2'(X)      L3'(X)" +
                          "      L4'(X)      L5'(X)");
        Console.WriteLine("");

        int ni = 2 * nd - 1;
        double[] xi = typeMethods.r8vec_linspace_new(ni, xlo, xhi);
        double[] lpi = FEM_1D_Lagrange.lagrange_derivative(nd, xd, ni, xi);

        for (i = 0; i < ni; i++)
        {
            string cout = "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                               + "  " +xi[i].ToString(CultureInfo.InvariantCulture).PadLeft(10);
            int j;
            for (j = 0; j < nd; j++)
            {
                cout += lpi[i + j * ni].ToString(CultureInfo.InvariantCulture).PadLeft(10);
            }

            Console.WriteLine(cout);
        }
    }

    private static void fem1d_lagrange_stiffness_test(int x_num, int q_num)
//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_LAGRANGE_STIFFNESS_TEST tests FEM1D_LAGRANGE_STIFFNESS.
//
//  Discussion:
//
//    The results are very sensitive to the quadrature rule.
//
//    In particular, if X_NUM points are used, the mass matrix will
//    involve integrals of polynomials of degree 2*(X_NUM-1), so the
//    quadrature rule should use at least Q_NUM = X_NUM - 1 points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 November 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X_NUM, the number of nodes.
//
//    Input, int Q_NUM, the number of quadrature points.
//
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("FEM1D_LAGRANGE_STIFFNESS_TEST");
        Console.WriteLine("  FEM1D_LAGRANGE_STIFFNESS computes the stiffness matrix,");
        Console.WriteLine("  the mass matrix, and right hand side vector for a");
        Console.WriteLine("  finite element problem using Lagrange interpolation");
        Console.WriteLine("  basis polynomials.");
        Console.WriteLine("");
        Console.WriteLine("  Solving:");
        Console.WriteLine("    -u''+u=x on 0 < x < 1");
        Console.WriteLine("    u(0) = u(1) = 0");
        Console.WriteLine("  Exact solution:");
        Console.WriteLine("    u(x) = x - sinh(x)/sinh(1)");
        Console.WriteLine("");
        Console.WriteLine("  Number of mesh points = " + x_num + "");
        Console.WriteLine("  Number of quadrature points = " + q_num + "");

        const double x_lo = 0.0;
        const double x_hi = 1.0;
        double[] x = typeMethods.r8vec_linspace_new(x_num, x_lo, x_hi);

        double[] a = new double[x_num * x_num];
        double[] m = new double[x_num * x_num];
        double[] b = new double[x_num];

        FEM_1D_Lagrange.fem1d_lagrange_stiffness(x_num, x, q_num, FEM_Test_Methods.f, ref a, ref m, ref b);

        double[] k = new double[x_num * x_num];

        for (j = 0; j < x_num; j++)
        {
            for (i = 0; i < x_num; i++)
            {
                k[i + j * x_num] = a[i + j * x_num] + m[i + j * x_num];
            }
        }

        for (j = 0; j < x_num; j++)
        {
            k[0 + j * x_num] = 0.0;
        }

        k[0 + 0 * x_num] = 1.0;
        b[0] = 0.0;

        for (j = 0; j < x_num; j++)
        {
            k[x_num - 1 + j * x_num] = 0.0;
        }

        k[x_num - 1 + (x_num - 1) * x_num] = 1.0;
        b[x_num - 1] = 0.0;

        double[] u = typeMethods.r8mat_fs_new(x_num, k, b);

        double[] u_e = FEM_Test_Methods.exact(x_num, x);

        Console.WriteLine("");
        Console.WriteLine("   I      X             U              U(exact)         Error");
        Console.WriteLine("");

        for (i = 0; i < x_num; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + u_e[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(u[i] - u_e[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }
}