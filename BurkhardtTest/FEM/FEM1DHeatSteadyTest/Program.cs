using System;
using System.Globalization;
using Burkardt.FEM;
using Burkardt.Types;

namespace FEM1DHeatSteadyTest;

internal static class Program
{
    private static void Main()
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM1D_HEAT_STEADY_TEST.
//
//  Discussion:
//
//    FEM1D_HEAT_STEADY_TEST tests the FEM1D_HEAT_STEADY library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("FEM1D_BVP_LINEAR_TEST");
        Console.WriteLine("  Test the FEM1D_HEAT_STEADY library.");

        fem1d_heat_steady_test01();

        Console.WriteLine("");
        Console.WriteLine("FEM1D_HEAT_STEADY_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void fem1d_heat_steady_test01()

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_HEAT_STEADY_TEST01 carries out test case #1.
//
//  Discussion:
//
//    Use K1, F1, EXACT1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 2011
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        const int n = 11;

        Console.WriteLine("");
        Console.WriteLine("FEM1D_HEAT_STEADY_TEST01");
        Console.WriteLine("  K1(X)  = 1.0");
        Console.WriteLine("  F1(X)  = X * ( X + 3 ) * exp ( X )");
        Console.WriteLine("  U1(X)  = X * ( 1 - X ) * exp ( X )");
//
//  Geometry definitions.
//
        const double a = 0.0;
        const double b = 1.0;
        const double ua = 0.0;
        const double ub = 0.0;
        double[] x = typeMethods.r8vec_even_new(n, a, b);

        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        Console.WriteLine("  Left endpoint A = " + a + "");
        Console.WriteLine("  Right endpoint B = " + b + "");
        Console.WriteLine("  Prescribed U(A) = " + ua + "");
        Console.WriteLine("  Prescribed U(B) = " + ub + "");

        double[] u = FEM_1D_Heat_Steady.fem1d_heat_steady(n, a, b, ua, ub, FEM_Test_Methods.k1, FEM_Test_Methods.f1, x);

        Console.WriteLine("");
        Console.WriteLine("     I         X          U                Uexact      Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            double uexact = FEM_Test_Methods.exact1(x[i]);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }
}