using System;
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
        double a;
        double b;
        int i;
        int n = 11;
        double[] u;
        double ua;
        double ub;
        double uexact;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("FEM1D_HEAT_STEADY_TEST01");
        Console.WriteLine("  K1(X)  = 1.0");
        Console.WriteLine("  F1(X)  = X * ( X + 3 ) * exp ( X )");
        Console.WriteLine("  U1(X)  = X * ( 1 - X ) * exp ( X )");
//
//  Geometry definitions.
//
        a = 0.0;
        b = 1.0;
        ua = 0.0;
        ub = 0.0;
        x = typeMethods.r8vec_even_new(n, a, b);

        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        Console.WriteLine("  Left endpoint A = " + a + "");
        Console.WriteLine("  Right endpoint B = " + b + "");
        Console.WriteLine("  Prescribed U(A) = " + ua + "");
        Console.WriteLine("  Prescribed U(B) = " + ub + "");

        u = FEM_1D_Heat_Steady.fem1d_heat_steady(n, a, b, ua, ub, FEM_Test_Methods.k1, FEM_Test_Methods.f1, x);

        Console.WriteLine("");
        Console.WriteLine("     I         X          U                Uexact      Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            uexact = FEM_Test_Methods.exact1(x[i]);
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i].ToString().PadLeft(8)
                                   + "  " + u[i].ToString().PadLeft(14)
                                   + "  " + uexact.ToString().PadLeft(14)
                                   + "  " + Math.Abs(u[i] - uexact).ToString().PadLeft(14) + "");
        }
    }
}