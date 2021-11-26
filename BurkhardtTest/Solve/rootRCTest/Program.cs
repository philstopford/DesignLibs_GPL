using System;
using System.Globalization;
using Burkardt.SolveNS;

namespace rootRCTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ROOT_RC_TEST.
        //
        //  Discussion:
        //
        //    ROOT_RC_TEST tests the ROOT_RC library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 December 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double ferr = 0;
        int i;
        double[] q = new double[9];
        double xerr = 0;

        Console.WriteLine("");
        Console.WriteLine("ROOT_RC_TEST:");
        Console.WriteLine("  ROOT_RC searches for an");
        Console.WriteLine("  approximate solution of F(X) = 0, using reverse communication.");
        Console.WriteLine("");
        Console.WriteLine("       X              XERR            FX              FERR");
        Console.WriteLine("");
        //
        //  Initialization.
        //
        int it = 0;
        const int it_max = 30;
        for ( i = 0; i < 9; i++ )
        {
            q[i] = 0.0;
        }
        double x = - 2.1;
        //
        //  Each call takes one more step of improvement.
        //
        for ( ; ; )
        {
            double fx = Math.Cos ( x ) - x;

            switch (it)
            {
                case 0:
                    Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                           + "  " + "              "
                                           + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    break;
                default:
                    Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                           + "  " + xerr.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                           + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                           + "  " + ferr.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    break;
            }

            x = RootRC.root_rc ( x, fx, ref ferr, ref xerr, ref q );

            if ( ferr < 1.0E-08 )
            {
                Console.WriteLine("");
                Console.WriteLine("  Uncertainty in F(X) less than tolerance");
                break;
            }

            if ( xerr < 1.0E-08 )
            {
                Console.WriteLine("");
                Console.WriteLine("  Width of X interal less than tolerance.");
                break;
            }

            if ( it_max < it )
            {
                Console.WriteLine("");
                Console.WriteLine("  Too many iterations!'");
                break;
            }
            it += 1;     
        }

        Console.WriteLine("");
        Console.WriteLine("ROOT_RC_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}