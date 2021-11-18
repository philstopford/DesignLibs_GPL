using System;
using Burkardt.SolveNS;

namespace rootsRCTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ROOTS_RC_TEST.
        //
        //  Discussion:
        //
        //    ROOTS_RC_TEST tests the ROOTS_RC library.
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
        double[] fx;
        int i;
        int it;
        int it_max = 30;
        int j;
        int n = 4;
        double[] q;
        double[] x;
        double[] xnew;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("ROOTS_RC_TEST:");

        Console.WriteLine("  ROOTS_RC seeks a solution of");
        Console.WriteLine("  the N-dimensional nonlinear system F(X) = 0.");

        fx = new double[n];
        q = new double[(2 * n + 2) * (n + 2)];
        x = new double[n];
        xnew = new double[n];

        Console.WriteLine("");
        Console.WriteLine("       FERR          X");
        Console.WriteLine("");
        //
        //  Initialization.
        //
        for (j = 0; j < n + 2; j++)
        {
            for (i = 0; i < 2 * n + 2; i++)
            {
                q[i + j * (2 * n + 2)] = 0.0;
            }
        }

        xnew[0] = 1.2;
        for (i = 1; i < n; i++)
        {
            xnew[i] = 1.0;
        }

        it = 0;

        for (;;)
        {
            for (i = 0; i < n; i++)
            {
                x[i] = xnew[i];
            }

            fx[0] = 1.0 - x[0];
            for (i = 1; i < n; i++)
            {
                fx[i] = 10.0 * (x[i] - x[i - 1] * x[i - 1]);
            }

            cout = it switch
            {
                0 => "                ",
                _ => "  " + ferr.ToString(CultureInfo.InvariantCulture).PadLeft(14)
            };

            for (i = 0; i < n; i++)
            {
                cout += "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
            cout = "";

            RootsRC.roots_rc(n, x, fx, ref ferr, ref xnew, ref q);

            if (ferr < 1.0E-07)
            {
                Console.WriteLine("");
                Console.WriteLine("  Sum of |f(x)| less than tolerance.");
                break;
            }

            if (it_max < it)
            {
                Console.WriteLine("");
                Console.WriteLine("  Too many iterations!");
                break;
            }

            it += 1;
        }

        Console.WriteLine("");
        Console.WriteLine("ROOTS_RC_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}