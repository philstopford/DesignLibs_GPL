using System;
using Burkardt.Types;

namespace Burkardt.Function;

public static class CompassSearch
{
    public static double[] compass_search(Func<int, double[], double> function_handle, int m,
            double[] x0, double delta_tol, double delta_init, int k_max, ref double fx,
            ref int k)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMPASS_SEARCH carries out a direct search minimization algorithm.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Tamara Kolda, Robert Michael Lewis, Virginia Torczon,
        //    Optimization by Direct Search: New Perspectives on Some Classical 
        //    and Modern Methods,
        //    SIAM Review,
        //    Volume 45, Number 3, 2003, pages 385-482. 
        //
        //  Parameters:
        //
        //    Input, double FUNCTION_HANDLE ( int m, double x[] ), the name of
        //    a function which evaluates the function to be minimized.
        //
        //    Input, int M, the number of variables.
        //
        //    Input, double X0[M], a starting estimate for the minimizer.
        //
        //    Input, double DELTA_TOL, the smallest step size that is allowed.
        //
        //    Input, double DELTA_INIT, the starting stepsize.  
        //
        //    Input, int K_MAX, the maximum number of steps allowed.
        //
        //    Output, double COMPASS_SEARCH[M], the estimated minimizer.
        //
        //    Output, double &FX, the function value at X.
        //
        //    Output, int &K, the number of steps taken.
        //
    {
        k = 0;
        double[] x = new double[m];
        double[] xd = new double[m];
        typeMethods.r8vec_copy(m, x0, ref x);
        fx = function_handle(m, x);

        switch (delta_tol)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("COMPASS_SEARCH - Fatal error!");
                Console.WriteLine("  DELTA_TOL <= 0.0.");
                Console.WriteLine("  DELTA_TOL = " + delta_tol + "");
                return null;
        }

        if (delta_init <= delta_tol)
        {
            Console.WriteLine("");
            Console.WriteLine("COMPASS_SEARCH - Fatal error!");
            Console.WriteLine("  DELTA_INIT < DELTA_TOL.");
            Console.WriteLine("  DELTA_INIT = " + delta_init + "");
            Console.WriteLine("  DELTA_TOL = " + delta_tol + "");
            return null;
        }

        double delta = delta_init;

        while (k < k_max)
        {
            k += 1;
            //
            //  For each coordinate direction I, seek a lower function value
            //  by increasing or decreasing X(I) by DELTA.
            //
            bool decrease = false;
            double s = +1.0;
            int i = 0;

            int ii;
            for (ii = 1; ii <= 2 * m; ii++)
            {
                typeMethods.r8vec_copy(m, x, ref xd);
                xd[i] += s * delta;
                double fxd = function_handle(m, xd);
                //
                //  As soon as a decrease is noticed, accept the new point.
                //
                if (fxd < fx)
                {
                    typeMethods.r8vec_copy(m, xd, ref x);
                    fx = fxd;
                    decrease = true;
                    break;
                }

                s = -s;
                switch (s)
                {
                    case +1.0:
                        i += 1;
                        break;
                }
            }

            //
            //  If no decrease occurred, reduce DELTA.
            //
            if (!decrease)
            {
                delta /= 2.0;
                if (delta < delta_tol)
                {
                    break;
                }
            }
        }
        return x;
    }
}