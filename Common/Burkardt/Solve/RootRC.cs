using System;
using Burkardt.Types;

namespace Burkardt.SolveNS;

public static class RootRC
{
    public static double root_rc(double x, double fx, ref double ferr, ref double xerr, ref double[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROOT_RC solves a single nonlinear equation using reverse communication.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Gaston Gonnet.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        // 
        //    Gaston Gonnet,
        //    On the Structure of Zero Finders,
        //    BIT Numerical Mathematics,
        //    Volume 17, Number 2, June 1977, pages 170-183.
        //
        //  Parameters:
        //
        //    Input, double X, an estimate for the root.  On the first
        //    call, this must be a value chosen by the user.  Thereafter, it may
        //    be a value chosen by the user, or the value of ROOT returned on the
        //    previous call to the function.
        //
        //    Input, double FX, the value of the function at X.
        //
        //    Output, double &FERR, the smallest value of F encountered.
        //
        //    Output, double &XERR, the width of the change-in-sign interval,
        //    if one was encountered.
        //
        //    Input/output, double Q[9], storage needed by the function.
        //    Before the first call, the user must set Q(1) to 0.
        //
        //    Output, double ROOT_RC, an improved estimate for the root.
        //
    {
        int i;
        double xnew;
        switch (fx)
        {
            //
            //  If we found an exact zero, there is nothing more to do.
            //
            case 0.0:
                ferr = 0.0;
                xerr = 0.0;
                xnew = x;
                return xnew;
        }

        ferr = Math.Abs(fx);
        switch (q[0])
        {
            //
            //  If this is the first time, initialize, estimate the first root, and exit.
            //
            case 0.0:
            {
                q[0] = fx;
                q[1] = x;
                for (i = 2; i < 9; i++)
                {
                    q[i] = 0.0;
                }

                xnew = x + fx;
                xerr = typeMethods.r8_huge();
                return xnew;
            }
        }

        //
        //  This is not the first call.
        //
        q[8] += 1.0;
        switch (q[8])
        {
            //
            //  Check for too many iterations.
            //
            case > 80.0:
                Console.WriteLine("");
                Console.WriteLine("ROOT_RC - Fatal error!");
                Console.WriteLine("  Number of iterations = " + (int) q[8] + "");
                return 1;
        }

        //
        //  Check for a repeated X value.
        //
        if (2.0 <= q[8] && Math.Abs(x - q[3]) <= double.Epsilon || Math.Abs(x - q[1]) <= double.Epsilon)
        {
            Console.WriteLine("");
            Console.WriteLine("ROOT_RC - Fatal error!");
            Console.WriteLine("  Value of X has been input before.");
            return 1;
        }

        //
        //  Push X -> A -> B -> C
        //
        for (i = 5; 2 <= i; i--)
        {
            q[i] = q[i - 2];
        }

        q[0] = fx;
        q[1] = x;
        //
        //  If we have a change-in-sign interval, store the opposite value.
        //
        if (Math.Abs(typeMethods.r8_sign(q[0]) - typeMethods.r8_sign(q[2])) > double.Epsilon)
        {
            q[6] = q[2];
            q[7] = q[3];
        }

        //
        //  Calculate XERR.
        //
        xerr = q[6] != 0.0 ? Math.Abs(q[7] - q[1]) : typeMethods.r8_huge();

        switch (q[8])
        {
            //
            //  If more than 30 iterations, and we have change-in-sign interval, bisect.
            //
            case > 30.0 when q[6] != 0.0:
                xnew = q[1] + (q[7] - q[1]) / 2.0;
                return xnew;
        }

        double v = (q[2] - q[0]) / (q[3] - q[1]);
        //
        //  If 3 or more points, try Muller.
        //
        if (q[4] != 0.0)
        {
            double u = (q[4] - q[2]) / (q[5] - q[3]);
            double w = q[3] - q[1];
            double z = (q[5] - q[1]) / w;
            double r = (z + 1.0) * v - u;

            if (r != 0.0)
            {
                double p = 2.0 * z * q[0] / r;
                double d = 2.0 * p / (w * r) * (v - u);
                switch (d)
                {
                    case >= -1.0:
                    {
                        xnew = q[1] - p / (1.0 + Math.Sqrt(1.0 + d));
                        if (q[6] == 0.0 ||
                            q[1] < xnew && xnew < q[7] ||
                            q[7] < xnew && xnew < q[1])
                        {
                            return xnew;
                        }

                        break;
                    }
                }
            }
        }

        //
        //  Try the secant step.
        //
        if (Math.Abs(q[0] - q[2]) > double.Epsilon || q[6] == 0.0)
        {
            if (Math.Abs(q[0] - q[2]) <= double.Epsilon)
            {
                Console.WriteLine("");
                Console.WriteLine("ROOT_RC - Fatal error!");
                Console.WriteLine("  Cannot apply any method.");
                return 1;
            }

            double decr = q[0] / v;
            if (Math.Abs(decr) * 4.6E+18 < Math.Abs(q[1]))
            {
                decr = 1.74E-18 * Math.Abs(q[1]) * typeMethods.r8_sign(decr);
            }

            xnew = q[1] - decr;
            if (q[6] == 0.0 ||
                q[1] < xnew && xnew < q[7] ||
                q[7] < xnew && xnew < q[1])
            {
                return xnew;
            }
        }

        //
        //  Apply bisection.
        //
        xnew = q[1] + (q[7] - q[1]) / 2.0;

        return xnew;
    }
}