using System;
using Burkardt.Function;

namespace Burkardt.WFunction
{
    public class WAPRData
    {
        public double an3 = 0;
        public double an4 = 0;
        public double an5 = 0;
        public double an6 = 0;
        public double c13 = 0;
        public double c23 = 0;
        public double d12 = 0;
        public double em = 0;
        public double em2 = 0;
        public double em9 = 0;
        public int init = 0;
        public int nbits = 0;
        public int niter = 1;
        public double s2 = 0;
        public double s21 = 0;
        public double s22 = 0;
        public double s23 = 0;
        public double tb = 0;
        public double x0 = 0;
        public double x1 = 0;
    }

    public static class WAPR
    {
        public static double wapr(ref WAPRData data, double x, int nb, ref int nerror, int l )

        //****************************************************************************80
        //
        // WAPR approximates the W function.
        //
        //  Discussion:
        //
        //    The call will fail if the input value X is out of range.
        //    The range requirement for the upper branch is:
        //      -Math.Exp(-1) <= X.
        //    The range requirement for the lower branch is:
        //      -Math.Exp(-1) < X < 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 June 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
        //    Patricia Culligan-Hensley.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
        //    Algorithm 743: WAPR - A Fortran routine for calculating real 
        //    values of the W-function,
        //    ACM Transactions on Mathematical Software,
        //    Volume 21, Number 2, June 1995, pages 172-181.
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Input, int NB, indicates the desired branch.
        //    * 0, the upper branch;
        //    * nonzero, the lower branch.
        //
        //    Output, int &NERROR, the error flag.
        //    * 0, successful call.
        //    * 1, failure, the input X is out of range.
        //
        //    Input, int L, indicates the interpretation of X.
        //    * 1, X is actually the offset from -(Math.Exp-1), so compute W(X-Math.Exp(-1)).
        //    * not 1, X is the argument; compute W(X);
        //
        //    Output, double WAPR, the approximate value of W(X).
        //
        {
            double an2;
            double delx;
            double eta;
            int i;
            double reta;
            double t;
            double temp;
            double temp2;
            double ts;
            double value;
            double xx;
            double zl;
            double zn;

            value = 0.0;
            nerror = 0;

            if (data.init == 0)
            {
                data.init = 1;

                data.nbits = NBITS.nbits_compute();

                if (56 <= data.nbits)
                {
                    data.niter = 2;
                }

                //
                //  Various mathematical constants.
                //
                data.em = -Math.Exp(-1.0);
                data.em9 = -Math.Exp(-9.0);
                data.c13 = 1.0 / 3.0;
                data.c23 = 2.0 * data.c13;
                data.em2 = 2.0 / data.em;
                data.d12 = -data.em2;
                data.tb = Math.Pow(0.5, data.nbits);
                data.x0 = Math.Pow(data.tb, 1.0 / 6.0) * 0.5;
                data.x1 = (1.0 - 17.0 * Math.Pow(data.tb, 2.0 / 7.0)) * data.em;
                data.an3 = 8.0 / 3.0;
                data.an4 = 135.0 / 83.0;
                data.an5 = 166.0 / 39.0;
                data.an6 = 3167.0 / 3549.0;
                data.s2 = Math.Sqrt(2.0);
                data.s21 = 2.0 * data.s2 - 3.0;
                data.s22 = 4.0 - 3.0 * data.s2;
                data.s23 = data.s2 - 2.0;
            }

            if (l == 1)
            {
                delx = x;

                if (delx < 0.0)
                {
                    nerror = 1;
                    Console.WriteLine("");
                    Console.WriteLine("WAPR - Fatal error!");
                    Console.WriteLine("  The offset X is negative.");
                    Console.WriteLine("  It must be nonnegative.");
                    return 1;
                }

                xx = x + data.em;
            }
            else
            {
                if (x < data.em)
                {
                    nerror = 1;
                    return value;
                }
                else if (x == data.em)
                {
                    value = -1.0;
                    return value;
                }

                xx = x;
                delx = xx - data.em;
            }

            //
            //  Calculations for Wp.
            //
            if (nb == 0)
            {
                if (Math.Abs(xx) <= data.x0)
                {
                    value = xx / (1.0 + xx / (1.0 + xx
                        / (2.0 + xx / (0.6 + 0.34 * xx))));
                    return value;
                }
                else if (xx <= data.x1)
                {
                    reta = Math.Sqrt(data.d12 * delx);
                    value = reta / (1.0 + reta / (3.0 + reta / (reta
                                / (data.an4 + reta / (reta * data.an6 + data.an5)) + data.an3)))
                            - 1.0;
                    return value;
                }
                else if (xx <= 20.0)
                {
                    reta = data.s2 * Math.Sqrt(1.0 - xx / data.em);
                    an2 = 4.612634277343749 * Math.Sqrt(Math.Sqrt(reta +
                                                        1.09556884765625));
                    value = reta / (1.0 + reta / (3.0 + (data.s21 * an2
                                                         + data.s22) * reta / (data.s23 * (an2 + reta)))) - 1.0;
                }
                else
                {
                    zl = Math.Log(xx);
                    value = Math.Log(xx / Math.Log(xx
                                         / Math.Pow(zl, Math.Exp(-1.124491989777808 /
                                                       (0.4225028202459761 + zl)))));
                }
            }
            //
            //  Calculations for Wm.
            //
            else
            {
                if (0.0 <= xx)
                {
                    nerror = 1;
                    return value;
                }
                else if (xx <= data.x1)
                {
                    reta = Math.Sqrt(data.d12 * delx);
                    value = reta / (reta / (3.0 + reta / (reta / (data.an4
                                                                  + reta / (reta * data.an6 - data.an5)) - data.an3)) - 1.0) - 1.0;
                    return value;
                }
                else if (xx <= data.em9)
                {
                    zl = Math.Log(-xx);
                    t = -1.0 - zl;
                    ts = Math.Sqrt(t);
                    value = zl - (2.0 * ts) / (data.s2 + (data.c13 - t
                        / (270.0 + ts * 127.0471381349219)) * ts);
                }
                else
                {
                    zl = Math.Log(-xx);
                    eta = 2.0 - data.em2 * xx;
                    value = Math.Log(xx / Math.Log(-xx / ((1.0
                                                 - 0.5043921323068457 * (zl + 1.0))
                        * (Math.Sqrt(eta) + eta / 3.0) + 1.0)));
                }

            }

            for (i = 1; i <= data.niter; i++)
            {
                zn = Math.Log(xx / value) - value;
                temp = 1.0 + value;
                temp2 = temp + data.c23 * zn;
                temp2 = 2.0 * temp * temp2;
                value = value * (1.0 + (zn / temp) * (temp2 - zn)
                    / (temp2 - 2.0 * zn));
            }

            return value;
        }
    }
}