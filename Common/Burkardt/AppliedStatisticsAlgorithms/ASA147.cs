using System;

namespace Burkardt.AppliedStatistics
{
    public static partial class Algorithms
    {
        public static double gammds(double x, double p, ref int ifault)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMDS computes the incomplete Gamma integral.
        //
        //  Discussion:
        //
        //    The parameters must be positive.  An infinite series is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 January 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Chi Leung Lau.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Chi Leung Lau,
        //    Algorithm AS 147:
        //    A Simple Series for the Incomplete Gamma Integral,
        //    Applied Statistics,
        //    Volume 29, Number 1, 1980, pages 113-114.
        //
        //  Parameters:
        //
        //    Input, double X, P, the arguments of the incomplete
        //    Gamma integral.  X and P must be greater than 0.
        //
        //    Output, int *IFAULT, error flag.
        //    0, no errors.
        //    1, X <= 0 or P <= 0.
        //    2, underflow during the computation.
        //
        //    Output, double GAMMDS, the value of the incomplete
        //    Gamma integral.
        //
        {
            double e = 1.0E-09;
            double uflo = 1.0E-37;
            double value;
            //
            //  Check the input.
            //
            if (x <= 0.0)
            {
                ifault = 1;
                value = 0.0;
                return value;
            }

            if (p <= 0.0)
            {
                ifault = 1;
                value = 0.0;
                return value;
            }

            //
            //  LGAMMA is the natural logarithm of the gamma function.
            //
            double arg = p * Math.Log(x) - Helpers.LogGamma(p + 1.0) - x;

            if (arg < Math.Log(uflo))
            {
                value = 0.0;
                ifault = 2;
                return value;
            }

            double f = Math.Exp(arg);

            if (f == 0.0)
            {
                value = 0.0;
                ifault = 2;
                return value;
            }

            ifault = 0;
            //
            //  Series begins.
            //
            double c = 1.0;
            value = 1.0;
            double a = p;

            for (;;)
            {
                a = a + 1.0;
                c = c * x / a;
                value = value + c;

                if (c <= e * value)
                {
                    break;
                }
            }

            value = value * f;

            return value;
        }
    }
}