using System;

namespace Burkardt.AppliedStatistics
{
    public static partial class Algorithms
    {
        public static double gamain(double x, double p, ref int ifault)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMAIN computes the incomplete gamma ratio.
        //
        //  Discussion:
        //
        //    A series expansion is used if P > X or X <= 1.  Otherwise, a
        //    continued fraction approximation is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 June 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by G Bhattacharjee.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    G Bhattacharjee,
        //    Algorithm AS 32:
        //    The Incomplete Gamma Integral,
        //    Applied Statistics,
        //    Volume 19, Number 3, 1970, pages 285-287.
        //
        //  Parameters:
        //
        //    Input, double X, P, the parameters of the incomplete 
        //    gamma ratio.  0 <= X, and 0 < P.
        //
        //    Output, int *IFAULT, error flag.
        //    0, no errors.
        //    1, P <= 0.
        //    2, X < 0.
        //    3, underflow.
        //    4, error return from the Log Gamma routine.
        //
        //    Output, double GAMAIN, the value of the incomplete gamma ratio.
        //
        {
            double acu = 1.0E-08;
            double gin;
            double oflo = 1.0E+37;
            double[] pn = new double[6];
            double rn;
            double term;
            double uflo = 1.0E-37;
            double value;

            ifault = 0;
            //
            //  Check the input.
            //
            if (p <= 0.0)
            {
                ifault = 1;
                value = 0.0;
                return value;
            }

            if (x < 0.0)
            {
                ifault = 2;
                value = 0.0;
                return value;
            }

            if (x == 0.0)
            {
                ifault = 0;
                value = 0.0;
                return value;
            }

            double g = Helpers.LogGamma(p);

            double arg = p * Math.Log(x) - x - g;

            if (arg < Math.Log(uflo))
            {
                ifault = 3;
                value = 0.0;
                return value;
            }

            ifault = 0;
            double factor = Math.Exp(arg);
            //
            //  Calculation by series expansion.
            //
            if (x <= 1.0 || x < p)
            {
                gin = 1.0;
                term = 1.0;
                rn = p;

                for (;;)
                {
                    rn = rn + 1.0;
                    term = term * x / rn;
                    gin = gin + term;

                    if (term <= acu)
                    {
                        break;
                    }
                }

                value = gin * factor / p;
                return value;
            }

            //
            //  Calculation by continued fraction.
            //
            double a = 1.0 - p;
            double b = a + x + 1.0;
            term = 0.0;

            pn[0] = 1.0;
            pn[1] = x;
            pn[2] = x + 1.0;
            pn[3] = x * b;

            gin = pn[2] / pn[3];

            for (;;)
            {
                a = a + 1.0;
                b = b + 2.0;
                term = term + 1.0;
                double an = a * term;
                for (int i = 0; i <= 1; i++)
                {
                    pn[i + 4] = b * pn[i + 2] - an * pn[i];
                }

                if (pn[5] != 0.0)
                {
                    rn = pn[4] / pn[5];
                    double dif = Math.Abs(gin - rn);
                    //
                    //  Absolute error tolerance satisfied?
                    //
                    if (dif <= acu)
                    {
                        //
                        //  Relative error tolerance satisfied?
                        //
                        if (dif <= acu * rn)
                        {
                            value = 1.0 - factor * gin;
                            break;
                        }
                    }

                    gin = rn;
                }

                for (int i = 0; i < 4; i++)
                {
                    pn[i] = pn[i + 2];
                }

                if (oflo <= Math.Abs(pn[4]))
                {
                    for (int i = 0; i < 4; i++)
                    {
                        pn[i] = pn[i] / oflo;
                    }
                }
            }

            return value;
        }

    }
}