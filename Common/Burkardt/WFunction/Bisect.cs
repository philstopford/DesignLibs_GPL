﻿using System;
using Burkardt.Function;

namespace Burkardt.WFunction;

public static class Bisect
{
    public class BisectData
    {
        public int nbits;
        public Crude.CrudeData data = new();
    }
        
    public static double bisect(ref BisectData data, double xx, int nb, ref int ner, int l )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BISECT approximates the W function using bisection.
        //
        //  Discussion:
        //
        //    The parameter TOL, which determines the accuracy of the bisection
        //    method, is calculated using NBITS (assuming the final bit is lost
        //    due to rounding error).
        //
        //    N0 is the maximum number of iterations used in the bisection
        //    method.
        //
        //    For XX close to 0 for Wp, the exponential approximation is used.
        //    The approximation is exact to O(XX^8) so, depending on the value
        //    of NBITS, the range of application of this formula varies. Outside
        //    this range, the usual bisection method is used.
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
        //    Input, double XX, the argument.
        //
        //    Input, int NB, indicates the branch of the W function.
        //    0, the upper branch;
        //    nonzero, the lower branch.
        //
        //    Output, int &NER, the error flag.
        //    0, success;
        //    1, the routine did not converge.  Perhaps reduce NBITS and try again.
        //
        //    Input, int L, the offset indicator.
        //    1, XX represents the offset of the argument from -Math.Exp(-1).
        //    not 1, XX is the actual argument.
        //
        //    Output, double BISECT, the value of W(X), as determined
        //
    {
        double d;
        double f;
        double fd;
        int i;
        const int n0 = 500;
        double r;
        double tol;
        double u;

        double value = 0.0;
        ner = 0;

        data.nbits = data.nbits switch
        {
            0 => NBITS.nbits_compute(),
            _ => data.nbits
        };

        double x = l switch
        {
            1 => xx - Math.Exp(-1.0),
            _ => xx
        };

        switch (nb)
        {
            case 0:
            {
                double test = 1.0 / Math.Pow(Math.Pow(2.0, data.nbits), 1.0 / 7.0);

                if (Math.Abs(x) < test)
                {
                    value = x
                            * Math.Exp(-x
                                       * Math.Exp(-x
                                                  * Math.Exp(-x
                                                             * Math.Exp(-x
                                                                        * Math.Exp(-x
                                                                            * Math.Exp(-x))))));

                    return value;
                }

                u = Crude.crude(ref data.data, x, nb) + 1.0E-03;
                tol = Math.Abs(u) / Math.Pow(2.0, data.nbits);
                d = Math.Max(u - 2.0E-03, -1.0);

                for (i = 1; i <= n0; i++)
                {
                    r = 0.5 * (u - d);
                    value = d + r;
                    //
                    //  Find root using w*Math.Exp(w)-x to avoid ln(0) error.
                    //
                    if (x < Math.Exp(1.0))
                    {
                        f = value * Math.Exp(value) - x;
                        fd = d * Math.Exp(d) - x;
                    }
                    //
                    //  Find root using ln(w/x)+w to avoid overflow error.
                    //
                    else
                    {
                        f = Math.Log(value / x) + value;
                        fd = Math.Log(d / x) + d;
                    }

                    switch (f)
                    {
                        case 0.0:
                            return value;
                    }

                    if (Math.Abs(r) <= tol)
                    {
                        return value;
                    }

                    switch (fd * f)
                    {
                        case > 0.0:
                            d = value;
                            break;
                        default:
                            u = value;
                            break;
                    }
                }

                break;
            }
            default:
            {
                d = Crude.crude(ref data.data, x, nb) - 1.0E-03;
                u = Math.Min(d + 2.0E-03, -1.0);
                tol = Math.Abs(u) / Math.Pow(2.0, data.nbits);

                for (i = 1; i <= n0; i++)
                {
                    r = 0.5 * (u - d);
                    value = d + r;
                    f = value * Math.Exp(value) - x;

                    switch (f)
                    {
                        case 0.0:
                            return value;
                    }

                    if (Math.Abs(r) <= tol)
                    {
                        return value;
                    }

                    fd = d * Math.Exp(d) - x;

                    switch (fd * f)
                    {
                        case > 0.0:
                            d = value;
                            break;
                        default:
                            u = value;
                            break;
                    }
                }

                break;
            }
        }

        //
        //  The iteration did not converge.
        //
        ner = 1;

        return value;
    }
}