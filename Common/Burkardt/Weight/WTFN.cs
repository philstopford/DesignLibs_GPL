using System;

namespace Burkardt.Weight;

public static class WTFN
{
    public static double[] wtfn(double[] t, int nt, int kind, double alpha, double beta )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WTFN evaluates the classical weight functions at given points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 February 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Sylvan Elhay, Jaroslav Kautsky,
        //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
        //    Interpolatory Quadrature,
        //    ACM Transactions on Mathematical Software,
        //    Volume 13, Number 4, December 1987, pages 399-415.
        //
        //  Parameters:
        //
        //    Input, double T[NT], the points where the weight function
        //    is to be evaluated.
        //
        //    Input, int NT, the number of evaluation points.
        //
        //    Input, int KIND, the rule.
        //    1, Legendre,             (a,b)       1.0
        //    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
        //    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
        //    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
        //    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
        //    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
        //    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
        //    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
        //    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
        //
        //    Input, double ALPHA, the value of Alpha, if needed.
        //
        //    Input, double BETA, the value of Beta, if needed.
        //
        //    Output, double WTFN[NT], the value of the weight function.
        //
    {
        int i;

        PARCHK.parchk(kind, 1, alpha, beta);

        double[] w = new double[nt];

        switch (kind)
        {
            case 1:
            {
                for (i = 0; i < nt; i++)
                {
                    w[i] = 1.0;
                }

                break;
            }
            case 2:
            {
                for (i = 0; i < nt; i++)
                {
                    w[i] = 1.0 / Math.Sqrt((1.0 - t[i]) * (1.0 + t[i]));
                }

                break;
            }
            case 3 when alpha == 0.0:
            {
                for (i = 0; i < nt; i++)
                {
                    w[i] = 1.0;
                }

                break;
            }
            case 3:
            {
                for (i = 0; i < nt; i++)
                {
                    w[i] = Math.Pow((1.0 - t[i]) * (1.0 + t[i]), alpha);
                }

                break;
            }
            case 4:
            {
                switch (alpha)
                {
                    case 0.0:
                    {
                        for (i = 0; i < nt; i++)
                        {
                            w[i] = 1.0;
                        }

                        break;
                    }
                    default:
                    {
                        for (i = 0; i < nt; i++)
                        {
                            w[i] = Math.Pow(1.0 - t[i], alpha);
                        }

                        break;
                    }
                }

                if (beta != 0.0)
                {
                    for (i = 0; i < nt; i++)
                    {
                        w[i] *= Math.Pow(1.0 + t[i], beta);
                    }
                }

                break;
            }
            case 5 when alpha == 0.0:
            {
                for (i = 0; i < nt; i++)
                {
                    w[i] = Math.Exp(-t[i]);
                }

                break;
            }
            case 5:
            {
                for (i = 0; i < nt; i++)
                {
                    w[i] = Math.Exp(-t[i]) * Math.Pow(t[i], alpha);
                }

                break;
            }
            case 6 when alpha == 0.0:
            {
                for (i = 0; i < nt; i++)
                {
                    w[i] = Math.Exp(-t[i] * t[i]);
                }

                break;
            }
            case 6:
            {
                for (i = 0; i < nt; i++)
                {
                    w[i] = Math.Exp(-t[i] * t[i]) * Math.Pow(Math.Abs(t[i]), alpha);
                }

                break;
            }
            case 7 when alpha != 0.0:
            {
                for (i = 0; i < nt; i++)
                {
                    w[i] = Math.Pow(Math.Abs(t[i]), alpha);
                }

                break;
            }
            case 7:
            {
                for (i = 0; i < nt; i++)
                {
                    w[i] = 1.0;
                }

                break;
            }
            case 8:
            {
                switch (alpha)
                {
                    case 0.0:
                    {
                        for (i = 0; i < nt; i++)
                        {
                            w[i] = 1.0;
                        }

                        break;
                    }
                    default:
                    {
                        for (i = 0; i < nt; i++)
                        {
                            w[i] = Math.Pow(t[i], alpha);
                        }

                        break;
                    }
                }

                if (beta != 0.0)
                {
                    for (i = 0; i < nt; i++)
                    {
                        w[i] *= Math.Pow(1.0 + t[i], beta);
                    }
                }

                break;
            }
            case 9:
            {
                for (i = 0; i < nt; i++)
                {
                    w[i] = Math.Sqrt((1.0 - t[i]) * (1.0 + t[i]));
                }

                break;
            }
        }

        return w;
    }


}