using System;
using Burkardt.Types;

namespace Burkardt.SolveNS;

public static class Theta
{
    public static double[] theta_solve(double a, double cl, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    THETA_SOLVE solves a pair of transcendental equations.
        //
        //  Discussion:
        //
        //    The vector THETA returned by this function is needed in order to define
        //    the terms in a Karhunen-Loeve expansion of a diffusion coefficient.
        //
        //    The two equations are:
        //
        //      1/CL - THETA * Math.Tan ( A * THETA ) = 0
        //      THETA - 1/CL * Math.Tan ( A * THETA ) = 0
        //
        //    A and CL are taken to be positive.  Over each open interval 
        //
        //      ( n - 1/2 pi, n + 1/2 Math.PI ) / A, for N = 0, 1, ...
        //
        //    the function Math.Tan ( A * THETA ) monotonically rises from -oo to +00; 
        //    therefore, it can be shown that there is one root of each equation 
        //    in every interval of this form.  Moreover, because of the positivity
        //    of A and CL, we can restrict our search to the interval 
        //
        //      [ n pi, n + 1/2 Math.PI ) / A, for N = 0, 1, ...
        //
        //    This function computes K such roots, starting in the first interval,
        //    finding those two roots, moving to the next interval, and so on, until
        //    the requested number of roots have been found.  Odd index roots will
        //    correspond to the first equation, and even index roots to the second.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Howard Elman, Darran Furnival,
        //    Solving the Stochastic Steady-State Diffusion Problem Using Multigrid,
        //    University of Maryland Department of Computer Science,
        //    Technical Report TR-4786.
        //
        //  Parameters:
        //
        //    Input, double A, the "radius" of the domain, D = (-A,A)x(-A,A).
        //    0 < A.
        //
        //    Input, double CL, the correlation length.
        //    0 < CL.
        //
        //    Input, int M, the number of values to compute.
        //
        //    Output, double THETA_SOLVE[M], the values of Theta.
        //
    {
        double bmatol;
        double eps;
        double fa;
        double fc;
        double ftol;
        int k;
        double[] theta;
        double xa;
        double xa_init;
        double xb;
        double xb_init;
        double xc = 0;

        theta = new double[m];
        for (k = 0; k < m; k++)
        {
            theta[k] = 0.0;
        }

        //
        //  [ XA_INIT, XB_INIT] = [ n * pi, n+1/2 Math.PI ] / a, n = 0, 1, 2, ...
        //
        xa_init = 0.0;
        xb_init = Math.PI / 2.0 / a;
        eps = typeMethods.r8_epsilon();

        k = 0;
        for (;;)
        {
            //
            //  Seek root of equation 1 in interval.
            //
            if (m <= k)
            {
                break;
            }

            k += 1;
            xa = xa_init;
            fa = 1.0 / cl - xa * Math.Tan(a * xa);
            ftol = eps * (Math.Abs(fa) + 1.0);
            xb = xb_init;
            fc = fa;
            bmatol = 100.0 * eps * (Math.Abs(xa) + Math.Abs(xb));

            while (bmatol < xb - xa)
            {
                xc = (xa + xb) / 2.0;
                fc = 1.0 / cl - xc * Math.Tan(a * xc);

                if (Math.Abs(fc) <= ftol)
                {
                    break;
                }

                switch (fc)
                {
                    case > 0.0:
                        xa = xc;
                        break;
                    default:
                        xb = xc;
                        break;
                }
            }

            theta[k - 1] = xc;
            //
            //  Seek root of equation 2 in interval.
            //
            if (m <= k)
            {
                break;
            }

            k += 1;
            switch (k)
            {
                //
                //  In the first interval, we need to skip the zero root of equation 2.
                //
                case 2:
                    k -= 1;
                    break;
                default:
                {
                    xa = xa_init;
                    fa = xa - Math.Tan(a * xa) / cl;
                    ftol = eps * (Math.Abs(fa) + 1.0);
                    xb = xb_init;

                    while (bmatol < xb - xa)
                    {
                        xc = (xa + xb) / 2.0;
                        fc = xc - Math.Tan(a * xc) / cl;

                        if (Math.Abs(fc) <= ftol)
                        {
                            break;
                        }

                        switch (fc)
                        {
                            case > 0.0:
                                xa = xc;
                                break;
                            default:
                                xb = xc;
                                break;
                        }
                    }

                    theta[k - 1] = xc;
                    break;
                }
            }

            //
            //  Advance the interval.
            //
            xa_init += Math.PI / a;
            xb_init += Math.PI / a;
        }

        return theta;
    }
}