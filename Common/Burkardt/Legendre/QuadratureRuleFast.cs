﻿using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.ODENS;
using Burkardt.PolynomialNS;
using Burkardt.Types;

namespace Burkardt.Legendre;

using Polynomial = Polynomial;
public static class QuadratureRuleFast
{
    public static void legendre_handle(int n, double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_HANDLE computes the requested Gauss-Legendre rule and outputs it.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 October 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the rule.
        //
        //    Input, double A, B, the left and right endpoints of the integration
        //    interval.
        // 
    {
        double[] r = new double[2];
        double[] w = new double[n];
        double[] x = new double[n];

        r[0] = a;
        r[1] = b;
        //
        //  Compute the rule.
        //
        DateTime dt = DateTime.Now;
        legendre_compute_glr(n, ref x, ref w);
        double t = (DateTime.Now - dt).TotalSeconds;
        Console.WriteLine("");
        Console.WriteLine("  Elapsed time during computation was " + t + " seconds.");
        //
        //  Rescale the rule.
        //
        ClenshawCurtis.rescale(a, b, n, ref x, ref w);
        //
        //  Write the data to files.
        //
        string output_w = "leg_o" + n + "_w.txt";
        string output_x = "leg_o" + n + "_x.txt";
        string output_r = "leg_o" + n + "_r.txt";

        Console.WriteLine("");
        Console.WriteLine("  Weight file will be   \"" + output_w + "\".");
        Console.WriteLine("  Abscissa file will be \"" + output_x + "\".");
        Console.WriteLine("  Region file will be   \"" + output_r + "\".");

        typeMethods.r8mat_write(output_w, 1, n, w);
        typeMethods.r8mat_write(output_x, 1, n, x);
        typeMethods.r8mat_write(output_r, 1, 2, r);
    }

    public static void legendre_compute_glr(int n, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_COMPUTE_GLR: Legendre quadrature by the Glaser-Liu-Rokhlin method.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 October 2009
        //
        //  Author:
        //
        //    Original C++ version by Nick Hale.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
        //    A fast algorithm for the calculation of the roots of special functions, 
        //    SIAM Journal on Scientific Computing,
        //    Volume 29, Number 4, pages 1420-1438, 2007.
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        int i;
        double p = 0;
        double pp = 0;
        //
        //  Get the value and derivative of the N-th Legendre polynomial at 0.
        //
        legendre_compute_glr0(n, ref p, ref pp);
        switch (n % 2)
        {
            //
            //  If N is odd, then zero is a root.
            //  
            case 1:
                x[(n - 1) / 2] = p;
                w[(n - 1) / 2] = pp;
                break;
            //
            default:
                legendre_compute_glr2(p, n, ref x[n / 2], ref w[n / 2]);
                break;
        }

        //
        //  Get the complete set of roots and derivatives.
        //
        legendre_compute_glr1(n, ref x, ref w);
        //
        //  Compute the W.
        //
        for (i = 0; i < n; i++)
        {
            w[i] = 2.0 / (1.0 - x[i]) / (1.0 + x[i]) / w[i] / w[i];
        }

        double w_sum = 0.0;
        for (i = 0; i < n; i++)
        {
            w_sum += w[i];
        }

        for (i = 0; i < n; i++)
        {
            w[i] = 2.0 * w[i] / w_sum;
        }
    }

    public static void legendre_compute_glr0(int n, ref double p, ref double pp)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_COMPUTE_GLR0 gets a starting value for the fast algorithm.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 October 2009
        //
        //  Author:
        //
        //    Original C++ version by Nick Hale.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
        //    A fast algorithm for the calculation of the roots of special functions, 
        //    SIAM Journal on Scientific Computing,
        //    Volume 29, Number 4, pages 1420-1438, 2007.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the Legendre polynomial.
        //
        //    Output, double *P, *PP, the value of the N-th Legendre polynomial
        //    and its derivative at 0.
        //
    {
        int k;

        double pm2 = 0.0;
        double pm1 = 1.0;
        double ppm2 = 0.0;
        double ppm1 = 0.0;

        for (k = 0; k < n; k++)
        {
            double dk = k;
            p = -dk * pm2 / (dk + 1.0);
            pp = ((2.0 * dk + 1.0) * pm1 - dk * ppm2) / (dk + 1.0);
            pm2 = pm1;
            pm1 = p;
            ppm2 = ppm1;
            ppm1 = pp;
        }
    }

    public static void legendre_compute_glr1(int n, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_COMPUTE_GLR1 gets the complete set of Legendre points and weights.
        //
        //  Discussion:
        //
        //    This routine requires that a starting estimate be provided for one
        //    root and its derivative.  This information will be stored in entry
        //    (N+1)/2 if N is odd, or N/2 if N is even, of X and W.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 October 2009
        //
        //  Author:
        //
        //    Original C++ version by Nick Hale.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
        //    A fast algorithm for the calculation of the roots of special functions, 
        //    SIAM Journal on Scientific Computing,
        //    Volume 29, Number 4, pages 1420-1438, 2007.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the Legendre polynomial.
        //
        //    Input/output, double X[N].  On input, a starting value
        //    has been set in one entry.  On output, the roots of the Legendre 
        //    polynomial.
        //
        //    Input/output, double W[N].  On input, a starting value
        //    has been set in one entry.  On output, the derivatives of the Legendre 
        //    polynomial at the zeros.
        //
        //  Local Parameters:
        //
        //    Local, int M, the number of terms in the Taylor expansion.
        //
    {
        int j;
        int k;
        const int m = 30;
        int n2;
        int s;

        switch (n % 2)
        {
            case 1:
                n2 = (n - 1) / 2 - 1;
                s = 1;
                break;
            default:
                n2 = n / 2 - 1;
                s = 0;
                break;
        }

        double[] u = new double[m + 2];
        double[] up = new double[m + 1];

        double dn = n;

        for (j = n2 + 1; j < n - 1; j++)
        {
            double xp = x[j];

            double h = RungeKutta.rk2_leg(Math.PI / 2.0, -Math.PI / 2.0, xp, n) - xp;

            u[0] = 0.0;
            u[1] = 0.0;
            u[2] = w[j];

            up[0] = 0.0;
            up[1] = u[2];

            for (k = 0; k <= m - 2; k++)
            {
                double dk = k;

                u[k + 3] =
                (
                    2.0 * xp * (dk + 1.0) * u[k + 2]
                    + (dk * (dk + 1.0) - dn * (dn + 1.0)) * u[k + 1] / (dk + 1.0)
                ) / (1.0 - xp) / (1.0 + xp) / (dk + 2.0);

                up[k + 2] = (dk + 2.0) * u[k + 3];
            }

            int l;
            for (l = 0; l < 5; l++)
            {
                h -= Polynomial.ts_mult(u, h, m) / Polynomial.ts_mult(up, h, m - 1);
            }

            x[j + 1] = xp + h;
            w[j + 1] = Polynomial.ts_mult(up, h, m - 1);
        }

        for (k = 0; k <= n2 + s; k++)
        {
            x[k] = -x[n - 1 - k];
            w[k] = w[n - 1 - k];
        }
    }

    public static void legendre_compute_glr2(double pn0, int n, ref double x1, ref double d1)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_COMPUTE_GLR2 finds the first real root.
        //
        //  Discussion:
        //
        //    This function is only called if N is even.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 October 2009
        //
        //  Author:
        //
        //    Original C++ version by Nick Hale.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
        //    A fast algorithm for the calculation of the roots of special functions, 
        //    SIAM Journal on Scientific Computing,
        //    Volume 29, Number 4, pages 1420-1438, 2007.
        //
        //  Parameters:
        //
        //    Input, double PN0, the value of the N-th Legendre polynomial
        //    at 0.
        //
        //    Input, int N, the order of the Legendre polynomial.
        //
        //    Output, double *X1, the first real root.
        //
        //    Output, double *D1, the derivative at X1.
        //
        //  Local Parameters:
        //
        //    Local, int M, the number of terms in the Taylor expansion.
        //
    {
        int k;
        int l;
        const int m = 30;

        const double t = 0.0;
        x1 = RungeKutta.rk2_leg(t, -Math.PI / 2.0, 0.0, n);

        double[] u = new double[m + 2];
        double[] up = new double[m + 1];

        double dn = n;
        //
        //  U[0] and UP[0] are never used.
        //  U[M+1] is set, but not used, and UP[M] is set and not used.
        //  What gives?
        //
        u[0] = 0.0;
        u[1] = pn0;

        up[0] = 0.0;

        for (k = 0; k <= m - 2; k += 2)
        {
            double dk = k;

            u[k + 2] = 0.0;
            u[k + 3] = (dk * (dk + 1.0) - dn * (dn + 1.0)) * u[k + 1]
                       / (dk + 1.0) / (dk + 2.0);

            up[k + 1] = 0.0;
            up[k + 2] = (dk + 2.0) * u[k + 3];
        }

        for (l = 0; l < 5; l++)
        {
            x1 -= Polynomial.ts_mult(u, x1, m) / Polynomial.ts_mult(up, x1, m - 1);
        }

        d1 = Polynomial.ts_mult(up, x1, m - 1);

    }
}