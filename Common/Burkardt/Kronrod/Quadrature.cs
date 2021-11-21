using System;

namespace Burkardt.Kronrod;

public static class Quadrature
{
    public static void kronrod(int n, double eps, ref double[] x, ref double[] w1, ref double[] w2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KRONROD adds N+1 points to an N-point Gaussian rule.
        //
        //  Discussion:
        //
        //    This subroutine calculates the abscissas and weights of the 2N+1
        //    point Gauss Kronrod quadrature formula which is obtained from the 
        //    N point Gauss quadrature formula by the optimal addition of N+1 points.
        //
        //    The optimally added points are called Kronrod abscissas.  The 
        //    abscissas and weights for both the Gauss and Gauss Kronrod rules
        //    are calculated for integration over the interval [-1,+1].
        //
        //    Since the quadrature formula is symmetric with respect to the origin,
        //    only the nonnegative abscissas are calculated.
        //
        //    Note that the code published in Mathematics of Computation 
        //    omitted the definition of the variable which is here called COEF2.
        //
        //  Storage:
        //
        //    Given N, let M = ( N + 1 ) / 2.  
        //
        //    The Gauss-Kronrod rule will include 2*N+1 points.  However, by symmetry,
        //    only N + 1 of them need to be listed.
        //
        //    The arrays X, W1 and W2 contain the nonnegative abscissas in decreasing
        //    order, and the weights of each abscissa in the Gauss-Kronrod and
        //    Gauss rules respectively.  This means that about half the entries
        //    in W2 are zero.
        //
        //    For instance, if N = 3, the output is:
        //
        //    I      X               W1              W2
        //
        //    1    0.960491        0.104656         0.000000   
        //    2    0.774597        0.268488         0.555556    
        //    3    0.434244        0.401397         0.000000
        //    4    0.000000        0.450917         0.888889
        //
        //    and if N = 4, (notice that 0 is now a Kronrod abscissa)
        //    the output is
        //
        //    I      X               W1              W2
        //
        //    1    0.976560        0.062977        0.000000   
        //    2    0.861136        0.170054        0.347855    
        //    3    0.640286        0.266798        0.000000   
        //    4    0.339981        0.326949        0.652145    
        //    5    0.000000        0.346443        0.000000
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 September 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Robert Piessens, Maria Branders.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Robert Piessens, Maria Branders,
        //    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
        //    of Gauss and Lobatto,
        //    Mathematics of Computation,
        //    Volume 28, Number 125, January 1974, pages 135-139.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the Gauss rule.
        //
        //    Input, double EPS, the requested absolute accuracy of the
        //    abscissas.
        //
        //    Output, double X[N+1], the abscissas.
        //
        //    Output, double W1[N+1], the weights for the Gauss-Kronrod rule.
        //
        //    Output, double W2[N+1], the weights for 
        //    the Gauss rule.
        //
    {
        int i;
        int k;
        int l;

        double[] b = new double[(n + 1) / 2 + 1];
        double[] tau = new double[(n + 1) / 2];

        int m = (n + 1) / 2;
        bool even = 2 * m == n;

        double d = 2.0;
        double an = 0.0;
        for (k = 1; k <= n; k++)
        {
            an += 1.0;
            d = d * an / (an + 0.5);
        }

        //
        //  Calculation of the Chebyshev coefficients of the orthogonal polynomial.
        //
        tau[0] = (an + 2.0) / (an + an + 3.0);
        b[m - 1] = tau[0] - 1.0;
        double ak = an;

        for (l = 1; l < m; l++)
        {
            ak += 2.0;
            tau[l] = ((ak - 1.0) * ak
                      - an * (an + 1.0)) * (ak + 2.0) * tau[l - 1]
                     / (ak * ((ak + 3.0) * (ak + 2.0)
                              - an * (an + 1.0)));
            b[m - l - 1] = tau[l];

            int ll;
            for (ll = 1; ll <= l; ll++)
            {
                b[m - l - 1] += tau[ll - 1] * b[m - l + ll - 1];
            }
        }

        b[m] = 1.0;
        //
        //  Calculation of approximate values for the abscissas.
        //
        double bb = Math.Sin(1.570796 / (an + an + 1.0));
        double x1 = Math.Sqrt(1.0 - bb * bb);
        double s = 2.0 * bb * x1;
        double c = Math.Sqrt(1.0 - s * s);
        double coef = 1.0 - (1.0 - 1.0 / an) / (8.0 * an * an);
        double xx = coef * x1;
        //
        //  Coefficient needed for weights.
        //
        //  COEF2 = 2^(2*n+1) * n! * n! / (2n+1)! 
        //        = 2 * 4^n * n! / product( (n+1)*...*(2*n+1))
        //
        double coef2 = 2.0 / (2 * n + 1);
        for (i = 1; i <= n; i++)
        {
            coef2 = coef2 * 4.0 * i / (n + i);
        }

        //
        //  Calculation of the K-th abscissa (a Kronrod abscissa) and the
        //  corresponding weight.
        //
        for (k = 1; k <= n; k += 2)
        {
            Abscissa_Weight.abwe1(n, m, eps, coef2, even, b, ref xx, ref w1[+k - 1]);
            w2[k - 1] = 0.0;

            x[k - 1] = xx;
            double y = x1;
            x1 = y * c - bb * s;
            bb = y * s + bb * c;

            if (k == n)
            {
                xx = 0.0;
            }
            else
            {
                xx = coef * x1;
            }

            //
            //  Calculation of the K+1 abscissa (a Gaussian abscissa) and the
            //  corresponding weights.
            //
            Abscissa_Weight.abwe2(n, m, eps, coef2, even, b, ref xx, ref w1[k], ref w2[k]);

            x[k] = xx;
            y = x1;
            x1 = y * c - bb * s;
            bb = y * s + bb * c;
            xx = coef * x1;
        }

        switch (even)
        {
            //
            //  If N is even, we have one more Kronrod abscissa to compute,
            //  namely the origin.
            //
            case true:
                xx = 0.0;
                Abscissa_Weight.abwe1(n, m, eps, coef2, even, b, ref xx, ref w1[+n]);
                w2[n] = 0.0;
                x[n] = xx;
                break;
        }
    }

    public static void kronrod_adjust(double a, double b, int n, ref double[] x, ref double[] w1,
            ref double[] w2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KRONROD_ADJUST adjusts a Gauss-Kronrod rule from [-1,+1] to [A,B].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the endpoints of the new interval.
        //
        //    Input, int N, the order of the rule.
        //
        //    Input/output, double X[N+1], W1[N+1], W2[N+1], the abscissas
        //    and weights.
        //
    {
        int i;

        for (i = 0; i < n + 1; i++)
        {
            x[i] = ((1.0 - x[i]) * a
                    + (1.0 + x[i]) * b)
                   / 2.0;

            w1[i] = (b - a) / 2.0 * w1[i];
            w2[i] = (b - a) / 2.0 * w2[i];
        }

    }
}