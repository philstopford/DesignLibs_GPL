using System;

namespace Burkardt.Interpolation;

public static class LagrangeInterpolation
{
    public static double[] lagrange_rule(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_RULE computes the weights of a Lagrange interpolation rule.
        //
        //  Discussion:
        //
        //    Given N abscissas X, an arbitrary function F(X) can be
        //    interpolated by a polynomial P(X) of order N (and degree N-1)
        //    using weights that depend only on X.
        //
        //    Standard Lagrange interpolation can be rewritten into this form,
        //    which is more economical than evaluating the individual Lagrange
        //    basis polynomials.
        //
        //    If we define
        //
        //      W(I) = 1 / product ( 1 <= J <= N, J /= I ) ( X(J) - X(I) )
        //
        //    then
        //
        //      P(XV) = sum ( 1 <= I <= N ) W(I) * F( X(I) ) / ( XV - X(I) )
        //            / sum ( 1 <= I <= N ) W(I)             / ( XV - X(I) )
        //
        //    except when XV = X(J), for some J, when we set:
        //
        //      P(X(J)) = F(X(J))
        //
        //  Modified:
        //
        //    24 May 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jean-Paul Berrut, Lloyd Trefethen,
        //    Barycentric Lagrange Interpolation,
        //    SIAM Review,
        //    Volume 46, Number 3, September 2004, pages 501-517.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the rule.
        //
        //    Input, double X[N], the abscissas of the rule.
        //
        //    Output, double LAGRANGE_RULE[N], the weights of the rule.
        //
    {
        int i;

        double[] w = new double[n];

        for (i = 0; i < n; i++)
        {
            w[i] = 1.0;
        }

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                if (i != j)
                {
                    w[j] /= x[i] - x[j];
                }
            }
        }

        return w;
    }

    public static double lagrange_sum(int n, double[] x, double[] w, double[] y, double xv )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_SUM carries out a Lagrange interpolation rule.
        //
        //  Discussion:
        //
        //    It is assumed that LAGRANGE_RULE has already been called to compute
        //    the appropriate weights for the given set of abscissas.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 May 2001
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jean-Paul Berrut, Lloyd Trefethen,
        //    Barycentric Lagrange Interpolation,
        //    SIAM Review,
        //    Volume 46, Number 3, September 2004, pages 501-517.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the rule.
        //
        //    Input, double X[N], the abscissas of the rule.
        //
        //    Input, double W[N], the weights of the rule.
        //
        //    Input, double Y[N], the function values at the abscissas.
        //
        //    Input, double XV, a point where an interpolated value is
        //    needed.
        //
        //    Output, double LAGRANGE_SUM, the interpolated function value.
        //
    {
        int i;
        double yv;

        for (i = 0; i < n; i++)
        {
            if (!(Math.Abs(xv - x[i]) <= double.Epsilon))
            {
                continue;
            }

            yv = y[i];
            return yv;
        }

        double top = 0.0;
        double bot = 0.0;

        for (i = 0; i < n; i++)
        {
            top += w[i] * y[i] / (xv - x[i]);
            bot += w[i] / (xv - x[i]);
        }

        yv = top / bot;

        return yv;
    }

    public static double lagrange_val(int n, double[] x, double[] y, double xv )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_VAL applies a naive form of Lagrange interpolation.
        //
        //  Discussion:
        //
        //    Given N abscissas X, an arbitrary function Y(X) can be
        //    interpolated by a polynomial P(X) of order N (and degree N-1)
        //    using Lagrange basis polynomials of degree N-1.
        //
        //    Standard Lagrange interpolation can be rewritten into this form,
        //    which is more economical than evaluating the individual Lagrange
        //    basis polynomials.
        //
        //    If we define
        //
        //      L(I)(XV) = product ( 1 <= J <= N, J /= I )
        //        ( XV - X(J) ) / ( X(I) - X(J) )
        //
        //    then
        //
        //      P(XV) = sum ( 1 <= I <= N ) Y( X(I) ) * L(I)(XV)
        //
        //    Applying this form of the interpolation rule directly involves 
        //    about N^2 work.  There are more efficient forms of the rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 May 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data points.
        //
        //    Input, double X[N], the abscissas.
        //
        //    Input, double Y[N], the function values at the abscissas.
        //
        //    Input, double XV, a point where an interpolated value is
        //    needed.
        //
        //    Output, double LAGRANGE_VAL, the interpolated function value.
        //
    {
        int i;

        double yv = 0.0;

        for (i = 0; i < n; i++)
        {
            double poly = 1.0;
            int j;
            for (j = 0; j < n; j++)
            {
                if (j != i)
                {
                    poly = poly * (xv - x[j]) / (x[i] - x[j]);
                }
            }

            yv += y[i] * poly;
        }

        return yv;
    }
}