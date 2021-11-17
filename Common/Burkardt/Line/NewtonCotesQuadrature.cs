using Burkardt.Types;

namespace Burkardt.LineNS;

public static class NewtonCotesQuadrature
{
    public static void line_ncc_rule ( int n, double a, double b, double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_NCC_RULE computes a Newton-Cotes Closed (NCC) quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      Integral ( A <= X <= B ) F(X) dx
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //
        //    Input, double A, B, the endpoints of the interval.
        //
        //    Input, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        double[] d;
        int i;
        int j;
        int k;
        double y_a;
        double y_b;
        //
        //  Define the points X.
        //
        typeMethods.r8vec_linspace ( n, a, b, ref x );

        d = new double[n];

        for ( i = 0; i < n; i++ )
        {
            //
            //  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
            //  and zero at the other nodes.
            //
            for ( j = 0; j < n; j++ )
            {
                d[j] = 0.0;
            }
            d[i] = 1.0;

            for ( j = 2; j <= n; j++ )
            {
                for ( k = j; k <= n; k++ )
                {
                    d[n+j-k-1] = ( d[n+j-k-2] - d[n+j-k-1] ) / ( x[n-k] - x[n+j-k-1] );
                }
            }
            for ( j = 1; j <= n - 1; j++ )
            {
                for ( k = 1; k <= n - j; k++ )
                {
                    d[n-k-1] -= x[n-k-j] * d[n-k];
                }
            }
            //
            //  Evaluate the antiderivative of the polynomial at the endpoints.
            //
            y_a = d[n-1] / n;
            for ( j = n - 2; 0 <= j; j-- )
            {
                y_a = y_a * a + d[j] / (j + 1);
            }
            y_a *= a;

            y_b = d[n-1] / n;
            for ( j = n - 2; 0 <= j; j-- )
            {
                y_b = y_b * b + d[j] / (j + 1);
            }
            y_b *= b;

            w[i] = y_b - y_a;
        }
    }
        
    public static void line_nco_rule ( int n, double a, double b, double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_NCO_RULE computes a Newton-Cotes Open (NCO) quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      Integral ( A <= X <= B ) F(X) dx
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //
        //    Input, double A, B, the endpoints of the interval.
        //
        //    Input, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        double[] d;
        int i;
        int j;
        int k;
        double y_a;
        double y_b;
        //
        //  Define the points X.
        //
        typeMethods.r8vec_linspace2 ( n, a, b, ref x );

        d = new double[n];

        for ( i = 0; i < n; i++ )
        {
            //
            //  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
            //  and zero at the other nodes.
            //
            for ( j = 0; j < n; j++ )
            {
                d[j] = 0.0;
            }
            d[i] = 1.0;

            for ( j = 2; j <= n; j++ )
            {
                for ( k = j; k <= n; k++ )
                {
                    d[n+j-k-1] = ( d[n+j-k-2] - d[n+j-k-1] ) / ( x[n-k] - x[n+j-k-1] );
                }
            }
            for ( j = 1; j <= n - 1; j++ )
            {
                for ( k = 1; k <= n - j; k++ )
                {
                    d[n-k-1] -= x[n-k-j] * d[n-k];
                }
            }
            //
            //  Evaluate the antiderivative of the polynomial at the endpoints.
            //
            y_a = d[n-1] / n;
            for ( j = n - 2; 0 <= j; j-- )
            {
                y_a = y_a * a + d[j] / (j + 1);
            }
            y_a *= a;

            y_b = d[n-1] / n;
            for ( j = n - 2; 0 <= j; j-- )
            {
                y_b = y_b * b + d[j] / (j + 1);
            }
            y_b *= b;

            w[i] = y_b - y_a;
        }
    }
}