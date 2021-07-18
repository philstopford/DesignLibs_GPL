namespace Burkardt
{
    public static class Vandermonde
    {
        public static double[] cheby_van1(int m, double a, double b, int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_VAN1 returns the CHEBY_VAN1 matrix.
            //
            //  Discussion:
            //
            //    Normally, the Chebyshev polynomials are defined on -1 <= XI <= +1.
            //    Here, we assume the Chebyshev polynomials have been defined on the
            //    interval A <= X <= B, using the mapping
            //      XI = ( - ( B - X ) + ( X - A ) ) / ( B - A )
            //    so that
            //      ChebyAB(A,B;X) = Cheby(XI).
            //
            //    if ( I == 1 ) then
            //      V(1,1:N) = 1;
            //    elseif ( I == 2 ) then
            //      V(2,1:N) = XI(1:N);
            //    else
            //      V(I,1:N) = 2.0 * XI(1:N) * V(I-1,1:N) - V(I-2,1:N);
            //
            //  Example:
            //
            //    N = 5, A = -1, B = +1, X = ( 1, 2, 3, 4, 5 )
            //
            //    1  1   1    1    1
            //    1  2   3    4    5
            //    1  7  17   31   49
            //    1 26  99  244  485
            //    1 97 577 1921 4801
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Nicholas Higham,
            //    Stability analysis of algorithms for solving confluent
            //    Vandermonde-like systems,
            //    SIAM Journal on Matrix Analysis and Applications,
            //    Volume 11, 1990, pages 23-41.
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows of the matrix.
            //
            //    Input, double A, B, the interval.
            //
            //    Input, int N, the number of columns of the matrix.
            //
            //    Input, double X[N], the vector that defines the matrix.
            //
            //    Output, double CHEBY_VAN1[M*N], the matrix.
            //
        {
            int i;
            int j;
            double[] v;
            double xi;

            v = new double[m * n];

            for (j = 0; j < n; j++)
            {
                xi = (-(b - x[j]) + (x[j] - a)) / (b - a);
                for (i = 0; i < m; i++)
                {
                    if (i == 0)
                    {
                        v[i + j * m] = 1.0;
                    }
                    else if (i == 1)
                    {
                        v[i + j * m] = xi;
                    }
                    else
                    {
                        v[i + j * m] = 2.0 * xi * v[i - 1 + j * m] - v[i - 2 + j * m];
                    }
                }
            }

            return v;
        }

        public static double[] legendre_van(int m, double a, double b, int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_VAN returns the LEGENDRE_VAN matrix.
            //
            //  Discussion:
            //
            //    The LEGENDRE_VAN matrix is the Legendre Vandermonde-like matrix.
            //
            //    Normally, the Legendre polynomials are defined on -1 <= XI <= +1.
            //    Here, we assume the Legendre polynomials have been defined on the
            //    interval A <= X <= B, using the mapping
            //      XI = ( - ( B - X ) + ( X - A ) ) / ( B - A )
            //    so that
            //      Lab(A,B;X) = L(XI).
            //
            //    if ( I = 1 ) then
            //      V(1,1:N) = 1
            //    else if ( I = 2 ) then
            //      V(2,1:N) = XI(1:N)
            //    else
            //      V(I,1:N) = ( (2*I-1) * XI(1:N) * V(I-1,1:N) - (I-1)*V(I-2,1:N) ) / I
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows of the matrix.
            //
            //    Input, double A, B, the interval.
            //
            //    Input, int N, the number of columns of the matrix.
            //
            //    Input, double X[N], the vector that defines the matrix.
            //
            //    Output, double LEGENDRE_VAN[M*N], the matrix.
            //
        {
            int i;
            int j;
            double[] v;
            double xi;

            v = new double[m * n];

            for (j = 0; j < n; j++)
            {
                xi = (-(b - x[j]) + (x[j] - a)) / (b - a);
                for (i = 0; i < m; i++)
                {
                    if (i == 0)
                    {
                        v[i + j * m] = 1.0;
                    }
                    else if (i == 1)
                    {
                        v[i + j * m] = xi;
                    }
                    else
                    {
                        v[i + j * m] = ((double) (2 * i - 1) * xi * v[i - 1 + j * m] +
                                        (double) (-i + 1) * v[i - 2 + j * m])
                                       / (double) (i);
                    }
                }
            }

            return v;
        }
    }
}