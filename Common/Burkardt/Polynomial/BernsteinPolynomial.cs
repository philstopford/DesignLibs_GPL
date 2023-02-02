using System;
using Burkardt.Types;

namespace Burkardt.PolynomialNS;

public static class BernsteinPolynomial
{
    public static double[] bernstein_matrix(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_MATRIX returns the Bernstein matrix.
        //
        //  Discussion:
        //
        //    The Bernstein matrix of order N is an NxN matrix A which can be used to
        //    transform a vector of power basis coefficients C representing a polynomial 
        //    P(X) to a corresponding Bernstein basis coefficient vector B:
        //
        //      B = A * C
        //
        //    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
        //    Bernstein basis vectors as ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
        //
        //    For N = 5, the matrix has the form:
        //
        //      1 -4   6  -4  1
        //      0  4 -12  12 -4
        //      0  0   6 -12  6
        //      0  0   0   4 -4
        //      0  0   0   0  1
        //
        //    and the numbers in each column represent the coefficients in the power
        //    series expansion of a Bernstein polynomial, so that 
        //
        //      B(5,4) = - 4 x^4 + 12 x^3 - 12 x^2 + 4 x
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double BERNSTEIN_MATRIX[N*N], the Bernstein matrix.
        //
    {
        int j;

        double[] a = new double[n * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i <= j; i++)
            {
                a[i + j * n] = typeMethods.r8_mop(j - i) * typeMethods.r8_choose(n - 1 - i, j - i)
                                                         * typeMethods.r8_choose(n - 1, i);
            }

            for (i = j + 1; i < n; i++)
            {
                a[i + j * n] = 0.0;
            }
        }

        return a;
    }

    public static double bernstein_matrix_determinant(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_MATRIX_DETERMINANT returns the determinant of the BERNSTEIN matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double BERNSTEIN_MATRIX_DETERMINANT, the determinant.
        //
    {
        int i;

        double value = 1.0;
        for (i = 0; i < n; i++)
        {
            value *= typeMethods.r8_choose(n - 1, i);
        }

        return value;
    }

    public static double[] bernstein_matrix_inverse(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_MATRIX_INVERSE returns the inverse Bernstein matrix.
        //
        //  Discussion:
        //
        //    The inverse Bernstein matrix of order N is an NxN matrix A which can 
        //    be used to transform a vector of Bernstein basis coefficients B
        //    representing a polynomial P(X) to a corresponding power basis 
        //    coefficient vector C:
        //
        //      C = A * B
        //
        //    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
        //    Bernstein basis vectors as ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
        //
        //    For N = 5, the matrix has the form:
        //
        //      1   1    1    1   1
        //      0  1/4  1/2  3/4  1
        //      0   0   1/6  1/2  1
        //      0   0    0   1/4  1
        //      0   0    0    0   1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double BERNSTEIN_MATRIX_INVERSE[N*N], the inverse Bernstein matrix.
        //
    {
        int j;

        double[] a = new double[n * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i <= j; i++)
            {
                a[i + j * n] = typeMethods.r8_choose(j, i) / typeMethods.r8_choose(n - 1, i);
            }

            for (i = j + 1; i < n; i++)
            {
                a[i + j * n] = 0.0;
            }
        }

        return a;
    }

    public static void bernstein_poly(int n, double x, ref double[] bern)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_POLY evaluates the Bernstein polynomials at a point X.
        //
        //  Discussion:
        //
        //    The Bernstein polynomials are assumed to be based on [0,1].
        //
        //    The formula is:
        //
        //      B(N,I,X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
        //
        //  First values:
        //
        //    B(0,0,X) = 1
        //
        //    B(1,0,X) =      1-X
        //    B(1,1,X) =               X
        //
        //    B(2,0,X) =     (1-X)^2
        //    B(2,1,X) = 2 * (1-X)   * X
        //    B(2,2,X) =               X^2
        //
        //    B(3,0,X) =     (1-X)^3
        //    B(3,1,X) = 3 * (1-X)^2 * X
        //    B(3,2,X) = 3 * (1-X)   * X^2
        //    B(3,3,X) =               X^3
        //
        //    B(4,0,X) =     (1-X)^4
        //    B(4,1,X) = 4 * (1-X)^3 * X
        //    B(4,2,X) = 6 * (1-X)^2 * X^2
        //    B(4,3,X) = 4 * (1-X)   * X^3
        //    B(4,4,X) =               X^4
        //
        //  Special values:
        //
        //    B(N,I,X) has a unique maximum value at X = I/N.
        //
        //    B(N,I,X) has an I-fold zero at 0 and and N-I fold zero at 1.
        //
        //    B(N,I,1/2) = C(N,K) / 2^N
        //
        //    For a fixed X and N, the polynomials add up to 1:
        //
        //      Sum ( 0 <= I <= N ) B(N,I,X) = 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2003
        //
        //  Parameters:
        //
        //    Input, int N, the degree of the Bernstein polynomials to be
        //    used.  For any N, there is a set of N+1 Bernstein polynomials,
        //    each of degree N, which form a basis for polynomials on [0,1].
        //
        //    Input, double X, the point at which the polynomials are to be evaluated.
        //
        //    Output, double BERN[N+1], the values of the Bernstein polynomials 
        //    of orders 0 through N at X.
        //
    {
        switch (n)
        {
            case 0:
                bern[0] = 1.0;
                break;
            case > 0:
            {
                bern[0] = 1.0 - x;
                bern[1] = x;

                int i;
                for (i = 2; i <= n; i++)
                {
                    bern[i] = x * bern[i - 1];
                    int j;
                    for (j = i - 1; 1 <= j; j--)
                    {
                        bern[j] = x * bern[j - 1] + (1.0 - x) * bern[j];
                    }

                    bern[0] = (1.0 - x) * bern[0];
                }

                break;
            }
        }
    }

    public static double[] bernstein_poly_01(int n, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_POLY_01 evaluates the Bernstein polynomials based in [0,1].
        //
        //  Discussion:
        //
        //    The Bernstein polynomials are assumed to be based on [0,1].
        //
        //    The formula is:
        //
        //      B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
        //
        //  First values:
        //
        //    B(0,0)(X) = 1
        //
        //    B(1,0)(X) =      1-X
        //    B(1,1)(X) =                X
        //
        //    B(2,0)(X) =     (1-X)^2
        //    B(2,1)(X) = 2 * (1-X)    * X
        //    B(2,2)(X) =                X^2
        //
        //    B(3,0)(X) =     (1-X)^3
        //    B(3,1)(X) = 3 * (1-X)^2 * X
        //    B(3,2)(X) = 3 * (1-X)   * X^2
        //    B(3,3)(X) =               X^3
        //
        //    B(4,0)(X) =     (1-X)^4
        //    B(4,1)(X) = 4 * (1-X)^3 * X
        //    B(4,2)(X) = 6 * (1-X)^2 * X^2
        //    B(4,3)(X) = 4 * (1-X)   * X^3
        //    B(4,4)(X) =               X^4
        //
        //  Special values:
        //
        //    B(N,I)(X) has a unique maximum value at X = I/N.
        //
        //    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
        //
        //    B(N,I)(1/2) = C(N,K) / 2^N
        //
        //    For a fixed X and N, the polynomials add up to 1:
        //
        //      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the degree of the Bernstein polynomials 
        //    to be used.  For any N, there is a set of N+1 Bernstein polynomials,
        //    each of degree N, which form a basis for polynomials on [0,1].
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double BERNSTEIN_POLY[N+1], the values of the N+1 
        //    Bernstein polynomials at X.
        //
    {
        double[] bern = new double[n + 1];

        switch (n)
        {
            case 0:
                bern[0] = 1.0;
                break;
            case > 0:
            {
                bern[0] = 1.0 - x;
                bern[1] = x;

                int i;
                for (i = 2; i <= n; i++)
                {
                    bern[i] = x * bern[i - 1];
                    int j;
                    for (j = i - 1; 1 <= j; j--)
                    {
                        bern[j] = x * bern[j - 1]
                                  + (1.0 - x) * bern[j];
                    }

                    bern[0] = (1.0 - x) * bern[0];
                }

                break;
            }
        }

        return bern;
    }

    public static double[] bernstein_poly_01_matrix(int m, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_POLY_01_MATRIX evaluates the Bernstein polynomials based in [0,1].
        //
        //  Discussion:
        //
        //    The Bernstein polynomials are assumed to be based on [0,1].
        //
        //    The formula is:
        //
        //      B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
        //
        //  First values:
        //
        //    B(0,0)(X) = 1
        //
        //    B(1,0)(X) =      1-X
        //    B(1,1)(X) =                X
        //
        //    B(2,0)(X) =     (1-X)^2
        //    B(2,1)(X) = 2 * (1-X)    * X
        //    B(2,2)(X) =                X^2
        //
        //    B(3,0)(X) =     (1-X)^3
        //    B(3,1)(X) = 3 * (1-X)^2 * X
        //    B(3,2)(X) = 3 * (1-X)   * X^2
        //    B(3,3)(X) =               X^3
        //
        //    B(4,0)(X) =     (1-X)^4
        //    B(4,1)(X) = 4 * (1-X)^3 * X
        //    B(4,2)(X) = 6 * (1-X)^2 * X^2
        //    B(4,3)(X) = 4 * (1-X)   * X^3
        //    B(4,4)(X) =               X^4
        //
        //  Special values:
        //
        //    B(N,I)(X) has a unique maximum value at X = I/N.
        //
        //    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
        //
        //    B(N,I)(1/2) = C(N,K) / 2^N
        //
        //    For a fixed X and N, the polynomials add up to 1:
        //
        //      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of evaluation points.
        //
        //    Input, int N, the degree of the Bernstein polynomials 
        //    to be used.  For any N, there is a set of N+1 Bernstein polynomials,
        //    each of degree N, which form a basis for polynomials on [0,1].
        //
        //    Input, double X[M], the evaluation points.
        //
        //    Output, double BERNSTEIN_POLY_01_MATRIX[M*(N+1)], the values of the N+1 
        //    Bernstein polynomials at the evaluation points.
        //
    {
        int i;

        double[] b = new double[m * (n + 1)];

        for (i = 0; i < m; i++)
        {
            switch (n)
            {
                case 0:
                    b[i + 0 * m] = 1.0;
                    break;
                case > 0:
                {
                    b[i + 0 * m] = 1.0 - x[i];
                    b[i + 1 * m] = x[i];

                    int j;
                    for (j = 2; j <= n; j++)
                    {
                        b[i + j * m] = x[i] * b[i + (j - 1) * m];
                        int k;
                        for (k = j - 1; 1 <= k; k--)
                        {
                            b[i + k * m] = x[i] * b[i + (k - 1) * m]
                                           + (1.0 - x[i]) * b[i + k * m];
                        }

                        b[i + 0 * m] = (1.0 - x[i]) * b[i + 0 * m];
                    }

                    break;
                }
            }
        }

        return b;
    }

    public static double[] bernstein_poly_ab(int n, double a, double b, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_POLY_AB evaluates at X the Bernstein polynomials based in [A,B].
        //
        //  Discussion:
        //
        //    The formula is:
        //
        //      BERN(N,I)(X) = [N!/(I!*(N-I)!)] * (B-X)^(N-I) * (X-A)^I / (B-A)^N
        //
        //  First values:
        //
        //    B(0,0)(X) =   1
        //
        //    B(1,0)(X) = (      B-X                ) / (B-A)
        //    B(1,1)(X) = (                 X-A     ) / (B-A)
        //
        //    B(2,0)(X) = (     (B-X)^2             ) / (B-A)^2
        //    B(2,1)(X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)^2
        //    B(2,2)(X) = (                (X-A)^2  ) / (B-A)^2
        //
        //    B(3,0)(X) = (     (B-X)^3             ) / (B-A)^3
        //    B(3,1)(X) = ( 3 * (B-X)^2  * (X-A)    ) / (B-A)^3
        //    B(3,2)(X) = ( 3 * (B-X)    * (X-A)^2  ) / (B-A)^3
        //    B(3,3)(X) = (                (X-A)^3  ) / (B-A)^3
        //
        //    B(4,0)(X) = (     (B-X)^4             ) / (B-A)^4
        //    B(4,1)(X) = ( 4 * (B-X)^3  * (X-A)    ) / (B-A)^4
        //    B(4,2)(X) = ( 6 * (B-X)^2  * (X-A)^2  ) / (B-A)^4
        //    B(4,3)(X) = ( 4 * (B-X)    * (X-A)^3  ) / (B-A)^4
        //    B(4,4)(X) = (                (X-A)^4  ) / (B-A)^4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the degree of the Bernstein polynomials 
        //    to be used.  For any N, there is a set of N+1 Bernstein polynomials, 
        //    each of degree N, which form a basis for polynomials on [A,B].
        //
        //    Input, double A, B, the endpoints of the interval on which the
        //    polynomials are to be based.  A and B should not be equal.
        //
        //    Input, double X, the point at which the polynomials 
        //    are to be evaluated.
        //
        //    Output, double BERNSTEIN_POLY_AB[N+1], the values of the N+1
        //    Bernstein polynomials at X.
        //
    {
        if (Math.Abs(b - a) <= typeMethods.r8_epsilon())
        {
            Console.WriteLine("");
            Console.WriteLine("BERNSTEIN_POLY_AB - Fatal error!");
            Console.WriteLine("  A = B = " + a + "");
            return new double[1];
        }

        double[] bern = new double[n + 1];

        switch (n)
        {
            case 0:
                bern[0] = 1.0;
                break;
            case > 0:
            {
                bern[0] = (b - x) / (b - a);
                bern[1] = (x - a) / (b - a);

                int i;
                for (i = 2; i <= n; i++)
                {
                    bern[i] = (x - a) * bern[i - 1] / (b - a);
                    int j;
                    for (j = i - 1; 1 <= j; j--)
                    {
                        bern[j] = ((b - x) * bern[j]
                                   + (x - a) * bern[j - 1])
                                  / (b - a);
                    }

                    bern[0] = (b - x) * bern[0] / (b - a);
                }

                break;
            }
        }

        return bern;
    }

    public static void bpab(int n, double x, double a, double b, ref double[] bern)
    {
        bern = bernstein_poly_ab(n, a, b, x);
    }

    public static double bpab_approx(int n, double a, double b, double[] ydata, double xval )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BPAB_APPROX evaluates the Bernstein polynomial for F(X) on [A,B].
        //
        //  Formula:
        //
        //    BERN(F)(X) = sum ( 0 <= I <= N ) F(X(I)) * B_BASE(I,X)
        //
        //    where
        //
        //      X(I) = ( ( N - I ) * A + I * B ) / N
        //      B_BASE(I,X) is the value of the I-th Bernstein basis polynomial at X.
        //
        //  Discussion:
        //
        //    The Bernstein polynomial BERN(F) for F(X) is an approximant, not an
        //    interpolant; in other words, its value is not guaranteed to equal
        //    that of F at any particular point.  However, for a fixed interval
        //    [A,B], if we let N increase, the Bernstein polynomial converges
        //    uniformly to F everywhere in [A,B], provided only that F is continuous.
        //    Even if F is not continuous, but is bounded, the polynomial converges
        //    pointwise to F(X) at all points of continuity.  On the other hand,
        //    the convergence is quite slow compared to other interpolation
        //    and approximation schemes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input, int N, the degree of the Bernstein polynomial to be used.
        //
        //    Input, double A, B, the endpoints of the interval on which the
        //    approximant is based.  A and B should not be equal.
        //
        //    Input, double YDATA[0:N], the data values at N+1 equally spaced points
        //    in [A,B].  If N = 0, then the evaluation point should be 0.5 * ( A + B).
        //    Otherwise, evaluation point I should be ( (N-I)*A + I*B ) / N ).
        //
        //    Input, double XVAL, the point at which the Bernstein polynomial
        //    approximant is to be evaluated.  XVAL does not have to lie in the
        //    interval [A,B].
        //
        //    Output, double BPAB_APPROX, the value of the Bernstein polynomial approximant
        //    for F, based in [A,B], evaluated at XVAL.
        //
    {
        int i;
        //
        //  Evaluate the Bernstein basis polynomials at XVAL.
        //
        double[] bvec = bpab(n, a, b, xval);
        //
        //  Now compute the sum of YDATA(I) * BVEC(I).
        //
        double yval = 0.0;

        for (i = 0; i <= n; i++)
        {
            yval += ydata[i] * bvec[i];
        }
            
        return yval;
    }

    public static double[] bpab(int n, double a, double b, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BPAB evaluates the Bernstein basis polynomials for [A,B] at a point.
        //
        //  Formula:
        //
        //    BERN(N,I,X) = [N!/(I!*(N-I)!)] * (B-X)^(N-I) * (X-A)^I / (B-A)^N
        //
        //  First values:
        //
        //    B(0,0,X) =   1
        //
        //    B(1,0,X) = (      B-X              ) / (B-A)
        //    B(1,1,X) = (                X-A    ) / (B-A)
        //
        //    B(2,0,X) = (     (B-X)^2           ) / (B-A)^2
        //    B(2,1,X) = ( 2 * (B-X)   * (X-A)   ) / (B-A)^2
        //    B(2,2,X) = (               (X-A)^2 ) / (B-A)^2
        //
        //    B(3,0,X) = (     (B-X)^3           ) / (B-A)^3
        //    B(3,1,X) = ( 3 * (B-X)^2 * (X-A)   ) / (B-A)^3
        //    B(3,2,X) = ( 3 * (B-X)   * (X-A)^2 ) / (B-A)^3
        //    B(3,3,X) = (               (X-A)^3 ) / (B-A)^3
        //
        //    B(4,0,X) = (     (B-X)^4           ) / (B-A)^4
        //    B(4,1,X) = ( 4 * (B-X)^3 * (X-A)   ) / (B-A)^4
        //    B(4,2,X) = ( 6 * (B-X)^2 * (X-A)^2 ) / (B-A)^4
        //    B(4,3,X) = ( 4 * (B-X)   * (X-A)^3 ) / (B-A)^4
        //    B(4,4,X) = (               (X-A)^4 ) / (B-A)^4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input, integer N, the degree of the Bernstein basis polynomials.
        //    For any N greater than or equal to 0, there is a set of N+1
        //    Bernstein basis polynomials, each of degree N, which form a basis
        //    for polynomials on [A,B].
        //
        //    Input, double A, B, the endpoints of the interval on which the
        //    polynomials are to be based.  A and B should not be equal.
        //
        //    Input, double X, the point at which the polynomials are to be
        //    evaluated.  X need not lie in the interval [A,B].
        //
        //    Output, double BERN[0:N], the values of the N+1 Bernstein basis
        //    polynomials at X.
        //
    {
        if (Math.Abs(b - a) <= typeMethods.r8_epsilon())
        {
            Console.WriteLine("");
            Console.WriteLine("BPAB - Fatal error!");
            Console.WriteLine("  A = B = " + a + "");
            return null;
        }

        double[] bern = new double[n + 1];

        switch (n)
        {
            case 0:
                bern[0] = 1.0;
                break;
            case > 0:
            {
                bern[0] = (b - x) / (b - a);
                bern[1] = (x - a) / (b - a);

                int i;
                for (i = 2; i <= n; i++)
                {
                    bern[i] = (x - a) * bern[i - 1] / (b - a);
                    int j;
                    for (j = i - 1; 1 <= j; j--)
                    {
                        bern[j] = ((b - x) * bern[j] + (x - a) * bern[j - 1]) / (b - a);
                    }

                    bern[0] = (b - x) * bern[0] / (b - a);
                }

                break;
            }
        }

        return bern;
    }

    public static double[] bernstein_poly_ab_approx(int n, double a, double b, double[] ydata,
            int nval, double[] xval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_POLY_AB_APPROX: Bernstein approximant to F(X) on [A,B].
        //
        //  Formula:
        //
        //    BPAB(F)(X) = sum ( 0 <= I <= N ) F(X(I)) * B_BASE(I,X)
        //
        //    where
        //
        //      X(I) = ( ( N - I ) * A + I * B ) / N
        //      B_BASE(I,X) is the value of the I-th Bernstein basis polynomial at X.
        //
        //  Discussion:
        //
        //    The Bernstein polynomial BPAB(F) for F(X) over [A,B] is an approximant, 
        //    not an interpolant; in other words, its value is not guaranteed to equal
        //    that of F at any particular point.  However, for a fixed interval
        //    [A,B], if we let N increase, the Bernstein polynomial converges
        //    uniformly to F everywhere in [A,B], provided only that F is continuous.
        //    Even if F is not continuous, but is bounded, the polynomial converges
        //    pointwise to F(X) at all points of continuity.  On the other hand,
        //    the convergence is quite slow compared to other interpolation
        //    and approximation schemes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input, int N, the degree of the Bernstein polynomial
        //    to be used.  N must be at least 0.
        //
        //    Input, double A, B, the endpoints of the interval on which the
        //    approximant is based.  A and B should not be equal.
        //
        //    Input, double YDATA[N+1], the data values at N+1 equally
        //    spaced points in [A,B].  If N = 0, then the evaluation point should
        //    be 0.5 * ( A + B).  Otherwise, evaluation point I should be
        //    ( (N-I)*A + I*B ) / N ).
        //
        //    Input, int NVAL, the number of points at which the
        //    approximant is to be evaluated.
        //
        //    Input, double XVAL[NVAL], the point at which the Bernstein 
        //    polynomial approximant is to be evaluated.  The entries of XVAL do not 
        //    have to lie in the interval [A,B].
        //
        //    Output, double BPAB_APPROX[NVAL], the values of the Bernstein 
        //    polynomial approximant for F, based in [A,B], evaluated at XVAL.
        //
    {
        int i;

        double[] yval = new double[nval];

        for (i = 0; i < nval; i++)
        {
            //
            //  Evaluate the Bernstein basis polynomials at XVAL.
            //
            double[] bvec = bernstein_poly_ab(n, a, b, xval[i]);
            //
            //  Now compute the sum of YDATA(I) * BVEC(I).
            //
            yval[i] = typeMethods.r8vec_dot_product(n + 1, ydata, bvec);
        }

        return yval;
    }

    public static double[] bernstein_to_legendre(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_TO_LEGENDRE returns the Bernstein-to-Legendre matrix.
        //
        //  Discussion:
        //
        //    The Legendre polynomials are often defined on [-1,+1], while the
        //    Bernstein polynomials are defined on [0,1].  For this function,
        //    the Legendre polynomials have been shifted to share the [0,1]
        //    interval of definition.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the maximum degree of the polynomials.
        //
        //    Output, double BERNSTEIN_TO_LEGENDRE[(N+1)*(N+1)],
        //    the Bernstein-to-Legendre matrix.
        //
    {
        int i;

        double[] a = typeMethods.r8mat_zeros_new(n + 1, n + 1);

        for (i = 0; i <= n; i++)
        {
            int j;
            for (j = 0; j <= n; j++)
            {
                int k;
                for (k = 0; k <= i; k++)
                {
                    a[i + j * (n + 1)] += typeMethods.r8_mop(i + k) * Math.Pow(typeMethods.r8_choose(i, k), 2)
                                          / typeMethods.r8_choose(n + i, j + k);
                }

                a[i + j * (n + 1)] = a[i + j * (n + 1)] * typeMethods.r8_choose(n, j)
                                                        * (2 * i + 1) / (n + i + 1);
            }
        }

        return a;
    }

    public static double[] bernstein_to_power(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_TO_POWER returns the Bernstein-to-Power matrix.
        //
        //  Discussion:
        //
        //    The Bernstein-to-Power matrix of degree N is an N+1xN+1 matrix A which can 
        //    be used to transform the N+1 coefficients of a polynomial of degree N
        //    from a vector B of Bernstein basis polynomial coefficients ((1-x)^n,...,x^n).
        //    to a vector P of coefficients of the power basis (1,x,x^2,...,x^n).
        //
        //    If we are using N=4-th degree polynomials, the matrix has the form:
        //
        //      1   0   0   0  0
        //     -4   4   0   0  0
        //      6 -12   6   0  0
        //     -4  12 -12   4  0
        //      1  -4   6  -4  1
        //
        //   and a polynomial with the Bernstein basis representation
        //     p(x) = 3/4 * b(4,1) + 1/2 b(4,2)
        //   whose Bernstein coefficient vector is
        //     B = ( 0, 3/4, 1/2, 0, 0 )
        //   will have the Bernstein basis coefficients 
        //     P = A * B = ( 0, 3, -6, 3, 0 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the maximum degree of the polynomials.
        //
        //    Output, double BERNSTEIN_TO_POWER[(N+1)*(N+1)],
        //    the Bernstein-to-Power matrix.
        //
    {
        int j;

        double[] a = typeMethods.r8mat_zeros_new(n + 1, n + 1);

        for (j = 0; j <= n; j++)
        {
            int i;
            for (i = 0; i <= j; i++)
            {
                a[i + j * (n + 1)] = typeMethods.r8_mop(j - i) * typeMethods.r8_choose(n - i, j - i)
                                     / typeMethods.r8_choose(n, i);
            }
        }

        return a;
    }

    public static double[] bernstein_vandermonde(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_VANDERMONDE returns the Bernstein Vandermonde matrix.
        //
        //  Discussion:
        //
        //    The Bernstein Vandermonde matrix of order N is constructed by
        //    evaluating the N Bernstein polynomials of degree N-1 at N equally
        //    spaced points between 0 and 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double BERNSTEIN_VANDERMONDE[N*N], the Bernstein Vandermonde matrix.
        //
    {
        int i;

        double[] v = new double[n * n];

        switch (n)
        {
            case 1:
                v[0 + 0 * 1] = 1.0;
                return v;
        }

        for (i = 0; i < n; i++)
        {
            double x = i / (double)(n - 1);
            double[] b = bernstein_poly_01(n - 1, x);
            int j;
            for (j = 0; j < n; j++)
            {
                v[i + j * n] = b[j];
            }
        }

        return v;
    }

    public static double[] legendre_to_bernstein(int n)

        //*****************************************************************************/
        //
        //  Purpose:
        //
        //    LEGENDRE_TO_BERNSTEIN returns the Legendre-to-Bernstein matrix.
        //
        //  Discussion:
        //
        //    The Legendre polynomials are often defined on [-1,+1], while the
        //    Bernstein polynomials are defined on [0,1].  For this function,
        //    the Legendre polynomials have been shifted to share the [0,1]
        //    interval of definition.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the maximum degree of the polynomials.
        //
        //    Output, double LEGENDRE_TO_BERNSTEIN[(N+1)*(N+1)], the 
        //    Legendre-to-Bernstein matrix.
        //
    {
        int i;

        double[] a = typeMethods.r8mat_zeros_new(n + 1, n + 1);

        for (i = 0; i <= n; i++)
        {
            int j;
            for (j = 0; j <= n; j++)
            {
                int k;
                for (k = Math.Max(0, i + j - n); k <= Math.Min(i, j); k++)
                {
                    a[i + j * (n + 1)] += typeMethods.r8_mop(j + k) * Math.Pow(typeMethods.r8_choose(j, k), 2)
                                                                    * typeMethods.r8_choose(n - j, i - k);
                }

                a[i + j * (n + 1)] /= typeMethods.r8_choose(n, i);
            }
        }

        return a;
    }

    public static double[] power_to_bernstein(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_TO_BERNSTEIN returns the Power-to-Bernstein matrix.
        //
        //  Discussion:
        //
        //    The Power-to-Bernstein matrix of degree N is an N+1xN+1 matrix A which can 
        //    be used to transform the N+1 coefficients of a polynomial of degree N
        //    from a vector P of coefficients of the power basis (1,x,x^2,...,x^n)
        //    to a vector B of Bernstein basis polynomial coefficients ((1-x)^n,...,x^n).
        //
        //    If we are using N=4-th degree polynomials, the matrix has the form:
        //
        //          1   0    0    0   0
        //          1  1/4   0    0   0
        //      A = 1  1/2  1/6   0   0
        //          1  3/4  1/2  1/4  1
        //          1   1    1    1   1
        //
        //   and a polynomial 
        //     p(x) = 3x - 6x^2 + 3x^3
        //   whose power coefficient vector is
        //     P = ( 0, 3, -6, 3, 0 )
        //   will have the Bernstein basis coefficients 
        //     B = A * P = ( 0, 3/4, 1/2, 0, 0 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the maximum degree of the polynomials.
        //
        //    Output, double POWER_TO_BERNSTEIN[(N+1)*(N+1)], the 
        //    Power-to-Bernstein matrix.
        //
    {
        int j;

        double[] a = typeMethods.r8mat_zeros_new(n + 1, n + 1);

        for (j = 0; j <= n; j++)
        {
            int i;
            for (i = 0; i <= j; i++)
            {
                a[n - 1 + (n - j) * (n + 1)] = typeMethods.r8_choose(j, i) / typeMethods.r8_choose(n, i);
            }
        }

        return a;
    }


}