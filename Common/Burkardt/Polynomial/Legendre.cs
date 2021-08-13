using System;
using Burkardt.MatrixNS;
using Burkardt.MonomialNS;
using Burkardt.PolynomialNS;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace Burkardt.PolynomialNS
{
    public static class Legendre
    {
        public static void legendre_associated(int n, int m, double x, ref double[] cx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_ASSOCIATED evaluates the associated Legendre functions.
            //
            //  Differential equation:
            //
            //    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
            //
            //  First terms:
            //
            //    M = 0  ( = Legendre polynomials of first kind P(N,X) )
            //
            //    P00 =    1
            //    P10 =    1 X
            //    P20 = (  3 X^2 -   1)/2
            //    P30 = (  5 X^3 -   3 X)/2
            //    P40 = ( 35 X^4 -  30 X^2 +   3)/8
            //    P50 = ( 63 X^5 -  70 X^3 +  15 X)/8
            //    P60 = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
            //    P70 = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
            //
            //    M = 1
            //
            //    P01 =   0
            //    P11 =   1 * SQRT(1-X*X)
            //    P21 =   3 * SQRT(1-X*X) * X
            //    P31 = 1.5 * SQRT(1-X*X) * (5*X*X-1)
            //    P41 = 2.5 * SQRT(1-X*X) * (7*X*X*X-3*X)
            //
            //    M = 2
            //
            //    P02 =   0
            //    P12 =   0
            //    P22 =   3 * (1-X*X)
            //    P32 =  15 * (1-X*X) * X
            //    P42 = 7.5 * (1-X*X) * (7*X*X-1)
            //
            //    M = 3
            //
            //    P03 =   0
            //    P13 =   0
            //    P23 =   0
            //    P33 =  15 * (1-X*X)^1.5
            //    P43 = 105 * (1-X*X)^1.5 * X
            //
            //    M = 4
            //
            //    P04 =   0
            //    P14 =   0
            //    P24 =   0
            //    P34 =   0
            //    P44 = 105 * (1-X*X)^2
            //
            //  Recursion:
            //
            //    if N < M:
            //      P(N,M) = 0
            //    if N = M:
            //      P(N,M) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
            //      all the odd integers less than or equal to N.
            //    if N = M+1:
            //      P(N,M) = X*(2*M+1)*P(M,M)
            //    if M+1 < N:
            //      P(N,M) = ( X*(2*N-1)*P(N-1,M) - (N+M-1)*P(N-2,M) )/(N-M)
            //
            //  Special values:
            //
            //    P(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
            //    function of the first kind equals the Legendre polynomial of the
            //    first kind.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 March 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //  Parameters:
            //
            //    Input, int N, the maximum first index of the Legendre
            //    function, which must be at least 0.
            //
            //    Input, int M, the second index of the Legendre function,
            //    which must be at least 0, and no greater than N.
            //
            //    Input, double X, the point at which the function is to be
            //    evaluated.  X must satisfy -1 <= X <= 1.
            //
            //    Output, double CX[N+1], the values of the first N+1 function.
            //
        {
            double factor;
            int i;
            double somx2;

            if (m < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_ASSOCIATED - Fatal error!");
                Console.WriteLine("  Input value of M is " + m + "");
                Console.WriteLine("  but M must be nonnegative.");
                return;
            }

            if (n < m)
            {
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_ASSOCIATED - Fatal error!");
                Console.WriteLine("  Input value of M = " + m + "");
                Console.WriteLine("  Input value of N = " + n + "");
                Console.WriteLine("  but M must be less than or equal to N.");
                return;
            }

            if (x < -1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_ASSOCIATED - Fatal error!");
                Console.WriteLine("  Input value of X = " + x + "");
                Console.WriteLine("  but X must be no less than -1.");
                return;
            }

            if (1.0 < x)
            {
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_ASSOCIATED - Fatal error!");
                Console.WriteLine("  Input value of X = " + x + "");
                Console.WriteLine("  but X must be no more than 1.");
                return;
            }

            for (i = 0; i <= m - 1; i++)
            {
                cx[i] = 0.0;
            }

            cx[m] = 1.0;

            somx2 = Math.Sqrt(1.0 - x * x);

            factor = 1.0;
            for (i = 1; i <= m; i++)
            {
                cx[m] = -cx[m] * factor * somx2;
                factor = factor + 2.0;
            }

            if (m == n)
            {
                return;
            }

            cx[m + 1] = x * (double)(2 * m + 1) * cx[m];

            for (i = m + 2; i <= n; i++)
            {
                cx[i] = ((double)(2 * i - 1) * x * cx[i - 1]
                         + (double)(-i - m + 1) * cx[i - 2])
                        / (double)(i - m);
            }
        }

        public static void legendre_associated_normalized(int n, int m, double x, ref double[] cx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_ASSOCIATED_NORMALIZED: normalized associated Legendre functions.
            //
            //  Discussion:
            //
            //    The unnormalized associated Legendre functions P_N^M(X) have
            //    the property that
            //
            //      Integral ( -1 <= X <= 1 ) ( P_N^M(X) )^2 dX
            //      = 2 * ( N + M )! / ( ( 2 * N + 1 ) * ( N - M )! )
            //
            //    By dividing the function by the square root of this term,
            //    the normalized associated Legendre functions have norm 1.
            //
            //    However, we plan to use these functions to build spherical
            //    harmonics, so we use a slightly different normalization factor of
            //
            //      sqrt ( ( ( 2 * N + 1 ) * ( N - M )! ) / ( 4 * pi * ( N + M )! ) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 March 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //  Parameters:
            //
            //    Input, int N, the maximum first index of the Legendre
            //    function, which must be at least 0.
            //
            //    Input, int M, the second index of the Legendre function,
            //    which must be at least 0, and no greater than N.
            //
            //    Input, double X, the point at which the function is to be
            //    evaluated.  X must satisfy -1 <= X <= 1.
            //
            //    Output, double CX[N+1], the values of the first N+1 function.
            //
        {
            double factor;
            int i;
            int mm;
            const double r8_pi = 3.141592653589793;
            double somx2;

            if (m < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!");
                Console.WriteLine("  Input value of M is " + m + "");
                Console.WriteLine("  but M must be nonnegative.");
                return;
            }

            if (n < m)
            {
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!");
                Console.WriteLine("  Input value of M = " + m + "");
                Console.WriteLine("  Input value of N = " + n + "");
                Console.WriteLine("  but M must be less than or equal to N.");
                return;
            }

            if (x < -1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!");
                Console.WriteLine("  Input value of X = " + x + "");
                Console.WriteLine("  but X must be no less than -1.");
                return;
            }

            if (1.0 < x)
            {
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!");
                Console.WriteLine("  Input value of X = " + x + "");
                Console.WriteLine("  but X must be no more than 1.");
                return;
            }

            for (i = 0; i <= m - 1; i++)
            {
                cx[i] = 0.0;
            }

            cx[m] = 1.0;

            somx2 = Math.Sqrt(1.0 - x * x);

            factor = 1.0;
            for (i = 1; i <= m; i++)
            {
                cx[m] = -cx[m] * factor * somx2;
                factor = factor + 2.0;
            }

            if (m + 1 <= n)
            {
                cx[m + 1] = x * (double)(2 * m + 1) * cx[m];
            }

            for (i = m + 2; i <= n; i++)
            {
                cx[i] = ((double)(2 * i - 1) * x * cx[i - 1]
                         + (double)(-i - m + 1) * cx[i - 2])
                        / (double)(i - m);
            }

            //
            //  Normalization.
            //
            for (mm = m; mm <= n; mm++)
            {
                factor = Math.Sqrt(((double)(2 * mm + 1) * typeMethods.r8_factorial(mm - m))
                              / (4.0 * r8_pi * typeMethods.r8_factorial(mm + m)));
                cx[mm] = cx[mm] * factor;
            }
        }

        public static double[] legendre_zeros(int order)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_ZEROS returns the zeros of the Legendre polynomial of degree N.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 June 2011
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Philip Davis, Philip Rabinowitz,
            //    Methods of Numerical Integration,
            //    Second Edition,
            //    Dover, 2007,
            //    ISBN: 0486453391,
            //    LC: QA299.3.D28.
            //
            //  Parameters:
            //
            //    Input, int ORDER, the order.
            //    ORDER must be greater than 0.
            //
            //    Output, double LEGENDRE_ZEROS[ORDER], the zeros.
            //
        {
            double d1;
            double d2pn;
            double d3pn;
            double d4pn;
            double dp;
            double dpn;
            double e1;
            //double fx;
            double h;
            int i;
            int iback;
            int k;
            int m;
            int mp1mi;
            int ncopy;
            int nmove;
            double p;
            double pk;
            double pkm1;
            double pkp1;
            const double r8_pi = 3.141592653589793;
            double t;
            double u;
            double v;
            double x0;
            double[] xtab;
            double xtemp;

            xtab = new double[order];

            e1 = (double) (order * (order + 1));

            m = (order + 1) / 2;

            for (i = 1; i <= m; i++)
            {
                mp1mi = m + 1 - i;

                t = (double) (4 * i - 1) * r8_pi / (double) (4 * order + 2);

                x0 = Math.Cos(t) * (1.0 - (1.0 - 1.0 / (double) (order))
                    / (double) (8 * order * order));

                pkm1 = 1.0;
                pk = x0;

                for (k = 2; k <= order; k++)
                {
                    pkp1 = 2.0 * x0 * pk - pkm1 - (x0 * pk - pkm1) / (double) (k);
                    pkm1 = pk;
                    pk = pkp1;
                }

                d1 = (double) (order) * (pkm1 - x0 * pk);

                dpn = d1 / (1.0 - x0 * x0);

                d2pn = (2.0 * x0 * dpn - e1 * pk) / (1.0 - x0 * x0);

                d3pn = (4.0 * x0 * d2pn + (2.0 - e1) * dpn) / (1.0 - x0 * x0);

                d4pn = (6.0 * x0 * d3pn + (6.0 - e1) * d2pn) / (1.0 - x0 * x0);

                u = pk / dpn;
                v = d2pn / dpn;
                //
                //  Initial approximation H:
                //
                h = -u * (1.0 + 0.5 * u * (v + u * (v * v - d3pn / (3.0 * dpn))));
                //
                //  Refine H using one step of Newton's method:
                //
                p = pk + h * (dpn + 0.5 * h * (d2pn + h / 3.0
                    * (d3pn + 0.25 * h * d4pn)));

                dp = dpn + h * (d2pn + 0.5 * h * (d3pn + h * d4pn / 3.0));

                h = h - p / dp;

                xtemp = x0 + h;

                xtab[mp1mi - 1] = xtemp;

                //  fx = d1 - h * e1 * ( pk + 0.5 * h * ( dpn + h / 3.0 
                //    * ( d2pn + 0.25 * h * ( d3pn + 0.2 * h * d4pn ) ) ) );
            }

            if ((order % 2) == 1)
            {
                xtab[0] = 0.0;
            }

            //
            //  Shift the data up.
            //
            nmove = (order + 1) / 2;
            ncopy = order - nmove;

            for (i = 1; i <= nmove; i++)
            {
                iback = order + 1 - i;
                xtab[iback - 1] = xtab[iback - ncopy - 1];
            }

            //
            //  Reflect values for the negative abscissas.
            //
            for (i = 1; i <= order - nmove; i++)
            {
                xtab[i - 1] = -xtab[order - i];
            }

            return xtab;
        }

        public static void legendre_poly(int n, double x, ref double[] cx, ref double[] cpx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_POLY evaluates the Legendre polynomials.
            //
            //  Discussion:
            //
            //    P(N,1) = 1.
            //    P(N,-1) = (-1)^N.
            //    | P(N,X) | <= 1 in [-1,1].
            //
            //    P(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
            //    function of the first kind and order N equals the Legendre polynomial
            //    of the first kind and order N.
            //
            //    The N zeroes of P(N,X) are the abscissas used for Gauss-Legendre
            //    quadrature of the integral of a function F(X) with weight function 1
            //    over the interval [-1,1].
            //
            //    The Legendre polynomials are orthogonal under the inner product defined
            //    as integration from -1 to 1:
            //
            //      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX
            //        = 0 if I =/= J
            //        = 2 / ( 2*I+1 ) if I = J.
            //
            //    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
            //
            //    A function F(X) defined on [-1,1] may be approximated by the series
            //      C0*P(0,X) + C1*P(1,X) + ... + CN*P(N,X)
            //    where
            //      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,X) dx.
            //
            //    The formula is:
            //
            //      P(N,X) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
            //
            //  Differential equation:
            //
            //    (1-X*X) * P(N,X)'' - 2 * X * P(N,X)' + N * (N+1) = 0
            //
            //  First terms:
            //
            //    P( 0,X) =       1
            //    P( 1,X) =       1 X
            //    P( 2,X) =  (    3 X^2 -       1)/2
            //    P( 3,X) =  (    5 X^3 -     3 X)/2
            //    P( 4,X) =  (   35 X^4 -    30 X^2 +     3)/8
            //    P( 5,X) =  (   63 X^5 -    70 X^3 +    15 X)/8
            //    P( 6,X) =  (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
            //    P( 7,X) =  (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
            //    P( 8,X) =  ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
            //    P( 9,X) =  (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
            //    P(10,X) =  (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2
            //                 -63 ) /256
            //
            //  Recursion:
            //
            //    P(0,X) = 1
            //    P(1,X) = X
            //    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
            //
            //    P'(0,X) = 0
            //    P'(1,X) = 1
            //    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Daniel Zwillinger, editor,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition,
            //    CRC Press, 1996.
            //
            //  Parameters:
            //
            //    Input, int N, the highest order polynomial to evaluate.
            //    Note that polynomials 0 through N will be evaluated.
            //
            //    Input, double X, the point at which the polynomials are to be evaluated.
            //
            //    Output, double CX[N+1], the values of the Legendre polynomials
            //    of order 0 through N at the point X.
            //
            //    Output, double CPX[N+1], the values of the derivatives of the
            //    Legendre polynomials of order 0 through N at the point X.
            //
        {
            int i;

            if (n < 0)
            {
                return;
            }

            cx[0] = 1.0;
            cpx[0] = 0.0;

            if (n < 1)
            {
                return;
            }

            cx[1] = x;
            cpx[1] = 1.0;

            for (i = 2; i <= n; i++)
            {
                cx[i] = ((double)(2 * i - 1) * x * cx[i - 1]
                         + (double)(-i + 1) * cx[i - 2])
                        / (double)(i);

                cpx[i] = ((double)(2 * i - 1) * (cx[i - 1] + x * cpx[i - 1])
                          + (double)(-i + 1) * cpx[i - 2])
                         / (double)(i);

            }
        }

        public static void legendre_poly_coef(int n, ref double[] c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_POLY_COEF evaluates the Legendre polynomial coefficients.
            //
            //  First terms:
            //
            //     1
            //     0     1
            //    -1/2   0      3/2
            //     0    -3/2    0     5/2
            //     3/8   0    -30/8   0     35/8
            //     0    15/8    0   -70/8    0     63/8
            //    -5/16  0    105/16  0   -315/16   0    231/16
            //     0   -35/16   0   315/16   0   -693/16   0    429/16
            //
            //     1.00000
            //     0.00000  1.00000
            //    -0.50000  0.00000  1.50000
            //     0.00000 -1.50000  0.00000  2.5000
            //     0.37500  0.00000 -3.75000  0.00000  4.37500
            //     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
            //    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
            //     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 February 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Daniel Zwillinger, editor,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition,
            //    CRC Press, 1996.
            //
            //  Parameters:
            //
            //    Input, int N, the highest order polynomial to evaluate.
            //    Note that polynomials 0 through N will be evaluated.
            //
            //    Output, double C[(N+1)*(N+1)], the coefficients of the Legendre polynomials
            //    of degree 0 through N.  Each polynomial is stored as a row.
            //
        {
            int i;
            int j;

            if (n < 0)
            {
                return;
            }

            for (i = 0; i <= n; i++)
            {
                for (j = 0; j <= n; j++)
                {
                    c[i + j * (n + 1)] = 0.0;
                }
            }

            c[0 + 0 * (n + 1)] = 1.0;

            if (n <= 0)
            {
                return;
            }

            c[1 + 1 * (n + 1)] = 1.0;

            for (i = 2; i <= n; i++)
            {
                for (j = 0; j <= i - 2; j++)
                {
                    c[i + j * (n + 1)] =
                        (double)(-i + 1) * c[i - 2 + j * (n + 1)] / (double)i;
                }

                for (j = 1; j <= i; j++)
                {
                    c[i + j * (n + 1)] = c[i + j * (n + 1)]
                                         + (double)(i + i - 1) * c[i - 1 + (j - 1) * (n + 1)] / (double)i;
                }
            }
        }

        public static void lp_coefficients(int n, ref int o, ref double[] c, ref int[] f)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LP_COEFFICIENTS: coefficients of Legendre polynomials P(n,x).
            //
            //  First terms:
            //
            //     1
            //     0     1
            //    -1/2   0      3/2
            //     0    -3/2    0     5/2
            //     3/8   0    -30/8   0     35/8
            //     0    15/8    0   -70/8    0     63/8
            //    -5/16  0    105/16  0   -315/16   0    231/16
            //     0   -35/16   0   315/16   0   -693/16   0    429/16
            //
            //     1.00000
            //     0.00000  1.00000
            //    -0.50000  0.00000  1.50000
            //     0.00000 -1.50000  0.00000  2.5000
            //     0.37500  0.00000 -3.75000  0.00000  4.37500
            //     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
            //    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
            //     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Daniel Zwillinger, editor,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition,
            //    CRC Press, 1996.
            //
            //  Parameters:
            //
            //    Input, int N, the highest order polynomial to evaluate.
            //    Note that polynomials 0 through N will be evaluated.
            //
            //    Output, int &O, the number of coefficients.
            //
            //    Output, double C[(N+2)/2], the coefficients of the Legendre
            //    polynomial of degree N.
            //
            //    Output, int F[(N+2)/2], the exponents.
            //
        {
            double[] ctable;
            int i;
            int j;
            int k;

            ctable = new double[(n + 1) * (n + 1)];

            for (i = 0; i <= n; i++)
            {
                for (j = 0; j <= n; j++)
                {
                    ctable[i + j * (n + 1)] = 0.0;
                }
            }

            ctable[0 + 0 * (n + 1)] = 1.0;

            if (0 < n)
            {
                ctable[1 + 1 * (n + 1)] = 1.0;

                for (i = 2; i <= n; i++)
                {
                    for (j = 0; j <= i - 2; j++)
                    {
                        ctable[i + j * (n + 1)] =
                            (double) (-i + 1) * ctable[i - 2 + j * (n + 1)] / (double) i;
                    }

                    for (j = 1; j <= i; j++)
                    {
                        ctable[i + j * (n + 1)] = ctable[i + j * (n + 1)]
                                                  + (double) (i + i - 1) * ctable[i - 1 + (j - 1) * (n + 1)] /
                                                  (double) i;
                    }
                }
            }

            //
            //  Extract the nonzero data from the alternating columns of the last row.
            //
            o = (n + 2) / 2;

            k = o;
            for (j = n; 0 <= j; j = j - 2)
            {
                k = k - 1;
                c[k] = ctable[n + j * (n + 1)];
                f[k] = j;
            }
        }

        public static double[] lp_value(int n, int o, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LP_VALUE evaluates the Legendre polynomials P(n,x).
            //
            //  Discussion:
            //
            //    P(n,1) = 1.
            //    P(n,-1) = (-1)^N.
            //    | P(n,x) | <= 1 in [-1,1].
            //
            //    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
            //    quadrature of the integral of a function F(X) with weight function 1
            //    over the interval [-1,1].
            //
            //    The Legendre polynomials are orthogonal under the inner product defined
            //    as integration from -1 to 1:
            //
            //      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX
            //        = 0 if I =/= J
            //        = 2 / ( 2*I+1 ) if I = J.
            //
            //    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
            //
            //    A function F(X) defined on [-1,1] may be approximated by the series
            //      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
            //    where
            //      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,x) dx.
            //
            //    The formula is:
            //
            //      P(n,x) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
            //
            //  Differential equation:
            //
            //    (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0
            //
            //  First terms:
            //
            //    P( 0,x) =      1
            //    P( 1,x) =      1 X
            //    P( 2,x) = (    3 X^2 -       1)/2
            //    P( 3,x) = (    5 X^3 -     3 X)/2
            //    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
            //    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
            //    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
            //    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
            //    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
            //    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
            //    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
            //
            //  Recursion:
            //
            //    P(0,x) = 1
            //    P(1,x) = x
            //    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n
            //
            //    P'(0,x) = 0
            //    P'(1,x) = 1
            //    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Daniel Zwillinger, editor,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition,
            //    CRC Press, 1996.
            //
            //  Parameters:
            //
            //    Input, int N, the number of evaluation points.
            //
            //    Input, int O, the degree of the polynomial.
            //
            //    Input, double X[N], the evaluation points.
            //
            //    Output, double LP_VALUE[N], the value of the Legendre polynomial 
            //    of degree N at the points X.
            //
        {
            int i;
            int j;
            double[] v;
            double[] vtable;

            vtable = new double[n * (o + 1)];

            for (i = 0; i < n; i++)
            {
                vtable[i + 0 * n] = 1.0;
            }

            if (1 <= o)
            {
                for (i = 0; i < n; i++)
                {
                    vtable[i + 1 * n] = x[i];
                }

                for (j = 2; j <= o; j++)
                {
                    for (i = 0; i < n; i++)
                    {
                        vtable[i + j * n] =
                            ((double) (2 * j - 1) * x[i] * vtable[i + (j - 1) * n]
                             - (double) (j - 1) * vtable[i + (j - 2) * n])
                            / (double) (j);
                    }
                }
            }

            v = new double[n];

            for (i = 0; i < n; i++)
            {
                v[i] = vtable[i + o * n];
            }

            return v;
        }


        public static void lpp_to_polynomial(int m, int[] l, int o_max, ref int o, ref double[] c, ref int[] e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LPP_TO_POLYNOMIAL writes a Legendre Product Polynomial as a polynomial.
            //
            //  Discussion:
            //
            //    For example, if 
            //      M = 3,
            //      L = ( 1, 0, 2 ),
            //    then
            //      L(1,0,2)(X,Y,Z) 
            //      = L(1)(X) * L(0)(Y) * L(2)(Z)
            //      = X * 1 * ( 3Z^2-1)/2
            //      = - 1/2 X + (3/2) X Z^2
            //    so
            //      O = 2 (2 nonzero terms)
            //      C = -0.5
            //           1.5
            //      E = 4    <-- index in 3-space of exponent (1,0,0)
            //          15   <-- index in 3-space of exponent (1,0,2)
            //
            //    The output value of O is no greater than
            //      O_MAX = product ( 1 <= I <= M ) (L(I)+2)/2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int L[M], the index of each Legendre product 
            //    polynomial factor.  0 <= L(*).
            //
            //    Input, int O_MAX, an upper limit on the size of the 
            //    output arrays.
            //      O_MAX = product ( 1 <= I <= M ) (L(I)+2)/2.
            //
            //    Output, int &O, the "order" of the polynomial product.
            //
            //    Output, double C[O], the coefficients of the polynomial product.
            //
            //    Output, int E[O], the indices of the exponents of the 
            //    polynomial product.
            //
        {
            double[] c1;
            double[] c2;
            int[] e1;
            int[] e2;
            int[] f2;
            int i;
            int i1;
            int i2;
            int j1;
            int j2;
            int o1;
            int o2 = 0;
            int[] p = new int[1];
            int[] pp;

            c1 = new double[o_max];
            c2 = new double[o_max];
            e1 = new int[o_max];
            e2 = new int[o_max];
            f2 = new int[o_max];
            pp = new int[m];

            o1 = 1;
            c1[0] = 1.0;
            e1[0] = 1;
            //
            //  Implicate one factor at a time.
            //
            for (i = 0; i < m; i++)
            {
                lp_coefficients(l[i], ref o2, ref c2, ref f2);

                o = 0;

                for (j2 = 0; j2 < o2; j2++)
                {
                    for (j1 = 0; j1 < o1; j1++)
                    {
                        c[o] = c1[j1] * c2[j2];
                        if (0 < i)
                        {
                            p = Monomial.mono_unrank_grlex(i, e1[j1]);
                        }

                        for (i2 = 0; i2 < i; i2++)
                        {
                            pp[i2] = p[i2];
                        }

                        pp[i] = f2[j2];
                        e[o] = Monomial.mono_rank_grlex(i + 1, pp);
                        o = o + 1;
                        if (0 < i)
                        {
                            p = null;
                        }
                    }
                }

                Polynomial.polynomial_sort(o, ref c, ref e);
                Polynomial.polynomial_compress(o, c, e, ref o, ref c, ref e);

                o1 = o;
                for (i1 = 0; i1 < o; i1++)
                {
                    c1[i1] = c[i1];
                    e1[i1] = e[i1];
                }
            }
        }
        
        public static double[] lpp_value ( int m, int n, int[] o, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LPP_VALUE evaluates a Legendre Product Polynomial at several points X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of evaluation points.
        //
        //    Input, int O[M], the degree of the polynomial factors.
        //    0 <= O(*).
        //
        //    Input, double X[M*N], the evaluation points.
        //
        //    Output, double LPP_VALUE[N], the value of the Legendre Product 
        //    Polynomial of degree O at the points X.
        //
        {
            int i;
            int j;
            double[] v;
            double[] vi;
            double[] xi;

            v = new double[n];

            for ( j = 0; j < n; j++ )
            {
                v[j] = 1.0;
            }

            xi = new double[n];

            for ( i = 0; i < m; i++ )
            {
                for ( j = 0; j < n; j++ )
                {
                    xi[j] = x[i+j*m];
                }
                vi = lp_value ( n, o[i], xi );
                for ( j = 0; j < n; j++ )
                {
                    v[j] = v[j] * vi[j];
                }
            }

            return v;
        }

        public static double[] p_exponential_product(int p, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P_EXPONENTIAL_PRODUCT: exponential products for P(n,x).
            //
            //  Discussion:
            //
            //    Let P(n,x) represent the Legendre polynomial of degree n.  
            //
            //    For polynomial chaos applications, it is of interest to know the
            //    value of the integrals of products of exp(B*X) with every possible pair
            //    of basis functions.  That is, we'd like to form
            //
            //      Tij = Integral ( -1.0 <= X <= +1.0 ) exp(B*X) * P(I,X) * P(J,X) dx
            //
            //    We will estimate these integrals using Gauss-Legendre quadrature.
            //    Because of the exponential factor exp(B*X), the quadrature will not 
            //    be exact.
            //
            //    However, when B = 0, the quadrature is exact, and moreoever, the
            //    table will be the identity matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int P, the maximum degree of the polyonomial 
            //    factors.  0 <= P.
            //
            //    Input, double B, the coefficient of X in the exponential factor.
            //
            //    Output, double P_EXPONENTIAL_PRODUCT[(P+1)*(P+1)], the table of integrals.  
            //
        {
            double[] h_table;
            int i;
            int j;
            int k;
            int order;
            double[] table;
            double[] w_table;
            double x;
            double[] x_table;

            table = new double[(p + 1) * (p + 1)];

            for (j = 0; j <= p; j++)
            {
                for (i = 0; i <= p; i++)
                {
                    table[i + j * (p + 1)] = 0.0;
                }
            }

            order = (3 * p + 4) / 2;

            x_table = new double[order];
            w_table = new double[order];

            LegendreQuadrature.p_quadrature_rule(order, ref x_table, ref w_table);

            for (k = 0; k < order; k++)
            {
                x = x_table[k];
                h_table = p_polynomial_value(1, p, x_table, xIndex: +k);
                //
                //  The following formula is an outer product in H_TABLE.
                //
                for (j = 0; j <= p; j++)
                {
                    for (i = 0; i <= p; i++)
                    {
                        table[i + j * (p + 1)] = table[i + j * (p + 1)]
                                                 + w_table[k] * Math.Exp(b * x) * h_table[i] * h_table[j];
                    }
                }
            }

            return table;
        }

        public static double p_integral(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P_INTEGRAL evaluates a monomial integral associated with P(n,x).
            //
            //  Discussion:
            //
            //    The integral:
            //
            //      integral ( -1 <= x < +1 ) x^n dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the exponent.
            //    0 <= N.
            //
            //    Output, double P_INTEGRAL, the value of the integral.
            //
        {
            double value;

            if ((n % 2) == 1)
            {
                value = 0.0;
            }
            else
            {
                value = 2.0 / (double) (n + 1);
            }

            return value;
        }

        public static double[] p_polynomial_coefficients(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P_POLYNOMIAL_COEFFICIENTS: coefficients of Legendre polynomial P(n,x).
            //
            //  Discussion:
            //
            //     1
            //     0     1
            //    -1/2   0      3/2
            //     0    -3/2    0     5/2
            //     3/8   0    -30/8   0     35/8
            //     0    15/8    0   -70/8    0     63/8
            //    -5/16  0    105/16  0   -315/16   0    231/16
            //     0   -35/16   0   315/16   0   -693/16   0    429/16
            //
            //     1.00000
            //     0.00000  1.00000
            //    -0.50000  0.00000  1.50000
            //     0.00000 -1.50000  0.00000  2.5000
            //     0.37500  0.00000 -3.75000  0.00000  4.37500
            //     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
            //    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
            //     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Daniel Zwillinger, editor,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition,
            //    CRC Press, 1996.
            //
            //  Parameters:
            //
            //    Input, int N, the highest order polynomial to evaluate.
            //    Note that polynomials 0 through N will be evaluated.
            //
            //    Output, double P_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients of 
            //    the Legendre polynomials of degree 0 through N.
            //
        {
            double[] c;
            int i;
            int j;

            if (n < 0)
            {
                return null;
            }

            c = new double[(n + 1) * (n + 1)];

            for (i = 0; i <= n; i++)
            {
                for (j = 0; j <= n; j++)
                {
                    c[i + j * (n + 1)] = 0.0;
                }
            }

            c[0 + 0 * (n + 1)] = 1.0;

            if (0 < n)
            {
                c[1 + 1 * (n + 1)] = 1.0;
            }

            for (i = 2; i <= n; i++)
            {
                for (j = 0; j <= i - 2; j++)
                {
                    c[i + j * (n + 1)] =
                        (double) (-i + 1) * c[i - 2 + j * (n + 1)] / (double) i;
                }

                for (j = 1; j <= i; j++)
                {
                    c[i + j * (n + 1)] = c[i + j * (n + 1)]
                                         + (double) (i + i - 1) * c[i - 1 + (j - 1) * (n + 1)] / (double) i;
                }
            }

            return c;
        }

        public static double[] p_polynomial_prime(int m, int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P_POLYNOMIAL_PRIME evaluates the derivative of Legendre polynomials P(n,x).
            //
            //  Discussion:
            //
            //    P(0,X) = 1
            //    P(1,X) = X
            //    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
            //
            //    P'(0,X) = 0
            //    P'(1,X) = 1
            //    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
            //
            //    Thanks to Dimitriy Morozov for pointing out a memory leak caused by
            //    not deleting the work array V before return, 19 March 2013.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 March 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Daniel Zwillinger, editor,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition,
            //    CRC Press, 1996.
            //
            //  Parameters:
            //
            //    Input, int M, the number of evaluation points.
            //
            //    Input, int N, the highest order polynomial to evaluate.
            //    Note that polynomials 0 through N will be evaluated.
            //
            //    Input, double X[M], the evaluation points.
            //
            //    Output, double P_POLYNOMIAL_PRIME[M*(N+1)], the values of the derivatives 
            //    of the Legendre polynomials of order 0 through N at the points.
            //
        {
            int i;
            int j;
            double[] v;
            double[] vp;

            if (n < 0)
            {
                return null;
            }

            vp = new double[m * (n + 1)];

            for (i = 0; i < m; i++)
            {
                vp[i + 0 * m] = 0.0;
            }

            if (n < 1)
            {
                return vp;
            }

            v = new double[m * (n + 1)];

            for (i = 0; i < m; i++)
            {
                v[i + 0 * m] = 1.0;
            }

            for (i = 0; i < m; i++)
            {
                v[i + 1 * m] = x[i];
                vp[i + 1 * m] = 1.0;
            }

            for (j = 2; j <= n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    v[i + j * m] = ((double) (2 * j - 1) * x[i] * v[i + (j - 1) * m]
                                    - (double) (j - 1) * v[i + (j - 2) * m])
                                   / (double) (j);

                    vp[i + j * m] = ((double) (2 * j - 1) * (v[i + (j - 1) * m] + x[i] * vp[i + (j - 1) * m])
                                     - (double) (j - 1) * vp[i + (j - 2) * m])
                                    / (double) (j);
                }
            }

            return vp;
        }

        public static double[] p_polynomial_prime2(int m, int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P_POLYNOMIAL_PRIME2: second derivative of Legendre polynomials P(n,x).
            //
            //  Discussion:
            //
            //    P(0,X) = 1
            //    P(1,X) = X
            //    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
            //
            //    P'(0,X) = 0
            //    P'(1,X) = 1
            //    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
            //
            //    P"(0,X) = 0
            //    P"(1,X) = 0
            //    P"(N,X) = ( (2*N-1)*(2*P'(N-1,X)+X*P"(N-1,X)-(N-1)*P"(N-2,X) ) / N
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 May 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Daniel Zwillinger, editor,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition,
            //    CRC Press, 1996.
            //
            //  Parameters:
            //
            //    Input, int M, the number of evaluation points.
            //
            //    Input, int N, the highest order polynomial to evaluate.
            //    Note that polynomials 0 through N will be evaluated.
            //
            //    Input, double X[M], the evaluation points.
            //
            //    Output, double P_POLYNOMIAL_PRIME2[M*(N+1)], the second derivative
            //    of the Legendre polynomials of order 0 through N at the points.
            //
        {
            int i;
            int j;
            double[] v;
            double[] vp;
            double[] vpp;

            if (n < 0)
            {
                return null;
            }

            vpp = new double[m * (n + 1)];

            for (i = 0; i < m; i++)
            {
                vpp[i + 0 * m] = 0.0;
            }

            if (n < 1)
            {
                return vpp;
            }

            v = new double[m * (n + 1)];
            vp = new double[m * (n + 1)];

            for (i = 0; i < m; i++)
            {
                v[i + 0 * m] = 1.0;
                vp[i + 0 * m] = 0.0;
            }

            for (i = 0; i < m; i++)
            {
                v[i + 1 * m] = x[i];
                vp[i + 1 * m] = 1.0;
                vpp[i + 1 * m] = 0.0;
            }

            for (j = 2; j <= n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    v[i + j * m] =
                        ((double) (2 * j - 1) * x[i] * v[i + (j - 1) * m]
                         - (double) (j - 1) * v[i + (j - 2) * m])
                        / (double) (j);

                    vp[i + j * m] =
                        ((double) (2 * j - 1) * (v[i + (j - 1) * m] + x[i] * vp[i + (j - 1) * m])
                         - (double) (j - 1) * vp[i + (j - 2) * m])
                        / (double) (j);

                    vpp[i + j * m] =
                        ((double) (2 * j - 1) * (2.0 * vp[i + (j - 1) * m] + x[i] * vpp[i + (j - 1) * m])
                         - (double) (j - 1) * vpp[i + (j - 2) * m])
                        / (double) (j);
                }
            }

            return vpp;
        }

        public static double[] p_polynomial_value(int m, int n, double[] x, int xIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P_POLYNOMIAL_VALUE evaluates the Legendre polynomials P(n,x).
            //
            //  Discussion:
            //
            //    P(n,1) = 1.
            //    P(n,-1) = (-1)^N.
            //    | P(n,x) | <= 1 in [-1,1].
            //
            //    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
            //    quadrature of the integral of a function F(X) with weight function 1
            //    over the interval [-1,1].
            //
            //    The Legendre polynomials are orthogonal under the inner product defined
            //    as integration from -1 to 1:
            //
            //      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX 
            //        = 0 if I =/= J
            //        = 2 / ( 2*I+1 ) if I = J.
            //
            //    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
            //
            //    A function F(X) defined on [-1,1] may be approximated by the series
            //      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
            //    where
            //      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,x) dx.
            //
            //    The formula is:
            //
            //      P(n,x) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
            //
            //  Differential equation:
            //
            //    (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0
            //
            //  First terms:
            //
            //    P( 0,x) =      1
            //    P( 1,x) =      1 X
            //    P( 2,x) = (    3 X^2 -       1)/2
            //    P( 3,x) = (    5 X^3 -     3 X)/2
            //    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
            //    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
            //    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
            //    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
            //    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
            //    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
            //    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
            //
            //  Recursion:
            //
            //    P(0,x) = 1
            //    P(1,x) = x
            //    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n
            //
            //    P'(0,x) = 0
            //    P'(1,x) = 1
            //    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Daniel Zwillinger, editor,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition,
            //    CRC Press, 1996.
            //
            //  Parameters:
            //
            //    Input, int M, the number of evaluation points.
            //
            //    Input, int N, the highest order polynomial to evaluate.
            //    Note that polynomials 0 through N will be evaluated.
            //
            //    Input, double X[M], the evaluation points.
            //
            //    Output, double P_POLYNOMIAL_VALUE[M*(N+1)], the values of the Legendre
            //    polynomials of order 0 through N.
            //
        {
            int i;
            int j;
            double[] v;

            if (n < 0)
            {
                return null;
            }

            v = new double[m * (n + 1)];

            for (i = 0; i < m; i++)
            {
                v[i + 0 * m] = 1.0;
            }

            if (n < 1)
            {
                return v;
            }

            for (i = 0; i < m; i++)
            {
                v[i + 1 * m] = x[xIndex + i];
            }

            for (j = 2; j <= n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    v[i + j * m] = ((double) (2 * j - 1) * x[xIndex + i] * v[i + (j - 1) * m]
                                    - (double) (j - 1) * v[i + (j - 2) * m])
                                   / (double) (j);
                }
            }

            return v;
        }


        public static double[] p_polynomial_zeros(int nt)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P_POLYNOMIAL_ZEROS: zeros of Legendre function P(n,x).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NT, the order of the rule.
            //
            //    Output, double P_POLYNOMIAL_ZEROS[NT], the zeros.
            //
        {
            double[] bj;
            int i;
            double[] t;
            double[] wts;

            t = new double[nt];

            for (i = 0; i < nt; i++)
            {
                t[i] = 0.0;
            }

            bj = new double[nt];
            for (i = 0; i < nt; i++)
            {
                bj[i] = (double) ((i + 1) * (i + 1))
                        / (double) (4 * (i + 1) * (i + 1) - 1);
                bj[i] = Math.Sqrt(bj[i]);
            }

            wts = new double[nt];
            wts[0] = Math.Sqrt(2.0);
            for (i = 1; i < nt; i++)
            {
                wts[i] = 0.0;
            }

            IMTQLX.imtqlx(nt, ref t, ref bj, ref wts);

            return t;
        }

        public static double[] p_power_product(int p, int e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P_POWER_PRODUCT: power products for Legendre polynomial P(n,x).
            //
            //  Discussion:
            //
            //    Let P(n,x) represent the Legendre polynomial of degree n.  
            //
            //    For polynomial chaos applications, it is of interest to know the
            //    value of the integrals of products of X with every possible pair
            //    of basis functions.  That is, we'd like to form
            //
            //      Tij = Integral ( -1.0 <= X <= +1.0 ) X^E * P(I,x) * P(J,x) dx
            //
            //    We will estimate these integrals using Gauss-Legendre quadrature.
            //
            //    When E is 0, the computed table should be the identity matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int P, the maximum degree of the polyonomial 
            //    factors.  0 <= P.
            //
            //    Input, int E, the exponent of X in the integrand.
            //    0 <= E.
            //
            //    Output, double P_POWER_PRODUCT[(P+1)*(P+1)], the table of integrals.  
            //
        {
            double[] h_table;
            int i;
            int j;
            int k;
            int order;
            double[] table;
            double[] w_table;
            double x;
            double[] x_table;

            table = new double[(p + 1) * (p + 1)];

            for (j = 0; j <= p; j++)
            {
                for (i = 0; i <= p; i++)
                {
                    table[i + j * (p + 1)] = 0.0;
                }
            }

            order = p + 1 + ((e + 1) / 2);

            x_table = new double[order];
            w_table = new double[order];

            LegendreQuadrature.p_quadrature_rule(order, ref x_table, ref w_table);

            for (k = 0; k < order; k++)
            {
                x = x_table[k];
                h_table = p_polynomial_value(1, p, x_table, xIndex: +k);
                //
                //  The following formula is an outer product in H_TABLE.
                //
                if (e == 0)
                {
                    for (i = 0; i <= p; i++)
                    {
                        for (j = 0; j <= p; j++)
                        {
                            table[i + j * (p + 1)] = table[i + j * (p + 1)] + w_table[k] * h_table[i] * h_table[j];
                        }
                    }
                }
                else
                {
                    for (i = 0; i <= p; i++)
                    {
                        for (j = 0; j <= p; j++)
                        {
                            table[i + j * (p + 1)] = table[i + j * (p + 1)]
                                                     + w_table[k] * Math.Pow(x, e) * h_table[i] * h_table[j];
                        }
                    }
                }
            }

            return table;
        }

        public static double[] pm_polynomial_value(int mm, int n, int m, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PM_POLYNOMIAL_VALUE evaluates the Legendre polynomials Pm(n,m,x).
            //
            //  Differential equation:
            //
            //    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
            //
            //  First terms:
            //
            //    M = 0  ( = Legendre polynomials of first kind P(N,X) )
            //
            //    Pm(0,0,x) =    1
            //    Pm(1,0,x) =    1 X
            //    Pm(2,0,x) = (  3 X^2 -   1)/2
            //    Pm(3,0,x) = (  5 X^3 -   3 X)/2
            //    Pm(4,0,x) = ( 35 X^4 -  30 X^2 +   3)/8
            //    Pm(5,0,x) = ( 63 X^5 -  70 X^3 +  15 X)/8
            //    Pm(6,0,x) = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
            //    Pm(7,0,x) = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
            //
            //    M = 1
            //
            //    Pm(0,1,x) =   0
            //    Pm(1,1,x) =   1 * SQRT(1-X^2)
            //    Pm(2,1,x) =   3 * SQRT(1-X^2) * X
            //    Pm(3,1,x) = 1.5 * SQRT(1-X^2) * (5*X^2-1)
            //    Pm(4,1,x) = 2.5 * SQRT(1-X^2) * (7*X^3-3*X)
            //
            //    M = 2
            //
            //    Pm(0,2,x) =   0
            //    Pm(1,2,x) =   0
            //    Pm(2,2,x) =   3 * (1-X^2)
            //    Pm(3,2,x) =  15 * (1-X^2) * X
            //    Pm(4,2,x) = 7.5 * (1-X^2) * (7*X^2-1)
            //
            //    M = 3
            //
            //    Pm(0,3,x) =   0
            //    Pm(1,3,x) =   0
            //    Pm(2,3,x) =   0
            //    Pm(3,3,x) =  15 * (1-X^2)^1.5
            //    Pm(4,3,x) = 105 * (1-X^2)^1.5 * X
            //
            //    M = 4
            //
            //    Pm(0,4,x) =   0
            //    Pm(1,4,x) =   0
            //    Pm(2,4,x) =   0
            //    Pm(3,4,x) =   0
            //    Pm(4,4,x) = 105 * (1-X^2)^2
            //
            //  Recursion:
            //
            //    if N < M:
            //      Pm(N,M,x) = 0
            //    if N = M:
            //      Pm(N,M,x) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
            //      all the odd integers less than or equal to N.
            //    if N = M+1:
            //      Pm(N,M,x) = X*(2*M+1)*Pm(M,M,x)
            //    if M+1 < N:
            //      Pm(N,M,x) = ( X*(2*N-1)*Pm(N-1,M,x) - (N+M-1)*Pm(N-2,M,x) )/(N-M)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //  Parameters:
            //
            //    Input, int MM, the number of evaluation points.
            //
            //    Input, int N, the maximum first index of the Legendre
            //    function, which must be at least 0.
            //
            //    Input, int M, the second index of the Legendre function,
            //    which must be at least 0, and no greater than N.
            //
            //    Input, double X[MM], the point at which the function is to be
            //    evaluated.
            //
            //    Output, double PM_POLYNOMIAL_VALUE[MM*(N+1)], the function values.
            //
        {
            double fact;
            int i;
            int j;
            int k;
            double[] v;

            v = new double[mm * (n + 1)];

            for (j = 0; j < n + 1; j++)
            {
                for (i = 0; i < mm; i++)
                {
                    v[i + j * mm] = 0.0;
                }
            }

            //
            //  J = M is the first nonzero function.
            //
            if (m <= n)
            {
                for (i = 0; i < mm; i++)
                {
                    v[i + m * mm] = 1.0;
                }

                fact = 1.0;
                for (k = 0; k < m; k++)
                {
                    for (i = 0; i < mm; i++)
                    {
                        v[i + m * mm] = -v[i + m * mm] * fact * Math.Sqrt(1.0 - x[i] * x[i]);
                    }

                    fact = fact + 2.0;
                }
            }

            //
            //  J = M + 1 is the second nonzero function.
            //
            if (m + 1 <= n)
            {
                for (i = 0; i < mm; i++)
                {
                    v[i + (m + 1) * mm] = x[i] * (double) (2 * m + 1) * v[i + m * mm];
                }
            }

            //
            //  Now we use a three term recurrence.
            //
            for (j = m + 2; j <= n; j++)
            {
                for (i = 0; i < mm; i++)
                {
                    v[i + j * mm] = ((double) (2 * j - 1) * x[i] * v[i + (j - 1) * mm]
                                     + (double) (-j - m + 1) * v[i + (j - 2) * mm])
                                    / (double) (j - m);
                }
            }

            return v;
        }


        public static double[] pmn_polynomial_value(int mm, int n, int m, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PMN_POLYNOMIAL_VALUE: normalized Legendre polynomial Pmn(n,m,x).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //  Parameters:
            //
            //    Input, int MM, the number of evaluation points.
            //
            //    Input, int N, the maximum first index of the Legendre
            //    function, which must be at least 0.
            //
            //    Input, int M, the second index of the Legendre function,
            //    which must be at least 0, and no greater than N.
            //
            //    Input, double X[MM], the evaluation points.
            //
            //    Output, double PMN_POLYNOMIAL_VALUE[MM*(N+1)], the function values.
            //
        {
            double factor;
            int i;
            int j;
            double[] v;

            v = pm_polynomial_value(mm, n, m, x);
            //
            //  Normalization.
            //
            for (j = m; j <= n; j++)
            {
                factor = Math.Sqrt(((double) (2 * j + 1) * typeMethods.r8_factorial(j - m))
                              / (2.0 * typeMethods.r8_factorial(j + m)));
                for (i = 0; i < mm; i++)
                {
                    v[i + j * mm] = v[i + j * mm] * factor;
                }
            }

            return v;
        }


        public static double[] pmns_polynomial_value(int mm, int n, int m, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PMNS_POLYNOMIAL_VALUE: sphere-normalized Legendre polynomial Pmn(n,m,x).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //  Parameters:
            //
            //    Input, int MM, the number of evaluation points.
            //
            //    Input, int N, the maximum first index of the Legendre
            //    function, which must be at least 0.
            //
            //    Input, int M, the second index of the Legendre function,
            //    which must be at least 0, and no greater than N.
            //
            //    Input, double X[MM], the evaluation points.
            //
            //    Output, double PMNS_POLYNOMIAL_VALUE[MM*(N+1)], the function values.
            //
        {
            double factor;
            int i;
            int j;
            const double pi = 3.141592653589793;
            double[] v;

            v = pm_polynomial_value(mm, n, m, x);
            //
            //  Normalization.
            //
            for (j = m; j <= n; j++)
            {
                factor = Math.Sqrt(((double) (2 * j + 1) * typeMethods.r8_factorial(j - m))
                              / (4.0 * pi * typeMethods.r8_factorial(j + m)));
                for (i = 0; i < mm; i++)
                {
                    v[i + j * mm] = v[i + j * mm] * factor;
                }
            }

            return v;
        }


        public static double[] pn_pair_product(int p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PN_PAIR_PRODUCT: pair products for normalized Legendre polynomial Pn(n,x).
            //
            //  Discussion:
            //
            //    Let P(n,x) represent the Legendre polynomial of degree n.  
            //
            //    To check orthonormality, we compute
            //
            //      Tij = Integral ( -1.0 <= X <= +1.0 ) Pn(i,x) * Pn(j,x) dx
            //
            //    We will estimate these integrals using Gauss-Legendre quadrature.
            //
            //    The computed table should be the identity matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int P, the maximum degree of the polyonomial 
            //    factors.  0 <= P.
            //
            //    Input, int E, the exponent of X in the integrand.
            //    0 <= E.
            //
            //    Output, double PN_PAIR_PRODUCT[(P+1)*(P+1)], the table of integrals.  
            //
        {
            double[] h_table;
            int i;
            int j;
            int k;
            int order;
            double[] table;
            double[] w_table;
            double[] x_table;

            table = new double[(p + 1) * (p + 1)];

            for (j = 0; j <= p; j++)
            {
                for (i = 0; i <= p; i++)
                {
                    table[i + j * (p + 1)] = 0.0;
                }
            }

            order = p + 1;

            x_table = new double[order];
            w_table = new double[order];

            LegendreQuadrature.p_quadrature_rule(order, ref x_table, ref w_table);

            for (k = 0; k < order; k++)
            {
                h_table = pn_polynomial_value(1, p, x_table, xIndex: + k);

                for (i = 0; i <= p; i++)
                {
                    for (j = 0; j <= p; j++)
                    {
                        table[i + j * (p + 1)] = table[i + j * (p + 1)]
                                                 + w_table[k] * h_table[i] * h_table[j];
                    }
                }
            }

            return table;
        }

        public static double[] pn_polynomial_coefficients(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PN_POLYNOMIAL_COEFFICIENTS: coefficients of normalized Legendre Pn(n,x).
            //
            //  Discussion:
            //
            //    Pn(n,x) = P(n,x) * sqrt ( (2n+1)/2 )
            //
            //          1       x       x^2     x^3     x^4      x^5    x^6     x^7
            //
            //    0   0.707
            //    1   0.000   1.224
            //    2  -0.790   0.000   2.371
            //    3   0.000  -2.806   0.000   4.677
            //    4   0.795   0.000  -7.954   0.000   9.280
            //    5   0.000   4.397   0.000 -20.520   0.000   18.468
            //    6  -0.796   0.000  16.731   0.000 -50.193    0.000  36.808
            //    7   0.000  -5.990   0.000  53.916   0.000 -118.616   0.000  73.429 
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Daniel Zwillinger, editor,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition,
            //    CRC Press, 1996.
            //
            //  Parameters:
            //
            //    Input, int N, the highest order polynomial to evaluate.
            //    Note that polynomials 0 through N will be evaluated.
            //
            //    Output, double PN_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients of 
            //    the normalized Legendre polynomials of degree 0 through N.
            //
        {
            double[] c;
            int i;
            int j;
            double t;

            if (n < 0)
            {
                return null;
            }

            //
            //  Compute P(i,x) coefficients.
            //
            c = new double[(n + 1) * (n + 1)];

            for (i = 0; i <= n; i++)
            {
                for (j = 0; j <= n; j++)
                {
                    c[i + j * (n + 1)] = 0.0;
                }
            }

            c[0 + 0 * (n + 1)] = 1.0;

            if (0 < n)
            {
                c[1 + 1 * (n + 1)] = 1.0;
            }

            for (i = 2; i <= n; i++)
            {
                for (j = 0; j <= i - 2; j++)
                {
                    c[i + j * (n + 1)] =
                        (double) (-i + 1) * c[i - 2 + j * (n + 1)] / (double) i;
                }

                for (j = 1; j <= i; j++)
                {
                    c[i + j * (n + 1)] = c[i + j * (n + 1)]
                                         + (double) (i + i - 1) * c[i - 1 + (j - 1) * (n + 1)] / (double) i;
                }
            }

            //
            //  Normalize them.
            //
            for (i = 0; i <= n; i++)
            {
                t = Math.Sqrt((double) (2 * i + 1) / 2.0);
                for (j = 0; j <= i; j++)
                {
                    c[i + j * (n + 1)] = c[i + j * (n + 1)] * t;
                }
            }

            return c;
        }

        public static double[] pn_polynomial_value(int m, int n, double[] x, int xIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PN_POLYNOMIAL_VALUE evaluates the normalized Legendre polynomials Pn(n,x).
            //
            //  Discussion:
            //
            //    The normalized Legendre polynomials are orthonormal under the inner product 
            //    defined as integration from -1 to 1:
            //
            //      Integral ( -1 <= x <= +1 ) Pn(i,x) * Pn(j,x) dx = delta(i,j)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Daniel Zwillinger, editor,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition,
            //    CRC Press, 1996.
            //
            //  Parameters:
            //
            //    Input, int M, the number of evaluation points.
            //
            //    Input, int N, the highest order polynomial to evaluate.
            //    Note that polynomials 0 through N will be evaluated.
            //
            //    Input, double X[M], the evaluation points.
            //
            //    Output, double PN_POLYNOMIAL_VALUE[M*(N+1)], the values of the Legendre
            //    polynomials of order 0 through N.
            //
        {
            int i;
            int j;
            double norm;
            double[] v;

            v = p_polynomial_value(m, n, x, xIndex:xIndex);

            for (j = 0; j <= n; j++)
            {
                norm = Math.Sqrt(2 / (double) (2 * j + 1));
                for (i = 0; i < m; i++)
                {
                    v[i + j * m] = v[i + j * m] / norm;
                }
            }

            return v;
        }
    }
}