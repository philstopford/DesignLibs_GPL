using System;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.Laguerre;

public static partial class Functions
{
    public static double lm_integral(int n, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LM_INTEGRAL evaluates a monomial integral associated with Lm(n,m,x).
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( 0 <= x < +oo ) x^n * x^m * exp ( -x ) dx
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
        //  Parameters:
        //
        //    Input, int N, the exponent.
        //    0 <= N.
        //
        //    Input, int M, the parameter.
        //    0 <= M.
        //
        //    Output, double LM_INTEGRAL, the value of the integral.
        //
    {
        double value = 0;

        value = typeMethods.r8_factorial(n + m);

        return value;
    }

    public static double[] lm_polynomial(int mm, int n, int m, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LM_POLYNOMIAL evaluates Laguerre polynomials Lm(n,m,x).
        //
        //  First terms:
        //
        //    M = 0
        //
        //    Lm(0,0,X) =   1
        //    Lm(1,0,X) =  -X   +  1
        //    Lm(2,0,X) =   X^2 -  4 X   +  2
        //    Lm(3,0,X) =  -X^3 +  9 X^2 -  18 X   +    6
        //    Lm(4,0,X) =   X^4 - 16 X^3 +  72 X^2 -   96 X +     24
        //    Lm(5,0,X) =  -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x   +  120
        //    Lm(6,0,X) =   X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720
        //
        //    M = 1
        //
        //    Lm(0,1,X) =    0
        //    Lm(1,1,X) =   -1,
        //    Lm(2,1,X) =    2 X - 4,
        //    Lm(3,1,X) =   -3 X^2 + 18 X - 18,
        //    Lm(4,1,X) =    4 X^3 - 48 X^2 + 144 X - 96
        //
        //    M = 2
        //
        //    Lm(0,2,X) =    0
        //    Lm(1,2,X) =    0,
        //    Lm(2,2,X) =    2,
        //    Lm(3,2,X) =   -6 X + 18,
        //    Lm(4,2,X) =   12 X^2 - 96 X + 144
        //
        //    M = 3
        //
        //    Lm(0,3,X) =    0
        //    Lm(1,3,X) =    0,
        //    Lm(2,3,X) =    0,
        //    Lm(3,3,X) =   -6,
        //    Lm(4,3,X) =   24 X - 96
        //
        //    M = 4
        //
        //    Lm(0,4,X) =    0
        //    Lm(1,4,X) =    0
        //    Lm(2,4,X) =    0
        //    Lm(3,4,X) =    0
        //    Lm(4,4,X) =   24
        //
        //  Recursion:
        //
        //    Lm(0,M,X)   = 1 
        //    Lm(1,M,X)   = (M+1-X)
        //
        //    if 2 <= N:
        //
        //      Lm(N,M,X)   = ( (M+2*N-1-X) * Lm(N-1,M,X) 
        //                   +   (1-M-N)    * Lm(N-2,M,X) ) / N
        //
        //  Special values:
        //
        //    For M = 0, the associated Laguerre polynomials Lm(N,M,X) are equal 
        //    to the Laguerre polynomials L(N,X).
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
        //  Parameters:
        //
        //    Input, int MM, the number of evaluation points.
        //
        //    Input, int N, the highest order polynomial to compute.
        //
        //    Input, int M, the parameter.  M must be nonnegative.
        //
        //    Input, double X[MM], the evaluation points.
        //
        //    Output, double LM_POLYNOMIAL[MM*(N+1)], the function values.
        //
    {
        int i;
        int j;
        double[] v;

        switch (n)
        {
            case < 0:
                return null;
        }

        v = new double[mm * (n + 1)];

        for (j = 0; j <= n; j++)
        {
            for (i = 0; i < mm; i++)
            {
                v[i + j * mm] = 0.0;
            }
        }

        for (i = 0; i < mm; i++)
        {
            v[i + 0 * mm] = 1.0;
        }

        switch (n)
        {
            case 0:
                return v;
        }

        for (i = 0; i < mm; i++)
        {
            v[i + 1 * mm] = m + 1 - x[i];
        }

        for (j = 2; j <= n; j++)
        {
            for (i = 0; i < mm; i++)
            {
                v[i + j * mm] = ((m + 2 * j - 1 - x[i]) * v[i + (j - 1) * mm]
                                 + (-m - j + 1) * v[i + (j - 2) * mm])
                                / j;
            }
        }

        return v;
    }

    public static double[] lm_polynomial_coefficients(int n, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LM_POLYNOMIAL_COEFFICIENTS: coefficients of Laguerre polynomial Lm(n,m,x).
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
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to compute.
        //    Note that polynomials 0 through N will be computed.
        //
        //    Input, int M, the parameter.
        //
        //    Output, double LM_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients
        //    of the Laguerre polynomials of degree 0 through N. 
        //
    {
        double[] c;
        int i;
        int j;

        switch (n)
        {
            case < 0:
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

        switch (n)
        {
            case 0:
                return c;
        }

        c[1 + 0 * (n + 1)] = m + 1;
        c[1 + 1 * (n + 1)] = -1.0;

        for (i = 2; i <= n; i++)
        {
            for (j = 0; j <= i; j++)
            {
                c[i + j * (n + 1)] = (
                                         (m + 2 * i - 1) * c[i - 1 + j * (n + 1)]
                                         + (-m - i + 1) * c[i - 2 + j * (n + 1)])
                                     / i;
            }

            for (j = 1; j <= i; j++)
            {
                c[i + j * (n + 1)] -= c[i - 1 + (j - 1) * (n + 1)] / i;
            }
        }

        return c;
    }

    public static void lm_polynomial_values(ref int n_data, ref int n, ref int m, ref double x,
            ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LM_POLYNOMIAL_VALUES: some values of the Laguerre polynomial Lm(n,m,x).
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      LaguerreL[n,m,x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &N, the order of the function.
        //
        //    Output, int &M, the parameter.
        //
        //    Output, double &X, the point where the function is evaluated.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 20;

        double[] fx_vec =
            {
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1500000000000000E+01,
                0.1625000000000000E+01,
                0.1479166666666667E+01,
                0.1148437500000000E+01,
                0.4586666666666667E+00,
                0.2878666666666667E+01,
                0.8098666666666667E+01,
                0.1711866666666667E+02,
                0.1045328776041667E+02,
                0.1329019368489583E+02,
                0.5622453647189670E+02,
                0.7484729341779436E+02,
                0.3238912982762806E+03,
                0.4426100000097533E+03,
                0.1936876572288250E+04
            }
            ;

        int[] m_vec =
            {
                0, 0, 0, 0,
                0, 1, 1, 1,
                1, 0, 1, 2,
                3, 2, 2, 3,
                3, 4, 4, 5
            }
            ;

        int[] n_vec =
            {
                1, 2, 3, 4,
                5, 1, 2, 3,
                4, 3, 3, 3,
                3, 4, 5, 6,
                7, 8, 9, 10
            }
            ;

        double[] x_vec =
            {
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00
            }
            ;

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            n = 0;
            m = 0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            n = n_vec[n_data - 1];
            m = m_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static double[] lm_polynomial_zeros(int n, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LM_POLYNOMIAL_ZEROS returns the zeros for Lm(n,m,x).
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
        //    John Burkardt.
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
        //    Input, int N, the order.
        //
        //    Input, int M, the parameter.
        //    0 <= M.
        //
        //    Output, double X[N], the zeros.
        //
    {
        double[] bj;
        int i;
        double[] w;
        double[] x;
        double zemu;
        //
        //  Define the zero-th moment.
        //
        zemu = typeMethods.r8_factorial(m);
        //
        //  Define the Jacobi matrix.
        //
        bj = new double[n];
        for (i = 0; i < n; i++)
        {
            bj[i] = (double) (i + 1) * (i + 1 + m);
        }

        x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = 2 * i + 1 + m;
        }

        w = new double[n];
        w[0] = Math.Sqrt(zemu);
        for (i = 1; i < n; i++)
        {
            w[i] = 0.0;
        }

        //
        //  Diagonalize the Jacobi matrix.
        //
        IMTQLX.imtqlx(n, ref x, ref bj, ref w);

        return x;
    }

    public static void lm_quadrature_rule(int n, int m, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LM_QUADRATURE_RULE: Gauss-Laguerre quadrature rule for Lm(n,m,x);
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
        //    John Burkardt.
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
        //    Input, int N, the order.
        //
        //    Input, int M, the parameter.
        //    0 <= M.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        double[] bj;
        int i;
        double zemu;
        //
        //  Define the zero-th moment.
        //
        zemu = typeMethods.r8_factorial(m);
        //
        //  Define the Jacobi matrix.
        //
        bj = new double[n];
        for (i = 0; i < n; i++)
        {
            bj[i] = (double) (i + 1) * (i + 1 + m);
        }

        for (i = 0; i < n; i++)
        {
            x[i] = 2 * i + 1 + m;
        }

        w[0] = Math.Sqrt(zemu);
        for (i = 1; i < n; i++)
        {
            w[i] = 0.0;
        }

        //
        //  Diagonalize the Jacobi matrix.
        //
        IMTQLX.imtqlx(n, ref x, ref bj, ref w);

        for (i = 0; i < n; i++)
        {
            w[i] *= w[i];
        }
    }
}