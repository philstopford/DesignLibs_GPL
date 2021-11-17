using System;

namespace Burkardt.Function;

public static class Fraction
{

    public static void cfrac_to_rat(int n, int[] a, ref int[] p, ref int[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CFRAC_TO_RAT converts a monic continued fraction to an ordinary fraction.
        //
        //  Discussion:
        //
        //    The routine is given the monic or "simple" continued fraction with
        //    integer coefficients:
        //
        //      A(1) + 1 / ( A(2) + 1 / ( A(3) ... + 1 / A(N) ) )
        //
        //    and returns the N successive approximants P(I)/Q(I)
        //    to the value of the rational number represented by the continued
        //    fraction, with the value exactly equal to the final ratio P(N)/Q(N).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 June 2004
        //
        //  Author:
        //
        //    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
        //    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher, 
        //    Christoph Witzgall.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
        //    John Rice, Henry Thatcher, Christoph Witzgall,
        //    Computer Approximations,
        //    Wiley, 1968.
        //
        //  Parameters:
        //
        //    Input, int N, the number of continued fraction coefficients.
        //
        //    Input, int A[N], the continued fraction coefficients.
        //
        //    Output, int P[N], Q[N], the N successive approximations
        //    to the value of the continued fraction.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            switch (i)
            {
                case 0:
                    p[i] = a[i] * 1 + 0;
                    q[i] = a[i] * 0 + 1;
                    break;
                case 1:
                    p[i] = a[i] * p[i - 1] + 1;
                    q[i] = a[i] * q[i - 1] + 0;
                    break;
                default:
                    p[i] = a[i] * p[i - 1] + p[i - 2];
                    q[i] = a[i] * q[i - 1] + q[i - 2];
                    break;
            }
        }
    }

    public static void cfrac_to_rfrac(int m, double[] g, double[] h, ref double[] p, ref double[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CFRAC_TO_RFRAC converts a polynomial fraction from continued to rational form.
        //
        //  Discussion:
        //
        //    The routine accepts a continued polynomial fraction:
        //
        //      G(1)     / ( H(1) +
        //      G(2) * X / ( H(2) +
        //      G(3) * X / ( H(3) + ...
        //      G(M) * X / ( H(M) )...) ) )
        //
        //    and returns the equivalent rational polynomial fraction:
        //
        //      P(1) + P(2) * X + ... + P(L1) * X^(L1)
        //      -------------------------------------------------------
        //      Q(1) + Q(2) * X + ... + Q(L2) * X^(L2-1)
        //
        //    where
        //
        //      L1 = (M+1)/2
        //      L2 = (M+2)/2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 June 2004
        //
        //  Author:
        //
        //    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
        //    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher, 
        //    Christoph Witzgall.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
        //    John Rice, Henry Thatcher, Christoph Witzgall,
        //    Computer Approximations,
        //    Wiley, 1968.
        //
        //  Parameters:
        //
        //    Input, int M, the number of continued fraction polynomial coefficients.
        //
        //    Input, double G[M], H[M], the continued polynomial fraction coefficients.
        //
        //    Output, double P[(M+1)/2], Q[(M+2)/2], the rational polynomial fraction
        //    coefficients.
        //
    {
        double[] a;
        int i;
        int j;

        switch (m)
        {
            case 1:
                p[0] = g[0];
                q[0] = h[0];
                return;
        }

        a = new double[m * ((m + 2) / 2)];

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < (m + 2) / 2; j++)
            {
                a[i + j * m] = 0.0;
            }
        }

        //
        //  Solve for P's
        //
        a[0 + 0 * m] = g[0];
        a[1 + 0 * m] = g[0] * h[1];

        for (i = 2; i < m; i++)
        {
            a[i + 0 * m] = h[i] * a[i - 1 + 0 * m];
            for (j = 1; j < (i + 2) / 2; j++)
            {
                a[i + j * m] = h[i] * a[i - 1 + j * m] + g[i] * a[i - 2 + (j - 1) * m];
            }
        }

        for (j = 0; j < (m + 1) / 2; j++)
        {
            p[j] = a[m - 1 + j * m];
        }

        //
        //  Solve for Q's.
        //
        a[0 + 0 * m] = h[0];
        a[1 + 0 * m] = h[0] * h[1];
        a[1 + 1 * m] = g[1];

        for (i = 2; i < m; i++)
        {
            a[i + 0 * m] = h[i] * a[i - 1 + 0 * m];
            for (j = 1; j < (i + 3) / 2; j++)
            {
                a[i + j * m] = h[i] * a[i - 1 + j * m] + g[i] * a[i - 2 + (j - 1) * m];
            }
        }

        for (j = 0; j < (m + 2) / 2; j++)
        {
            q[j] = a[m - 1 + j * m];
        }
    }

    public static void jfrac_to_rfrac(int m, double[] r, double[] s, ref double[] p, ref double[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JFRAC_TO_RFRAC converts a J-fraction into a rational polynomial fraction.
        //
        //  Discussion:
        //
        //    The routine accepts a J-fraction:
        //
        //        R(1) / ( X + S(1)
        //      + R(2) / ( X + S(2)
        //      + R(3) / ...
        //      + R(M) / ( X + S(M) )... ))
        //
        //    and returns the equivalent rational polynomial fraction:
        //
        //      P(1) + P(2) * X + ... + P(M) * X^(M-1)
        //      -------------------------------------------------------
        //      Q(1) + Q(2) * X + ... + Q(M) * X^(M-1) + Q(M+1) * X^M
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
        //    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher,
        //    Christop Witzgall.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
        //    John Rice, Henry Thatcher, Christoph Witzgall,
        //    Computer Approximations,
        //    Wiley, 1968.
        //
        //  Parameters:
        //
        //    Input, int M, defines the number of P, R, and S
        //    coefficients, and is one less than the number of Q
        //    coefficients.
        //
        //    Input, double R[M], S[M], the coefficients defining the J-fraction.
        //
        //    Output, double P[M], Q[M+1], the coefficients defining the rational
        //    polynomial fraction.  The algorithm used normalizes the coefficients
        //    so that Q[M+1] = 1.0.
        //
    {
        double[] a;
        double[] b;
        int i;
        int k;

        a = new double[m * m];
        b = new double[m * m];

        a[0 + 0 * m] = r[0];
        b[0 + 0 * m] = s[0];

        switch (m)
        {
            case > 1:
            {
                for (k = 1; k < m; k++)
                {
                    a[k + k * m] = r[0];
                    b[k + k * m] = b[k - 1 + (k - 1) * m] + s[k];
                }

                a[0 + 1 * m] = r[0] * s[1];
                b[0 + 1 * m] = r[1] + s[0] * s[1];

                for (k = 2; k < m; k++)
                {
                    a[0 + k * m] = s[k] * a[0 + (k - 1) * m] + r[k] * a[0 + (k - 2) * m];
                    a[k - 1 + k * m] = a[k - 2 + (k - 1) * m] + s[k] * r[0];
                    b[0 + k * m] = s[k] * b[0 + (k - 1) * m] + r[k] * b[0 + (k - 2) * m];
                    b[k - 1 + k * m] = b[k - 2 + (k - 1) * m] + s[k] * b[k - 1 + (k - 1) * m] + r[k];
                }

                for (k = 2; k < m; k++)
                {
                    for (i = 1; i < k - 1; i++)
                    {
                        a[i + k * m] = a[i - 1 + (k - 1) * m] + s[k] * a[i + (k - 1) * m]
                                                              + r[k] * a[i + (k - 2) * m];
                        b[i + k * m] = b[i - 1 + (k - 1) * m] + s[k] * b[i + (k - 1) * m]
                                                              + r[k] * b[i + (k - 2) * m];
                    }
                }

                break;
            }
        }

        for (i = 0; i < m; i++)
        {
            p[i] = a[i + (m - 1) * m];
        }

        for (i = 0; i < m; i++)
        {
            q[i] = b[i + (m - 1) * m];
        }

        q[m] = 1.0;
    }

    public static void rfrac_to_cfrac(int m, double[] p, double[] q, ref double[] t, ref bool error)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RFRAC_TO_CFRAC converts a rational polynomial fraction to a continued fraction.
        //
        //  Discussion:
        //
        //    That is, it accepts
        //
        //      P(1) + P(2) * X + ... + P(M) * X^(M-1)
        //      -------------------------------------------------------
        //      Q(1) + Q(2) * X + ... + Q(M) * X^(M-1) + Q(M+1) * X^M
        //
        //    and returns the equivalent continued fraction:
        //
        //      1 / (T(1) + X/(T(2) + X/(...T(2*M-1)+X/(T(2*M) ... )))
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 June 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
        //    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher,
        //    Christoph Witzgall.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
        //    John Rice, Henry Thatcher, Christoph Witzgall,
        //    Computer Approximations,
        //    Wiley, 1968.
        //
        //  Parameters:
        //
        //    Input, int M, defines the number of P coefficients,
        //    and is one less than the number of Q coefficients, and one
        //    half the number of T coefficients.
        //
        //    Input, double P[M], Q[M+1], the coefficients defining the rational
        //    polynomial fraction.
        //
        //    Output, double T[2*M], the coefficients defining the continued fraction.
        //
        //    Output, bool &ERROR, is TRUE if an error occurred.
        //
    {
        double[] a;
        int i;
        int k;
        double ta;

        a = new double[(m + 1) * (2 * m + 1)];

        error = false;

        for (i = 0; i <= m; i++)
        {
            a[i + 0 * (m + 1)] = q[i];
        }

        for (i = 0; i < m; i++)
        {
            a[i + 1 * (m + 1)] = p[i];
        }

        t[0] = a[0 + 0 * (m + 1)] / a[0 + 1 * (m + 1)];
        ta = a[m + 0 * (m + 1)];

        for (i = 1; i <= m; i++)
        {
            a[m - i + 2 * i * (m + 1)] = ta;
        }

        for (k = 1; k <= 2 * m - 2; k++)
        {
            for (i = 1; i <= (2 * m - k) / 2; i++)
            {
                a[i - 1 + (k + 1) * (m + 1)] = a[i + (k - 1) * (m + 1)] - t[k - 1] * a[i + k * (m + 1)];
            }

            switch (a[0 + (k + 1) * (m + 1)])
            {
                case 0.0:
                    error = true;
                    Console.WriteLine("");
                    Console.WriteLine("RFRAC_TO_CFRAC - Fatal error!");
                    Console.WriteLine("  A[0,K+1] is zero for K = " + k + "");
                    return;
                default:
                    t[k] = a[0 + k * (m + 1)] / a[0 + (k + 1) * (m + 1)];
                    break;
            }
        }

        t[2 * m - 1] = a[0 + (2 * m - 1) * (m + 1)] / a[0 + 2 * m * (m + 1)];
    }

    public static void rfrac_to_jfrac(int m, double[] p, double[] q, ref double[] r, ref double[] s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RFRAC_TO_JFRAC converts a rational polynomial fraction to a J fraction.
        //
        //  Discussion:
        //
        //    The routine accepts
        //
        //    P(1) + P(2) * X + ... + P(M) * X^(M-1)
        //    -------------------------------------------------------
        //    Q(1) + Q(2) * X + ... + Q(M) * X^(M-1) + Q(M+1) * X^M
        //
        //    and returns the equivalent J-fraction:
        //
        //    R(1) / ( X + S(1) + 
        //    R(2) / ( X + S(2) + 
        //    R(3) / ...        +
        //    R(M) / ( X + S(M) )... ))
        //
        //    Thanks to Henry Amuasi for noticing and correcting an error in a
        //    previous formulation of this routine, 02 October 2010.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 October 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
        //    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher,
        //    Christoph Witzgall.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
        //    John Rice, Henry Thatcher, Christoph Witzgall,
        //    Computer Approximations,
        //    Wiley, 1968.
        //
        //  Parameters:
        //
        //    Input, int M, defines the number of P, R, and S coefficients,
        //    and is one less than the number of Q coefficients.
        //    1 <= M.
        //
        //    Input, double P[M], Q[M+1], the coefficients defining the rational
        //    polynomial fraction.
        //
        //    Output, double R[M], S[M], the coefficients defining the
        //    J-fraction.
        //
    {
        double[] a;
        int i;
        int k;

        switch (m)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("RFRAC_TO_JFRAC - Fatal error!");
                Console.WriteLine("  Input M < 1.");
                return;
        }

        a = new double[(m + 1) * (m + 1)];

        for (i = 0; i <= m; i++)
        {
            a[i + 0 * (m + 1)] = q[i];
        }

        for (i = 0; i < m; i++)
        {
            a[i + 1 * (m + 1)] = p[i];
        }

        switch (m)
        {
            case > 1:
            {
                r[0] = a[m - 1 + 1 * (m + 1)] / a[m + 0 * (m + 1)];
                s[0] = (r[0] * a[m - 1 + 0 * (m + 1)] - a[m - 2 + 1 * (m + 1)]) / a[m - 1 + 1 * (m + 1)];

                for (k = 0; k < m - 2; k++)
                {
                    a[0 + (k + 2) * (m + 1)] = r[k] * a[0 + k * (m + 1)] - s[k] * a[0 + (k + 1) * (m + 1)];

                    for (i = 1; i < m - k; i++)
                    {
                        a[i + (k + 2) * (m + 1)] = r[k] * a[i + k * (m + 1)]
                                                   - a[i - 1 + (k + 1) * (m + 1)] - s[k] * a[i + (k + 1) * (m + 1)];
                    }

                    switch (a[m - 2 - k + (k + 2) * (m + 1)])
                    {
                        case 0.0:
                            Console.WriteLine("");
                            Console.WriteLine("RFRAC_TO_JFRAC - Fatal error!");
                            Console.WriteLine("  A[M-2-K,K+2] = 0 for K = " + k + "");
                            return;
                        default:
                            r[k + 1] = a[m - k - 2 + (k + 2) * (m + 1)] / a[m - 2 - k + 1 + (k + 1) * (m + 1)];
                            s[k + 1] = (r[k + 1] * a[m - 2 - k + (k + 1) * (m + 1)]
                                        - a[m - 2 - k - 1 + (k + 2) * (m + 1)]) / a[m - 2 - k + (k + 2) * (m + 1)];
                            break;
                    }
                }

                a[0 + m * (m + 1)] = r[m - 2] * a[0 + (m - 2) * (m + 1)] - s[m - 2] * a[0 + (m - 1) * (m + 1)];
                break;
            }
        }

        r[m - 1] = a[0 + m * (m + 1)] / a[1 + (m - 1) * (m + 1)];
        s[m - 1] = a[0 + (m - 1) * (m + 1)] / a[1 + (m - 1) * (m + 1)];

    }
}