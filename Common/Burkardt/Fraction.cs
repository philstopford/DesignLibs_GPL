namespace Burkardt
{
    public static class Fraction
    {

        public static void cfrac_to_rat(int n, int[] a, ref int[] p, ref int[] q )

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
                if (i == 0)
                {
                    p[i] = a[i] * 1 + 0;
                    q[i] = a[i] * 0 + 1;
                }
                else if (i == 1)
                {
                    p[i] = a[i] * p[i - 1] + 1;
                    q[i] = a[i] * q[i - 1] + 0;
                }
                else
                {
                    p[i] = a[i] * p[i - 1] + p[i - 2];
                    q[i] = a[i] * q[i - 1] + q[i - 2];
                }

            }
        }

        public static void cfrac_to_rfrac(int m, double[] g, double[] h, ref double[] p, ref double[] q )

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

            if (m == 1)
            {
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
        
        public static void jfrac_to_rfrac ( int m, double[] r, double[] s, ref double[] p, ref double[] q )

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

            if (1 < m)
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
    }
}