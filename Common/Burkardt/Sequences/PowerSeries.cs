using System;

namespace Burkardt
{
    public static class PowerSeries
    {
        public static void power_series1(int n, double alpha, double[] a, ref double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POWER_SERIES1 computes a power series for a function G(Z) = (1+F(Z))^ALPHA.
            //
            //  Discussion:
            //
            //    The power series for F(Z) is given.
            //
            //    The form of the power series are:
            //
            //      F(Z) = A1*Z + A2*Z^2 + A3*Z^3 + ... + AN*Z^N
            //
            //      G(Z) = B1*Z + B2*Z^2 + B3*Z^3 + ... + BN*Z^N
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
            //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms for Computers and Calculators,
            //    Second Edition,
            //    Academic Press, 1978,
            //    ISBN: 0-12-519260-6,
            //    LC: QA164.N54.
            //
            //  Parameters:
            //
            //    Input, int N, the number of terms in the power series.
            //
            //    Input, double ALPHA, the exponent of 1+F(Z) in the definition of G(Z).
            //
            //    Input, double A[N], the power series coefficients for F(Z).
            //
            //    Output, double B[N], the power series coefficients for G(Z).
            //
        {
            int i;
            int j;
            double v;

            for (j = 1; j <= n; j++)
            {
                v = 0.0;
                for (i = 1; i <= j - 1; i++)
                {
                    v = v + b[i - 1] * a[j - i - 1] * (alpha * (j - i) - i);
                }

                b[j - 1] = alpha * a[j - 1] + v / ((double)j);
            }

        }

        public static void power_series2(int n, double[] a, ref double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POWER_SERIES2 computes the power series for a function G(Z) = EXP(F(Z)) - 1.
            //
            //  Discussion:
            //
            //    The power series for F(Z) is given.
            //
            //    The power series have the form:
            //
            //      F(Z) = A1*Z + A2*Z^2 + A3*Z^3 + ... + AN*Z^N
            //
            //      G(Z) = B1*Z + B2*Z^2 + B3*Z^3 + ... + BN*Z^N
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 June 2003
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms for Computers and Calculators,
            //    Second Edition,
            //    Academic Press, 1978,
            //    ISBN: 0-12-519260-6,
            //    LC: QA164.N54.
            //
            //  Parameters:
            //
            //    Input, int N, the number of terms in the power series.
            //
            //    Input, double A[N], the power series coefficients for F(Z).
            //
            //    Output, double B[N], the power series coefficients for G(Z).
            //
        {
            int i;
            int j;
            double v;

            for (j = 1; j <= n; j++)
            {
                v = 0.0;

                for (i = 1; i <= j - 1; i++)
                {
                    v = v + b[i - 1] * a[j - i - 1] * (double)(j - i);
                }

                b[j - 1] = a[j - 1] + v / (double)j;
            }
        }

        public static void power_series3(int n, double[] a, double[] b, ref double[] c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POWER_SERIES3 computes the power series for a function H(Z) = G(F(Z)).
            //
            //  Discussion:
            //
            //    The power series for F and G are given.
            //
            //    We assume that
            //
            //      F(Z) = A1*Z + A2*Z^2 + A3*Z^3 + ... + AN*Z^N
            //      G(Z) = B1*Z + B2*Z^2 + B3*Z^3 + ... + BN*Z^N
            //      H(Z) = C1*Z + C2*Z^2 + C3*Z^3 + ... + CN*Z^N
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 June 2003
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms for Computers and Calculators,
            //    Second Edition,
            //    Academic Press, 1978,
            //    ISBN: 0-12-519260-6,
            //    LC: QA164.N54.
            //
            //  Parameters:
            //
            //    Input, int N, the number of terms in the power series.
            //
            //    Input, double A[N], the power series for F
            //
            //    Input, double B[N], the power series for G.
            //
            //    Output, double C[N], the power series for H.
            //
        {
            int i;
            int iq;
            int j;
            int m;
            double r;
            double v;
            double[] work;

            work = new double[n];

            for (i = 0; i < n; i++)
            {
                work[i] = b[0] * a[i];
            }

            //
            //  Search for IQ, the index of the first nonzero entry in A.
            //
            iq = 0;

            for (i = 1; i <= n; i++)
            {
                if (a[i - 1] != 0.0)
                {
                    iq = i;
                    break;
                }
            }

            if (iq != 0)
            {
                m = 1;

                for (;;)
                {
                    m = m + 1;

                    if (n < m * iq)
                    {
                        break;
                    }

                    if (b[m - 1] == 0.0)
                    {
                        continue;
                    }

                    r = b[m - 1] * Math.Pow(a[iq - 1], m);
                    work[m * iq - 1] = work[m * iq - 1] + r;

                    for (j = 1; j <= n - m * iq; j++)
                    {
                        v = 0.0;
                        for (i = 1; i <= j - 1; i++)
                        {
                            v = v + c[i - 1] * a[j - i + iq - 1] * (double)(m * (j - i) - i);
                        }

                        c[j - 1] = ((double)m * a[j - 1] + v / (double)j) / a[iq - 1];

                    }

                    for (i = 1; i <= n - m * iq; i++)
                    {
                        work[i + m * iq - 1] = work[i + m * iq - 1] + c[i - 1] * r;
                    }
                }
            }

            for (i = 0; i < n; i++)
            {
                c[i] = work[i];
            }
        }

        public static void power_series4(int n, double[] a, double[] b, ref double[] c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POWER_SERIES4 computes the power series for a function H(Z) = G ( 1/F(Z) ).
            //
            //  Discussion:
            //
            //    POWER_SERIES4 is given the power series for the functions F and G.
            //
            //    We assume that
            //
            //      F(Z) = A1*Z + A2*Z^2 + A3*Z^3 + ... + AN*Z^N
            //      G(Z) = B1*Z + B2*Z^2 + B3*Z^3 + ... + BN*Z^N
            //      H(Z) = C1*Z + C2*Z^2 + C3*Z^3 + ... + CN*Z^N
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 June 2003
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms for Computers and Calculators,
            //    Second Edition,
            //    Academic Press, 1978,
            //    ISBN: 0-12-519260-6,
            //    LC: QA164.N54.
            //
            //  Parameters:
            //
            //    Input, int N, the number of terms in the power series.
            //
            //    Input, double A[N], the power series for F.  For this problem, A(1)
            //    may not be 0.0.
            //
            //    Input, double B(N), the power series for G.
            //
            //    Output, double C(N), the power series for H.
            //
        {
            int i;
            int l;
            int m;
            double s;
            double t;
            double[] work;

            if (a[0] == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("POWER_SERIES4 - Fatal error!");
                Console.WriteLine("  First entry of A is zero.");
                return;
            }

            work = new double[n];

            t = 1.0;

            for (i = 0; i < n; i++)
            {
                t = t / a[0];
                c[i] = b[i] * t;
                work[i] = a[i] * t;
            }

            for (m = 2; m <= n; m++)
            {
                s = -work[m - 1];
                for (i = m; i <= n; i++)
                {
                    for (l = i; l <= n; l++)
                    {
                        c[l - 1] = c[l - 1] + s * c[l - m];
                        work[l - 1] = work[l - 1] + s * work[l - m];
                    }
                }
            }
        }
    }
}