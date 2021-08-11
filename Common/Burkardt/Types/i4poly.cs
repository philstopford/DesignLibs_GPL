using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void i4poly(int n, ref int[] a, int x0, int iopt, ref int val )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY performs operations on I4POLY's in power or factorial form.
        //
        //  Discussion:
        //
        //    The power sum form of a polynomial is
        //
        //      P(X) = A1 + A2*X + A3*X^2 + ... + (AN+1)*X^N
        //
        //    The Taylor expansion at C has the form
        //
        //      P(X) = A1 + A2*(X-C) + A3*(X-C)^2 + ... + (AN+1)*(X-C)^N
        //
        //    The factorial form of a polynomial is
        //
        //      P(X) = A1 + A2*X + A3*(X)*(X-1) + A4*(X)*(X-1)*(X-2)+...
        //        + (AN+1)*(X)*(X-1)*...*(X-N+1)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
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
        //    Input, int N, the number of coefficients in the polynomial
        //    (in other words, the polynomial degree + 1)
        //
        //    Input/output, int A[N], the coefficients of the polynomial.  Depending
        //    on the option chosen, these coefficients may be overwritten by those
        //    of a different form of the polynomial.
        //
        //    Input, int X0, for IOPT = -1, 0, or positive, the value of the
        //    argument at which the polynomial is to be evaluated, or the
        //    Taylor expansion is to be carried out.
        //
        //    Input, int IOPT, a flag describing which algorithm is to
        //    be carried out:
        //    -3: Reverse Stirling.  Input the coefficients of the polynomial in
        //    factorial form, output them in power sum form.
        //    -2: Stirling.  Input the coefficients in power sum form, output them
        //    in factorial form.
        //    -1: Evaluate a polynomial which has been input in factorial form.
        //    0:  Evaluate a polynomial input in power sum form.
        //    1 or more:  Given the coefficients of a polynomial in
        //    power sum form, compute the first IOPT coefficients of
        //    the polynomial in Taylor expansion form.
        //
        //    Output, int &VAL, for IOPT = -1 or 0, the value of the
        //    polynomial at the point X0.
        //
        {
            int eps;
            int i;
            int m;
            int n1;
            int w;
            int z;

            n1 = Math.Min(n, iopt);
            n1 = Math.Max(1, n1);

            if (iopt < -1)
            {
                n1 = n;
            }

            eps = Math.Max(-iopt, 0) % 2;

            w = -n * eps;

            if (-2 < iopt)
            {
                w = w + x0;
            }

            for (m = 1; m <= n1; m++)
            {
                val = 0;
                z = w;

                for (i = m; i <= n; i++)
                {
                    z = z + eps;
                    val = a[n + m - i - 1] + z * val;
                    if (iopt != 0 && iopt != -1)
                    {
                        a[n + m - i - 1] = val;
                    }
                }

                if (iopt < 0)
                {
                    w = w + 1;
                }
            }
        }

        public static int[] i4poly_add(int na, int[] a, int nb, int[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_ADD adds two I4POLY's.
        //
        //  Discussion:
        //
        //    The polynomials are in power sum form.
        //
        //    The power sum form is:
        //
        //      p(x) = a(0) + a(1)*x + ... + a(n-1)*x^(n-1) + a(n)*x^(n)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NA, the degree of polynomial A.
        //
        //    Input, int A[NA+1], the coefficients of the first
        //    polynomial factor.
        //
        //    Input, int NB, the degree of polynomial B.
        //
        //    Input, int B[NB+1], the coefficients of the
        //    second polynomial factor.
        //
        //    Output, int C[max(NA,NB)+1], the coefficients of A + B.
        //
        {
            int i;
            int[] c;

            c = new int[Math.Max(na, nb) + 1];

            if (nb == na)
            {
                for (i = 0; i <= na; i++)
                {
                    c[i] = a[i] + b[i];
                }
            }
            else if (nb < na)
            {
                for (i = 0; i <= nb; i++)
                {
                    c[i] = a[i] + b[i];
                }

                for (i = nb + 1; i <= na; i++)
                {
                    c[i] = a[i];
                }
            }
            else if (na < nb)
            {
                for (i = 0; i <= na; i++)
                {
                    c[i] = a[i] + b[i];
                }

                for (i = na + 1; i <= nb; i++)
                {
                    c[i] = b[i];
                }
            }

            return c;
        }

        public static void i4poly_cyclo(int n, int[] phi)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4POLY_CYCLO computes a cyclotomic I4POLY.
            //
            //  Discussion:
            //
            //    For 1 <= N, let
            //
            //      I = SQRT ( - 1 )
            //      L = EXP ( 2 * PI * I / N )
            //
            //    Then the N-th cyclotomic polynomial is defined by
            //
            //      PHI(N;X) = Product ( 1 <= K <= N and GCD(K,N) = 1 ) ( X - L^K )
            //
            //    We can use the Moebius MU function to write
            //
            //      PHI(N;X) = Product ( mod ( D, N ) = 0 ) ( X^D - 1 )^MU(N/D)
            //
            //    There is a sort of inversion formula:
            //
            //      X^N - 1 = Product ( mod ( D, N ) = 0 ) PHI(D;X)
            //
            //  Example:
            //
            //     N  PHI
            //
            //     0  1
            //     1  X - 1
            //     2  X + 1
            //     3  X^2 + X + 1
            //     4  X^2 + 1
            //     5  X^4 + X^3 + X^2 + X + 1
            //     6  X^2 - X + 1
            //     7  X^6 + X^5 + X^4 + X^3 + X^2 + X + 1
            //     8  X^4 + 1
            //     9  X^6 + X^3 + 1
            //    10  X^4 - X^3 + X^2 - X + 1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Raymond Seroul,
            //    Programming for Mathematicians,
            //    Springer Verlag, 2000, page 269.
            //
            //  Parameters:
            //
            //    Input, int N, the index of the cyclotomic polynomial desired.
            //
            //    Output, int PHI[N+1], the N-th cyclotomic polynomial.
            //
        {
            int POLY_MAX = 100;

            int d = 0;
            int[] den = new int[POLY_MAX + 1];
            int den_n;
            int[] factor;
            int i = 0;
            int j = 0;
            int mu = 0;
            int nq = 0;
            int nr = 0;
            int[] num = new int[POLY_MAX + 1];
            int num_n = 0;
            int[] rem;

            factor = new int[n + 1];
            rem = new int[n + 1];

            num[0] = 1;
            for (i = 1; i <= POLY_MAX; i++)
            {
                num[i] = 0;
            }

            num_n = 0;

            den[0] = 1;
            for (i = 1; i <= POLY_MAX; i++)
            {
                den[i] = 0;
            }

            den_n = 0;

            for (i = 0; i <= n; i++)
            {
                phi[i] = 0;
            }

            for (d = 1; d <= n; d++)
            {
                //
                //  For each divisor D of N, ...
                //
                if ((n % d) == 0)
                {
                    mu = i4_moebius(n / d);
                    //
                    //  ...multiply the numerator or denominator by (X^D-1).
                    //
                    factor[0] = -1;
                    for (j = 1; j <= d - 1; j++)
                    {
                        factor[j] = 0;
                    }

                    factor[d] = 1;

                    if (mu == +1)
                    {
                        if (POLY_MAX < num_n + d)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("I4POLY_CYCLO - Fatal error!");
                            Console.WriteLine("  Numerator polynomial degree too high.");
                            return;
                        }

                        i4poly_mul(num_n, num, d, factor, ref num);

                        num_n = num_n + d;
                    }
                    else if (mu == -1)
                    {
                        if (POLY_MAX < den_n + d)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("I4POLY_CYCLO - Fatal error!");
                            Console.WriteLine("  Denominator polynomial degree too high.");
                            return;
                        }

                        i4poly_mul(den_n, den, d, factor, ref den);

                        den_n = den_n + d;
                    }
                }
            }

            //
            //  PHI = NUM / DEN
            //
            i4poly_div(num_n, num, den_n, den, ref nq, ref phi, ref nr, ref rem);
        }

        public static int i4poly_degree(int na, int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4POLY_DEGREE returns the degree of an I4POLY.
            //
            //  Discussion:
            //
            //    The degree of a polynomial is the index of the highest power
            //    of X with a nonzero coefficient.
            //
            //    The degree of a constant polynomial is 0.  The degree of the
            //    zero polynomial is debatable, but this routine returns the
            //    degree as 0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NA, the dimension of A.
            //
            //    Input, int A[NA+1], the coefficients of the polynomials.
            //
            //    Output, int I4POLY_DEGREE, the degree of the polynomial.
            //
        {
            int degree;

            degree = na;

            while (0 < degree)
            {
                if (a[degree] != 0)
                {
                    return degree;
                }

                degree = degree - 1;
            }

            return degree;
        }

        public static int[] i4poly_dif(int na, int[] a, int d )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_DIF differentiates an I4POLY.
        //
        //  Discussion:
        //
        //    The polynomials are in power sum form.
        //
        //    The power sum form is:
        //
        //      p(x) = a(0) + a(1)*x + ... + a(n-1)*x^(n-1) + a(n)*x^(n)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NA, the degree of polynomial A.
        //
        //    Input, int A[NA+1], the coefficients of a polynomial.
        //
        //    Input, int D, the number of times the polynomial
        //    is to be differentiated.
        //
        //    Output, int I4POLY_DIF[NA-D+1], the coefficients of the
        //    differentiated polynomial.
        //
        {
            int[] b;
            int i;

            if (na < d)
            {
                b = new int[1];
                b[0] = 0;
                return b;
            }

            b = new int[na - d + 1];
            for (i = 0; i <= na - d; i++)
            {
                b[i] = a[i + d] * i4_fall(i + d, d);
            }

            return b;
        }

        public static void i4poly_div(int na, int[] a, int nb, int[] b, ref int nq, ref int[] q,
        ref int nr, ref int[] r )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_DIV computes the quotient and remainder of two I4POLY's.
        //
        //  Discussion:
        //
        //    Normally, the quotient and remainder would have rational coefficients.
        //    This routine assumes that the special case applies that the quotient
        //    and remainder are known beforehand to be integral.
        //
        //    The polynomials are assumed to be stored in power sum form.
        //
        //    The power sum form is:
        //
        //      p(x) = a(0) + a(1)*x + ... + a(n-1)*x^(n-1) + a(n)*x^(n)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NA, the degree of polynomial A.
        //
        //    Input, int A[NA+1], the coefficients of the polynomial to be divided.
        //
        //    Input, int NB, the degree of polynomial B.
        //
        //    Input, int B[NB+1], the coefficients of the divisor polynomial.
        //
        //    Output, int &NQ, the degree of polynomial Q.
        //    If the divisor polynomial is zero, NQ is returned as -1.
        //
        //    Output, int Q[NA-NB+1], contains the quotient of A/B.
        //    If A and B have full degree, Q should be dimensioned Q(0:NA-NB).
        //    In any case, Q(0:NA) should be enough.
        //
        //    Output, int &NR, the degree of polynomial R.
        //    If the divisor polynomial is zero, NR is returned as -1.
        //
        //    Output, int R[NB], contains the remainder of A/B.
        //    If B has full degree, R should be dimensioned R(0:NB-1).
        //    Otherwise, R will actually require less space.
        //
        {
            int[] a2;
            int i;
            int j;
            int na2;
            int nb2;

            na2 = i4poly_degree(na, a);

            nb2 = i4poly_degree(nb, b);

            if (b[nb2] == 0)
            {
                nq = -1;
                nr = -1;
                return;
            }

            a2 = new int[na + 1];

            for (i = 0; i <= na2; i++)
            {
                a2[i] = a[i];
            }

            nq = na2 - nb2;
            nr = nb2 - 1;

            for (i = nq; 0 <= i; i--)
            {
                q[i] = a2[i + nb2] / b[nb2];
                a2[i + nb2] = 0;
                for (j = 0; j < nb2; j++)
                {
                    a2[i + j] = a2[i + j] - q[i] * b[j];
                }
            }

            for (i = 0; i <= nr; i++)
            {
                r[i] = a2[i];
            }
        }

        public static void i4poly_mul(int na, int[] a, int nb, int[] b, ref int[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_MUL computes the product of two I4POLY's.
        //
        //  Discussion:
        //
        //    The polynomials are in power sum form.
        //
        //    The power sum form is:
        //
        //      p(x) = a(0) + a(1)*x + ... + a(n-1)*x^(n-1) + a(n)*x^(n)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NA, the degree of polynomial A.
        //
        //    Input, int A[NA+1], the coefficients of the first polynomial factor.
        //
        //    Input, int NB, the degree of polynomial B.
        //
        //    Input, int B[NB+1], the coefficients of the second polynomial factor.
        //
        //    Output, int C[NA+NB+1], the coefficients of A * B.
        //
        {
            int[] d;
            int i;
            int j;

            d = new int[na + nb + 1];

            for (i = 0; i <= na + nb; i++)
            {
                d[i] = 0;
            }

            for (i = 0; i <= na; i++)
            {
                for (j = 0; j <= nb; j++)
                {
                    d[i + j] = d[i + j] + a[i] * b[j];
                }
            }

            for (i = 0; i <= na + nb; i++)
            {
                c[i] = d[i];
            }
        }

        public static void i4poly_print(int n, int[] a, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_PRINT prints out an I4POLY.
        //
        //  Discussion:
        //
        //    The power sum form is:
        //
        //      p(x) = a(0) + a(1)*x + ... + a(n-1)*x^(n-1) + a(n)*x^(n)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the degree of polynomial A.
        //
        //    Input, int A[N+1], the polynomial coefficients.
        //    A(0) is the constant term and
        //    A(N) is the coefficient of X^N.
        //
        //    Input, string TITLE, a title.
        //
        {
            int i;
            int mag;
            int n2;
            char plus_minus;

            if (0 < title.Length)
            {
                Console.WriteLine("");
                Console.WriteLine(title + "");
            }

            n2 = i4poly_degree(n, a);

            if (a[n2] < 0)
            {
                plus_minus = '-';
            }
            else
            {
                plus_minus = ' ';
            }

            mag = Math.Abs(a[n2]);

            if (2 <= n2)
            {
                Console.WriteLine("p(x) = " + plus_minus + mag + " * x^" + n2 + "");
            }
            else if (n2 == 1)
            {
                Console.WriteLine("p(x) = " + plus_minus + mag + " * x" + "");
            }
            else if (n2 == 0)
            {
                Console.WriteLine("p(x) = " + plus_minus + mag + "");
            }

            for (i = n2 - 1; 0 <= i; i--)
            {
                if (a[i] < 0.0)
                {
                    plus_minus = '-';
                }
                else
                {
                    plus_minus = '+';
                }

                mag = Math.Abs(a[i]);

                if (mag != 0)
                {
                    if (2 <= i)
                    {
                        Console.WriteLine("       " + plus_minus + mag + " * x^" + i + "");
                    }
                    else if (i == 1)
                    {
                        Console.WriteLine("       " + plus_minus + mag + " * x" + "");
                    }
                    else if (i == 0)
                    {
                        Console.WriteLine("       " + plus_minus + mag + "");
                    }
                }
            }
        }

        public static int i4poly_to_i4(int n, int[] a, int x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_TO_I4 evaluates an I4POLY.
        //
        //  Discussion:
        //
        //    The power sum form is:
        //
        //      p(x) = a(0) + a(1)*x + ... + a(n-1)*x^(n-1) + a(n)*x^(n)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the degree of the polynomial.
        //
        //    Input, int A[N+1], the polynomial coefficients.
        //    A[0] is the constant term and
        //    A[N] is the coefficient of X^N.
        //
        //    Input, int X, the point at which the polynomial is to be evaluated.
        //
        //    Output, int I4POLY_TO_I4, the value of the polynomial.
        //
        {
            int i;
            int value;

            value = 0;

            for (i = n; 0 <= i; i--)
            {
                value = value * x + a[i];
            }

            return value;
        }
    }
}