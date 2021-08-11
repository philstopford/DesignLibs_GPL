using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8poly_value(int n, double[] a, double x)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_VALUE evaluates a double precision polynomial.
            //
            //  Discussion:
            //
            //    For sanity's sake, the value of N indicates the NUMBER of 
            //    coefficients, or more precisely, the ORDER of the polynomial,
            //    rather than the DEGREE of the polynomial.  The two quantities
            //    differ by 1, but cause a great deal of confusion.
            //
            //    Given N and A, the form of the polynomial is:
            //
            //      p(x) = a[0] + a[1] * x + ... + a[n-2] * x^(n-2) + a[n-1] * x^(n-1)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the polynomial.
            //
            //    Input, double A[N], the coefficients of the polynomial.
            //    A[0] is the constant term.
            //
            //    Input, double X, the point at which the polynomial is to be evaluated.
            //
            //    Output, double R8POLY_VALUE, the value of the polynomial at X.
            //
        {
            double value = 0.0;

            for (int i = n - 1; 0 <= i; i--)
            {
                value = value * x + a[i];
            }

            return value;
        }

        public static double r8poly_value_horner(int m, double[] c, double x)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_VALUE_HORNER evaluates a polynomial using Horner's method.
            //
            //  Discussion:
            //
            //    The polynomial 
            //
            //      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
            //
            //    is to be evaluated at the value X.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the degree of the polynomial.
            //
            //    Input, double C[M+1], the coefficients of the polynomial.
            //    A[0] is the constant term.
            //
            //    Input, double X, the point at which the polynomial is to be evaluated.
            //
            //    Output, double R8POLY_VALUE_HORNER, the value of the polynomial at X.
            //
        {
            double value = c[m];

            for (int i = m - 1; 0 <= i; i--)
            {
                value = value * x + c[i];
            }

            return value;
        }

        public static void r8poly_ant_cof(int n, double[] poly_cof, ref double[] poly_cof2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_ANT_COF integrates an R8POLY in standard form.
            //
            //  Discussion:
            //
            //    The antiderivative of a polynomial P(X) is any polynomial Q(X)
            //    with the property that d/dX Q(X) = P(X).
            //
            //    This routine chooses the antiderivative whose constant term is zero.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the polynomial.
            //
            //    Input, double POLY_COF[N], the polynomial coefficients.
            //    POLY_COF[0] is the constant term, and POLY_COF[N-1] is the
            //    coefficient of X**(N-1).
            //
            //    Output, double POLY_COF2[N+1], the coefficients of the antiderivative
            //    polynomial, in standard form.  The constant term is set to zero.
            //
        {
            int i;
            //
            //  Set the constant term.
            //
            poly_cof2[0] = 0.0;
            //
            //  Integrate the polynomial.
            //
            for (i = 1; i <= n; i++)
            {
                poly_cof2[i] = poly_cof[i - 1] / (double)i;
            }
        }

        public static double r8poly_ant_val(int n, double[] poly_cof, double xval)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_ANT_VAL evaluates the antiderivative of an R8POLY in standard form.
            //
            //  Discussion:
            //
            //    The constant term of the antiderivative is taken to be zero.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the polynomial.
            //
            //    Input, double POLY_COF[N], the polynomial coefficients.  POLY_COF[0]
            //    is the constant term, and POLY_COF[N-1] is the coefficient of X**(N-1).
            //
            //    Input, double XVAL, the point where the antiderivative is to be
            //    evaluated.
            //
            //    Output, double R8POLY_ANT_VAL, the value of the antiderivative of the polynomial
            //    at XVAL.
            //
        {
            int i;
            double value;

            value = 0.0;

            for (i = n - 1; 0 <= i; i--)
            {
                value = (value + poly_cof[i] / (double)(i + 1)) * xval;
            }

            return value;
        }

        public static void r8poly_basis(int ntab, double[] xtab, ref double[] poly_cof)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_BASIS computes all Lagrange basis polynomials in standard form.
            //
            //  Discussion:
            //
            //    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
            //    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
            //    XTAB(J) for J not equal to I, and 1 when J is equal to I.
            //
            //    The Lagrange basis polynomials have the property that the interpolating
            //    polynomial through a set of NTAB data points (XTAB,YTAB) may be
            //    represented as
            //
            //      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
            //
            //    Higher order interpolation at selected points may be accomplished
            //    using repeated X values, and scaled derivative values.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NTAB, the number of data points XTAB.
            //
            //    Input, double XTAB[NTAB], the X values upon which the Lagrange basis
            //    polynomial is to be based.
            //
            //    Output, double *POLY_COF, points to NTAB * NTAB values.  The polynomial
            //    coefficients for the I-th Lagrange basis polynomial are stored in
            //    (logical) row I.  POLY_COF[0,*] is the constant term, and POLY_COF[NTAB-1,*] is
            //    the coefficient of X**(NTAB-1).
            //
        {
            int i;
            int j;
            double[] pointer1;
            double[] pointer2;

            pointer1 = poly_cof;

            for (i = 0; i < ntab; i++)
            {
                pointer2 = pointer1;

                for (j = 0; j < ntab; j++)
                {
                    if (j == i)
                    {
                        pointer1[i] = 1.0;
                    }
                    else
                    {
                        pointer1[i] = 0.0;
                    }
                }

                //
                //  Compute the divided difference table for the IVAL-th Lagrange basis
                //  polynomial.
                //
                Data.data_to_dif(ntab, xtab, pointer2, ref pointer2);
                //
                //  Convert the divided difference table coefficients to standard polynomial
                //  coefficients.
                //
                Dif.dif_to_r8poly(ntab, xtab, pointer2, ref pointer2);
            }
        }

        public static void r8poly_basis_1(int ival, int ntab, double[] xtab, ref double[] poly_cof)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_BASIS_1 computes the I-th Lagrange basis polynomial in standard form.
            //
            //  Discussion:
            //
            //    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
            //    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
            //    XTAB(J) for J not equal to I, and 1 when J is equal to I.
            //
            //    The Lagrange basis polynomials have the property that the interpolating
            //    polynomial through a set of NTAB data points (XTAB,YTAB) may be
            //    represented as
            //
            //      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
            //
            //    Higher order interpolation at selected points may be accomplished
            //    using repeated X values, and scaled derivative values.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int IVAL, the index of the desired Lagrange basis polynomial.
            //    IVAL should be between 1 and NTAB.
            //
            //    Input, int NTAB, the number of data points XTAB.
            //
            //    Input, double XTAB[NTAB], the X values upon which the Lagrange basis
            //    polynomial is to be based.
            //
            //    Output, double POLY_COF[NTAB], the polynomial coefficients for the
            //    IVAL-th Lagrange basis polynomial.
            //
        {
            int i;
            //
            //  Check IVAL.
            //
            if (ival < 1 || ntab < ival)
            {
                Console.WriteLine("");
                Console.WriteLine("R8POLY_BASIS_1 - Fatal error!");
                Console.WriteLine("  IVAL must be between 1 and " + ntab + ".");
                Console.WriteLine("  but your value is " + ival + ".");
                return;
            }

            //
            //  Initialize POLY_COF to the IVAL-th column of the identity matrix.
            //
            for (i = 0; i <= ntab - 1; i++)
            {
                poly_cof[i] = 0.0;
            }

            poly_cof[ival - 1] = 1.0;
            //
            //  Compute the divided difference table for the IVAL-th Lagrange basis
            //  polynomial.
            //
            Data.data_to_dif(ntab, xtab, poly_cof, ref poly_cof);
            //
            //  Convert the divided difference table coefficients to standard polynomial
            //  coefficients.
            //
            Dif.dif_to_r8poly(ntab, xtab, poly_cof, ref poly_cof);
        }

        public static void r8poly_der_cof(int n, double[] poly_cof, ref double[] poly_cof2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_DER_COF computes the coefficients of the derivative of a real polynomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the polynomial.
            //
            //    Input, double POLY_COF[N], the coefficients of the polynomial to
            //    be differentiated.  POLY_COF[0] is the constant term, and
            //    POLY_COF[N-1] is the coefficient of X**(N-1).
            //
            //    Output, double POLY_COF2[N-1], the coefficients of the derivative of
            //    the polynomial.
            //
        {
            // Safety due to some oddball tests.
            n = Math.Min(poly_cof.Length - 1, n);

            int i;

            for (i = 0; i < n; i++)
            {
                poly_cof2[i] = (double)(i + 1) * poly_cof[i + 1];
            }

            return;
        }

        public static double r8poly_der_val(int n, double[] poly_cof, double xval)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_DER_VAL evaluates the derivative of a real polynomial in standard form.
            //
            //  Discussion:
            //
            //    A polynomial in standard form, with coefficients POLY_COF(*),
            //    may be written:
            //
            //    P(X) = POLY_COF[0]
            //         + POLY_COF[1] * X
            //         ...
            //         + POLY_COF[N-1] * X**(N-1)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the polynomial.
            //
            //    Input, double POLY_COF[N], the polynomial coefficients.  POLY_COF[0]
            //    is the constant term, and POLY_COF[N-1] is the coefficient of
            //    X**(N-1).
            //
            //    Input, double XVAL, a value where the derivative of the polynomial
            //    is to be evaluated.
            //
            //    Output, double R8POLY_DER_VAL, the value of the derivative of the polynomial
            //    at XVAL.
            //
        {
            int i;
            double value;

            value = (double)(n - 1) * poly_cof[n - 1];

            for (i = n - 2; 1 <= i; i--)
            {
                value = value * xval + (double)i * poly_cof[i];
            }

            return value;
        }

        public static int r8poly_order(int na, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_ORDER returns the order of a polynomial.
            //
            //  Discussion:
            //
            //    The order of a polynomial is the degree plus 1.
            //
            //    The order of a constant polynomial is 1.  The order of the
            //    zero polynomial is debatable, but this routine returns the
            //    order as 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NA, the order of the polynomial (ignoring zero coefficients).
            //
            //    Input, double A[NA], the coefficients of the polynomial.
            //
            //    Output, int R8POLY_ORDER, the degree of the polynomial.
            //
        {
            int value;

            value = na;

            while (1 < value)
            {
                if (a[value - 1] != 0.0)
                {
                    return value;
                }

                value = value - 1;
            }

            return value;
        }

        public static void r8poly_print(int n, double[] a, string title, int aIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_PRINT prints out a polynomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of A.
            //
            //    Input, double A[N+1], the polynomial coefficients.
            //    A(0) is the constant term and
            //    A(N) is the coefficient of X^N.
            //
            //    Input, string TITLE, a title.
            //
        {
            // Safety measure due to some odd test configs.
            n = Math.Min(a.Length - 1, n);

            int i;
            double mag;
            char plus_minus;

            if (0 < title.Length)
            {
                Console.WriteLine("");
                Console.WriteLine(title + "");
            }

            Console.WriteLine("");

            if (n < 0)
            {
                Console.WriteLine("  p(x) = 0");
                return;
            }

            if (a[n] < 0.0)
            {
                plus_minus = '-';
            }
            else
            {
                plus_minus = ' ';
            }

            mag = Math.Abs(a[n]);

            if (2 <= n)
            {
                Console.WriteLine("  p(x) = " + plus_minus
                                              + mag.ToString().PadLeft(14) + " * x ^ " + n + "");
            }
            else if (n == 1)
            {
                Console.WriteLine("  p(x) = " + plus_minus
                                              + mag.ToString().PadLeft(14) + " * x");
            }
            else if (n == 0)
            {
                Console.WriteLine("  p(x) = " + plus_minus
                                              + mag.ToString().PadLeft(14) + "");
            }

            for (i = n - 1; 0 <= i; i--)
            {
                if (a[aIndex + i] < 0.0)
                {
                    plus_minus = '-';
                }
                else
                {
                    plus_minus = '+';
                }

                mag = Math.Abs(a[aIndex + i]);

                if (mag != 0.0)
                {
                    if (2 <= i)
                    {
                        Console.WriteLine("         " + plus_minus
                                                      + mag.ToString().PadLeft(14) + " * x ^ " + i + "");
                    }
                    else if (i == 1)
                    {
                        Console.WriteLine("         " + plus_minus
                                                      + mag.ToString().PadLeft(14) + " * x");
                    }
                    else if (i == 0)
                    {
                        Console.WriteLine("         " + plus_minus
                                                      + mag.ToString().PadLeft(14) + "");
                    }
                }
            }
        }

        public static double r8poly_pval ( int n, double[] a, double x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_PVAL evaluates a real polynomial in power sum form.
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
        //    Input, int N, the dimension of A.
        //
        //    Input, double A[N+1], the coefficients of the polynomial.
        //    A(0) is the constant term.
        //
        //    Input, double X, the point at which the polynomial is to be evaluated.
        //
        //    Output, double R8POLY_VAL, the value of the polynomial at X.
        //
        {
            int i;
            double value;

            value = 0.0;
            for ( i = n; 0 <= i; i-- )
            {
                value = value * x + a[i];
            }

            return value;
        }
        public static void r8poly_shift(double scale, double shift, int n, ref double[] poly_cof)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_SHIFT adjusts the coefficients of a polynomial for a new argument.
            //
            //  Discussion:
            //
            //    Assuming P(X) is a polynomial in the argument X, of the form:
            //
            //      P(X) =
            //          C(N-1) * X^(N-1)
            //        + ...
            //        + C(1) * X
            //        + C(0),
            //
            //    and that Z is related to X by the formula:
            //
            //      Z = SCALE * X + SHIFT
            //
            //    then this routine computes coefficients C for the polynomial Q(Z):
            //
            //      Q(Z) =
            //          C(N-1) * Z^(N-1)
            //        + ...
            //        + C(1) * Z
            //        + C(0)
            //
            //    so that:
            //
            //      Q(Z(X)) = P(X)
            //
            //  Example:
            //
            //    P(X) = 2 * X^2 - X + 6
            //
            //    Z = 2.0 * X + 3.0
            //
            //    Q(Z) = 0.5 *         Z^2 -  3.5 * Z + 12
            //
            //    Q(Z(X)) = 0.5 * ( 4.0 * X^2 + 12.0 * X +  9 )
            //            - 3.5 * (               2.0 * X +  3 )
            //                                            + 12
            //
            //            = 2.0         * X^2 -  1.0 * X +  6
            //
            //            = P(X)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double SHIFT, SCALE, the shift and scale applied to X,
            //    so that Z = SCALE * X + SHIFT.
            //
            //    Input, int N, the order of the polynomial.
            //
            //    Input/output, double POLY_COF[N].
            //    On input, the coefficient array in terms of the X variable.
            //    On output, the coefficient array in terms of the Z variable.
            //
        {
            int i;
            int j;

            for (i = 1; i <= n; i++)
            {
                for (j = i + 1; j <= n; j++)
                {
                    poly_cof[j - 1] = poly_cof[j - 1] / scale;
                }
            }

            for (i = 1; i <= n; i++)
            {
                for (j = n - 1; i <= j; j--)
                {
                    poly_cof[j - 1] = poly_cof[j - 1] - shift * poly_cof[j];
                }
            }

            return;
        }

        public static void r8poly_t2p ( int n, ref double[] a, double x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_T2P converts a real polynomial from Taylor form to power sum form
        //
        //  Discussion:
        //
        //    The Taylor form is
        //
        //      p(x) =   a(1)
        //             + a(2) * (x-x0)
        //             + a(3) * (x-x0)^2
        //             ...
        //             + a(n) * (x-x0)^(n-1)
        //
        //    The power sum form is
        //
        //      p(x) = a(1) + a(2)*x + a(3)*x^2 + ... + a(n)*x^(n-1)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of A.
        //
        //    Input/output, double A[N].  On input, the coefficients in Taylor form,
        //    and on output, the coefficients in power sum form.
        //
        //    Input, double X, the point at which the Taylor form polynomial is based.
        //
        {
            int i;
            int j;

            for ( i = n; 1 <= i; i-- )
            {
                for ( j = i; j <= n-1; j++ )
                {
                    a[j-1] = a[j-1] - a[j] * x;
                }
            }
        }
        public static double r8poly_val_horner(int n, double[] poly_cof, double xval, int polyCofIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_VAL_HORNER evaluates a real polynomial in standard form.
            //
            //  Discussion:
            //
            //    A polynomial in standard form, with coefficients POLY_COF(*),
            //    may be written:
            //
            //    P(X) = POLY_COF[0]
            //         + POLY_COF[1] * X
            //         ...
            //         + POLY_COF[N-1] * X^(N-1)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the polynomial.
            //
            //    Input, double POLY_COF[N], the polynomial coefficients.  POLY_COF[0]
            //    is the constant term, and POLY_COF[N-1] is the coefficient of
            //    X^(N-1).
            //
            //    Input, double XVAL, a value where the polynomial is to be evaluated.
            //
            //    Output, double R8POLY_VAL_HORNER, the value of the polynomial at XVAL.
            //
        {
            int i;
            double value;

            value = poly_cof[polyCofIndex + (n - 1)];

            for (i = n - 2; 0 <= i; i--)
            {
                value = value * xval + poly_cof[polyCofIndex + i];
            }

            return value;
        }

        public static void r8poly(int n, ref double[] a, double x0, int iopt, ref double val)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY performs operations on R8POLY's in power or factorial form.
            //
            //  Discussion:
            //
            //    The power sum form of a polynomial is
            //
            //      P(X) = A1 + A2*X + A3*X^2 + ... + (AN+1)*X^N
            //
            //    The Taylor expansion at C has the form
            //
            //      P(X) = A1 + A2*(X-C) + A3*(X-C)^2+... + (AN+1)*(X-C)^N
            //
            //    The factorial form of a polynomial is
            //
            //      P(X) = A1 + A2*X + A3*(X)*(X-1) + A4*(X)*(X-1)*(X-2) + ...
            //        + (AN+1)*(X)*(X-1)*...*(X-N+1)
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
            //    Input, int N, the number of coefficients in the polynomial
            //    (in other words, the polynomial degree + 1)
            //
            //    Input/output, double A[N], the coefficients of the polynomial.  Depending
            //    on the option chosen, these coefficients may be overwritten by those
            //    of a different form of the polynomial.
            //
            //    Input, double X0, for IOPT = -1, 0, or positive, the value of the
            //    argument at which the polynomial is to be evaluated, or the
            //    Taylor expansion is to be carried out.
            //
            //    Input, int IOPT, a flag describing which algorithm is to
            //    be carried out:
            //    -3: Reverse Stirling.  Input the coefficients of the polynomial in
            //    factorial form, output them in power sum form.
            //    -2: Stirling.  Input the coefficients in power sum
            //    form, output them in factorial form.
            //    -1: Evaluate a polynomial which has been input
            //    in factorial form.
            //    0:  Evaluate a polynomial input in power sum form.
            //    1 or more:  Given the coefficients of a polynomial in
            //    power sum form, compute the first IOPT coefficients of
            //    the polynomial in Taylor expansion form.
            //
            //    Output, double &VAL, for IOPT = -1 or 0, the value of the
            //    polynomial at the point X0.
            //
        {
            double eps;
            int i;
            int m;
            int n1;
            double w;
            double z;

            n1 = Math.Min(n, iopt);
            n1 = Math.Max(1, n1);

            if (iopt < -1)
            {
                n1 = n;
            }

            eps = (double)(Math.Max(-iopt, 0) % 2);

            w = -(double)n * eps;

            if (-2 < iopt)
            {
                w = w + x0;
            }

            for (m = 1; m <= n1; m++)
            {
                val = 0.0;
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
                    w = w + 1.0;
                }
            }

        }

        public static void r8poly_f2p(int n, ref double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_F2P converts a polynomial from factorial form to power sum form.
            //
            //  Discussion:
            //
            //    The (falling) factorial form is
            //
            //      p(x) =   a(1)
            //             + a(2) * x
            //             + a(3) * x*(x-1)
            //             ...
            //             + a(n) * x*(x-1)*...*(x-(n-2))
            //
            //    The power sum form is
            //
            //      p(x) = a(1) + a(2)*x + a(3)*x^2 + ... + a(n)*x^(n-1)
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
            //    Input, int N, the dimension of A.
            //
            //    Input/output, double A[N], on input, the polynomial
            //    coefficients in factorial form.  On output, the polynomial
            //    coefficients in power sum form.
            //
        {
            int i;
            int m;
            double val;
            double w;
            double z;

            w = -(double)n;

            for (m = 1; m <= n; m++)
            {
                val = 0.0;
                z = w;

                for (i = m; i <= n; i++)
                {
                    z = z + 1.0;
                    val = a[n + m - i - 1] + z * val;
                    a[n + m - i - 1] = val;
                }

                w = w + 1.0;
            }
        }

        public static double r8poly_fval(int n, double[] a, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_FVAL evaluates a double polynomial in factorial form.
            //
            //  Discussion:
            //
            //    The (falling) factorial form of a polynomial is:
            //
            //      p(x) = a(1)
            //           + a(2)  *x
            //           + a(3)  *x*(x-1)
            //           +...
            //           + a(n-1)*x*(x-1)*(x-2)...*(x-(n-3))
            //           + a(n)  *x*(x-1)*(x-2)...*(x-(n-3))*(x-(n-2))
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
            //    Input, int N, the dimension of A.
            //
            //    Input, double A[N], the coefficients of the polynomial.
            //    A(1) is the constant term.
            //
            //    Input, double X, the point at which the polynomial is to be evaluated.
            //
            //    Output, double R8POLY_FVAL, the value of the polynomial at X.
            //
        {
            int i;
            double value;

            value = 0.0;

            for (i = 1; i <= n; i++)
            {
                value = a[n - i] + (x - n + i) * value;
            }

            return value;
        }

        public static void r8poly_n2p(int n, ref double[] a, ref double[] xarray)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_N2P converts a polynomial from Newton form to power sum form.
            //
            //  Discussion:
            //
            //    This is done by shifting all the Newton abscissas to zero.
            //
            //    Actually, what happens is that the abscissas of the Newton form
            //    are all shifted to zero, which means that A is the power sum
            //    polynomial description and A, XARRAY is the Newton polynomial
            //    description.  It is only because all the abscissas are shifted to
            //    zero that A can be used as both a power sum and Newton polynomial
            //    coefficient array.
            //
            //    The Newton form of a polynomial is described by an array of N coefficients
            //    A and N abscissas X:
            //
            //      p(x) =   a(1)
            //             + a(2) * (x-x(1))
            //             + a(3) * (x-x(1)) * (x-x(2))
            //             ...
            //             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
            //
            //    X(N) does not occur explicitly in the formula for the evaluation of p(x),
            //    although it is used in deriving the coefficients A.
            //
            //    The power sum form of a polynomial is:
            //
            //      p(x) = a(1) + a(2)*x + ... + a(n-1)*x^(n-2) + a(n)*x^(n-1)
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
            //    Input, int N, the dimension of A.
            //
            //    Input/output, double A[N].  On input, the coefficients
            //    of the polynomial in Newton form, and on output, the coefficients
            //    in power sum form.
            //
            //    Input/output, double XARRAY[N].  On input, the abscissas of
            //    the Newton form of the polynomial.  On output, these values
            //    have all been set to zero.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                r8poly_nx(n, ref a, ref xarray, 0.0);
            }

        }

        public static double r8poly_nval(int n, double[] a, double[] xarray, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_NVAL evaluates a double polynomial in Newton form.
            //
            //  Definition:
            //
            //    The Newton form of a polynomial is;
            //
            //      p(x) = a(1)
            //           + a(2)  *(x-x1)
            //           + a(3)  *(x-x1)*(x-x2)
            //           +...
            //           + a(n-1)*(x-x1)*(x-x2)*(x-x3)...*(x-x(n-2))
            //           + a(n)  *(x-x1)*(x-x2)*(x-x3)...*(x-x(n-2))*(x-x(n-1))
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 July 2003
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
            //    Input, int N, the dimension of A.
            //
            //    Input, double A[N], the coefficients of the polynomial.
            //    A(1) is the constant term.
            //
            //    Input, double XARRAY[N-1], the N-1 points X which are part
            //    of the definition of the polynomial.
            //
            //    Input, double X, the point at which the polynomial is to be evaluated.
            //
            //    Output, double R8POLY_NVAL, the value of the polynomial at X.
            //
        {
            int i;
            double value;

            value = a[n - 1];

            for (i = n - 2; 0 <= i; i--)
            {
                value = a[i] + (x - xarray[i]) * value;
            }

            return value;
        }

        public static void r8poly_nx(int n, ref double[] a, ref double[] xarray, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_NX replaces one of the base points in a polynomial in Newton form.
            //
            //  Discussion:
            //
            //    The Newton form of a polynomial is described by an array of N coefficients
            //    A and N abscissas X:
            //
            //      p(x) =   a(1)
            //             + a(2) * (x-x(1))
            //             + a(3) * (x-x(1)) * (x-x(2))
            //             ...
            //             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
            //
            //    X(N) does not occur explicitly in the formula for the evaluation of p(x),
            //    although it is used in deriving the coefficients A.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 July 2003
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
            //   Input, int N, the dimension of A.
            //
            //   Input/output, double A[N], the polynomial coefficients of the Newton form.
            //
            //   Input/output, double XARRAY[N], the set of abscissas that
            //   are part of the Newton form of the polynomial.  On output,
            //   the abscissas have been shifted up one index, so that
            //   the first location now holds X, and the original value
            //   of the last entry is discarded.
            //
            //   Input, double X, the new point to be shifted into XARRAY.
            //
        {
            int i;

            for (i = n - 2; 0 <= i; i--)
            {
                a[i] = a[i] + (x - xarray[i]) * a[i + 1];
            }

            for (i = n - 1; 0 < i; i--)
            {
                xarray[i] = xarray[i - 1];
            }

            xarray[0] = x;

        }

        public static void r8poly_p2f(int n, ref double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_P2F converts a double polynomial from power sum form to factorial form.
            //
            //  Discussion:
            //
            //    The power sum form is
            //
            //      p(x) = a(1) + a(2)*x + a(3)*x^2 + ... + a(n)*x^(n-1)
            //
            //    The (falling) factorial form is
            //
            //      p(x) =   a(1)
            //             + a(2) * x
            //             + a(3) * x*(x-1)
            //             ...
            //             + a(n) * x*(x-1)*...*(x-(n-2))
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
            //    Input, int N, the dimension of A.
            //
            //    Input/output, double A[N], on input, the polynomial
            //    coefficients in the power sum form, on output, the polynomial
            //    coefficients in factorial form.
            //
        {
            int i;
            int m;
            double val;

            for (m = 1; m <= n; m++)
            {
                val = 0.0;
                for (i = m; i <= n; i++)
                {
                    val = a[n + m - i - 1] + (double)(m - 1) * val;
                    a[n + m - i - 1] = val;
                }
            }

        }

        public static void r8poly_p2n(int n, ref double[] a, double[] xarray)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_P2N converts a polynomial from power sum form to Newton form.
            //
            //  Discussion:
            //
            //    This is done by shifting all the Newton abscissas from zero.
            //
            //    The power sum form of a polynomial is:
            //
            //      p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
            //
            //    The Newton form of a polynomial is described by an array of N coefficients
            //    A and N abscissas X:
            //
            //      p(x) =   a(1)
            //             + a(2) * (x-x(1))
            //             + a(3) * (x-x(1)) * (x-x(2))
            //             ...
            //             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
            //
            //    X(N) does not occur explicitly in the formula for the evaluation of p(x),
            //    although it is used in deriving the coefficients A.
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
            //    Input, int N, the dimension of A.
            //
            //    Input/output, double A[N].  On input, the coefficients
            //    of the polynomial in power sum form, and on output, the
            //    coefficients in Newton form.
            //
            //    Input, double XARRAY[N].  On input, the desired abscissas of
            //    the Newton form of the polynomial.
            //
        {
            int i;
            double[] work;
            double x;

            work = new double[n];

            for (i = 0; i < n; i++)
            {
                work[i] = 0.0;
            }

            for (i = n - 1; 0 <= i; i--)
            {
                x = xarray[i];
                r8poly_nx(n, ref a, ref work, x);
            }
        }

        public static void r8poly_p2t(int n, ref double[] a, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_P2T converts a real polynomial from power sum form to Taylor form.
            //
            //  Discussion:
            //
            //    The power sum form is
            //
            //      p(x) = a(1) + a(2)*x + a(3)*x^2 + ... + a(n)*x^(n-1)
            //
            //    The Taylor form is
            //
            //      p(x) =   a(1)
            //             + a(2) * (x-x0)
            //             + a(3) * (x-x0)^2
            //             ...
            //             + a(n) * (x-x0)^(n-1)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of A.
            //
            //    Input/output, double A[N], on input, the coefficients in
            //    power sum form, and on output, the coefficients in Taylor form.
            //
            //    Input, double X, the point at which the Taylor form of the
            //    polynomial is to be based.
            //
        {
            int i;
            int m;
            double val;

            for (m = 1; m <= n; m++)
            {
                val = 0.0;
                for (i = m; i <= n; i++)
                {
                    val = a[n + m - i - 1] + x * val;
                    a[n + m - i - 1] = val;
                }
            }
        }

    }
}