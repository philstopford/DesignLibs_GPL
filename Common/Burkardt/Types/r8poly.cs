using System;
using System.Globalization;
using System.Numerics;
using Burkardt.MatrixNS;
using Burkardt.PolynomialNS;

namespace Burkardt.Types;

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
            poly_cof2[i] = poly_cof[i - 1] / i;
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
        double value = 0;

        value = 0.0;

        for (i = n - 1; 0 <= i; i--)
        {
            value = (value + poly_cof[i] / (i + 1)) * xval;
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

        double[] pointer1 = poly_cof;

        for (i = 0; i < ntab; i++)
        {
            double[] pointer2 = pointer1;

            int j;
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

    public static int r8poly_degree(int na, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_DEGREE returns the degree of a polynomial.
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
        //    Input, double A[NA+1], the coefficients of the polynomials.
        //
        //    Output, int R8POLY_DEGREE, the degree of the polynomial.
        //
    {
        int degree = na;

        while (0 < degree)
        {
            if (a[degree] != 0.0)
            {
                return degree;
            }

            degree -= 1;

        }

        return degree;
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
            poly_cof2[i] = (i + 1) * poly_cof[i + 1];
        }
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

        double value = (n - 1) * poly_cof[n - 1];

        for (i = n - 2; 1 <= i; i--)
        {
            value = value * xval + i * poly_cof[i];
        }

        return value;
    }

    public static double[] r8poly_deriv(int n, double[] c, int p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_DERIV returns the derivative of a polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the degree of the polynomial.
        //
        //    Input, double C[N+1], the polynomial coefficients.
        //    C[I] is the coefficient of X^I.
        //
        //    Input, int P, the order of the derivative.
        //    0 means no derivative is taken.
        //    1 means first derivative,
        //    2 means second derivative and so on.
        //    Values of P less than 0 are meaningless.  Values of P greater
        //    than N are meaningful, but the code will behave as though the
        //    value of P was N+1.
        //
        //    Output, double R8POLY_DERIV CP[N-P+1], the polynomial coefficients of
        //    the derivative.
        //
    {
        int d;

        if (n < p)
        {
            return null;
        }

        double[] cp_temp = r8vec_copy_new(n + 1, c);

        for (d = 1; d <= p; d++)
        {
            int i;
            for (i = 0; i <= n - d; i++)
            {
                cp_temp[i] = (i + 1) * cp_temp[i + 1];
            }

            cp_temp[n - d + 1] = 0.0;
        }

        double[] cp = r8vec_copy_new(n - p + 1, cp_temp);

        return cp;
    }

    public static double r8poly_lagrange_0(int npol, double[] xpol, double xval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_LAGRANGE_0 evaluates the Lagrange factor at a point.
        //
        //  Discussion:
        //
        //    W(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
        //
        //  Discussion:
        //
        //    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
        //    polynomial L(IPOL)(X), has the property:
        //
        //      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
        //
        //    and may be expressed as:
        //
        //      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NPOL, the number of abscissas.
        //    NPOL must be at least 1.
        //
        //    Input, double XPOL[NPOL], the abscissas, which should be distinct.
        //
        //    Input, double XVAL, the point at which the Lagrange factor is to be
        //    evaluated.
        //
        //    Output, double R8POLY_LAGRANGE_0, the value of the Lagrange factor at XVAL.
        //
    {
        int i;

        double wval = 1.0;
        for (i = 0; i < npol; i++)
        {
            wval *= xval - xpol[i];
        }

        return wval;
    }

    public static double r8poly_lagrange_1(int npol, double[] xpol, double xval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_LAGRANGE_1 evaluates the first derivative of the Lagrange factor.
        //
        //  Discussion:
        //
        //    W(XPOL(1:NPOL))(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
        //
        //    W'(XPOL(1:NPOL))(X)
        //      = Sum ( 1 <= J <= NPOL ) Product ( I /= J ) ( X - XPOL(I) )
        //
        //    We also have the recursion:
        //
        //      W'(XPOL(1:NPOL))(X) = d/dX ( ( X - XPOL(NPOL) ) * W(XPOL(1:NPOL-1))(X) )
        //                    = W(XPOL(1:NPOL-1))(X)
        //                    + ( X - XPOL(NPOL) ) * W'(XPOL(1:NPOL-1))(X)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NPOL, the number of abscissas.
        //
        //    Input, double XPOL[NPOL], the abscissas, which should be distinct.
        //
        //    Input, double XVAL, the point at which the Lagrange factor is to be
        //    evaluated.
        //
        //    Output, double R8POLY_LAGRANGE_1, the derivative of W with respect to XVAL.
        //
    {
        int i;

        double dwdx = 0.0;
        double w = 1.0;

        for (i = 0; i < npol; i++)
        {
            dwdx = w + (xval - xpol[i]) * dwdx;
            w *= xval - xpol[i];
        }

        return dwdx;
    }

    public static double r8poly_lagrange_2(int npol, double[] xpol, double xval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_LAGRANGE_2 evaluates the second derivative of the Lagrange factor.
        //
        //  Discussion:
        //
        //    W(X)  = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
        //
        //    W'(X) = Sum ( 1 <= J <= NPOL )
        //            Product ( I /= J ) ( X - XPOL(I) )
        //
        //    W"(X) = Sum ( 1 <= K <= NPOL )
        //            Sum ( J =/ K )
        //            Product ( I /= K, J ) ( X - XPOL(I) )
        //
        //    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
        //    polynomial L(IPOL)(X), has the property:
        //
        //      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
        //
        //    and may be expressed as:
        //
        //      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NPOL, the number of abscissas.
        //    NPOL must be at least 1.
        //
        //    Input, double XPOL[NPOL], the abscissas, which should be distinct.
        //
        //    Input, double XVAL, the point at which the Lagrange factor is to be
        //    evaluated.
        //
        //    Output, double R8POLY_LAGRANGE_2, the second derivative of W with respect to XVAL.
        //
    {
        int k;

        double dw2dx2 = 0.0;

        for (k = 0; k < npol; k++)
        {
            int j;
            for (j = 0; j < npol; j++)
            {
                if (j == k)
                {
                    continue;
                }

                double term = 1.0;
                int i;
                for (i = 0; i < npol; i++)
                {
                    if (i != j && i != k)
                    {
                        term *= xval - xpol[i];
                    }
                }

                dw2dx2 += term;
            }
        }

        return dw2dx2;
    }

    public static double[] r8poly_lagrange_coef(int npol, int ipol, double[] xpol)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_LAGRANGE_COEF returns the coefficients of a Lagrange polynomial.
        //
        //  Discussion:
        //
        //    Given NPOL distinct abscissas, XPOL(*), the IPOL-th Lagrange
        //    polynomial P(IPOL)(X) is defined as the polynomial of degree
        //    NPOL - 1 which is 1 at XPOL(IPOL) and 0 at the NPOL - 1 other
        //    abscissas.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NPOL, the number of abscissas.
        //    NPOL must be at least 1.
        //
        //    Input, int IPOL, the index of the polynomial to evaluate.
        //    IPOL must be between 1 and NPOL.
        //
        //    Input, double XPOL[NPOL], the abscissas of the Lagrange polynomials.
        //    The entries in XPOL must be distinct.
        //
        //    Output, double R8POLY_LAGRANGE_COEF[NPOL], the polynomial coefficients
        //    of the IPOL-th Lagrange polynomial.
        //
    {
        int i;
        //
        //  Make sure IPOL is legal.
        //
        if (ipol < 1 || npol < ipol)
        {
            Console.WriteLine("");
            Console.WriteLine("R8POLY_LAGRANGE_COEF - Fatal error!");
            Console.WriteLine("  1 <= IPOL <= NPOL is required.");
            Console.WriteLine("  but IPOL = " + ipol + "");
            Console.WriteLine("  and NPOL = " + npol + "");
            return null;
        }

        //
        //  Check that the abscissas are distinct.
        //
        if (!r8vec_is_distinct(npol, xpol))
        {
            Console.WriteLine("");
            Console.WriteLine("R8POLY_LAGRANGE_COEF - Fatal error!");
            Console.WriteLine("  Two entries of XPOL are equal:");
            return null;
        }

        double[] pcof = new double[npol];

        pcof[0] = 1.0;
        for (i = 1; i < npol; i++)
        {
            pcof[i] = 0.0;
        }

        int index = 0;

        for (i = 1; i <= npol; i++)
        {
            if (i == ipol)
            {
                continue;
            }

            index += 1;

            int j;
            for (j = index; 0 <= j; j--)
            {
                pcof[j] = -xpol[i - 1] * pcof[j] / (xpol[ipol - 1] - xpol[i - 1]);

                switch (j)
                {
                    case > 0:
                        pcof[j] += pcof[j - 1] / (xpol[ipol - 1] - xpol[i - 1]);
                        break;
                }
            }
        }

        return pcof;
    }

    public static void r8poly_lagrange_factor(int npol, double[] xpol, double xval,
            ref double wval, ref double dwdx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_LAGRANGE_FACTOR evaluates the polynomial Lagrange factor at a point.
        //
        //  Discussion:
        //
        //    Suppose F(X) is at least N times continuously differentiable in the
        //    interval [A,B].  Pick NPOL distinct points XPOL(I) in [A,B] and compute
        //    the interpolating polynomial P(X) of order NPOL ( and degree NPOL-1)
        //    which passes through all the points ( XPOL(I), F(XPOL(I)) ).
        //    Then in the interval [A,B], the maximum error
        //
        //      abs ( F(X) - P(X) )
        //
        //    is bounded by:
        //
        //      C * FNMAX * W(X)
        //
        //    where
        //
        //      C is a constant,
        //      FNMAX is the maximum value of the NPOL-th derivative of F in [A,B],
        //      W(X) is the Lagrange factor.
        //
        //    Thus, the value of W(X) is useful as part of an estimated bound
        //    for the interpolation error.
        //
        //    The formula is:
        //
        //      W(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
        //
        //    Note that the Chebyshev abscissas have the property that they minimize
        //    the value of W(X) over the interval [A,B].  Hence, if the abscissas may
        //    be chosen arbitrarily, the Chebyshev abscissas have this advantage over
        //    other choices.
        //
        //    For a set of points XPOL[I], 0 <= I <= NPOL-1, the IPOL-th Lagrange basis
        //    polynomial L(IPOL)(X), has the property:
        //
        //      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
        //
        //    and may be expressed as:
        //
        //      L(IPOL)(X) = W(X) / ( ( X - XPOL[IPOL] ) * W'(XPOL[IPOL]) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 May 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NPOL, the number of abscissas.
        //    NPOL must be at least 1.
        //
        //    Input, double XPOL[NPOL], the abscissas, which should be distinct.
        //
        //    Input, double XVAL, the point at which the Lagrange factor is to be evaluated.
        //
        //    Output, double &WVAL, the value of the Lagrange factor at XVAL.
        //
        //    Output, double &DWDX, the derivative of W with respect to XVAL.
        //
    {
        int i;

        wval = 1.0;
        for (i = 0; i < npol; i++)
        {
            wval *= xval - xpol[i];
        }

        dwdx = 0.0;

        for (i = 0; i < npol; i++)
        {
            double term = 1.0;

            int j;
            for (j = 0; j < npol; j++)
            {
                if (i != j)
                {
                    term *= xval - xpol[j];
                }
            }

            dwdx += term;
        }

    }

    public static int r8poly_lagrange_val(int npol, int ipol, double[] xpol, double xval,
            ref double pval, ref double dpdx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_LAGRANGE_VAL evaluates the IPOL-th Lagrange polynomial.
        //
        //  Discussion:
        //
        //    Given NPOL distinct abscissas, XPOL[*], the IPOL-th Lagrange
        //    polynomial P(IPOL)(X) is defined as the polynomial of degree
        //    NPOL - 1 which is 1 at XPOL[IPOL] and 0 at the NPOL - 1 other
        //    abscissas.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 May 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NPOL, the number of abscissas.
        //    NPOL must be at least 1.
        //
        //    Input, int IPOL, the index of the polynomial to evaluate.
        //    IPOL must be between 0 and NPOL-1.
        //
        //    Input, double XPOL[NPOL], the abscissas of the Lagrange polynomials.
        //    The entries in XPOL must be distinct.
        //
        //    Input, double XVAL, the point at which the IPOL-th Lagrange polynomial
        //    is to be evaluated.
        //
        //    Output, double &PVAL, the value of the IPOL-th Lagrange polynomial at XVAL.
        //
        //    Output, double &DPDX, the derivative of the IPOL-th Lagrange polynomial at XVAL.
        //
        //    Output, int R8POLY_LAGRANGE_VAL, 0 if no error.
        //
    {
        int i;
        int j;
        //
        //  Make sure IPOL is legal.
        //
        if (ipol < 0 || npol - 1 < ipol)
        {
            Console.WriteLine("");
            Console.WriteLine("R8POLY_LAGRANGE_VAL - Fatal error!");
            Console.WriteLine("  0 <= IPOL <= NPOL-1 is required.");
            return 1;
        }

        //
        //  Check that the abscissas are distinct.
        //
        for (i = 1; i < npol; i++)
        {
            for (j = 0; j < i; j++)
            {
                if (!(Math.Abs(xpol[i] - xpol[j]) <= typeMethods.r8_epsilon()))
                {
                    continue;
                }

                Console.WriteLine("");
                Console.WriteLine("R8POLY_LAGRANGE_VAL - Fatal error!");
                Console.WriteLine("  Two entries of XPOL are equal:");
                Console.WriteLine("  XPOL(" + i + ") = " + xpol[i] + ".");
                Console.WriteLine("  XPOL(" + j + ") = " + xpol[j] + ".");
                return 1;
            }
        }

        //
        //  Evaluate the polynomial.
        //
        pval = 1.0;

        for (i = 0; i < npol; i++)
        {
            if (i != ipol)
            {
                pval = pval * (xval - xpol[i]) / (xpol[ipol] - xpol[i]);
            }
        }

        //
        //  Evaluate the derivative, which can be found by summing up the result
        //  of differentiating one factor at a time, successively.
        //
        dpdx = 0.0;

        for (i = 0; i < npol; i++)
        {
            if (i == ipol)
            {
                continue;
            }

            double p2 = 1.0;

            for (j = 0; j < npol; j++)
            {
                if (j == i)
                {
                    p2 /= xpol[ipol] - xpol[j];
                }
                else if (j != ipol)
                {
                    p2 = p2 * (xval - xpol[j]) / (xpol[ipol] - xpol[j]);
                }
            }

            dpdx += p2;
        }

        return 0;
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
        int value = na;

        while (1 < value)
        {
            if (a[value - 1] != 0.0)
            {
                return value;
            }

            value -= 1;
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

        switch (title.Length)
        {
            case > 0:
                Console.WriteLine("");
                Console.WriteLine(title + "");
                break;
        }

        Console.WriteLine("");

        switch (n)
        {
            case < 0:
                Console.WriteLine("  p(x) = 0");
                return;
        }

        char plus_minus = a[n] switch
        {
            < 0.0 => '-',
            _ => ' '
        };

        double mag = Math.Abs(a[n]);

        switch (n)
        {
            case >= 2:
                Console.WriteLine("  p(x) = " + plus_minus
                                              + mag.ToString(CultureInfo.InvariantCulture).PadLeft(14) + " * x ^ " + n + "");
                break;
            case 1:
                Console.WriteLine("  p(x) = " + plus_minus
                                              + mag.ToString(CultureInfo.InvariantCulture).PadLeft(14) + " * x");
                break;
            case 0:
                Console.WriteLine("  p(x) = " + plus_minus
                                              + mag.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                break;
        }

        for (i = n - 1; 0 <= i; i--)
        {
            plus_minus = a[aIndex + i] switch
            {
                < 0.0 => '-',
                _ => '+'
            };

            mag = Math.Abs(a[aIndex + i]);

            if (mag != 0.0)
            {
                switch (i)
                {
                    case >= 2:
                        Console.WriteLine("         " + plus_minus
                                                      + mag.ToString(CultureInfo.InvariantCulture).PadLeft(14) + " * x ^ " + i + "");
                        break;
                    case 1:
                        Console.WriteLine("         " + plus_minus
                                                      + mag.ToString(CultureInfo.InvariantCulture).PadLeft(14) + " * x");
                        break;
                    case 0:
                        Console.WriteLine("         " + plus_minus
                                                      + mag.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                        break;
                }
            }
        }
    }

    public static double r8poly_pval(int n, double[] a, double x)

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
        double value = 0;

        value = 0.0;
        for (i = n; 0 <= i; i--)
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
                poly_cof[j - 1] /= scale;
            }
        }

        for (i = 1; i <= n; i++)
        {
            for (j = n - 1; i <= j; j--)
            {
                poly_cof[j - 1] -= shift * poly_cof[j];
            }
        }
    }

    public static void r8poly_t2p(int n, ref double[] a, double x)

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

        for (i = n; 1 <= i; i--)
        {
            int j;
            for (j = i; j <= n - 1; j++)
            {
                a[j - 1] -= a[j] * x;
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
        double value = 0;

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
        int m;

        int n1 = Math.Min(n, iopt);

        n1 = iopt switch
        {
            < -1 => n,
            _ => Math.Max(1, n1)
        };

        double eps = Math.Max(-iopt, 0) % 2;

        double w = -(double) n * eps;

        switch (iopt)
        {
            case > -2:
                w += x0;
                break;
        }

        for (m = 1; m <= n1; m++)
        {
            val = 0.0;
            double z = w;

            int i;
            for (i = m; i <= n; i++)
            {
                z += eps;
                val = a[n + m - i - 1] + z * val;
                if (iopt != 0 && iopt != -1)
                {
                    a[n + m - i - 1] = val;
                }
            }

            switch (iopt)
            {
                case < 0:
                    w += 1.0;
                    break;
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
        int m;

        double w = -(double) n;

        for (m = 1; m <= n; m++)
        {
            double val = 0.0;
            double z = w;

            int i;
            for (i = m; i <= n; i++)
            {
                z += 1.0;
                val = a[n + m - i - 1] + z * val;
                a[n + m - i - 1] = val;
            }

            w += 1.0;
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
        double value = 0;

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

        double value = a[n - 1];

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
            a[i] += (x - xarray[i]) * a[i + 1];
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
        int m;

        for (m = 1; m <= n; m++)
        {
            double val = 0.0;
            int i;
            for (i = m; i <= n; i++)
            {
                val = a[n + m - i - 1] + (m - 1) * val;
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

        double[] work = new double[n];

        for (i = 0; i < n; i++)
        {
            work[i] = 0.0;
        }

        for (i = n - 1; 0 <= i; i--)
        {
            double x = xarray[i];
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
        int m;

        for (m = 1; m <= n; m++)
        {
            double val = 0.0;
            int i;
            for (i = m; i <= n; i++)
            {
                val = a[n + m - i - 1] + x * val;
                a[n + m - i - 1] = val;
            }
        }
    }

    public static double[] r8poly_values_horner(int m, double[] c, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_VALUES_HORNER evaluates a polynomial using Horner's method.
        //
        //  Discussion:
        //
        //    The polynomial 
        //
        //      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
        //
        //    is to be evaluated at the vector of values X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the degree.
        //
        //    Input, double C[M+1], the polynomial coefficients.  
        //    C[I] is the coefficient of X^I.
        //
        //    Input, int N, the number of evaluation points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double R8POLY_VALUES_HORNER[N], the polynomial values.
        //
    {
        int i;
        int j;

        double[] p = new double[n];

        for (j = 0; j < n; j++)
        {
            p[j] = c[m];
        }

        for (i = m - 1; 0 <= i; i--)
        {
            for (j = 0; j < n; j++)
            {
                p[j] = p[j] * x[j] + c[i];
            }
        }

        return p;
    }

    public static double[] r8poly_value_2d(int m, double[] c, int n, double[] x, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_VALUE_2D evaluates a polynomial in 2 variables, X and Y.
        //
        //  Discussion:
        //
        //    We assume the polynomial is of total degree M, and has the form:
        //
        //      p(x,y) = c00 
        //             + c10 * x                + c01 * y
        //             + c20 * x^2   + c11 * xy + c02 * y^2
        //             + ...
        //             + cm0 * x^(m) + ...      + c0m * y^m.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the degree of the polynomial.
        //
        //    Input, double C[T(M+1)], the polynomial coefficients.  
        //    C[0] is the constant term.  T(M+1) is the M+1-th triangular number.
        //    The coefficients are stored consistent with the following ordering
        //    of monomials: 1, X, Y, X^2, XY, Y^2, X^3, X^2Y, XY^2, Y^3, X^4, ...
        //
        //    Input, int N, the number of evaluation points.
        //
        //    Input, double X[N], Y[N], the evaluation points.
        //
        //    Output, double R8POLY_VALUE_2D[N], the value of the polynomial at the 
        //    evaluation points.
        //
    {
        int i;
        int s;

        double[] p = new double[n];

        for (i = 0; i < n; i++)
        {
            p[i] = 0.0;
        }

        int j = 0;
        for (s = 0; s <= m; s++)
        {
            int ex;
            for (ex = s; 0 <= ex; ex--)
            {
                int ey = s - ex;
                for (i = 0; i < n; i++)
                {
                    p[i] += c[j] * Math.Pow(x[i], ex) * Math.Pow(y[i], ey);
                }

                j += 1;
            }
        }

        return p;
    }

    public static int r8poly2_ex(double x1, double y1, double x2, double y2, double x3,
            double y3, ref double x, ref double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY2_EX finds the extremal point of a parabola determined by three points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 October 1998
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
        //    on the parabola.  X1, X2 and X3 must be distinct.
        //
        //    Output, double &X, &Y, the X coordinate of the extremal point of the
        //    parabola, and the value of the parabola at that point.
        //
        //    Output, int R8POLY2_EX, error flag.
        //    0, no error.
        //    1, two of the X values are equal.
        //    2, the data lies on a straight line; there is no finite extremal
        //    point.
        //    3, the data lies on a horizontal line; every point is "extremal".
        //
    {
        x = 0.0;
        y = 0.0;

        if (Math.Abs(x1 - x2) <= typeMethods.r8_epsilon() || Math.Abs(x2 - x3) <= typeMethods.r8_epsilon() || Math.Abs(x3 - x1) <= typeMethods.r8_epsilon())
        {
            return 1;
        }

        if (Math.Abs(y1 - y2) <= typeMethods.r8_epsilon() && Math.Abs(y2 - y3) <= typeMethods.r8_epsilon() && Math.Abs(y3 - y1) <= typeMethods.r8_epsilon())
        {
            x = x1;
            y = y1;
            return 3;
        }

        double bot = (x2 - x3) * y1 + (x3 - x1) * y2 + (x1 - x2) * y3;

        switch (bot)
        {
            case 0.0:
                return 2;
        }

        x = 0.5 * (
                x1 * x1 * (y3 - y2)
                + x2 * x2 * (y1 - y3)
                + x3 * x3 * (y2 - y1)) /
            ((x2 - x3) * y1 + (x3 - x1) * y2 + (x1 - x2) * y3);

        y = -(
                (x - x2) * (x - x3) * (x2 - x3) * y1
                + (x - x1) * (x - x3) * (x3 - x1) * y2
                + (x - x1) * (x - x2) * (x1 - x2) * y3) /
            ((x1 - x2) * (x2 - x3) * (x3 - x1));

        return 0;
    }

    public static int r8poly2_ex2(double x1, double y1, double x2, double y2, double x3,
            double y3, ref double x, ref double y, ref double a, ref double b, ref double c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY2_EX2 finds the extremal point of a parabola determined by three points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
        //    on the parabola.  X1, X2 and X3 must be distinct.
        //
        //    Output, double &X, &Y, the X coordinate of the extremal point of the
        //    parabola, and the value of the parabola at that point.
        //
        //    Output, double &A, &B, &C, the coefficients that define the parabola:
        //    P(X) = A * X^2 + B * X + C.
        //
        //    Output, int R8POLY2_EX2, error flag.
        //    0, no error.
        //    1, two of the X values are equal.
        //    2, the data lies on a straight line; there is no finite extremal
        //    point.
        //    3, the data lies on a horizontal line; any point is an "extremal point".
        //
    {
        double[] v = new double[3 * 3];

        a = 0.0;
        b = 0.0;
        c = 0.0;
        x = 0.0;
        y = 0.0;

        if (Math.Abs(x1 - x2) <= typeMethods.r8_epsilon() || Math.Abs(x2 - x3) <= typeMethods.r8_epsilon() || Math.Abs(x3 - x1) <= typeMethods.r8_epsilon())
        {
            return 1;
        }

        if (Math.Abs(y1 - y2) <= typeMethods.r8_epsilon() && Math.Abs(y2 - y3) <= typeMethods.r8_epsilon() && Math.Abs(y3 - y1) <= typeMethods.r8_epsilon())
        {
            x = x1;
            y = y1;
            return 3;
        }

        //
        //  Set up the Vandermonde matrix.
        //
        v[0 + 0 * 3] = 1.0;
        v[0 + 1 * 3] = x1;
        v[0 + 2 * 3] = x1 * x1;

        v[1 + 0 * 3] = 1.0;
        v[1 + 1 * 3] = x2;
        v[1 + 2 * 3] = x2 * x2;

        v[2 + 0 * 3] = 1.0;
        v[2 + 1 * 3] = x3;
        v[2 + 2 * 3] = x3 * x3;
        //
        //  Get the inverse.
        //
        double[] w = r8mat_inverse_3d(v);
        //
        //  Compute the parabolic coefficients.
        //
        c = w[0 + 0 * 3] * y1 + w[0 + 1 * 3] * y2 + w[0 + 2 * 3] * y3;
        b = w[1 + 0 * 3] * y1 + w[1 + 1 * 3] * y2 + w[1 + 2 * 3] * y3;
        a = w[2 + 0 * 3] * y1 + w[2 + 1 * 3] * y2 + w[2 + 2 * 3] * y3;
        switch (a)
        {
            //
            //  Determine the extremal point.
            //
            case 0.0:
                return 2;
        }

        x = -b / (2.0 * a);
        y = a * x * x + b * x + c;

        return 0;
    }

    public static void r8poly2_root(double a, double b, double c, ref Complex r1,
            ref Complex r2)

        //****************************************************************************//
        //
        //  Purpose:
        //
        //    R8POLY2_ROOT returns the two roots of a quadratic polynomial.
        //
        //  Discussion:
        //
        //    The polynomial has the form:
        //
        //      A * X * X + B * X + C = 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the coefficients of the polynomial.
        //    A must not be zero.
        //
        //    Output, Complex &R1, &R2, the roots of the polynomial, which
        //    might be real and distinct, real and equal, or complex conjugates.
        //
    {
        switch (a)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8POLY2_ROOT - Fatal error!");
                Console.WriteLine("  The coefficient A is zero.");
                return;
        }

        Complex disc = b * b - 4.0 * a * c;
        Complex q = -0.5 * (b + r8_sign(b) * Complex.Sqrt(disc));
        r1 = q / a;
        r2 = c / q;

    }

    public static void r8poly2_rroot(double a, double b, double c, ref double r1, ref double r2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY2_RROOT returns the real parts of the roots of a quadratic polynomial.
        //
        //  Example:
        //
        //    A    B    C       roots              R1   R2
        //   --   --   --     ------------------   --   --
        //    1   -4    3     1          3          1    3
        //    1    0    4     2*i      - 2*i        0    0
        //    1   -6   10     3 +   i    3 -   i    3    3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 December 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the coefficients of the quadratic
        //    polynomial A * X^2 + B * X + C = 0 whose roots are desired.
        //    A must not be zero.
        //
        //    Output, double &R1, &R2, the real parts of the roots
        //    of the polynomial.
        //
    {
        switch (a)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8POLY2_RROOT - Fatal error!");
                Console.WriteLine("  The coefficient A is zero.");
                return;
        }

        double disc = b * b - 4.0 * a * c;
        switch (disc)
        {
            case >= 0.0:
                double q = b + r8_sign(b) * Math.Sqrt(disc);
                r1 = -0.5 * q / a;
                r2 = -2.0 * c / q;
                break;
            default:
                r1 = b / 2.0 / a;
                r2 = b / 2.0 / a;
                break;
        }

    }

    public static void r8poly2_val(double x1, double y1, double x2, double y2,
            double x3, double y3, double x, ref double y, ref double yp, ref double ypp)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY2_VAL evaluates a parabola defined by three data values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X1, Y1, X2, Y2, X3, Y3, three pairs of data values.
        //    If the X values are distinct, then all the Y values represent
        //    actual values of the parabola.
        //
        //    Three special cases are allowed:
        //
        //      X1 = X2 =/= X3: Y2 is the derivative at X1;
        //      X1 =/= X2 = X3: Y3 is the derivative at X3;
        //      X1 = X2 = X3:   Y2 is the derivative at X1, and
        //                      Y3 is the second derivative at X1.
        //
        //    Input, double X, an abscissa at which the parabola is to be
        //    evaluated.
        //
        //    Output, double &Y, &YP, &YPP, the values of the parabola and
        //    its first and second derivatives at X.
        //
    {
        int distinct;
        double dif1;
        double dif2;
        switch (Math.Abs(x1 - x2))
        {
            //
            //  If any X's are equal, put them and the Y data first.
            //
            case <= typeMethods._r8_epsilon when Math.Abs(x2 - x3) <= typeMethods.r8_epsilon():
                distinct = 1;
                break;
            case <= typeMethods._r8_epsilon:
                distinct = 2;
                break;
            default:
            {
                if (Math.Abs(x1 - x3) <= typeMethods.r8_epsilon())
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8POLY2_VAL - Fatal error!");
                    Console.WriteLine("  X1 = X3 =/= X2.");
                    return;
                }

                if (Math.Abs(x2 - x3) <= typeMethods.r8_epsilon())
                {
                    distinct = 2;
                    (x1, x3) = (x3, x1);
                    y1 = y2;
                    y2 = y3;
                    y3 = y1;
                }
                else
                {
                    distinct = 3;
                }

                break;
            }
        }

        switch (distinct)
        {
            //
            //  Set up the coefficients.
            //
            case 1:
                dif1 = y2;
                dif2 = 0.5 * y3;
                break;
            case 2:
                dif1 = y2;
                dif2 = ((y3 - y1) / (x3 - x1)
                        - y2) / (x3 - x2);
                break;
            default:
                dif1 = (y2 - y1) / (x2 - x1);
                dif2 = ((y3 - y1) / (x3 - x1)
                        - (y2 - y1) / (x2 - x1)) / (x3 - x2);
                break;
        }

        //
        //  Evaluate.
        //
        y = y1 + (x - x1) * dif1 + (x - x1) * (x - x2) * dif2;
        yp = dif1 + (2.0 * x - x1 - x2) * dif2;
        ypp = 2.0 * dif2;

    }

    public static void r8poly2_val2(int ndata, double[] tdata,
            double[] ydata, int left, double tval, ref double yval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY2_VAL2 evaluates a parabolic function through 3 points in a table.
        //
        //  Discussion:
        //
        //    This routine is a utility routine used by OVERHAUSER_SPLINE_VAL.
        //    It constructs the parabolic interpolant through the data in
        //    3 consecutive entries of a table and evaluates this interpolant
        //    at a given abscissa value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 March 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NDATA, the number of data points.
        //    NDATA must be at least 3.
        //
        //    Input, double TDATA[NDATA], the abscissas of the data points.  The
        //    values in TDATA must be in strictly ascending order.
        //
        //    Input, double YDATA[NDATA], the data points corresponding to
        //    the abscissas.
        //
        //    Input, int LEFT, the location of the first of the three
        //    consecutive data points through which the parabolic interpolant
        //    must pass.  0 <= LEFT <= NDATA - 3.
        //
        //    Input, double TVAL, the value of T at which the parabolic interpolant
        //    is to be evaluated.  Normally, TDATA[0] <= TVAL <= T[NDATA-1], and
        //    the data will be interpolated.  For TVAL outside this range,
        //    extrapolation will be used.
        //
        //    Output, double &YVAL, the value of the parabolic interpolant
        //    at TVAL.
        //
    {
        //
        //  Check.
        //
        if (left < 0 || ndata - 3 < left)
        {
            Console.WriteLine("");
            Console.WriteLine("RPOLY2_VAL2 - Fatal error!");
            Console.WriteLine("  LEFT < 0 or NDATA-3 < LEFT.");
            return;
        }

        //
        //  Copy out the three abscissas.
        //
        double t1 = tdata[left];
        double t2 = tdata[left + 1];
        double t3 = tdata[left + 2];

        if (t2 <= t1 || t3 <= t2)
        {
            Console.WriteLine("");
            Console.WriteLine("RPOLY2_VAL2 - Fatal error!");
            Console.WriteLine("  T2 <= T1 or T3 <= T2.");
            Console.WriteLine("  T1 = " + t1 + "");
            Console.WriteLine("  T2 = " + t2 + "");
            Console.WriteLine("  T3 = " + t3 + "");
            return;
        }

        //
        //  Construct and evaluate a parabolic interpolant for the data.
        //
        double y1 = ydata[left];
        double y2 = ydata[left + 1];
        double y3 = ydata[left + 2];

        double dif1 = (y2 - y1) / (t2 - t1);
        double dif2 = ((y3 - y1) / (t3 - t1)
                       - (y2 - y1) / (t2 - t1)) / (t3 - t2);

        yval = y1 + (tval - t1) * (dif1 + (tval - t2) * dif2);

    }

    public static void r8poly3_root(double a, double b, double c, double d,
            ref Complex r1, ref Complex r2, ref Complex r3)

        //****************************************************************************//
        /*
        Purpose:
        
        R8POLY3_ROOT returns the three roots of a cubic polynomial.
        
        Discussion:
        
        The polynomial has the form
        
        A * X^3 + B * X^2 + C * X + D = 0
        
        Licensing:
        
        This code is distributed under the GNU LGPL license.
        
        Modified:
        
        09 August 2018
        
        Parameters:
        
        Input, double A, B, C, D, the coefficients of the polynomial.
        A must not be zero.
        
        Output, Complex &R1, &R2, &R3, the roots of the polynomial, which
        will include at least one real root.
        */
    {
        switch (a)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8POLY3_ROOT - Fatal error!");
                Console.WriteLine("  A must not be zero!");
                return;
        }

        Complex i = Complex.Sqrt(-1.0);

        double q = (b / a * (b / a) - 3.0 * (c / a)) / 9.0;

        double r = (2.0 * (b / a) * (b / a) * (b / a) - 9.0 * (b / a) * (c / a)
                    + 27.0 * (d / a)) / 54.0;

        if (r * r < q * q * q)
        {
            double theta = Math.Acos(r / Math.Sqrt(q * q * q));
            r1 = -2.0 * Complex.Sqrt(q) * Complex.Cos(theta / 3.0);
            r2 = -2.0 * Complex.Sqrt(q) * Complex.Cos((theta + 2.0 * Math.PI) / 3.0);
            r3 = -2.0 * Complex.Sqrt(q) * Complex.Cos((theta + 4.0 * Math.PI) / 3.0);
        }
        else if (q * q * q <= r * r)
        {
            double temp = -r + Math.Sqrt(r * r - q * q * q);
            double s1 = r8_sign(temp) * Math.Pow(Math.Abs(temp), 1.0 / 3.0);

            temp = -r - Math.Sqrt(r * r - q * q * q);
            double s2 = r8_sign(temp) * Math.Pow(Math.Abs(temp), 1.0 / 3.0);

            r1 = s1 + s2;
            r2 = -0.5 * (s1 + s2) + i * 0.5 * Complex.Sqrt(3.0) * (s1 - s2);
            r3 = -0.5 * (s1 + s2) - i * 0.5 * Complex.Sqrt(3.0) * (s1 - s2);
        }

        r1 -= b / (3.0 * a);
        r2 -= b / (3.0 * a);
        r3 -= b / (3.0 * a);

    }

    public static void r8poly4_root(double a, double b, double c, double d, double e,
            ref Complex r1, ref Complex r2, ref Complex r3,
            ref Complex r4)

        //****************************************************************************//
        //
        //  Purpose:
        //
        //    R8POLY4_ROOT returns the four roots of a quartic polynomial.
        //
        //  Discussion:
        //
        //    The polynomial has the form:
        //
        //      A * X^4 + B * X^3 + C * X^2 + D * X + E = 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 August 2018
        //
        //  Parameters:
        //
        //    Input, double A, B, C, D, E, the coefficients of the polynomial.
        //    A must not be zero.
        //
        //    Output, Complex &R1, &R2, &R3, &R4, the roots of the polynomial.
        //
    {
        Complex p;
        Complex q;

        Complex zero = 0.0;

        switch (a)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8POLY4_ROOT - Fatal error!");
                Console.WriteLine("  A must not be zero!");
                return;
        }

        double a4 = b / a;
        double b4 = c / a;
        double c4 = d / a;
        double d4 = e / a;
        //
        //  Set the coefficients of the resolvent cubic equation.
        //
        const double a3 = 1.0;
        double b3 = -b4;
        double c3 = a4 * c4 - 4.0 * d4;
        double d3 = -a4 * a4 * d4 + 4.0 * b4 * d4 - c4 * c4;
        //
        //  Find the roots of the resolvent cubic.
        //
        r8poly3_root(a3, b3, c3, d3, ref r1, ref r2, ref r3);
        //
        //  Choose one root of the cubic, here R1.
        //
        //  Set R = sqrt ( 0.25 * A4 ^ 2 - B4 + R1 )
        //
        Complex r = Complex.Sqrt(0.25 * a4 * a4 - b4 + r1);

        if (r != zero)
        {
            p = Complex.Sqrt(0.75 * a4 * a4 - r * r - 2.0 * b4
                             + 0.25 * (4.0 * a4 * b4 - 8.0 * c4 - a4 * a4 * a4) / r);

            q = Complex.Sqrt(0.75 * a4 * a4 - r * r - 2.0 * b4
                             - 0.25 * (4.0 * a4 * b4 - 8.0 * c4 - a4 * a4 * a4) / r);
        }
        else
        {
            p = Complex.Sqrt(0.75 * a4 * a4 - 2.0 * b4 + 2.0 * Complex.Sqrt(r1 * r1 - 4.0 * d4));

            q = Complex.Sqrt(0.75 * a4 * a4 - 2.0 * b4 - 2.0 * Complex.Sqrt(r1 * r1 - 4.0 * d4));
        }

        //
        //  Set the roots.
        //
        r1 = -0.25 * a4 + 0.5 * r + 0.5 * p;
        r2 = -0.25 * a4 + 0.5 * r - 0.5 * p;
        r3 = -0.25 * a4 - 0.5 * r + 0.5 * q;
        r4 = -0.25 * a4 - 0.5 * r - 0.5 * q;

    }


}