using System;
using Burkardt.Types;

namespace Burkardt.PolynomialNS;

public static class LeastSquares
{
    public static void least_set(int point_num, double[] x, double[] f, double[] w,
            int nterms, ref double[] b, ref double[] c, ref double[] d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEAST_SET defines a least squares polynomial for given data.
        //
        //  Discussion:
        //
        //    This routine is based on ORTPOL by Conte and deBoor.
        //
        //    The polynomial may be evaluated at any point X by calling LEAST_VAL.
        //
        //    Thanks to Andrew Telford for pointing out a mistake in the form of
        //    the check that there are enough unique data points, 25 June 2008.
        //
        //    Thanks to Thomas Beutlich for pointing out that S needs to be
        //    freed on return, 10 October 2012.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 October 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Samuel Conte, Carl deBoor.
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Samuel Conte, Carl deBoor,
        //    Elementary Numerical Analysis,
        //    Second Edition,
        //    McGraw Hill, 1972,
        //    ISBN: 07-012446-4,
        //    LC: QA297.C65.
        //
        //  Parameters:
        //
        //    Input, int POINT_NUM, the number of data values.
        //
        //    Input, double X[POINT_NUM], the abscissas of the data points.
        //    At least NTERMS of the values in X must be distinct.
        //
        //    Input, double F[POINT_NUM], the data values at the points X(*).
        //
        //    Input, double W[POINT_NUM], the weights associated with
        //    the data points.  Each entry of W should be positive.
        //
        //    Input, int NTERMS, the number of terms to use in the
        //    approximating polynomial.  NTERMS must be at least 1.
        //    The degree of the polynomial is NTERMS-1.
        //
        //    Output, double B[NTERMS], C[NTERMS], D[NTERMS], are quantities
        //    defining the least squares polynomial for the input data,
        //    which will be needed to evaluate the polynomial.
        //
    {
        int i;
        int j;
        const double tol = 0.0;
        //
        //  Make sure at least NTERMS X values are unique.
        //
        int unique_num = typeMethods.r8vec_unique_count(point_num, x, tol);

        if (unique_num < nterms)
        {
            Console.WriteLine("");
            Console.WriteLine("LEAST_SET - Fatal error!");
            Console.WriteLine("  The number of distinct X values must be");
            Console.WriteLine("  at least NTERMS = " + nterms + "");
            Console.WriteLine("  but the input data has only " + unique_num + "");
            Console.WriteLine("  distinct entries.");
            return;
        }

        //
        //  Make sure all W entries are positive.
        //
        for (i = 0; i < point_num; i++)
        {
            switch (w[i])
            {
                case <= 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("LEAST_SET - Fatal error!");
                    Console.WriteLine("  All weights W must be positive,");
                    Console.WriteLine("  but weight " + i + "");
                    Console.WriteLine("  is " + w[i] + "");
                    return;
            }
        }

        double[] s = new double[nterms];
        //
        //  Start inner product summations at zero.
        //
        typeMethods.r8vec_zero(nterms, ref b);
        typeMethods.r8vec_zero(nterms, ref c);
        typeMethods.r8vec_zero(nterms, ref d);
        typeMethods.r8vec_zero(nterms, ref s);
        //
        //  Set the values of P(-1,X) and P(0,X) at all data points.
        //
        double[] pjm1 = new double[point_num];
        double[] pj = new double[point_num];

        typeMethods.r8vec_zero(point_num, ref pjm1);

        for (i = 0; i < point_num; i++)
        {
            pj[i] = 1.0;
        }

        //
        //  Now compute the value of P(J,X(I)) as
        //
        //    P(J,X(I)) = ( X(I) - B(J) ) * P(J-1,X(I)) - C(J) * P(J-2,X(I))
        //
        //  where
        //
        //    S(J) = < P(J,X), P(J,X) >
        //    B(J) = < x*P(J,X), P(J,X) > / < P(J,X), P(J,X) >
        //    C(J) = S(J) / S(J-1)
        //
        //  and the least squares coefficients are
        //
        //    D(J) = < F(X), P(J,X) > / < P(J,X), P(J,X) >
        //
        for (j = 1; j <= nterms; j++)
        {
            int k;
            for (k = 0; k < point_num; k++)
            {
                d[j - 1] += w[k] * f[k] * pj[k];
                b[j - 1] += w[k] * x[k] * pj[k] * pj[k];
                s[j - 1] += w[k] * pj[k] * pj[k];
            }

            d[j - 1] /= s[j - 1];

            if (j == nterms)
            {
                c[j - 1] = 0.0;
                break;
            }

            b[j - 1] /= s[j - 1];

            c[j - 1] = j switch
            {
                1 => 0.0,
                _ => s[j - 1] / s[j - 2]
            };

            for (i = 1; i <= point_num; i++)
            {
                double p = pj[i - 1];
                pj[i - 1] = (x[i - 1] - b[j - 1]) * pj[i - 1] - c[j - 1] * pjm1[i - 1];
                pjm1[i - 1] = p;
            }
        }
    }

    public static double least_val(int nterms, double[] b, double[] c, double[] d,
            double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEAST_VAL evaluates a least squares polynomial defined by LEAST_SET.
        //
        //  Discussion:
        //
        //    The least squares polynomial is assumed to be defined as a sum
        //
        //      P(X) = SUM ( I = 1 to NTERMS ) D(I) * P(I-1,X)
        //
        //    where the orthogonal basis polynomials P(I,X) satisfy the following
        //    three term recurrence:
        //
        //      P(-1,X) = 0
        //      P(0,X) = 1
        //      P(I,X) = ( X - B(I-1) ) * P(I-1,X) - C(I-1) * P(I-2,X)
        //
        //    Therefore, the least squares polynomial can be evaluated as follows:
        //
        //    If NTERMS is 1, then the value of P(X) is D(1) * P(0,X) = D(1).
        //
        //    Otherwise, P(X) is defined as the sum of NTERMS > 1 terms.  We can
        //    reduce the number of terms by 1, because the polynomial P(NTERMS,X)
        //    can be rewritten as a sum of polynomials;  Therefore, P(NTERMS,X)
        //    can be eliminated from the sum, and its coefficient merged in with
        //    those of other polynomials.  Repeat this process for P(NTERMS-1,X)
        //    and so on until a single term remains.
        //    P(NTERMS,X) of P(NTERMS-1,X)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Samuel Conte, Carl deBoor,
        //    Elementary Numerical Analysis,
        //    Second Edition,
        //    McGraw Hill, 1972,
        //    ISBN: 07-012446-4,
        //    LC: QA297.C65.
        //
        //  Parameters:
        //
        //    Input, int NTERMS, the number of terms in the least squares
        //    polynomial.  NTERMS must be at least 1.  The input value of NTERMS
        //    may be reduced from the value given to R8POLY_LS_SET.  This will
        //    evaluate the least squares polynomial of the lower degree specified.
        //
        //    Input, double B[NTERMS], C[NTERMS], D[NTERMS], the information
        //    computed by R8POLY_LS_SET.
        //
        //    Input, double X, the point at which the least squares polynomial
        //    is to be evaluated.
        //
        //    Output, double LEAST_VAL, the value of the least squares 
        //    polynomial at X.
        //
    {
        int i;

        double px = d[nterms - 1];
        double prev = 0.0;

        for (i = nterms - 1; 1 <= i; i--)
        {
            double prev2 = prev;
            prev = px;

            if (i == nterms - 1)
            {
                px = d[i - 1] + (x - b[i - 1]) * prev;
            }
            else
            {
                px = d[i - 1] + (x - b[i - 1]) * prev - c[i] * prev2;
            }
        }

        return px;
    }

    public static void least_val2(int nterms, double[] b, double[] c, double[] d, double x,
            ref double px, ref double pxp)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEAST_VAL2 evaluates a least squares polynomial defined by LEAST_SET.
        //
        //  Discussion:
        //
        //    This routine also computes the derivative of the polynomial.
        //
        //    The least squares polynomial is assumed to be defined as a sum
        //
        //      P(X) = SUM ( I = 1 to NTERMS ) D(I) * P(I-1,X)
        //
        //    where the orthogonal basis polynomials P(I,X) satisfy the following
        //    three term recurrence:
        //
        //      P(-1,X) = 0
        //      P(0,X) = 1
        //      P(I,X) = ( X - B(I-1) ) * P(I-1,X) - C(I-1) * P(I-2,X)
        //
        //    Therefore, the least squares polynomial can be evaluated as follows:
        //
        //    If NTERMS is 1, then the value of P(X) is D(1) * P(0,X) = D(1).
        //
        //    Otherwise, P(X) is defined as the sum of NTERMS > 1 terms.  We can
        //    reduce the number of terms by 1, because the polynomial P(NTERMS,X)
        //    can be rewritten as a sum of polynomials;  Therefore, P(NTERMS,X)
        //    can be eliminated from the sum, and its coefficient merged in with
        //    those of other polynomials.  Repeat this process for P(NTERMS-1,X)
        //    and so on until a single term remains.
        //    P(NTERMS,X) of P(NTERMS-1,X)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NTERMS, the number of terms in the least squares
        //    polynomial.  NTERMS must be at least 1.  The value of NTERMS
        //    may be reduced from the value given to R8POLY_LS_SET.
        //    This will cause R8POLY_LS_VAL to evaluate the least squares polynomial
        //    of the lower degree specified.
        //
        //    Input, double B[NTERMS], C[NTERMS], D[NTERMS], the information
        //    computed by R8POLY_LS_SET.
        //
        //    Input, double X, the point at which the least squares polynomial
        //    is to be evaluated.
        //
        //    Output, double *PX, *PXP, the value and derivative of the least
        //    squares polynomial at X.
        //
    {
        int i;

        px = d[nterms - 1];
        pxp = 0.0;
        double pxm1 = 0.0;
        double pxpm1 = 0.0;

        for (i = nterms - 1; 1 <= i; i--)
        {
            double pxm2 = pxm1;
            double pxpm2 = pxpm1;
            pxm1 = px;
            pxpm1 = pxp;

            if (i == nterms - 1)
            {
                px = d[i - 1] + (x - b[i - 1]) * pxm1;
                pxp = pxm1 + (x - b[i - 1]) * pxpm1;
            }
            else
            {
                px = d[i - 1] + (x - b[i - 1]) * pxm1 - c[i] * pxm2;
                pxp = pxm1 + (x - b[i - 1]) * pxpm1 - c[i] * pxpm2;
            }
        }
    }

    public static void least_set_old(int ntab, double[] xtab, double[] ytab, int ndeg,
            double[] ptab, double[] b, double[] c, double[] d, ref double eps, ref int ierror)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEAST_SET_OLD constructs the least squares polynomial approximation to data.
        //
        //  Discussion:
        //
        //    The least squares polynomial is not returned directly as a simple
        //    polynomial.  Instead, it is represented in terms of a set of
        //    orthogonal polynomials appopriate for the given data.  This makes
        //    the computation more accurate, but means that the user can not
        //    easily evaluate the computed polynomial.  Instead, the routine 
        //    LEAST_EVAL should be used to evaluate the least squares polynomial
        //    at any point.  (However, the value of the least squares polynomial
        //    at each of the data points is returned as part of this computation.)
        //
        //
        //    A discrete unweighted inner product is used, so that
        //
        //      ( F(X), G(X) ) = sum ( 1 <= I <= NTAB ) F(XTAB(I)) * G(XTAB(I)).
        //
        //    The least squares polynomial is determined using a set of
        //    orthogonal polynomials PHI.  These polynomials can be defined
        //    recursively by:
        //
        //      PHI(0)(X) = 1
        //      PHI(1)(X) = X - B(1)
        //      PHI(I)(X) = ( X - B(I) ) * PHI(I-1)(X) - D(I) * PHI(I-2)(X)
        //
        //    The array B(1:NDEG) contains the values
        //
        //      B(I) = ( X*PHI(I-1), PHI(I-1) ) / ( PHI(I-1), PHI(I-1) )
        //
        //    The array D(2:NDEG) contains the values
        //
        //      D(I) = ( PHI(I-1), PHI(I-1) ) / ( PHI(I-2), PHI(I-2) )
        //
        //    Using this basis, the least squares polynomial can be represented as
        //
        //      P(X)(I) = sum ( 0 <= I <= NDEG ) C(I) * PHI(I)(X)
        //
        //    The array C(0:NDEG) contains the values
        //
        //      C(I) = ( YTAB(I), PHI(I) ) / ( PHI(I), PHI(I) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 May 2004
        //
        //  Reference:
        //
        //    Gisela Engeln-Muellges, Frank Uhlig,
        //    Numerical Algorithms with C,
        //    Springer, 1996,
        //    ISBN: 3-540-60530-4.
        //
        //  Parameters:
        //
        //    Input, int NTAB, the number of data points.
        //
        //    Input, double XTAB[NTAB], the X data.  The values in XTAB
        //    should be distinct, and in increasing order.
        //
        //    Input, double YTAB[NTAB], the Y data values corresponding
        //    to the X data in XTAB.
        //
        //    Input, int NDEG, the degree of the polynomial which the
        //    program is to use.  NDEG must be at least 1, and less than or 
        //    equal to NTAB-1.
        //
        //    Output, double PTAB[NTAB], the value of the least squares polynomial 
        //    at the points XTAB(1:NTAB).
        //
        //    Output, double B[1:NDEG], C[0:NDEG], D[2:NDEG], arrays containing 
        //    data about the polynomial.
        //
        //    Output, double *EPS, the root-mean-square discrepancy of the
        //    polynomial fit.
        //
        //    Output, int *IERROR, error flag.
        //    zero, no error occurred;
        //    nonzero, an error occurred, and the polynomial could not be computed.
        //
    {
        const int B_OFFSET = -1;
        const int D_OFFSET = -2;
        int i;

        ierror = 0;
        double[] ztab = new double[2 * ntab];
        switch (ndeg)
        {
            //
            //  Check NDEG.
            //
            case < 1:
                ierror = 1;
                Console.WriteLine("");
                Console.WriteLine("LEAST_SET_OLD - Fatal error!");
                Console.WriteLine("  NDEG < 1.");
                return;
        }

        if (ntab <= ndeg)
        {
            ierror = 1;
            Console.WriteLine("");
            Console.WriteLine("LEAST_SET_OLD - Fatal error!");
            Console.WriteLine("  NTAB <= NDEG.");
            return;
        }

        //
        //  Check that the abscissas are strictly increasing.
        //
        for (i = 1; i <= ntab - 1; i++)
        {
            if (!(xtab[i] <= xtab[i - 1]))
            {
                continue;
            }

            ierror = 1;
            Console.WriteLine("");
            Console.WriteLine("LEAST_SET_OLD - Fatal error!");
            Console.WriteLine("  XTAB must be strictly increasing, but");
            Console.WriteLine("  XTAB(" + (i - 1) + ") = " + xtab[i - 1] + "");
            Console.WriteLine("  XTAB(" + i + ") = " + xtab[i] + "");
            return;
        }

        int i0l1 = 0;
        int i1l1 = ntab;
        //
        //  The polynomial is of degree at least zero.
        //
        double y_sum = 0.0;
        for (i = 0; i < ntab; i++)
        {
            y_sum += ytab[i];
        }

        double rn0 = ntab;
        c[0] = y_sum / ntab;

        for (i = 0; i < ntab; i++)
        {
            ptab[i] = y_sum / ntab;
        }

        switch (ndeg)
        {
            case 0:
            {
                eps = 0.0;
                for (i = 0; i < ntab; i++)
                {
                    eps += Math.Pow(y_sum / ntab - ytab[i], 2);
                }

                eps = Math.Sqrt(eps / ntab);
                return;
            }
        }

        //
        //  The polynomial is of degree at least 1.
        //
        ztab[0] = 0.0;
        for (i = 0; i < ntab; i++)
        {
            ztab[0] += xtab[i];
        }

        b[1 + B_OFFSET] = ztab[0] / ntab;

        double s = 0.0;
        double sum2 = 0.0;
        for (i = 0; i < ntab; i++)
        {
            ztab[i1l1 + i] = xtab[i] - b[1 + B_OFFSET];
            s += ztab[i1l1 + i] * ztab[i1l1 + i];
            sum2 += ztab[i1l1 + i] * (ytab[i] - ptab[i]);
        }

        double rn1 = s;
        c[1] = sum2 / s;

        for (i = 0; i < ntab; i++)
        {
            ptab[i] += c[1] * ztab[i1l1 + i];
        }


        switch (ndeg)
        {
            case 1:
            {
                eps = 0.0;
                for (i = 0; i < ntab; i++)
                {
                    eps += Math.Pow(ptab[i] - ytab[i], 2);
                }

                eps = Math.Sqrt(eps / ntab);
                return;
            }
        }

        for (i = 0; i < ntab; i++)
        {
            ztab[i] = 1.0;
        }

        int mdeg = 2;
        int k = 2;

        for (;;)
        {
            d[k + D_OFFSET] = rn1 / rn0;

            sum2 = 0.0;
            for (i = 0; i < ntab; i++)
            {
                sum2 += xtab[i] * ztab[i1l1 + i] * ztab[i1l1 + i];
            }

            b[k + B_OFFSET] = sum2 / rn1;

            s = 0.0;
            sum2 = 0.0;

            for (i = 0; i < ntab; i++)
            {
                ztab[i0l1 + i] = (xtab[i] - b[k + B_OFFSET]) * ztab[i1l1 + i]
                                 - d[k + D_OFFSET] * ztab[i0l1 + i];
                s += ztab[i0l1 + i] * ztab[i0l1 + i];
                sum2 += ztab[i0l1 + i] * (ytab[i] - ptab[i]);
            }

            rn0 = rn1;
            rn1 = s;

            c[k] = sum2 / rn1;

            (i0l1, i1l1) = (i1l1, i0l1);

            for (i = 0; i < ntab; i++)
            {
                ptab[i] += c[k] * ztab[i1l1 + i];
            }

            if (ndeg <= mdeg)
            {
                break;
            }

            mdeg += 1;
            k += 1;

        }

        //
        //  Compute the RMS error.
        //
        eps = 0.0;
        for (i = 0; i < ntab; i++)
        {
            eps += Math.Pow(ptab[i] - ytab[i], 2);
        }

        eps = Math.Sqrt(eps / ntab);

    }

    public static double least_val_old(double x, int ndeg, double[] b, double[] c, double[] d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEAST_VAL_OLD evaluates a least squares polynomial defined by LEAST_SET_OLD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 December 2004
        //
        //  Reference:
        //
        //    Gisela Engeln-Muellges, Frank Uhlig,
        //    Numerical Algorithms with C,
        //    Springer, 1996,
        //    ISBN: 3-540-60530-4.
        //
        //  Parameters:
        //
        //    Input, double X, the point at which the polynomial is to be evaluated.
        //
        //    Input, int NDEG, the degree of the polynomial fit used.
        //    This is the value of NDEG as returned from LEAST_SET_OLD.
        //
        //    Input, double B[1:NDEG], C[0:NDEG], D[2:NDEG], arrays defined by
        //    LEAST_SET, and needed to evaluate the polynomial.
        //
        //    Output, double LEAST_VALPOLD, the value of the polynomial at X.
        //
    {
        const int B_OFFSET = -1;
        const int D_OFFSET = -2;
        double sk = 0;
        double value;

        switch (ndeg)
        {
            case <= 0:
                value = c[0];
                break;
            case 1:
                value = c[0] + c[1] * (x - b[1 + B_OFFSET]);
                break;
            default:
            {
                double skp2 = c[ndeg];
                double skp1 = c[ndeg - 1] + (x - b[ndeg + B_OFFSET]) * skp2;

                int k;
                for (k = ndeg - 2; 0 <= k; k--)
                {
                    sk = c[k] + (x - b[k + 1 + B_OFFSET]) * skp1 - d[k + 2 + D_OFFSET] * skp2;
                    skp2 = skp1;
                    skp1 = sk;
                }

                value = sk;
                break;
            }
        }

        return value;
    }
}