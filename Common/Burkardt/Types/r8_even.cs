using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8vec_even_new(int n, double alo, double ahi)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_EVEN_NEW returns an R8VEC of values evenly spaced between ALO and AHI.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 May 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of values.
            //
            //    Input, double ALO, AHI, the low and high values.
            //
            //    Output, double R8VEC_EVEN_NEW[N], N evenly spaced values.
            //    Normally, A[0] = ALO and A[N-1] = AHI.
            //    However, if N = 1, then A[0] = 0.5*(ALO+AHI).
            //
        {
            double[] a = new double[n];

            if (n == 1)
            {
                a[0] = 0.5 * (alo + ahi);
            }
            else
            {
                for (int i = 0; i < n; i++)
                {
                    a[i] = ((double) (n - i - 1) * alo
                            + (double) (i) * ahi)
                           / (double) (n - 1);
                }
            }

            return a;
        }

        public static double[] r8vec_even(int n, double alo, double ahi)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_EVEN returns N real values, evenly spaced between ALO and AHI.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 February 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of values.
            //
            //    Input, double ALO, AHI, the low and high values.
            //
            //    Output, double R8VEC_EVEN[N], N evenly spaced values.
            //    Normally, A(1) = ALO and A(N) = AHI.
            //    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
            //
        {
            double[] a = new double[n];

            if (n == 1)
            {
                a[0] = 0.5 * (alo + ahi);
            }
            else
            {
                for (int i = 1; i <= n; i++)
                {
                    a[i - 1] = ((double) (n - i) * alo
                                + (double) (i - 1) * ahi)
                               / (double) (n - 1);
                }
            }

            return a;
        }

        public static double r8vec_even_select(int n, double xlo, double xhi, int ival)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / ( N - 1 )
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
            //    Input, int N, the number of values.
            //
            //    Input, double XLO, XHI, the low and high values.
            //
            //    Input, int IVAL, the index of the desired point.
            //    IVAL is normally between 1 and N, but may be any integer value.
            //
            //    Output, double R8VEC_EVEN_SELECT, the IVAL-th of N evenly spaced values
            //    between XLO and XHI.
            //    Unless N = 1, X(1) = XLO and X(N) = XHI.
            //    If N = 1, then X(1) = 0.5*(XLO+XHI).
            //
        {
            double xval;

            if (n == 1)
            {
                xval = 0.5 * (xlo + xhi);
            }
            else
            {
                xval = ((double) (n - ival) * xlo
                        + (double) (ival - 1) * xhi)
                       / (double) (n - 1);
            }

            return xval;
        }
        //****************************************************************************80

        public static void r8vec_even2(int maxval, int[] nfill, int nold, double[] xold,
        ref int nval, ref double[] xval )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_EVEN2 linearly interpolates new numbers into an R8VECa.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The number of values created between two old values can vary from
        //    one pair of values to the next.
        //
        //    The interpolated values are evenly spaced.
        //
        //    This routine is a generalization of R8VEC_EVEN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int MAXVAL, the size of the XVAL array, as declared by the
        //    user.  MAXVAL must be large enough to hold the NVAL values computed by
        //    this routine.  In other words, MAXVAL must be at least equal to
        //    NOLD + SUM (1 <= I <= NOLD-1) NFILL(I).
        //
        //    Input, int NFILL[NOLD-1], the number of values
        //    to be interpolated between XOLD(I) and XOLD(I+1).
        //    NFILL(I) does not count the endpoints.  Thus, if
        //    NFILL(I) is 1, there will be one new point generated
        //    between XOLD(I) and XOLD(I+1).
        //    NFILL(I) must be nonnegative.
        //
        //    Input, int NOLD, the number of values XOLD,
        //    between which extra values are to be interpolated.
        //
        //    Input, double XOLD[NOLD], the original vector of numbers
        //    between which new values are to be interpolated.
        //
        //    Output, int &NVAL, the number of values computed
        //    in the XVAL array.
        //    NVAL = NOLD + SUM ( 1 <= I <= NOLD-1 ) NFILL(I)
        //
        //    Output, double XVAL[MAXVAL].  On output, XVAL contains the
        //    NOLD values of XOLD, as well as the interpolated
        //    values, making a total of NVAL values.
        //
        {
            int i;
            int j;
            int nadd;

            nval = 1;

            for (i = 1; i <= nold - 1; i++)
            {

                if (nfill[i - 1] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8VEC_EVEN2 - Fatal error!");
                    Console.WriteLine("  NFILL[I-1] is negative for I = " + i + "");
                    Console.WriteLine("  NFILL[I-1] = " + nfill[i - 1] + "");
                    return;
                }

                if (maxval < nval + nfill[i - 1] + 1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8VEC_EVEN2 - Fatal error!");
                    Console.WriteLine("  MAXVAL = " + maxval + " is not large enough.");
                    Console.WriteLine("  for the storage for interval I = " + i + "");
                    return;
                }

                nadd = nfill[i - 1] + 2;

                for (j = 1; j <= nadd; j++)
                {
                    xval[nval + j - 2] = ((double) (nadd - j) * xold[i - 1]
                                          + (double) (j - 1) * xold[i])
                                         / (double) (nadd - 1);
                }

                nval = nval + nfill[i - 1] + 1;
            }

            return;
        }

        public static double r8vec_even2_select(int n, double xlo, double xhi, int ival)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_EVEN2_SELECT returns the I-th of N evenly spaced midpoint values.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    This function returns the I-th of N evenly spaced midpoints of N
            //    equal subintervals of [XLO,XHI].
            //
            //    XVAL = ( ( 2 * N - 2 * IVAL + 1 ) * XLO 
            //           + (         2 * IVAL - 1 ) * XHI ) 
            //           / ( 2 * N                )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 July 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of values.
            //
            //    Input, double XLO, XHI, the low and high values.
            //
            //    Input, int IVAL, the index of the desired point.
            //    IVAL is normally between 1 and N, but may be any integer value.
            //
            //    Output, double R8VEC_EVEN2_SELECT, the IVAL-th of N evenly spaced midpoints
            //    between XLO and XHI.
            //
        {
            double xval;

            xval = ((double) (2 * n - 2 * ival + 1) * xlo
                    + (double) (2 * ival - 1) * xhi)
                   / (double) (2 * n);

            return xval;
        }

        public static void r8vec_even3(int nold, int nval, double[] xold, ref double[] xval )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_EVEN3 evenly interpolates new data into an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    This routine accepts a short vector of numbers, and returns a longer
        //    vector of numbers, created by interpolating new values between
        //    the given values.
        //
        //    Between any two original values, new values are evenly interpolated.
        //
        //    Over the whole vector, the new numbers are interpolated in
        //    such a way as to try to minimize the largest distance interval size.
        //
        //    The algorithm employed is not "perfect".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NOLD, the number of values XOLD, between which extra
        //    values are to be interpolated.
        //
        //    Input, int NVAL, the number of values to be computed
        //    in the XVAL array.  NVAL should be at least NOLD.
        //
        //    Input, double XOLD[NOLD], the original vector of numbers
        //    between which new values are to be interpolated.
        //
        //    Output, double XVAL[NVAL].  On output, XVAL contains the
        //    NOLD values of XOLD, as well as interpolated
        //    values, making a total of NVAL values.
        //
        {
            double density;
            int i;
            int ival;
            int j;
            int nmaybe;
            int npts;
            int ntemp;
            int ntot;
            double xlen;
            double xleni;
            double xlentot;

            xlen = 0.0;
            for (i = 1; i <= nold - 1; i++)
            {
                xlen = xlen + Math.Abs(xold[i] - xold[i - 1]);
            }

            ntemp = nval - nold;

            density = (double) (ntemp) / xlen;

            ival = 1;
            ntot = 0;
            xlentot = 0.0;

            for (i = 1; i <= nold - 1; i++)
            {
                xleni = Math.Abs(xold[i] - xold[i - 1]);
                npts = (int) (density * xleni);
                ntot = ntot + npts;
                //
                //  Determine if we have enough left-over density that it should
                //  be changed into a point.  A better algorithm would agonize
                //  more over where that point should go.
                //
                xlentot = xlentot + xleni;
                nmaybe = (int)r8_nint(xlentot * density);

                if (ntot < nmaybe)
                {
                    npts = npts + nmaybe - ntot;
                    ntot = nmaybe;
                }

                for (j = 1; j <= npts + 2; j++)
                {
                    xval[ival + j - 2] = ((double) (npts + 2 - j) * xold[i - 1]
                                          + (double) (j - 1) * xold[i])
                                         / (double) (npts + 2 - 1);
                }

                ival = ival + npts + 1;
            }

            return;
        }
    }
}