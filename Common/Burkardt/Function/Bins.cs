using System;
using Burkardt.Types;

namespace Burkardt.Function
{
    public static class Bins
    {

        public static int[] r82_to_bin_even2(int nbin, double[] a, double[] b, double[] c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R82_TO_BIN_EVEN2 determines the appropriate "bin" for an R82 value.
            //
            //  Discussion:
            //
            //    The intervals [A(1),B(1)] and [A(2),B(2)] are each divided into NBIN
            //    equal subintervals or bins.  Boundary bins take care of extreme values.
            //
            //  Example:
            //
            //    NBIN = 5, A(1) = 5,  A(2) = 0,
            //              B(1) = 15, B(2) = 20.
            //
            //   20 +    +    +    +    +    +
            //        15 | 25 | 35 | 45 | 55
            //   16 +----+----+----+----+----+
            //        14 | 24 | 34 | 44 | 54
            //   12 +----+----+----+----+----+
            //        13 | 23 | 33 | 43 | 53
            //    8 +----+----+----+----+----+
            //        12 | 22 | 32 | 42 | 52
            //    4 +----+----+----+----+----+
            //        11 | 21 | 31 | 41 | 51
            //    0 +    +    +    +    +    +
            //      5    7    9   11   13   15
            //
            //      C      BIN
            //   ------  ------
            //    8 -2    2  1
            //    0  1    1  1
            //    6  9    1  3
            //   10 11    3  3
            //   14 23    5  5
            //   25 13    5  4
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NBIN, the number of bins in each dimension.
            //    NBIN must be at least 1.
            //
            //    Input, double A[2], B[2], the lower and upper limits of the bin interval.
            //    While A(I) is expected to be less than B(I), the code should return useful
            //    results if A(I) is actually greater than B(I).
            //
            //    Input, double C[2], a value to be placed in a bin.
            //
            //    Output, int R82_TO_BIN_EVEN2[2], the index of the bin to which C is assigned.
            //
        {
            int[] bin = new int[2];

            for (int i = 0; i < 2; i++)
            {
                bin[i] = r8_to_bin_even2(nbin, a[i], b[i], c[i]);
            }

            return bin;
        }

        public static int[] r82_to_bin_even3(int[] nbin, double[] a, double[] b, double[] c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R82_TO_BIN_EVEN3 determines the appropriate "bin" for an R82 value.
            //
            //  Discussion:
            //
            //    The interval [A(I),B(I)] is divided into NBIN(I) equal subintervals
            //    or bins.
            //
            //  Example:
            //
            //    NBIN = (/ 4, 5 /),
            //
            //      A(1) = 1,  A(2) = 0,
            //      B(1) = 17, B(2) = 20.
            //
            //   20 +    +    +    +    +
            //        15 | 25 | 35 | 45
            //   16 +----+----+----+----+
            //        14 | 24 | 34 | 44
            //   12 +----+----+----+----+
            //        13 | 23 | 33 | 43
            //    8 +----+----+----+----+
            //        12 | 22 | 32 | 42
            //    4 +----+----+----+----+
            //        11 | 21 | 31 | 41
            //    0 +    +    +    +    +
            //      1    5    9   13   17
            //
            //      C      BIN
            //   ------  ------
            //    8 -2    2  1
            //    0  1    1  1
            //    6  9    2  3
            //   10 11    3  3
            //   14 23    4  5
            //   25 13    4  4
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NBIN[2], the number of bins in each dimension.
            //
            //    Input, double A[2], B[2], the lower and upper limits of the bin interval.
            //    While A(I) is expected to be less than B(I), the code should return useful
            //    results if A(I) is actually greater than B(I).
            //
            //    Input, double C[2], a value to be placed in a bin.
            //
            //    Output, int R82_TO_BIN_EVEN3[2], the bin assignment for the value.
            //
        {
            int[] bin = new int[2];

            for (int i = 0; i < 2; i++)
            {
                bin[i] = r8_to_bin_even2(nbin[i], a[i], b[i], c[i]);
            }

            return bin;
        }

        public static int[] r83_to_bin_even2(int nbin, double[] a, double[] b, double[] c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R83_TO_BIN_EVEN2 determines the appropriate "bin" for an R83.
            //
            //  Discussion:
            //
            //    The intervals [A(I),B(I)] are each divided into NBIN
            //    equal subintervals or bins.  Boundary bins take care of extreme values.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NBIN, the number of bins in each dimension.
            //    NBIN must be at least 1.
            //
            //    Input, double A[3], B[3], the lower and upper limits of the bin interval.
            //    While A(I) is expected to be less than B(I), the code should return useful
            //    results if A(I) is actually greater than B(I).
            //
            //    Input, double C[3], a value to be placed in a bin.
            //
            //    Output, int R83_TO_BIN_EVEN2[3], the index of the bin to which C is assigned.
            //
        {
            int[] bin;
            int i;

            bin = new int[3];

            for (i = 0; i < 3; i++)
            {
                bin[i] = r8_to_bin_even2(nbin, a[i], b[i], c[i]);
            }

            return bin;
        }

        public static int[] r83_to_bin_even3(int[] nbin, double[] a, double[] b, double[] c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R83_TO_BIN_EVEN3 determines the appropriate "bin" for an R83.
            //
            //  Discussion:
            //
            //    The interval [A(I),B(I)] is divided into NBIN(I) equal subintervals
            //    or bins.
            //
            //  Example:
            //
            //    NBIN = (/ 4, 5, 2 /),
            //
            //      A(1) = 1,  A(2) = 0,  A(3) = 8
            //      B(1) = 17, B(2) = 20, B(3) = 10
            //
            //
            //            8 < Z < 9                    9 < Z < 10
            //
            //   20 +     +     +     +     +     20 +     +     +     +     +
            //        151 | 251 | 351 | 451            152 | 252 | 352 | 452
            //   16 +-----+-----+-----+-----+     16 +-----+-----+-----+-----+
            //        141 | 241 | 341 | 441            142 | 242 | 342 | 442
            //   12 +-----+-----+-----+-----+     12 +-----+-----+-----+-----+
            //        131 | 231 | 331 | 431            132 | 232 | 332 | 432
            //    8 +-----+-----+-----+-----+      8 +-----+-----+-----+-----+
            //        121 | 221 | 321 | 421            122 | 222 | 322 | 422
            //    4 +-----+-----+-----+-----+      4 +-----+-----+-----+-----+
            //        111 | 211 | 311 | 411            112 | 212 | 312 | 412
            //    0 +     +     +     +     +      0 +     +     +     +     +
            //      1     5     9    13    17        1     5     9    13    17
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NBIN[3], the number of bins in each dimension.
            //
            //    Input, double A[3], B[3], the lower and upper limits of the bin interval.
            //    While A(I) is expected to be less than B(I), the code should return useful
            //    results if A(I) is actually greater than B(I).
            //
            //    Input, double C[3], a value to be placed in a bin.
            //
            //    Output, int R83_TO_BIN_EVEN3[3], the index of the bin to which C is assigned.
            //
        {
            int[] bin;
            int i;

            bin = new int[3];

            for (i = 0; i < 3; i++)
            {
                bin[i] = r8_to_bin_even2(nbin[i], a[i], b[i], c[i]);
            }

            return bin;
        }


        public static int r8_to_bin_even2(int nbin, double a, double b, double c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_TO_BIN_EVEN2 determines the appropriate "bin" for C in [A,B].
            //
            //  Discussion:
            //
            //    The interval from A to B is divided into NBIN equal subintervals or bins.
            //
            //  Example:
            //
            //    NBIN = 5, A = 5, B = 15
            //
            //    <-1-+-2-+-3-+-4-+-5->
            //    5   7   9  11  13  15
            //
            //    C   BIN
            //
            //    1    1
            //    3    1
            //    4.9  1
            //    5    1
            //    6    1
            //    7.1  2
            //    8    2
            //    9.5  3
            //   12    4
            //   14    5
            //   15    5
            //   15.1  5
            //   99    5
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NBIN, the number of bins.
            //
            //    Input, double A, B, the lower and upper limits of the bin interval.
            //    While A is expected to be less than B, the code should return useful
            //    results if A is actually greater than B.
            //
            //    Input, double C, a value to be placed in a bin.
            //
            //    Output, int R_TO_BIN_EVEN_2, the index of the bin to which C is assigned.
            //
        {
            double a2;
            double b2;
            int bin;
            bool reorder;
            //
            //  Take care of special cases.
            //
            if (nbin < 1)
            {
                bin = 0;
                return bin;
            }

            if (nbin == 1)
            {
                bin = 1;
                return bin;
            }

            if (b == a)
            {
                bin = 1;
                return bin;
            }

            //
            //  If the limits are descending, then we switch them now, and
            //  unswitch the results at the end.
            //
            if (a < b)
            {
                reorder = false;
                a2 = a;
                b2 = b;
            }
            else
            {
                reorder = true;
                a2 = b;
                b2 = a;
            }

            //
            //  Compute the bin.
            //
            if (c <= a2)
            {
                bin = 1;
            }
            else if (b2 <= c)
            {
                bin = nbin;
            }
            else
            {
                bin = 1 + (int) (((double) nbin) * (c - a2) / (b2 - a2));
                bin = Math.Max(bin, 1);
                bin = Math.Min(bin, nbin);
            }

            //
            //  Reverse the switching.
            //
            if (reorder)
            {
                bin = nbin + 1 - bin;
            }

            return bin;
        }

        public static void bin_search_one_2d(int[] bin, int nset, double[] pset, int[] nbin,
                int[] bin_start, int[] bin_next, double[] ptest, ref bool found_a_neighbor,
                ref int i_min, ref double d_min_sq, ref int compares)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BIN_SEARCH_ONE_2D searches one cell in a 2D array of bins.
            //
            //  Discussion:
            //
            //    The bins are presumed to have been set up by successive calls to:
            //
            //      R82VEC_BIN_EVEN2,
            //      R82VEC_BINNED_REORDER, and
            //      R82VEC_BINNED_SORT_A.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 October 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int BIN[2], the indices of the cell to be examined.
            //
            //    Input, int NSET, the number of points in the set.
            //
            //    Input, double PSET[2*NSET], the coordinates of the points in the set.
            //
            //    Input, int NBIN[2], the number of cells in the horizontal and
            //    vertical directions.
            //
            //    Input, int BIN_START[NBIN[1]*NBIN[2]], BIN_LAST(NBIN(1),NBIN(2)),
            //    indicates the index of the first and last element in the bin, or -1
            //    if none.
            //
            //    Input, int BIN_NEXT[NSET], the index of the next element of the bin
            //    containing this element.
            //
            //    Input, double PTEST[2], the coordinates of the test point.
            //
            //    Input/output, bool *FOUND_A_NEIGHBOR, is set to TRUE if at least
            //    one point of PSET is found in the current bin.  Otherwise, it retains its
            //    input value.
            //
            //    Input/output, int *I_MIN, the index of the nearest neighbor in
            //    PSET to PTEST, if at least one neighbor has been found.
            //
            //    Input/output, double *D_MIN_SQ, the square of the distance from the nearest
            //    neighbor in PSET to PTEST, if at least one neighbor has been found.
            //
            //    Input/output, int *COMPARES, the number of elements of PSET whose
            //    distance to PTEST has been computed.
            //
        {
            int NDIM = 2;

            double d_sq;
            int i;
            int node;

            node = bin_start[bin[0] + bin[1] * nbin[0]];

            while (0 < node)
            {
                found_a_neighbor = true;

                d_sq = 0.0;
                for (i = 0; i < NDIM; i++)
                {
                    d_sq = d_sq + (ptest[i] - pset[i + node * NDIM]) *
                        (ptest[i] - pset[i + node * NDIM]);
                }

                compares = compares + 1;

                if (d_sq < d_min_sq)
                {
                    d_min_sq = d_sq;
                    i_min = node;
                }

                node = bin_next[node];

            }
        }

        public static void bin_to_r8_even(int nbin, int bin, double a, double b, ref double cmin,
                ref double cmax)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BIN_TO_R8_EVEN returns the limits for a given "bin" in [A,B].
            //
            //  Discussion:
            //
            //    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
            //    An initial bin takes everything less than A, and a final bin takes
            //    everything greater than B.
            //
            //  Example:
            //
            //    NBIN = 7, A = 10, B = 20
            //
            //    BIN      CMIN  CMAX
            //
            //    1         -HUGE 10
            //    2         10    12
            //    3         12    14
            //    4         14    16
            //    5         16    18
            //    6         18    20
            //    7         20    HUGE
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NBIN, the number of bins.  NBIN is normally at least 3.
            //    If NBIN is 1 or 2, then everything is assigned to bin 1.
            //
            //    Input, int BIN, the index of the bin to be considered.
            //    If BIN is less than 1, or greater than NBIN, the user will get what
            //    the user deserves.
            //
            //    Input, double A, B, the lower and upper limits of the bin interval.
            //    While A is expected to be less than B, the code should return useful
            //    results if A is actually greater than B.
            //
            //    Output, double *CMIN, *CMAX, the minimum and maximum limits on the bin.
            //
        {
            //
            //  Take care of special cases.
            //
            if (nbin <= 2)
            {
                cmin = -typeMethods.r8_huge();
                cmax = typeMethods.r8_huge();
                return;
            }

            if (b == a)
            {
                cmin = -typeMethods.r8_huge();
                cmax = typeMethods.r8_huge();
                return;
            }

            //
            //  Compute the bin limits.
            //
            if (bin == 1)
            {
                cmin = -typeMethods.r8_huge();
                cmax = a;
            }
            else if (bin < nbin)
            {
                cmin = ((double) (nbin - bin) * a
                        + (double) (bin - 2) * b)
                       / (double) (nbin - 2);

                cmax = ((double) (nbin - bin - 1) * a
                        + (double) (bin - 1) * b)
                       / (double) (nbin - 2);
            }
            else if (bin == nbin)
            {
                cmin = b;
                cmax = typeMethods.r8_huge();
            }
            else
            {
                cmin = -typeMethods.r8_huge();
                cmax = typeMethods.r8_huge();
            }
        }

        public static void bin_to_r8_even2(int nbin, int bin, double a, double b,
                ref double cmin, ref double cmax)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BIN_TO_R8_EVEN2 returns the limits for a given "bin" in [A,B].
            //
            //  Discussion:
            //
            //    The interval from A to B is divided into NBIN equal subintervals or bins.
            //
            //  Example:
            //
            //    NBIN = 5, A = 10, B = 20
            //
            //    BIN      CMIN  CMAX
            //
            //    1         10    12
            //    2         12    14
            //    3         14    16
            //    4         16    18
            //    5         18    20
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NBIN, the number of bins.
            //
            //    Input, int BIN, the index of the bin to be considered.
            //    If BIN is less than 1, or greater than NBIN, the user will get what
            //    the user deserves.
            //
            //    Input, double A, B, the lower and upper limits of the bin interval.
            //    While A is expected to be less than B, the code should return useful
            //    results if A is actually greater than B.
            //
            //    Output, double *CMIN, *CMAX, the minimum and maximum limits on the bin.
            //
        {
            //
            //  Compute the bin limits.
            //
            if (bin < 1)
            {
                cmin = -typeMethods.r8_huge();
                cmax = a;
            }
            else if (bin <= nbin)
            {
                cmin = ((double) (nbin - bin + 1) * a + (double) (bin - 1) * b)
                       / (double) (nbin);
                cmax = ((double) (nbin - bin) * a + (double) (bin) * b)
                       / (double) (nbin);
            }
            else if (nbin < bin)
            {
                cmin = b;
                cmax = typeMethods.r8_huge();
            }
        }

        public static void bin_to_r82_even(int nbin, int[] bin, double[] a, double[] b,
                ref double[] cmin, ref double[] cmax)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BIN_TO_R82_EVEN returns the limits for a given R82 "bin" in [A,B].
            //
            //  Discussion:
            //
            //    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
            //    An initial bin takes everything less than A, and a final bin takes
            //    everything greater than B.
            //
            //  Example:
            //
            //    NBIN = 7, A(1) = 5, B(1) = 15
            //              A(2) = 0, B(2) = 20
            //
            //     BIN         CMIN      CMAX
            //    ------   -----------  --------
            //    1, 1     -HUGE -HUGE   5     0
            //    2, 2       5     0     7     4
            //    3, 3       7     4     9     8
            //    4, 4       9     8    11    12
            //    5, 5      11    12    13    16
            //    6, 6      13    16    15    20
            //    7, 7      15    20    HUGE HUGE
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 February 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NBIN, the number of bins.  NBIN is normally
            //    at least 3.  If NBIN is 1 or 2, then everything is assigned to bin 1.
            //
            //    Input, int BIN[2], the index of the bin to be considered.
            //    If BIN(I) is less than 1, or greater than NBIN, the user will get what
            //    the user deserves.
            //
            //    Input, double A[2], B[2], the lower and upper limits of the
            //    bin interval.  While A(I) is expected to be less than B(I), the code
            //    should return useful results if A(I) is actually greater than B(I).
            //
            //    Output, double CMIN[2], CMAX[2], the minimum and maximum
            //    limits on the bin.
            //
        {
            double cmax_1d = 0;
            double cmin_1d = 0;
            int i;

            for (i = 0; i < 2; i++)
            {
                bin_to_r8_even(nbin, bin[i], a[i], b[i], ref cmin_1d, ref cmax_1d);
                cmin[i] = cmin_1d;
                cmax[i] = cmax_1d;
            }

            return;
        }

        public static void bin_to_r82_even2(int nbin, int[] bin, double[] a, double[] b,
                ref double[] cmin, ref double[] cmax)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BIN_TO_R82_EVEN2 returns the limits for a given R82 "bin" in [A,B].
            //
            //
            //  Discussion:
            //
            //    The interval from A to B is divided into NBIN equal subintervals or bins.
            //
            //  Example:
            //
            //    NBIN = 5, A(1) = 5, B(1) = 15
            //              A[2] = 0, B[2] = 20
            //
            //     BIN         CMIN      CMAX
            //    ------   -----------  --------
            //    1, 1       5     0     7     4
            //    2, 2       7     4     9     8
            //    3, 3       9     8    11    12
            //    4, 4      11    12    13    16
            //    5, 5      13    16    15    20
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NBIN, the number of bins.
            //
            //    Input, int BIN[2], the index of the bin to be considered.
            //
            //    Input, double A[2], B[2], the lower and upper limits of the bin interval.
            //    While A(I) is expected to be less than B(I), the code should return useful
            //    results if A(I) is actually greater than B(I).
            //
            //    Output, double CMIN[2], CMAX[2], the minimum and maximum limits on the bin.
            //
        {
            int i;
            int ndim = 2;

            for (i = 0; i < ndim; i++)
            {
                bin_to_r8_even2(nbin, bin[i], a[i], b[i], ref cmin[i], ref cmax[i]);
            }
        }

        public static void bin_to_r82_even3(int[] nbin, int[] bin, double[] a, double[] b,
                ref double[] cmin, ref double[] cmax)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BIN_TO_R82_EVEN3 returns the limits for a given R82 "bin" in [A,B].
            //
            //  Discussion:
            //
            //    The interval from A(I) to B(I) is divided into NBIN(I) equal
            //    subintervals or bins.
            //
            //  Example:
            //
            //    NBIN = (/ 4, 5, /)
            //
            //    A(1) = 5, B(1) = 15
            //    A[2] = 0, B[2] = 20
            //
            //     BIN         CMIN      CMAX
            //    ------   -----------  --------
            //    1, 1       5     0     7     4
            //    2, 2       7     4     9     8
            //    3, 3       9     8    11    12
            //    4, 4      11    12    13    16
            //    5, 5      13    16    15    20
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NBIN[2], the number of bins in each dimension.
            //
            //    Input, int BIN[2], the index of the bin to be considered.
            //
            //    Input, double A[2], B[2], the lower and upper limits of the bin interval.
            //    While A(I) is expected to be less than B(I), the code should return useful
            //    results if A(I) is actually greater than B(I).
            //
            //    Output, double CMIN[2], CMAX[2], the minimum and maximum limits on the bin.
            //
        {
            int NDIM = 2;

            int i;

            for (i = 0; i < NDIM; i++)
            {
                bin_to_r8_even2(nbin[i], bin[i], a[i], b[i], ref cmin[i], ref cmax[i]);
            }
        }

        public static void bin_to_r83_even2(int nbin, int[] bin, double[] a, double[] b,
                ref double[] cmin, ref double[] cmax)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BIN_TO_R83_EVEN2 returns the limits for a given R83 "bin" in [A,B].
            //
            //  Discussion:
            //
            //    The interval from A to B is divided into NBIN equal subintervals or bins.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NBIN, the number of bins.
            //
            //    Input, int BIN[3], the index of the bin to be considered.
            //
            //    Input, double A[3], B[3], the lower and upper limits of the bin interval.
            //    While A(I) is expected to be less than B(I), the code should return useful
            //    results if A(I) is actually greater than B(I).
            //
            //    Output, double CMIN[3], CMAX[3], the minimum and maximum limits on the bin.
            //
        {
            int NDIM = 3;

            int i;

            for (i = 0; i < NDIM; i++)
            {
                bin_to_r8_even2(nbin, bin[i], a[i], b[i], ref cmin[i], ref cmax[i]);
            }
        }

        public static void bin_to_r83_even3(int[] nbin, int[] bin, double[] a, double[] b,
                ref double[] cmin, ref double[] cmax)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BIN_TO_R83_EVEN3 returns the limits for a given R83 "bin" in [A,B].
            //
            //  Discussion:
            //
            //    The interval from A(I) to B(I) is divided into NBIN(I) equal
            //    subintervals or bins.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NBIN[3], the number of bins in each dimension.
            //
            //    Input, int BIN[3], the index of the bin to be considered.
            //
            //    Input, double A[3], B[3], the lower and upper limits of the bin interval.
            //    While A(I) is expected to be less than B(I), the code should return useful
            //    results if A(I) is actually greater than B(I).
            //
            //    Output, double CMIN[3], CMAX[3], the minimum and maximum limits on the bin.
            //
        {
            int NDIM = 3;

            int i;

            for (i = 0; i < NDIM; i++)
            {
                bin_to_r8_even2(nbin[i], bin[i], a[i], b[i], ref cmin[i], ref cmax[i]);
            }
        }

        public static void index_box2_next_2d(int n1, int n2, int ic, int jc, ref int i, ref int j,
                ref bool more)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_BOX2_NEXT_2D produces indices on the surface of a box in 2D.
            //
            //  Discussion:
            //
            //    The box has center at (IC,JC), and has half-widths N1 and N2.
            //    The indices are exactly those which are between (IC-N1,JC-N2) and
            //    (IC+N1,JC+N2) with the property that at least one of I and J
            //    is an "extreme" value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N1, N2, the half-widths of the box, that is, the
            //    maximum distance allowed between (IC,JC) and (I,J).
            //
            //    Input, int IC, JC, the central cell of the box.
            //
            //    Input/output, int *I, *J.  On input, the previous index set.
            //    On output, the next index set.  On the first call, MORE should
            //    be set to FALSE, and the input values of I and J are ignored.
            //
            //    Input/output, bool *MORE.
            //    On the first call for a given box, the user should set MORE to FALSE.
            //    On return, the routine sets MORE to TRUE.
            //    When there are no more indices, the routine sets MORE to FALSE.
            //
        {
            if (!(more))
            {
                more = true;
                i = ic - n1;
                j = jc - n2;
                return;
            }

            if (i == ic + n1 &&
                j == jc + n2)
            {
                more = false;
                return;
            }

            //
            //  Increment J.
            //
            j = j + 1;
            //
            //  Check J.
            //
            if (jc + n2 < j)
            {
                j = jc - n2;
                i = i + 1;
            }
            else if (j < jc + n2 && (i == ic - n1 || i == ic + n1))
            {
                return;
            }
            else
            {
                j = jc + n2;
                return;
            }

            return;
        }

        public static void index_box2_next_3d(int n1, int n2, int n3, int ic, int jc, int kc,
                ref int i, ref int j, ref int k, ref bool more)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_BOX2_NEXT_3D produces indices on the surface of a box in 3D.
            //
            //  Discussion:
            //
            //    The box has a central cell of (IC,JC,KC), with a half widths of
            //    (N1,N2,N3).  The indices are exactly those between (IC-N1,JC-N2,KC-N3) and
            //    (IC+N1,JC+N2,KC+N3) with the property that at least one of I, J, and K
            //    is an "extreme" value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N1, N2, N3, the "half widths" of the box, that is, the
            //    maximum distances from the central cell allowed for I, J and K.
            //
            //    Input, int IC, JC, KC, the central cell of the box.
            //
            //    Input/output, int *I, *J, *K.  On input, the previous index set.
            //    On output, the next index set.  On the first call, MORE should
            //    be set to FALSE, and the input values of I, J, and K are ignored.
            //
            //    Input/output, bool *MORE.
            //    On the first call for a given box, the user should set MORE to FALSE.
            //    On return, the routine sets MORE to TRUE.
            //    When there are no more indices, the routine sets MORE to FALSE.
            //
        {
            if (!(more))
            {
                more = true;
                i = ic - n1;
                j = jc - n2;
                k = kc - n3;
                return;
            }

            if (i == ic + n1 &&
                j == jc + n2 &&
                k == kc + n3)
            {
                more = false;
                return;
            }

            //
            //  Increment K.
            //
            k = k + 1;
            //
            //  Check K.
            //
            if (kc + n3 < k)
            {
                k = kc - n3;
                j = j + 1;
            }
            else if (k < kc + n3 &&
                     (i == ic - n1 ||
                      i == ic + n1 ||
                      j == jc - n2 ||
                      j == jc + n2))
            {
                return;
            }
            else
            {
                k = kc + n3;
                return;
            }

            //
            //  Check J.
            //
            if (jc + n2 < j)
            {
                j = jc - n2;
                i = i + 1;
            }
            else if (j < jc + n2 &&
                     (i == ic - n1 ||
                      i == ic + n1 ||
                      k == kc - n3 ||
                      k == kc + n3))
            {
                return;
            }
            else
            {
                j = jc + n2;
                return;
            }

            return;
        }

        public static int lrline(double xu, double yu, double xv1, double yv1, double xv2,
                double yv2, double dv)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LRLINE determines where a point lies in relation to a directed line.
            //
            //  Discussion:
            //
            //    LRLINE determines whether a point is to the left of, right of,
            //    or on a directed line parallel to a line through given points.
            //
            //  Modified:
            //
            //    22 June 2009
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Barry Joe.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Barry Joe,
            //    GEOMPACK - a software package for the generation of meshes
            //    using geometric algorithms,
            //    Advances in Engineering Software,
            //    Volume 13, pages 325-331, 1991.
            //
            //  Parameters:
            //
            //    Input, double XU, YU, XV1, YV1, XV2, YV2, are vertex coordinates; the
            //    directed line is parallel to and at signed distance DV to the left of
            //    the directed line from (XV1,YV1) to (XV2,YV2); (XU,YU) is the vertex for
            //    which the position relative to the directed line is to be determined.
            //
            //    Input, double DV, the signed distance, positive for left.
            //
            //    Output, int LRLINE, is +1, 0, or -1 depending on whether (XU,YU) is
            //    to the right of, on, or left of the directed line.  LRLINE is 0 if
            //    the line degenerates to a point.
            //
        {
            double dx;
            double dxu;
            double dy;
            double dyu;
            double t;
            double tol;
            double tolabs;
            int value = 0;

            tol = 100.0 * double.Epsilon;

            dx = xv2 - xv1;
            dy = yv2 - yv1;
            dxu = xu - xv1;
            dyu = yu - yv1;

            tolabs = tol * Math.Max(Math.Abs(dx),
                Math.Max(Math.Abs(dy),
                    Math.Max(Math.Abs(dxu),
                        Math.Max(Math.Abs(dyu), Math.Abs(dv)))));

            t = dy * dxu - dx * dyu + dv * Math.Sqrt(dx * dx + dy * dy);

            if (tolabs < t)
            {
                value = 1;
            }
            else if (-tolabs <= t)
            {
                value = 0;
            }
            else if (t < -tolabs)
            {
                value = -1;
            }

            return value;
        }

        public static bool perm_check(int n, int[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_CHECK checks that a vector represents a permutation.
            //
            //  Discussion:
            //
            //    The routine verifies that each of the integers from 1
            //    to N occurs among the N entries of the permutation.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 January 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries.
            //
            //    Input, int P[N], the array to check.
            //
            //    Output, bool PERM_CHECK, is TRUE if the permutation is OK.
            //
        {
            bool found;
            int i;
            int seek;

            for (seek = 1; seek <= n; seek++)
            {
                found = false;

                for (i = 0; i < n; i++)
                {
                    if (p[i] == seek)
                    {
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    return false;
                }

            }

            return true;
        }

        public static void perm_inv(int n, ref int[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_INV inverts a permutation "in place".
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 January 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects being permuted.
            //
            //    Input/output, int P[N], the permutation, in standard index form.
            //    On output, P describes the inverse permutation
            //
        {
            int i;
            int i0;
            int i1;
            int i2;
            int is_;

            if (n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_INV - Fatal error!");
                Console.WriteLine("  Input value of N = " + n + "");
                return;

            }

            if (!perm_check(n, p))
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_INV - Fatal error!");
                Console.WriteLine("  The input array does not represent");
                Console.WriteLine("  a proper permutation.");
                return;
            }

            is_ = 1;

            for (i = 1; i <= n; i++)
            {
                i1 = p[i - 1];

                while (i < i1)
                {
                    i2 = p[i1 - 1];
                    p[i1 - 1] = -i2;
                    i1 = i2;
                }

                is_ = -typeMethods.i4_sign(p[i - 1]);
                p[i - 1] = typeMethods.i4_sign(is_) * Math.Abs(p[i - 1]);
            }

            for (i = 1; i <= n; i++)
            {
                i1 = -p[i - 1];

                if (0 <= i1)
                {
                    i0 = i;

                    for (;;)
                    {
                        i2 = p[i1 - 1];
                        p[i1 - 1] = i0;

                        if (i2 < 0)
                        {
                            break;
                        }

                        i0 = i1;
                        i1 = i2;
                    }
                }
            }

            return;
        }

        public static int points_nearest_point_naive_2d(int nset, double[] pset, double[] ptest,
                ref double d_min)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POINTS_NEAREST_POINT_NAIVE_2D finds the nearest point to a given point in 2D.
            //
            //  Discussion:
            //
            //    A naive algorithm is used.  The distance to every point is calculated,
            //    in order to determine the smallest.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 October 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NSET, the number of points in the set.
            //
            //    Input, double PSET[2*NSET], the coordinates of the points in the set.
            //
            //    Input, double PTEST[2], the point whose nearest neighbor is sought.
            //
            //    Output, double *D_MIN, the distance between P and PSET(*,I_MIN).
            //
            //    Output, int POINTS_NEAREST_POINT_NAIVE_2D, the index of the nearest
            //    point in PSET to P.
            //
        {
            int NDIM = 2;

            double d;
            int i;
            int j;
            int p_min;

            d_min = typeMethods.r8_huge();
            p_min = 0;

            for (j = 0; j < nset; j++)
            {
                d = 0.0;
                for (i = 0; i < NDIM; i++)
                {
                    d = d + (ptest[i] - pset[i + j * NDIM]) * (ptest[i] - pset[i + j * NDIM]);
                }

                if (d < d_min)
                {
                    d_min = d;
                    p_min = j;
                }
            }

            d_min = Math.Sqrt(d_min);

            return p_min;
        }

        public static int points_nearest_point_naive_3d(int nset, double[] pset, double[] ptest,
                ref double d_min)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POINTS_NEAREST_POINT_NAIVE_3D finds the nearest point to a given point in 3D.
            //
            //  Discussion:
            //
            //    A naive algorithm is used.  The distance to every point is calculated,
            //    in order to determine the smallest.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 October 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NSET, the number of points in the set.
            //
            //    Input, double PSET[3*NSET], the coordinates of the points in the set.
            //
            //    Input, double PTEST[3], the point whose nearest neighbor is sought.
            //
            //    Output, double *D_MIN, the distance between P and PSET(*,I_MIN).
            //
            //    Output, int POINTS_NEAREST_POINT_NAIVE_3D, the index of the nearest
            //    point in PSET to P.
            //
        {
            int NDIM = 3;

            double d;
            int i;
            int j;
            int p_min;

            d_min = typeMethods.r8_huge();
            p_min = 0;

            for (j = 0; j < nset; j++)
            {
                d = 0.0;
                for (i = 0; i < NDIM; i++)
                {
                    d = d + (ptest[i] - pset[i + j * NDIM]) * (ptest[i] - pset[i + j * NDIM]);
                }

                if (d < d_min)
                {
                    d_min = d;
                    p_min = j;
                }
            }

            d_min = Math.Sqrt(d_min);

            return p_min;

        }

        public static int points_nearest_point_naive_nd(int ndim, int nset, double[] pset,
                double[] ptest, ref double d_min)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POINTS_NEAREST_POINT_NAIVE_ND finds the nearest point to a given point in ND.
            //
            //  Discussion:
            //
            //    A naive algorithm is used.  The distance to every point is calculated,
            //    in order to determine the smallest.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 October 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NSET, the number of points in the set.
            //
            //    Input, double PSET[NDIM*NSET], the coordinates of the points in the set.
            //
            //    Input, double PTEST[NDIM], the point whose nearest neighbor is sought.
            //
            //    Output, double *D_MIN, the distance between P and PSET(*,I_MIN).
            //
            //    Output, int POINTS_NEAREST_POINT_NAIVE_ND, the index of the nearest
            //    point in PSET to P.
            //
        {
            double d;
            int i;
            int j;
            int p_min;

            d_min = typeMethods.r8_huge();
            p_min = 0;

            for (j = 0; j < nset; j++)
            {
                d = 0.0;
                for (i = 0; i < ndim; i++)
                {
                    d = d + (ptest[i] - pset[i + j * ndim]) * (ptest[i] - pset[i + j * ndim]);
                }

                if (d < d_min)
                {
                    d_min = d;
                    p_min = j;
                }
            }

            d_min = Math.Sqrt(d_min);

            return p_min;
        }

        public static int[] points_nearest_points_naive_2d(int nset, double[] pset, int ntest,
                double[] ptest)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POINTS_NEAREST_POINTS_NAIVE_2D finds the nearest point to given points in 2D.
            //
            //  Discussion:
            //
            //    A naive algorithm is used.  The distance to every point is calculated,
            //    in order to determine the smallest.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 November 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NSET, the number of points in the set.
            //
            //    Input, double PSET[2*NSET], the coordinates of the points in the set.
            //
            //    Input, int NTEST, the number of test points.
            //
            //    Input, double PTEST[2*NTEST], the coordinates of the test points.
            //
            //    Output, int POINTS_NEAREST_POINTS_NAIVE_2D[NTEST], the index of the
            //    nearest point in PSET to each point in PTEST.
            //
        {
            int NDIM = 2;

            double d;
            double d_min;
            int i;
            int[] nearest;
            int set;
            int test;

            nearest = new int[ntest];

            for (test = 0; test < ntest; test++)
            {
                d_min = typeMethods.r8_huge();
                nearest[test] = -1;

                for (set = 0; set < nset; set++)
                {
                    d = 0.0;
                    for (i = 0; i < NDIM; i++)
                    {
                        d = d + (ptest[i + test * NDIM] - pset[i + set * NDIM])
                            * (ptest[i + test * NDIM] - pset[i + set * NDIM]);
                    }

                    if (d < d_min)
                    {
                        d_min = d;
                        nearest[test] = set;
                    }
                }
            }

            return nearest;
        }

        public static int[] points_nearest_points_naive_3d(int nset, double[] pset, int ntest,
                double[] ptest)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POINTS_NEAREST_POINTS_NAIVE_3D finds the nearest point to given points in 3D.
            //
            //  Discussion:
            //
            //    A naive algorithm is used.  The distance to every point is calculated,
            //    in order to determine the smallest.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 November 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NSET, the number of points in the set.
            //
            //    Input, double PSET[3*NSET], the coordinates of the points in the set.
            //
            //    Input, int NTEST, the number of test points.
            //
            //    Input, double PTEST[3*NTEST], the coordinates of the test points.
            //
            //    Output, int POINTS_NEAREST_POINTS_NAIVE_3D[NTEST], the index of the
            //    nearest point in PSET to each point in PTEST.
            //
        {
            int NDIM = 3;

            double d;
            double d_min;
            int i;
            int[] nearest;
            int set;
            int test;

            nearest = new int[ntest];

            for (test = 0; test < ntest; test++)
            {
                d_min = typeMethods.r8_huge();
                nearest[test] = -1;

                for (set = 0; set < nset; set++)
                {
                    d = 0.0;
                    for (i = 0; i < NDIM; i++)
                    {
                        d = d + (ptest[i + test * NDIM] - pset[i + set * NDIM])
                            * (ptest[i + test * NDIM] - pset[i + set * NDIM]);
                    }

                    if (d < d_min)
                    {
                        d_min = d;
                        nearest[test] = set;
                    }
                }
            }

            return nearest;
        }



    }
}