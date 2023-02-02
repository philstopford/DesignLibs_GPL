using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8vec_bin(int n, double[] x, int bin_num, double bin_min, double bin_max,
            ref int[] bin, ref double[] bin_limit )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_BIN computes bins based on a given R8VEC.
        //
        //  Discussion:
        //
        //    The user specifies minimum and maximum bin values, BIN_MIN and
        //    BIN_MAX, and the number of bins, BIN_NUM.  This determines a
        //    "bin width":
        //
        //      H = ( BIN_MAX - BIN_MIN ) / BIN_NUM
        //
        //    so that bin I will count all entries X(J) such that
        //
        //      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
        //
        //    The array X does NOT have to be sorted.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries of X.
        //
        //    Input, double X[N], an (unsorted) array to be binned.
        //
        //    Input, int BIN_NUM, the number of bins.  Two extra bins,
        //    #0 and #BIN_NUM+1, count extreme values.
        //
        //    Input, double BIN_MIN, BIN_MAX, define the range and size
        //    of the bins.  BIN_MIN and BIN_MAX must be distinct.
        //    Normally, BIN_MIN < BIN_MAX, and the documentation will assume
        //    this, but proper results will be computed if BIN_MIN > BIN_MAX.
        //
        //    Output, int BIN[BIN_NUM+2].
        //    BIN(0) counts entries of X less than BIN_MIN.
        //    BIN(BIN_NUM+1) counts entries greater than or equal to BIN_MAX.
        //    For 1 <= I <= BIN_NUM, BIN(I) counts the entries X(J) such that
        //      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
        //    where H is the bin spacing.
        //
        //    Output, double BIN_LIMIT[BIN_NUM+1], the "limits" of the bins.
        //    BIN(I) counts the number of entries X(J) such that
        //      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
        //
    {
        int i;

        if (Math.Abs(bin_max - bin_min) <= typeMethods.r8_epsilon())
        {
            Console.WriteLine("");
            Console.WriteLine("R8VEC_BIN - Fatal error!");
            Console.WriteLine("  BIN_MIN = BIN_MAX = " + bin_max + ".");
            return;
        }

        for (i = 0; i <= bin_num + 1; i++)
        {
            bin[i] = 0;
        }

        for (i = 0; i < n; i++)
        {
            double t = (x[i] - bin_min) / (bin_max - bin_min);

            int j = t switch
            {
                < 0.0 => 0,
                >= 1.0 => bin_num + 1,
                _ => 1 + (int) (bin_num * t)
            };

            bin[j] += 1;
        }

        //
        //  Compute the bin limits.
        //
        for (i = 0; i <= bin_num; i++)
        {
            bin_limit[i] = ((bin_num - i) * bin_min
                            + i * bin_max)
                           / bin_num;
        }
    }
}