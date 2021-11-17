namespace Burkardt.Types;

public static partial class typeMethods
{
    public static int[] r8vec_histogram ( int n, double[] a, double a_lo, double a_hi,
            int histo_num )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_HISTOGRAM histograms an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    Values between A_LO and A_HI will be histogrammed into the bins
        //    1 through HISTO_NUM.  Values below A_LO are counted in bin 0,
        //    and values greater than A_HI are counted in bin HISTO_NUM+1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Input, double A[N], the array to examine.
        //
        //    Input, double A_LO, A_HI, the lowest and highest
        //    values to be histogrammed.  These values will also define the bins.
        //
        //    Input, int HISTO_NUM, the number of bins to use.
        //
        //    Output, int HISTO_GRAM[HISTO_NUM+2], contains the number of
        //    entries of A in each bin.
        //
    {
        int i;

        int[] histo_gram = new int[histo_num+2];

        i4vec_zeros ( histo_num+2, ref histo_gram );

        double delta = ( a_hi - a_lo ) / (2 * histo_num);

        for ( i = 0; i < n; i++ )
        {
            if ( a[i] < a_lo )
            {
                histo_gram[0] += 1;
            }
            else if ( a[i] <= a_hi )
            {
                int j = (int)r8_nint (
                    ( ( a_hi -       delta - a[i]        ) * 1
                      + (      -       delta + a[i] - a_lo ) * histo_num )
                    / ( a_hi - 2.0 * delta        - a_lo ) );

                histo_gram[j] += 1;
            }
            else if ( a_hi < a[i] )
            {
                histo_gram[histo_num+1] += 1;
            }
        }

        return histo_gram;
    }
}