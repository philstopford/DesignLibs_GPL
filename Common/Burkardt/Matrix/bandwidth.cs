namespace Burkardt.MatrixNS
{
    public static partial class Matrix
    {
        public static void bandwidth ( int m, int n, double[] a, ref int b, ref int l, ref int d, ref int u )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BANDWIDTH returns the bandwidth of a matrix.
        //
        //  Discussion:
        //
        //    If the nonzeros of a matrix only occur in entries that are "close"
        //    to the main diagonal, we say the matrix is banded.
        //
        //    Roughly speaking, the bandwidth B of a matrix is the number of 
        //    diagonals containing nonzeros.  More precisely, it is the minimum number
        //    of contiguous diagonals that contain all the nonzeros.  It is presumed
        //    that the main diagonal is nonzero.
        //
        //    We can also measure U and L, the upper and lower "half-bandwidths" which
        //    count the number of contiguous diagonals above or below the main
        //    diagonal.
        //
        //    We may write
        //      B = L + D + U
        //    where D is presumably 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the matrix.
        //
        //    Output, int &B, the total bandwidth.
        //
        //    Output, int &L, &D, &U, the lower, diagonal, and upper 
        //    bandwidths.
        //
        {
            int i;
            int j;

            l = 0;
            d = 0;
            u = 0;

            for ( i = 0; i < n; i++ )
            {
                j = 0;
                while ( l < i - j )
                {
                    if ( a[i+j*m] != 0.0 )
                    {
                        l = i - j;
                        break;
                    }
                    j = j + 1;
                }

                if ( a[i+i*m] != 0.0 )
                {
                    d = 1;
                }

                j = n - 1;
                while ( u < j - i )
                {
                    if ( a[i+j*m] != 0.0 )
                    {
                        u = j - i;
                        break;
                    }
                    j = j - 1;
                }
            }

            b = l + d + u;

            return;
        }
    }
}