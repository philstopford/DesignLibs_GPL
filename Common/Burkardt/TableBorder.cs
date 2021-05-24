namespace Burkardt.Table
{
    public class TableBorder
    {
        static int[] i4mat_border_add ( int m, int n, int[] table )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_BORDER_ADD adds a "border" to an I4MAT.
        //
        //  Discussion:
        //
        //    An I4MAT is an array of I4's.
        //
        //    We suppose the input data gives values of a quantity on nodes
        //    in the interior of a 2D grid, and we wish to create a new table
        //    with additional positions for the nodes that would be on the
        //    border of the 2D grid.
        //
        //                  0 0 0 0 0 0
        //      * * * *     0 * * * * 0
        //      * * * * --> 0 * * * * 0
        //      * * * *     0 * * * * 0
        //                  0 0 0 0 0 0
        //
        //    The illustration suggests the situation in which a 3 by 4 array
        //    is input, and a 5 by 6 array is to be output.
        //
        //    The old data is shifted to its correct positions in the new array.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, int TABLE[M*N], the data.
        //
        //    Output, int TABLE2[(M+2)*(N+2)], the augmented data.
        //
        {
            int[] table2 = new int[(m+2)*(n+2)];

            for (int j = 0; j < n+2; j++ )
            {
                for (int i = 0; i < m+2; i++ )
                {
                    if ( i == 0 || i == m+1 || j == 0 || j == n+1 )
                    {
                        table2[i+j*(m+2)] = 0;
                    }
                    else
                    {
                        table2[i+j*(m+2)] = table[(i-1)+(j-1)*m];
                    }
                }
            }

            return table2;
        }

        
       static  int[] i4mat_border_cut ( int m, int n, int[] table )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_BORDER_CUT cuts the "border" of an I4MAT.
        //
        //  Discussion:
        //
        //    An I4MAT is an array of I4's.
        //
        //    We suppose the input data gives values of a quantity on nodes
        //    on a 2D grid, and we wish to create a new table corresponding only
        //    to those nodes in the interior of the 2D grid.
        //
        //      0 0 0 0 0 0
        //      0 * * * * 0    * * * *
        //      0 * * * * 0 -> * * * *
        //      0 * * * * 0    * * * *
        //      0 0 0 0 0 0
        //
        //    The illustration suggests the situation in which a 5 by 6 array
        //    is input, and a 3 by 4 array is to be output.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, int TABLE[M*N], the data.
        //
        //    Output, int TABLE2[(M-2)*(N-2)], the "interior" data.
        //
        {
            if ( m <= 2 || n <= 2 )
            {
                return new int[0];
            }

            int[] table2 = new int[(m-2)*(n-2)];

            for (int j = 0; j < n-2; j++ )
            {
                for (int i = 0; i < m-2; i++ )
                {
                    table2[i+j*(m-2)] = table[(i+1)+(j+1)*m];
                }
            }

            return table2;
        }

        
    }
}