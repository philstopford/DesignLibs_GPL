namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8mat_border_add(int m, int n, double[] table)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_BORDER_ADD adds a "border" to an R8MAT.
            //
            //  Discussion:
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
            //    Input, double TABLE[M*N], the table data.
            //
            //    Output, double TABLE2[(M+2)*(N+2)], the augmented table data.
            //
        {
            int i;
            int j;
            double[] table2;

            table2 = new double[(m + 2) * (n + 2)];

            for (j = 0; j < n + 2; j++)
            {
                for (i = 0; i < m + 2; i++)
                {
                    if (i == 0 || i == m + 1 || j == 0 || j == n + 1)
                    {
                        table2[i + j * (m + 2)] = 0.0;
                    }
                    else
                    {
                        table2[i + j * (m + 2)] = table[(i - 1) + (j - 1) * m];
                    }
                }
            }

            return table2;
        }

        public static double[] r8mat_border_cut(int m, int n, double[] table)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_BORDER_CUT cuts the "border" of an R8MAT.
            //
            //  Discussion:
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
            //    Input, double TABLE[M*N], the table data.
            //
            //    Output, double TABLE2[(M-2)*(N-2)], the "interior" table data.
            //
        {
            int i;
            int j;
            double[] table2;

            if (m <= 2 || n <= 2)
            {
                return null;
            }

            table2 = new double[(m - 2) * (n - 2)];

            for (j = 0; j < n - 2; j++)
            {
                for (i = 0; i < m - 2; i++)
                {
                    table2[i + j * (m - 2)] = table[(i + 1) + (j + 1) * m];
                }
            }

            return table2;
        }

    }
}