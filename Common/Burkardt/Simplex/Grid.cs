using Burkardt.Composition;

namespace Burkardt.SimplexNS
{
    public static class Grid
    {
        public static int[] simplex_grid_index_all(int m, int n, int ng)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_GRID_INDEX_ALL returns all the simplex grid indices.
            //
            //  Discussion:
            //
            //    The number of grid indices can be determined by calling 
            //      ng = simplex_grid_size ( m, n )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 2020
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of subintervals.
            //
            //    Input, int NG, the number of values in the grid.
            //
            //    Output, int SIMPLEX_GRID_INDEX_ALL[(M+1)*NG], the current, 
            //    and then the next, grid index.
            //
        {
            int[] g;
            int[] grid;
            int i;
            int k;

            g = new int[m + 1];

            for (i = 0; i < m; i++)
            {
                g[i] = 0;
            }

            g[m] = n;

            grid = new int[(m + 1) * ng];

            k = 0;
            for (i = 0; i <= m; i++)
            {
                grid[i + k * (m + 1)] = g[i];
            }

            for (k = 1; k < ng; k++)
            {
                Comp.comp_next_grlex(m + 1, ref g);

                for (i = 0; i <= m; i++)
                {
                    grid[i + k * (m + 1)] = g[i];
                }
            }

            return grid;
        }

        public static void simplex_grid_index_next(int m, int n, int[] g)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_GRID_INDEX_NEXT returns the next simplex grid index.
            //
            //  Discussion:
            //
            //    The vector G has dimension M+1.  The first M entries may be regarded
            //    as grid coordinates.  These coordinates must have a sum between 0 and N.
            //    The M+1 entry contains the remainder, that is N minus the sum of the
            //    first M coordinates.
            //
            //    Each time the function is called, it is given a current grid index, and
            //    computes the next one.  The very first index is all zero except for a 
            //    final value of N, and the very last index has all zero except for an'
            //    intial value of N.
            //
            //    For example, here are the coordinates in order for M = 3, N = 3:
            //
            //     0  0 0 0 3
            //     1  0 0 1 2
            //     2  0 0 2 1
            //     3  0 0 3 0
            //     4  0 1 0 2
            //     5  0 1 1 1
            //     6  0 1 2 0
            //     7  0 2 0 1
            //     8  0 2 1 0
            //     9  0 3 0 0
            //    10  1 0 0 2
            //    11  1 0 1 1
            //    12  1 0 2 0
            //    13  1 1 0 1
            //    14  1 1 1 0
            //    15  1 2 0 0
            //    16  2 0 0 1
            //    17  2 0 1 0
            //    18  2 1 0 0
            //    19  3 0 0 0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of subintervals.
            //
            //    Input/output, int G[M+1], the current, and then the next,
            //    grid index.
            //
        {
            Comp.comp_next_grlex(m + 1, ref g);
        }

        public static int[] simplex_grid_index_sample(int m, int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_GRID_INDEX_SAMPLE returns a random simplex grid index.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of subintervals in
            //    each dimension.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, int SIMPLEX_GRID_INDEX_SAMPLE[M+1], a randomly selected index 
            //    in the simplex grid.
            //
        {
            int[] g;

            g = Comp.comp_random_new(n, m + 1, ref seed);

            return g;
        }

        public static double[] simplex_grid_index_to_point(int m, int n, int ng, int[] g,
                double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_GRID_INDEX_TO_POINT returns  points corresponding to simplex indices.
            //
            //  Discussion:
            //
            //    The M-dimensional simplex is defined by M+1 vertices.
            //
            //    Given a regular grid that uses N subintervals along the edge between
            //    each pair of vertices, a simplex grid index G is a set of M+1 values
            //    each between 0 and N, and summing to N. 
            //
            //    This function determines the coordinates X of the point corresponding
            //    to the index G.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of subintervals.
            //
            //    Input, int NG, the number of grid indices to be converted.
            //
            //    Input, int G[(M+1)*NG], the grid indices of 1 
            //    or more points.
            //
            //    Input, double V[M*(M+1)], the coordinates of the vertices 
            //    of the simplex.
            //
            //    Output, double SIMPLEX_GRID_INDEX_TO_POINT[M*NG], the coordinates of one 
            //    or more points.
            //
        {
            int i;
            int j;
            int k;
            double[] x;

            x = new double[m * ng];

            for (j = 0; j < ng; j++)
            {
                for (i = 0; i < m; i++)
                {
                    x[i + j * m] = 0.0;
                    for (k = 0; k < m + 1; k++)
                    {
                        x[i + j * m] = x[i + j * m] + v[i + k * m] * (double)(g[k + j * (m + 1)]);
                    }

                    x[i + j * m] = x[i + j * m] / (double)(n);
                }
            }

            return x;
        }

        public static int simplex_grid_size(int m, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_GRID_SIZE counts the grid points inside a simplex.
            //
            //  Discussion:
            //
            //    The size of a grid with parameters M, N is C(M+N,N).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of subintervals.
            //
            //    Output, int SIMPLEX_GRID_SIZE, the number of grid points.
            //
        {
            int i;
            int ng;

            ng = 1;

            for (i = 1; i <= m; i++)
            {
                ng = (ng * (n + i)) / i;
            }

            return ng;
        }
    }
}