using System;
using Burkardt.Types;

namespace Burkardt.Grid
{
    public static class Multigrid
    {
        public static int[] multigrid_index0(int dim_num, int[] order_1d, int order_nd )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTIGRID_INDEX0 returns an indexed multidimensional grid.
        //
        //  Discussion:
        //
        //    For dimension DIM, the second index of INDX may vary from 
        //    0 to ORDER_1D[DIM]-1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 November 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int ORDER_1D[DIM_NUM], the order of the
        //    rule in each dimension.
        //
        //    Input, int ORDER_ND, the product of the entries of ORDER_1D.
        //
        //    Output, int INDX[DIM_NUM*ORDER_ND], the indices of the points in
        //    the grid.  The second dimension of this array is equal to the
        //    product of the entries of ORDER_1D.
        //
        {
            int[] a;
            int dim;
            bool more;
            int p;
            int[] indx;

            indx = new int[dim_num * order_nd];
            a = new int[dim_num];
            more = false;
            p = 0;

            for (;;)
            {
                typeMethods.vec_colex_next2(dim_num, order_1d, ref a, ref more);

                if (!more)
                {
                    break;
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    indx[dim + p * dim_num] = a[dim];
                }

                p = p + 1;
            }
            
            return indx;
        }

        public static int[] multigrid_index_z(int dim_num, int[] order_1d, int order_nd)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MULTIGRID_INDEX_Z returns an indexed multidimensional grid.
            //
            //  Discussion:
            //
            //    For dimension DIM, the number of points is ORDER_1D[DIM].
            //
            //    We assume that ORDER_1D[DIM] is an odd number,
            //      ORDER_1D[DIM] = N = 2 * M + 1
            //    so that the points have coordinates
            //      -M/M, -(M-1)/M, ..., -1/M, 0/M, 1/M, 2/M, 3/M, ..., (M-1)/M, M/M.
            //    and we index them as
            //      -M,   -(M-1),        -1,   0,   1,   2,   3,   ...,  M-1,    M.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 September 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Fabio Nobile, Raul Tempone, Clayton Webster,
            //    A Sparse Grid Stochastic Collocation Method for Partial Differential
            //    Equations with Random Input Data,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 46, Number 5, 2008, pages 2309-2345.
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension of the points.
            //
            //    Input, int ORDER_1D[DIM_NUM], the order of the
            //    rule in each dimension.
            //
            //    Input, int ORDER_ND, the product of the entries of ORDER_1D.
            //
            //    Output, int INDX[DIM_NUM*ORDER_ND], the indices of the points in
            //    the grid.  The second dimension of this array is equal to the
            //    product of the entries of ORDER_1D.
            //
        {
            int[] a;
            int dim;
            bool more;
            int p;
            int[] indx;

            indx = new int[dim_num * order_nd];
            a = new int[dim_num];
            more = false;
            p = 0;

            for (;;)
            {
                typeMethods.vec_colex_next2(dim_num, order_1d, ref a, ref more);

                if (!more)
                {
                    break;
                }

                //
                //  The values of A(DIM) are between 0 and ORDER_1D(DIM)-1 = N - 1 = 2 * M.
                //  Subtracting M sets the range to -M to +M, as we wish.
                //
                for (dim = 0; dim < dim_num; dim++)
                {
                    indx[dim + p * dim_num] = a[dim] - (order_1d[dim] - 1) / 2;
                }

                p = p + 1;
            }

            return indx;
        }

        public static void multigrid_scale_closed(int dim_num, int order_nd, int level_max,
            int[] level_1d, int[] grid_index )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTIGRID_SCALE_CLOSED renumbers a grid as a subgrid on a higher level.
        //
        //  Discussion:
        //
        //    This routine takes a grid associated with a given value of
        //    LEVEL, and multiplies all the indices by a power of 2, so that
        //    the indices reflect the position of the same points, but in
        //    a grid of level LEVEL_MAX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 November 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int ORDER_ND, the number of points in the grid.
        //
        //    Input, int LEVEL_MAX, the maximum value of LEVEL.
        //
        //    Input, int LEVEL_1D[DIM_NUM], the level in each dimension.
        //
        //    Input/output, int GRID_INDEX[DIM_NUM*POINT_NUM], the index
        //    values for each grid point.  On input, these indices are based in
        //    the level for which the grid was generated; on output, the
        //    indices are appropriate for the grid as a subgrid of a grid
        //    of level LEVEL_MAX.
        //
        {
            int dim;
            int factor;
            int order;
            int order_max;

            for (dim = 0; dim < dim_num; dim++)
            {
                if (level_1d[dim] == 0)
                {
                    if (0 == level_max)
                    {
                        order_max = 1;
                    }
                    else
                    {
                        order_max = (int)Math.Pow(2, level_max) + 1;
                    }

                    for (order = 0; order < order_nd; order++)
                    {
                        grid_index[dim + order * dim_num] = (order_max - 1) / 2;
                    }
                }
                else
                {
                    factor = (int)Math.Pow(2, level_max - level_1d[dim]);
                    for (order = 0; order < order_nd; order++)
                    {
                        grid_index[dim + order * dim_num] = grid_index[dim + order * dim_num] * factor;
                    }
                }
            }
        }
    }
}