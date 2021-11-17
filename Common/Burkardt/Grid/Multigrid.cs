using System;
using Burkardt.Types;

namespace Burkardt.Grid;

public static class Multigrid
{
    public static int[] multigrid_index0(int dim_num, int[] order_1d, int order_nd)

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

            p += 1;
        }

        return indx;
    }
        
    public static int[] multigrid_index1 ( int dim_num, int[] order_1d, int order_nd )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTIGRID_INDEX1 returns an indexed multidimensional grid.
        //
        //  Discussion:
        //
        //    For dimension DIM, the second index of INDX may vary from 
        //    1 to ORDER_1D[DIM].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
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

        indx = new int[dim_num*order_nd];
        a = new int[dim_num];
        more = false;
        p = 0;

        for ( ; ; )
        {
            typeMethods.vec_colex_next2 ( dim_num, order_1d, ref a, ref more );

            if ( !more )
            {
                break;
            }

            for ( dim = 0; dim < dim_num; dim++ )
            {
                indx[dim+p*dim_num] = a[dim] + 1;
            }
            p += 1;
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

            p += 1;
        }

        return indx;
    }

    public static int[] multigrid_index_one(int dim_num, int[] order_1d, int order_nd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTIGRID_INDEX_ONE returns an indexed multidimensional grid.
        //
        //  Discussion:
        //
        //    For dimension DIM, the number of points is M (the order of the 1D rule).
        //
        //    We index the points as:
        //      1,   2,   3,   ...,  M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 October 2007
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
                indx[dim + p * dim_num] = a[dim] + 1;
            }

            p += 1;
        }

        return indx;
    }

    public static int[] multigrid_index_cfn(int dim_num, int[] order_1d, int order_nd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTIGRID_INDEX_CFN indexes a sparse grid based on CFN 1D rules.
        //
        //  Discussion:
        //
        //    The sparse grid is presumed to have been created from products
        //    of CLOSED FULLY NESTED 1D quadrature rules.
        //
        //    CFN rules include Clenshaw Curtis rules.
        //
        //    For dimension DIM, the second index of INDX may vary from 
        //    0 to ORDER_1D(DIM)-1.
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

            p += 1;
        }

        return indx;
    }

    public static int[] multigrid_index_ofn(int dim_num, int[] order_1d, int order_nd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTIGRID_INDEX_OFN indexes a sparse grid based on OFN 1D rules.
        //
        //  Discussion:
        //
        //    The sparse grid is presumed to have been created from products
        //    of OPEN FULLY NESTED 1D quadrature rules.
        //
        //    OFN rules include Fejer 1, Fejer 2, and Gauss Patterson rules.
        //
        //    For dimension DIM, the second index of INDX may vary from 
        //    1 to ORDER_1D(DIM).
        //
        //  Modified:
        //
        //    26 March 2008
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

            for (dim = 0; dim < dim_num; dim++)
            {
                indx[dim + p * dim_num] = a[dim] + 1;
            }

            p += 1;
        }

        return indx;
    }

    public static int[] multigrid_index_onn(int dim_num, int[] order_1d, int order_nd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTIGRID_INDEX_ONN indexes a sparse grid based on ONN 1D rules.
        //
        //  Discussion:
        //
        //    The sparse grid is presumed to have been created from products
        //    of OPEN NON-NESTED 1D quadrature rules.
        //
        //    ONN rules include Gauss Laguerre.
        //
        //    For dimension DIM, the number of points is ORDER_1D(DIM).
        //
        //    We index the points as
        //      1, 2, 3, ..., ORDER_1D(DIM).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 October 2007
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
                indx[dim + p * dim_num] = a[dim] + 1;
            }

            p += 1;
        }

        return indx;
    }

    public static int[] multigrid_index_own(int dim_num, int[] order_1d, int order_nd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTIGRID_INDEX_OWN returns an indexed multidimensional grid.
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

            p += 1;
        }

        return indx;
    }

    public static void multigrid_scale_open(int dim_num, int order_nd, int level_max,
            ref int[] level_1d, ref int[] grid_index)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTIGRID_SCALE_OPEN renumbers a grid as a subgrid on a higher level.
        //
        //  Discussion:
        //
        //    This routine takes a grid associated with a given value of
        //    LEVEL, and multiplies all the indices by a power of 2, so that
        //    the indices reflect the position of the same points, but in
        //    a grid of level LEVEL_MAX.
        //
        //    For an open grid, going from one level to the next, a set of indices
        //    will be rescaled by 2*INDEX-1.
        //
        //  Modified:
        //
        //    08 June 2007
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

        for (dim = 0; dim < dim_num; dim++)
        {
            factor = (int) Math.Pow(2, level_max - level_1d[dim]);

            for (order = 0; order < order_nd; order++)
            {
                grid_index[dim + order * dim_num] *= factor;
            }
        }
    }

    public static void multigrid_scale_closed(int dim_num, int order_nd, int level_max,
            int[] level_1d, ref int[] grid_index)

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
            switch (level_1d[dim])
            {
                case 0:
                {
                    order_max = level_max switch
                    {
                        0 => 1,
                        _ => (int) Math.Pow(2, level_max) + 1
                    };

                    for (order = 0; order < order_nd; order++)
                    {
                        grid_index[dim + order * dim_num] = (order_max - 1) / 2;
                    }

                    break;
                }
                default:
                {
                    factor = (int) Math.Pow(2, level_max - level_1d[dim]);
                    for (order = 0; order < order_nd; order++)
                    {
                        grid_index[dim + order * dim_num] *= factor;
                    }

                    break;
                }
            }
        }
    }
}