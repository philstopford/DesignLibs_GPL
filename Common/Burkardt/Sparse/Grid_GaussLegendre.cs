using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Composition;
using Burkardt.Grid;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace Burkardt.Sparse;

public static class Grid_GaussLegendre
{
    public static void sparse_grid_gl(int dim_num, int level_max, int point_num,
            ref double[] grid_weight, ref double[] grid_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPARSE_GRID_GL computes a sparse grid of Gauss-Legendre points.
        //
        //  Discussion:
        //
        //    The quadrature rule is associated with a sparse grid derived from 
        //    a Smolyak construction using a 1D Gauss-Legendre quadrature rule.  
        // 
        //    The user specifies: 
        //    * the spatial dimension of the quadrature region, 
        //    * the level that defines the Smolyak grid. 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 July 2008
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
        //    Input, int LEVEL_MAX, controls the size of the final sparse grid.
        //
        //    Input, int POINT_NUM, the number of points in the grid, as determined
        //    by SPARSE_GRID_GL_SIZE.
        //
        //    Output, double GRID_WEIGHT[POINT_NUM], the weights.
        //
        //    Output, double GRID_POINT[DIM_NUM*POINT_NUM], the points.
        //
    {
        int coeff;
        int dim;
        int[] grid_base2;
        int[] grid_index2;
        int[] grid_level;
        double[] grid_point_temp;
        double[] grid_weight2;
        int h;
        int level;
        int[] level_1d;
        int level_min;
        bool more;
        int[] order_1d;
        int order_nd;
        int point;
        int point_num2;
        int point2;
        int point3 = 0;
        int t;

        for (point = 0; point < point_num; point++)
        {
            grid_weight[point] = 0.0;
        }

        //
        //  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
        //
        point_num2 = 0;

        level_min = Math.Max(0, level_max + 1 - dim_num);

        grid_base2 = new int[dim_num];
        level_1d = new int[dim_num];
        order_1d = new int[dim_num];

        for (level = level_min; level <= level_max; level++)
        {
            //
            //  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
            //  that adds up to LEVEL.
            //
            more = false;
            h = 0;
            t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);
                //
                //  Transform each 1D level to a corresponding 1D order.
                //  The relationship is the same as for other OPEN rules.
                //  The GL rule differs from the other OPEN rules only in the nesting behavior.
                //
                ClenshawCurtis.level_to_order_open(dim_num, level_1d, ref order_1d);

                for (dim = 0; dim < dim_num; dim++)
                {
                    grid_base2[dim] = (order_1d[dim] - 1) / 2;
                }

                //
                //  The product of the 1D orders gives us the number of points in this grid.
                //
                order_nd = typeMethods.i4vec_product(dim_num, order_1d);
                //
                //  Compute the weights for this product grid.
                //
                grid_weight2 = GaussQuadrature.product_weight_gl(dim_num, order_1d, order_nd);
                //
                //  Now determine the coefficient of the weight.
                //
                coeff = (int)Math.Pow(-1, level_max - level)
                        * typeMethods.i4_choose(dim_num - 1, level_max - level);
                //
                //  The inner (hidden) loop generates all points corresponding to given grid.
                //  The grid indices will be between -M to +M, where 2*M + 1 = ORDER_1D(DIM).
                //
                grid_index2 = Multigrid.multigrid_index_z(dim_num, order_1d, order_nd);
                //
                //  Determine the first level of appearance of each of the points.
                //  This allows us to flag certain points as being repeats of points
                //  generated on a grid of lower level.  
                //
                //  This is SLIGHTLY tricky.
                //
                grid_level = index_level_gl(level, level_max, dim_num, order_nd,
                    grid_index2, grid_base2);
                //
                //  Only keep those points which first appear on this level.
                //
                for (point = 0; point < order_nd; point++)
                {
                    //
                    //  Either a "new" point (increase count, create point, create weight)
                    //
                    if (grid_level[point] == level)
                    {
                        GaussQuadrature.gl_abscissa(dim_num, 1, grid_index2,
                            grid_base2, ref grid_point, gridIndex: +point * dim_num,
                            gridPointIndex: +point_num2 * dim_num);

                        grid_weight[point_num2] = coeff * grid_weight2[point];

                        point_num2 += 1;
                    }
                    //
                    //  or an already existing point (create point temporarily, find match,
                    //  add weight to matched point's weight).
                    //
                    else
                    {
                        grid_point_temp = new double[dim_num];

                        GaussQuadrature.gl_abscissa(dim_num, 1, grid_index2,
                            grid_base2, ref grid_point_temp, gridIndex: +point * dim_num);

                        for (point2 = 0; point2 < point_num2; point2++)
                        {
                            point3 = point2;
                            for (dim = 0; dim < dim_num; dim++)
                            {
                                if (!(Math.Abs(grid_point[dim + point2 * dim_num] - grid_point_temp[dim]) >
                                      double.Epsilon))
                                {
                                    continue;
                                }

                                point3 = -1;
                                break;
                            }

                            if (point3 == point2)
                            {
                                break;
                            }
                        }

                        switch (point3)
                        {
                            case -1:
                                Console.WriteLine("");
                                Console.WriteLine("SPARSE_GRID_GL - Fatal error!");
                                Console.WriteLine("  Could not match point.");
                                return;
                            default:
                                grid_weight[point3] += coeff * grid_weight2[point];
                                break;
                        }
                    }
                }

                if (!more)
                {
                    break;
                }
            }
        }
    }

    public static void sparse_grid_gl_index(int dim_num, int level_max, int point_num,
            ref int[] grid_index, ref int[] grid_base)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPARSE_GRID_GL_INDEX indexes the points forming a sparse grid of GL points.
        //
        //  Discussion:
        //
        //    The sparse grid is assumed to be formed from 1D Gauss-Legendre rules
        //    of ODD order, which have the property that only the central abscissa,
        //    X = 0.0, is "nested".
        //
        //    The necessary dimensions of GRID_INDEX can be determined by 
        //    calling SPARSE_GRID_GL_SIZE first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 July 2008
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
        //    Input, int LEVEL_MAX, the maximum value of LEVEL.
        //
        //    Input, int POINT_NUM, the total number of points in the grids.
        //
        //    Output, int GRID_INDEX[DIM_NUM*POINT_NUM], a list of 
        //    point indices, representing a subset of the product grid of level 
        //    LEVEL_MAX, representing (exactly once) each point that will show up in a
        //    sparse grid of level LEVEL_MAX.
        //
        //    Output, int GRID_BASE[DIM_NUM*POINT_NUM], a list of 
        //    the orders of the Gauss-Legendre rules associated with each point and dimension.
        //
    {
        int dim;
        int[] grid_base2;
        int[] grid_index2;
        int[] grid_level;
        int h;
        int level;
        int[] level_1d;
        int level_min;
        bool more;
        int[] order_1d;
        int order_nd;
        int point;
        int point_num2;
        int t;
        //
        //  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
        //
        point_num2 = 0;

        level_min = Math.Max(0, level_max + 1 - dim_num);

        grid_base2 = new int[dim_num];
        level_1d = new int[dim_num];
        order_1d = new int[dim_num];

        for (level = level_min; level <= level_max; level++)
        {
            //
            //  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
            //  that adds up to LEVEL.
            //
            more = false;
            h = 0;
            t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);
                //
                //  Transform each 1D level to a corresponding 1D order.
                //
                ClenshawCurtis.level_to_order_open(dim_num, level_1d, ref order_1d);
                for (dim = 0; dim < dim_num; dim++)
                {
                    grid_base2[dim] = (order_1d[dim] - 1) / 2;
                }

                //
                //  The product of the 1D orders gives us the number of points in this grid.
                //
                order_nd = typeMethods.i4vec_product(dim_num, order_1d);
                //
                //  The inner (hidden) loop generates all points corresponding to given grid.
                //
                grid_index2 = Multigrid.multigrid_index_z(dim_num, order_1d, order_nd);
                //
                //  Determine the first level of appearance of each of the points.
                //  This allows us to flag certain points as being repeats of points
                //  generated on a grid of lower level.  
                //
                //  This is SLIGHTLY tricky.
                //
                grid_level = index_level_gl(level, level_max, dim_num, order_nd,
                    grid_index2, grid_base2);
                //
                //  Only keep those points which first appear on this level.
                //
                for (point = 0; point < order_nd; point++)
                {
                    if (grid_level[point] == level)
                    {
                        for (dim = 0; dim < dim_num; dim++)
                        {
                            grid_index[dim + point_num2 * dim_num] =
                                grid_index2[dim + point * dim_num];
                            grid_base[dim + point_num2 * dim_num] = grid_base2[dim];
                        }

                        point_num2 += 1;
                    }
                }

                if (!more)
                {
                    break;
                }
            }
        }
    }

    public static int sparse_grid_gl_size(int dim_num, int level_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPARSE_GRID_GL_SIZE sizes a sparse grid of Gauss-Legendre points.
        //
        //  Discussion:
        //
        //    The grid is defined as the sum of the product rules whose LEVEL
        //    satisfies:
        //
        //      LEVEL_MIN <= LEVEL <= LEVEL_MAX.
        //
        //    where LEVEL_MAX is user specified, and 
        //
        //      LEVEL_MIN = max ( 0, LEVEL_MAX + 1 - DIM_NUM ).
        //
        //    The grids are only very weakly nested, since Gauss-Legendre rules
        //    only have the origin in common.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 July 2008
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
        //    Input, int LEVEL_MAX, the maximum value of LEVEL.
        //
        //    Output, int SPARSE_GRID_GL_SIZE, the number of points in the grid.
        //
    {
        int dim;
        int h;
        int level;
        int[] level_1d;
        int level_min;
        bool more;
        int[] order_1d;
        int point_num;
        int t;
        switch (level_max)
        {
            //
            //  Special case.
            //
            case 0:
                point_num = 1;
                return point_num;
        }

        //
        //  The outer loop generates LEVELs from 0 to LEVEL_MAX.
        //
        point_num = 0;

        level_min = Math.Max(0, level_max + 1 - dim_num);

        level_1d = new int[dim_num];
        order_1d = new int[dim_num];

        for (level = level_min; level <= level_max; level++)
        {
            //
            //  The middle loop generates the next partition that adds up to LEVEL.
            //
            more = false;
            h = 0;
            t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);
                //
                //  Transform each 1D level to a corresponding 1D order.
                //
                ClenshawCurtis.level_to_order_open(dim_num, level_1d, ref order_1d);

                for (dim = 0; dim < dim_num; dim++)
                {
                    //
                    //  If we can reduce the level in this dimension by 1 and
                    //  still not go below LEVEL_MIN.
                    //
                    if (level_min < level && 1 < order_1d[dim])
                    {
                        order_1d[dim] -= 1;
                    }
                }

                point_num += typeMethods.i4vec_product(dim_num, order_1d);

                if (!more)
                {
                    break;
                }
            }
        }

        return point_num;
    }

    public static int[] index_level_gl(int level, int level_max, int dim_num, int point_num,
            int[] grid_index, int[] grid_base)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_LEVEL_GL: determine first level at which given index is generated.
        //
        //  Discussion:
        //
        //    We are constructing a sparse grid of Gauss-Legendre points.  The grid
        //    is built up of product grids, with a characteristic LEVEL.  
        //
        //    We are concerned with identifying points in this product grid which
        //    have actually been generated previously, on a lower value of LEVEL.
        //
        //    This routine determines the lowest value of LEVEL at which each of
        //    the input points would be generated.
        //
        //    In 1D, given LEVEL, the number of points is ORDER = 2**(LEVEL+1) + 1,
        //    (except that LEVEL = 0 implies ORDER = 1//), the BASE is (ORDER-1)/2, 
        //    and the point INDEX values range from -BASE to +BASE.
        //
        //    The values of INDEX and BASE allow us to determine the abstract
        //    properties of the point.  In particular, if INDEX is 0, the corresponding
        //    Gauss-Legendre abscissa is 0, the special "nested" value we need
        //    to take care of.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 September 2007
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
        //    Input, int LEVEL, the level at which these points were 
        //    generated.  LEVEL_MIN <= LEVEL <= LEVEL_MAX.
        //
        //    Input, int LEVEL_MAX, the maximum level.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int POINT_NUM, the number of points to be tested.
        //
        //    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], the indices of the 
        //    points to be tested.
        //
        //    Input, int GRID_BASE[DIM_NUM], the "base", which is essentially
        //    the denominator of the index.
        //
        //    Output, int INDEX_LEVEL_GL[POINT_NUM], the value of LEVEL at 
        //    which the point would first be generated.  This will be the same as
        //    the input value of LEVEL, unless the point has an INDEX of 0 and
        //    a corresponding BASE that is NOT zero.
        //
    {
        int dim;
        int[] grid_level;
        int level_min;
        int point;

        grid_level = new int[point_num];

        level_min = Math.Max(0, level_max + 1 - dim_num);
        //
        //  If a point has a DIM-th component whose INDEX is 0, then the 
        //  value of LEVEL at which this point would first be generated is
        //  less than LEVEL, unless the DIM-th component of GRID_BASE is 0.
        //
        for (point = 0; point < point_num; point++)
        {
            grid_level[point] = Math.Max(level, level_min);

            for (dim = 0; dim < dim_num; dim++)
            {
                grid_level[point] = grid_index[dim + point * dim_num] switch
                {
                    0 => Math.Max(grid_level[point] - grid_base[dim], level_min),
                    _ => grid_level[point]
                };
            }
        }

        return grid_level;
    }
}