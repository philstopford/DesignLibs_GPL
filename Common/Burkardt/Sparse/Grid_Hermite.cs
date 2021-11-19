using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Composition;
using Burkardt.Grid;
using Burkardt.Quadrature;
using Burkardt.Types;
using Burkardt.Values;

namespace Burkardt.Sparse;

public static class Grid_Hermite
{
    public static void sparse_grid_hermite(int dim_num, int level_max, int point_num,
            ref double[] grid_weight, ref double[] grid_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPARSE_GRID_HERMITE computes a sparse grid of Gauss-Hermite points.
        //
        //  Discussion:
        //
        //    The quadrature rule is associated with a sparse grid derived from 
        //    a Smolyak construction using a 1D Gauss-Hermite quadrature rule.  
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
        //    05 July 2008
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
        //    by SPARSE_GRID_HERM_SIZE.
        //
        //    Output, double GRID_WEIGHT[POINT_NUM], the weights.
        //
        //    Output, double GRID_POINT[DIM_NUM*POINT_NUM], the points.
        //
    {
        int level;
        int point;
        int point3 = 0;

        for (point = 0; point < point_num; point++)
        {
            grid_weight[point] = 0.0;
        }

        //
        //  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
        //
        int point_num2 = 0;

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        int[] grid_base2 = new int[dim_num];
        int[] level_1d = new int[dim_num];
        int[] order_1d = new int[dim_num];

        for (level = level_min; level <= level_max; level++)
        {
            //
            //  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
            //  that adds up to LEVEL.
            //
            bool more = false;
            int h = 0;
            int t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);
                //
                //  Transform each 1D level to a corresponding 1D order.
                //  The relationship is the same as for other OPEN rules.
                //  The GL rule differs from the other OPEN rules only in the nesting behavior.
                //
                ClenshawCurtis.level_to_order_open(dim_num, level_1d, ref order_1d);

                int dim;
                for (dim = 0; dim < dim_num; dim++)
                {
                    grid_base2[dim] = (order_1d[dim] - 1) / 2;
                }

                //
                //  The product of the 1D orders gives us the number of points in this grid.
                //
                int order_nd = typeMethods.i4vec_product(dim_num, order_1d);
                //
                //  Compute the weights for this product grid.
                //
                double[] grid_weight2 = HermiteQuadrature.product_weight_hermite(dim_num, order_1d, order_nd);
                //
                //  Now determine the coefficient of the weight.
                //
                int coeff = (int)(Math.Pow(-1, level_max - level)
                                  * Binomial.choose(dim_num - 1, level_max - level));
                //
                //  The inner (hidden) loop generates all points corresponding to given grid.
                //  The grid indices will be between -M to +M, where 2*M + 1 = ORDER_1D(DIM).
                //
                int[] grid_index2 = Multigrid.multigrid_index_z(dim_num, order_1d, order_nd);
                //
                //  Determine the first level of appearance of each of the points.
                //  This allows us to flag certain points as being repeats of points
                //  generated on a grid of lower level.  
                //
                //  This is SLIGHTLY tricky.
                //
                int[] grid_level = HermiteQuadrature.index_level_hermite(level, level_max, dim_num, order_nd,
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
                        HermiteQuadrature.hermite_abscissa(dim_num, 1, grid_index2,
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
                        double[] grid_point_temp = new double[dim_num];

                        HermiteQuadrature.hermite_abscissa(dim_num, 1, grid_index2,
                            grid_base2, ref grid_point_temp, gridIndex: +point * dim_num);

                        int point2;
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
                                Console.WriteLine("SPARSE_GRID_HERM - Fatal error!");
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

    public static void sparse_grid_hermite_index(int dim_num, int level_max, int point_num,
            ref int[] grid_index, ref int[] grid_base)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPARSE_GRID_HERMITE_INDEX indexes points in a Gauss-Hermite sparse grid.
        //
        //  Discussion:
        //
        //    The sparse grid is assumed to be formed from 1D Gauss-Hermite rules
        //    of ODD order, which have the property that only the central abscissa,
        //    X = 0.0, is "nested".
        //
        //    The necessary dimensions of GRID_INDEX can be determined by 
        //    calling SPARSE_GRID_HERMITE_SIZE first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2008
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
        //    the orders of the Gauss-Hermite rules associated with each point 
        //    and dimension.
        //
    {
        int level;
        //
        //  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
        //
        int point_num2 = 0;

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        int[] grid_base2 = new int[dim_num];
        int[] level_1d = new int[dim_num];
        int[] order_1d = new int[dim_num];

        for (level = level_min; level <= level_max; level++)
        {
            //
            //  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
            //  that adds up to LEVEL.
            //
            bool more = false;
            int h = 0;
            int t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);
                //
                //  Transform each 1D level to a corresponding 1D order.
                //
                ClenshawCurtis.level_to_order_open(dim_num, level_1d, ref order_1d);
                int dim;
                for (dim = 0; dim < dim_num; dim++)
                {
                    grid_base2[dim] = (order_1d[dim] - 1) / 2;
                }

                //
                //  The product of the 1D orders gives us the number of points in this grid.
                //
                int order_nd = typeMethods.i4vec_product(dim_num, order_1d);
                //
                //  The inner (hidden) loop generates all points corresponding to given grid.
                //
                int[] grid_index2 = Multigrid.multigrid_index_z(dim_num, order_1d, order_nd);
                //
                //  Determine the first level of appearance of each of the points.
                //  This allows us to flag certain points as being repeats of points
                //  generated on a grid of lower level.  
                //
                //  This is SLIGHTLY tricky.
                //
                int[] grid_level = HermiteQuadrature.index_level_hermite(level, level_max, dim_num, order_nd,
                    grid_index2, grid_base2);
                //
                //  Only keep those points which first appear on this level.
                //
                int point;
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

    public static int sparse_grid_hermite_size(int dim_num, int level_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPARSE_GRID_HERMITE_SIZE sizes a sparse grid of Gauss-Hermite points.
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
        //    The grids are only very weakly nested, since Gauss-Hermite rules
        //    only have the origin in common.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2008
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
        //    Output, int SPARSE_GRID_HERM_SIZE, the number of points in the grid.
        //
    {
        int level;
        int point_num;
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

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        int[] level_1d = new int[dim_num];
        int[] order_1d = new int[dim_num];

        for (level = level_min; level <= level_max; level++)
        {
            //
            //  The middle loop generates the next partition that adds up to LEVEL.
            //
            bool more = false;
            int h = 0;
            int t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);
                //
                //  Transform each 1D level to a corresponding 1D order.
                //
                ClenshawCurtis.level_to_order_open(dim_num, level_1d, ref order_1d);

                int dim;
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
}