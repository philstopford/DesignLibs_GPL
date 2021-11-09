using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Composition;
using Burkardt.Grid;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace Burkardt.Sparse
{
    public static class Grid
    {
        public static void levels_index(int dim_num, int level_max, int rule, int point_num,
                ref int[] grid_index, ref int[] grid_base)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEVELS_INDEX indexes a sparse grid.
            //
            //  Discussion:
            //
            //    The sparse grid is the logical sum of product grids with total LEVEL 
            //    between LEVEL_MIN and LEVEL_MAX.
            //
            //    The necessary dimensions of GRID_INDEX can be determined by 
            //    calling LEVELS_INDEX_SIZE first.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
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
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int LEVEL_MAX, the maximum value of LEVEL.
            //
            //    Input, int RULE, the index of the rule.
            //    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
            //    2, "F1", Fejer 1 Open Fully Nested rule.
            //    3, "F2", Fejer 2 Open Fully Nested rule.
            //    4, "GP", Gauss Patterson Open Fully Nested rule.
            //    5, "GL", Gauss Legendre Open Weakly Nested rule.
            //    6, "GH", Gauss Hermite Open Weakly Nested rule.
            //    7, "LG", Gauss Laguerre Open Non Nested rule.
            //
            //    Input, int POINT_NUM, the total number of points 
            //    in the grids.
            //
            //    Output, int GRID_INDEX[DIM_NUM*POINT_NUM], a list of 
            //    point indices, representing a subset of the product grid of level 
            //    LEVEL_MAX, representing (exactly once) each point that will show up in a
            //    sparse grid of level LEVEL_MAX.
            //
            //    Output, int GRID_BASE[DIM_NUM*POINT_NUM], a list of 
            //    the orders of the rules associated with each point and dimension.
            //
        {
            if (rule == 1)
            {
                levels_index_cfn(dim_num, level_max, point_num, ref grid_index, ref grid_base);
            }
            else if (2 <= rule && rule <= 4)
            {
                levels_index_ofn(dim_num, level_max, point_num, ref grid_index, ref grid_base);
            }
            else if (5 <= rule && rule <= 6)
            {
                levels_index_own(dim_num, level_max, point_num, ref grid_index, ref grid_base);
            }
            else if (7 == rule)
            {
                levels_index_onn(dim_num, level_max, point_num, ref grid_index, ref grid_base);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("LEVELS_INDEX - Fatal error!");
                Console.WriteLine("  Unrecognized rule number = " + rule + "");
            }
        }

        public static void levels_index_cfn(int dim_num, int level_max, int point_num,
                ref int[] grid_index, ref int[] grid_base)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEVELS_INDEX_CFN indexes a sparse grid made from CFN 1D rules.
            //
            //  Discussion:
            //
            //    The sparse grid is presumed to have been created from products
            //    of CLOSED FULLY NESTED 1D quadrature rules.
            //
            //    CFN rules include Clenshaw Curtis rules.
            //
            //    The sparse grid is the logical sum of product grids with total LEVEL 
            //    between LEVEL_MIN and LEVEL_MAX.
            //
            //    The necessary dimensions of GRID_INDEX can be determined by 
            //    calling LEVELS_INDEX_SIZE_CFN first.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 July 2008
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
            //    Output, int LEVELS_INDEX_CFN[DIM_NUM*POINT_NUM], a list of point 
            //    indices, representing a subset of the product grid of level LEVEL_MAX,
            //    representing (exactly once) each point that will show up in a
            //    sparse grid of level LEVEL_MAX.
            //
        {
            int dim;
            int[] grid_index2;
            int[] grid_level;
            int h;
            int level;
            int[] level_1d;
            bool more;
            int[] order_1d;
            int order_nd;
            int point;
            int point_num2;
            int t;
            //
            //  The outer loop generates LEVELs from 0 to LEVEL_MAX.
            //
            point_num2 = 0;

            level_1d = new int[dim_num];
            order_1d = new int[dim_num];

            for (level = 0; level <= level_max; level++)
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
                    LevelToOrder.level_to_order_closed(dim_num, level_1d, ref order_1d);
                    //
                    //  The product of the 1D orders gives us the number of points in this grid.
                    //
                    order_nd = typeMethods.i4vec_product(dim_num, order_1d);
                    //
                    //  The inner (hidden) loop generates all points corresponding to given grid.
                    //
                    grid_index2 = Multigrid.multigrid_index_cfn(dim_num, order_1d, order_nd);
                    //
                    //  Adjust these grid indices to reflect LEVEL_MAX.
                    //
                    Multigrid.multigrid_scale_closed(dim_num, order_nd, level_max, level_1d,
                        ref grid_index2);
                    //
                    //  Determine the first level of appearance of each of the points.
                    //
                    grid_level = Abscissa.abscissa_level_closed_nd(level_max, dim_num, order_nd,
                        grid_index2);
                    //
                    //  Only keep those points which first appear on this level.
                    //
                    for (point = 0; point < order_nd; point++)
                    {
                        if (grid_level[point] == level)
                        {
                            if (point_num <= point_num2)
                            {
                                Console.WriteLine("");
                                Console.WriteLine("LEVELS_INDEX_CFN - Fatal error!");
                                Console.WriteLine("  Exceeding maximum point index POINT_NUM = "
                                                  + point_num + "");
                                return;
                            }

                            for (dim = 0; dim < dim_num; dim++)
                            {
                                grid_base[dim + point_num2 * dim_num] = order_1d[dim];
                                grid_index[dim + point_num2 * dim_num] =
                                    grid_index2[dim + point * dim_num];
                            }

                            point_num2 = point_num2 + 1;
                        }
                    }

                    if (!more)
                    {
                        break;
                    }
                }
            }

            if (point_num2 < point_num)
            {
                Console.WriteLine("");
                Console.WriteLine("LEVELS_INDEX_CFN - Fatal error!");
                Console.WriteLine("  Set fewer points than POINT_NUM = " + point_num + "");
            }

        }

        public static void levels_index_ofn(int dim_num, int level_max, int point_num,
                ref int[] grid_index, ref int[] grid_base)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEVELS_INDEX_OFN indexes a sparse grid made from OFN 1D rules.
            //
            //  Discussion:
            //
            //    The sparse grid is presumed to have been created from products
            //    of OPEN FULLY NESTED 1D quadrature rules.
            //
            //    OFN rules include Fejer 1, Fejer 2, and Gauss Patterson rules.
            //
            //    The sparse grid is the logical sum of product grids with total LEVEL 
            //    between LEVEL_MIN and LEVEL_MAX.
            //
            //    The necessary dimensions of GRID_INDEX can be determined by 
            //    calling LEVELS_INDEX_SIZE_OFN first.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 July 2008
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
            //    the orders of the rules associated with each point and dimension.
            //
        {
            int dim;
            int[] grid_index2;
            int h;
            int level;
            int[] level_1d;
            bool more;
            int[] order_1d;
            int order_nd;
            int point;
            int point_num2;
            int t;
            bool test;
            //
            //  The outer loop generates LEVELs from 0 to LEVEL_MAX.
            //
            level_1d = new int[dim_num];
            order_1d = new int[dim_num];

            point_num2 = 0;

            for (level = 0; level <= level_max; level++)
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
                    LevelToOrder.level_to_order_open(dim_num, level_1d, ref order_1d);
                    //
                    //  The product of the 1D orders gives us the number of points in this grid.
                    //
                    order_nd = typeMethods.i4vec_product(dim_num, order_1d);
                    //
                    //  The inner (hidden) loop generates all points corresponding to given grid.
                    //
                    grid_index2 = Multigrid.multigrid_index_ofn(dim_num, order_1d, order_nd);
                    //
                    //  Only keep those points which first appear on this level.
                    //  If you keep a point, it is necessary to rescale each of its components
                    //  so that we save the coordinates as they apply on the final grid.
                    //
                    for (point = 0; point < order_nd; point++)
                    {
                        test = true;
                        for (dim = 0; dim < dim_num; dim++)
                        {
                            if (grid_index2[dim + point * dim_num] % 2 == 0)
                            {
                                test = false;
                            }
                        }

                        if (test)
                        {

                            if (point_num <= point_num2)
                            {
                                Console.WriteLine("LEVELS_INDEX_OFN - Fatal error!");
                                Console.WriteLine("  Exceeding maximum point index POINT_NUM = "
                                                  + point_num + "");
                                return;
                            }

                            for (dim = 0; dim < dim_num; dim++)
                            {
                                grid_base[dim + point_num2 * dim_num] = order_1d[dim];

                                grid_index[dim + point_num2 * dim_num] =
                                    (int) Math.Pow(2, level_max - level_1d[dim])
                                    * grid_index2[dim + point * dim_num];
                            }

                            point_num2 = point_num2 + 1;
                        }
                    }


                    if (!more)
                    {
                        break;
                    }
                }
            }

            if (point_num2 < point_num)
            {
                Console.WriteLine("");
                Console.WriteLine("LEVELS_INDEX_OFN - Fatal error!");
                Console.WriteLine("  Set fewer points than POINT_NUM = " + point_num + "");
            }

        }

        public static void levels_index_onn(int dim_num, int level_max, int point_num,
                ref int[] grid_index, ref int[] grid_base)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEVELS_INDEX_ONN indexes a sparse grid made from ONN 1D rules.
            //
            //  Discussion:
            //
            //    The sparse grid is presumed to have been created from products
            //    of OPEN NON NESTED 1D quadrature rules.
            //
            //    ONN rules include Gauss Laguerre.
            //
            //    The sparse grid is the logical sum of product grids with total LEVEL 
            //    between LEVEL_MIN and LEVEL_MAX.
            //
            //    The necessary dimensions of GRID_INDEX can be determined by 
            //    calling LEVELS_INDEX_SIZE_ONN first.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 July 2008
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
            //    the orders of the rules associated with each point and dimension.
            //
        {
            int dim;
            int[] grid_base2;
            int[] grid_index2;
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
                    LevelToOrder.level_to_order_open(dim_num, level_1d, ref order_1d);

                    for (dim = 0; dim < dim_num; dim++)
                    {
                        grid_base2[dim] = order_1d[dim];
                    }

                    //
                    //  The product of the 1D orders gives us the number of points in this grid.
                    //
                    order_nd = typeMethods.i4vec_product(dim_num, order_1d);
                    //
                    //  The inner (hidden) loop generates all points corresponding to given grid.
                    //
                    grid_index2 = Multigrid.multigrid_index_onn(dim_num, order_1d, order_nd);
                    //
                    //  Only keep those points which first appear on this level.
                    //
                    for (point = 0; point < order_nd; point++)
                    {
                        if (point_num <= point_num2)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("LEVELS_INDEX_ONN - Fatal error!");
                            Console.WriteLine("  Exceeding maximum point index POINT_NUM = "
                                              + point_num + "");
                            return;
                        }

                        for (dim = 0; dim < dim_num; dim++)
                        {
                            grid_index[dim + point_num2 * dim_num] = grid_index2[dim + point * dim_num];
                            grid_base[dim + point_num2 * dim_num] = grid_base2[dim];
                        }

                        point_num2 = point_num2 + 1;
                    }

                    if (!more)
                    {
                        break;
                    }
                }
            }

            if (point_num2 < point_num)
            {
                Console.WriteLine("");
                Console.WriteLine("LEVELS_INDEX_ONN - Fatal error!");
                Console.WriteLine("  Set fewer points than POINT_NUM = " + point_num + "");
            }

        }

        public static void levels_index_own(int dim_num, int level_max, int point_num,
                ref int[] grid_index, ref int[] grid_base)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEVELS_INDEX_OWN indexes a sparse grid made from OWN 1D rules.
            //
            //  Discussion:
            //
            //    The sparse grid is presumed to have been created from products
            //    of OPEN WEAKLY NESTED 1D quadrature rules.
            //
            //    OWN rules include Gauss Hermite and Gauss Legendre.
            //
            //    The sparse grid is the logical sum of product grids with total LEVEL 
            //    between LEVEL_MIN and LEVEL_MAX.
            //
            //    The necessary dimensions of GRID_INDEX can be determined by 
            //    calling LEVELS_INDEX_SIZE_OWN first.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 July 2008
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
            //    the orders of the rules associated with each point and dimension.
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

            if (dim_num == 1)
            {
                level_min = level_max;
            }
            else
            {
                level_min = 0;
            }

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
                    LevelToOrder.level_to_order_open(dim_num, level_1d, ref order_1d);

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
                    grid_index2 = Multigrid.multigrid_index_own(dim_num, order_1d, order_nd);
                    //
                    //  Determine the first level of appearance of each of the points.
                    //  This allows us to flag certain points as being repeats of points
                    //  generated on a grid of lower level.  
                    //
                    //  This is SLIGHTLY tricky.
                    //
                    grid_level = LevelToOrder.index_level_own(level, level_max, dim_num, order_nd,
                        grid_index2, grid_base2);
                    //
                    //  Only keep those points which first appear on this level.
                    //
                    for (point = 0; point < order_nd; point++)
                    {
                        if (grid_level[point] == level)
                        {
                            if (point_num <= point_num2)
                            {
                                Console.WriteLine("");
                                Console.WriteLine("LEVELS_INDEX_OWN - Fatal error!");
                                Console.WriteLine("  Exceeding maximum point index POINT_NUM = "
                                                  + point_num + "");
                                return;
                            }

                            for (dim = 0; dim < dim_num; dim++)
                            {
                                grid_index[dim + point_num2 * dim_num] =
                                    grid_index2[dim + point * dim_num];
                                grid_base[dim + point_num2 * dim_num] = grid_base2[dim];
                            }

                            point_num2 = point_num2 + 1;
                        }
                    }

                    if (!more)
                    {
                        break;
                    }
                }
            }

            if (point_num2 < point_num)
            {
                Console.WriteLine("");
                Console.WriteLine("LEVELS_INDEX_OWN - Fatal error!");
                Console.WriteLine("  Set fewer points than POINT_NUM = " + point_num + "");
            }

        }

        public static int levels_index_size(int dim_num, int level_max, int rule)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEVELS_INDEX_SIZE sizes a sparse grid.
            //
            //  Discussion:
            //
            //    The sparse grid is the logical sum of product grids with total LEVEL 
            //    between LEVEL_MIN and LEVEL_MAX.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 December 2009
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
            //    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
            //
            //    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
            //
            //    Input, integer ( kind = 4 ) RULE, the index of the rule.
            //    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
            //    2, "F1", Fejer 1 Open Fully Nested rule.
            //    3, "F2", Fejer 2 Open Fully Nested rule.
            //    4, "GP", Gauss Patterson Open Fully Nested rule.
            //    5, "GL", Gauss Legendre Open Weakly Nested rule.
            //    6, "GH", Gauss Hermite Open Weakly Nested rule.
            //    7, "LG", Gauss Laguerre Open Non Nested rule.
            //
            //    Output, int  LEVELS_INDEX_SIZE, the total number of unique 
            //    points in the grids.
            //
        {
            int point_num;

            if (rule == 1)
            {
                point_num = sparse_grid_cc_size(dim_num, level_max);
            }
            else if (2 <= rule && rule <= 4)
            {
                point_num = sparse_grid_ofn_size(dim_num, level_max);
            }
            else if (5 <= rule && rule <= 6)
            {
                point_num = levels_index_size_own(dim_num, level_max);
            }
            else if (7 == rule)
            {
                point_num = levels_index_size_onn(dim_num, level_max);
            }
            else
            {
                point_num = -1;
                Console.WriteLine("");
                Console.WriteLine("LEVELS_INDEX_SIZE - Fatal error!");
                Console.WriteLine("  Unrecognized value of RULE = " + rule + "");
            }

            return point_num;
        }

        public static int levels_index_size_cfn(int dim_num, int level_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEVELS_INDEX_SIZE_CFN sizes a sparse grid made from CFN 1D rules.
            //
            //  Discussion:
            //
            //    The sparse grid is presumed to have been created from products
            //    of CLOSED FULLY NESTED 1D quadrature rules.
            //
            //    CFN rules include Clenshaw Curtis rules.
            //
            //    The sparse grid is the logical sum of product grids with total LEVEL 
            //    between LEVEL_MIN and LEVEL_MAX.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 July 2008
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
            //    Output, int LEVELS_INDEX_SIZE_CFN, the number of points in the grid.
            //
        {
            int[] grid_index;
            int[] grid_level;
            int h;
            int level;
            int[] level_1d;
            bool more;
            int[] order_1d;
            int order_nd;
            int point;
            int point_num;
            int t;
            //
            //  Special case.
            //
            if (level_max == 0)
            {
                point_num = 1;
                return point_num;
            }

            //
            //  The outer loop generates LEVELs from 0 to LEVEL_MAX.
            //
            point_num = 0;

            level_1d = new int[dim_num];
            order_1d = new int[dim_num];

            for (level = 0; level <= level_max; level++)
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
                    LevelToOrder.level_to_order_closed(dim_num, level_1d, ref order_1d);
                    //
                    //  The product of the 1D orders gives us the number of points in this grid.
                    //
                    order_nd = typeMethods.i4vec_product(dim_num, order_1d);
                    //
                    //  The inner (hidden) loop generates all points corresponding to given grid.
                    //
                    grid_index = Multigrid.multigrid_index_cfn(dim_num, order_1d, order_nd);
                    //
                    //  Adjust these grid indices to reflect LEVEL_MAX.
                    //
                    Multigrid.multigrid_scale_closed(dim_num, order_nd, level_max, level_1d,
                        ref grid_index);
                    //
                    //  Determine the first level of appearance of each of the points.
                    //
                    grid_level = Abscissa.abscissa_level_closed_nd(level_max, dim_num, order_nd,
                        grid_index);
                    //
                    //  Only keep those points which first appear on this level.
                    //
                    for (point = 0; point < order_nd; point++)
                    {
                        if (grid_level[point] == level)
                        {
                            point_num = point_num + 1;
                        }
                    }

                    if (!more)
                    {
                        break;
                    }
                }
            }

            return point_num;
        }

        public static int levels_index_size_onn(int dim_num, int level_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEVELS_INDEX_SIZE_ONN sizes a sparse grid made from ONN 1D rules.
            //
            //  Discussion:
            //
            //    The sparse grid is presumed to have been created from products
            //    of OPEN NON-NESTED 1D quadrature rules.
            //
            //    ONN rules include Gauss Laguerre.
            //
            //    The sparse grid is the logical sum of product grids with total LEVEL 
            //    between LEVEL_MIN and LEVEL_MAX.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 July 2008
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
            //    Output, int LEVELS_INDEX_SIZE_ONN, the number of points in the grid.
            //
        {
            int h;
            int level;
            int[] level_1d;
            int level_min;
            bool more;
            int[] order_1d;
            int point_num;
            int t;
            //
            //  Special case.
            //
            if (level_max == 0)
            {
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
                    LevelToOrder.level_to_order_open(dim_num, level_1d, ref order_1d);

                    point_num = point_num + typeMethods.i4vec_product(dim_num, order_1d);

                    if (!more)
                    {
                        break;
                    }
                }
            }

            return point_num;
        }

        public static int levels_index_size_own(int dim_num, int level_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEVELS_INDEX_SIZE_OWN sizes a sparse grid made from OWN 1D rules.
            //
            //  Discussion:
            //
            //    The sparse grid is presumed to have been created from products
            //    of OPEN WEAKLY NESTED 1D quadrature rules.
            //
            //    OWN rules include Gauss Hermite and Gauss Legendre.
            //
            //    The sparse grid is the logical sum of product grids with total LEVEL 
            //    between LEVEL_MIN and LEVEL_MAX.
            //
            //    Oddly enough, in order to count the number of points, we will
            //    behave as though LEVEL_MIN was zero.  This is because our computation
            //    concentrates on throwing away all points generated at lower levels,
            //    but, in fact, if we start at a nonzero level, we need to include
            //    on that level all the points that would have been generated on lower
            //    levels.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 July 2008
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
            //    Output, int LEVELS_INDEX_SIZE_OWN, the number of points in the grid.
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
            //
            //  Special case.
            //
            if (level_max == 0)
            {
                point_num = 1;
                return point_num;
            }

            //
            //  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
            //
            //  The normal definition of LEVEL_MIN:
            //
            //   level_min = max ( 0, level_max + 1 - dim_num )
            //
            //  Our somewhat artificial temporary local definition of LEVEL_MIN:
            //
            if (dim_num == 1)
            {
                level_min = level_max;
                point_num = 1;
            }
            else
            {
                level_min = 0;
                point_num = 0;
            }

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
                    LevelToOrder.level_to_order_open(dim_num, level_1d, ref order_1d);

                    for (dim = 0; dim < dim_num; dim++)
                    {
                        //
                        //  If we can reduce the level in this dimension by 1 and
                        //  still not go below LEVEL_MIN.
                        //
                        if (1 < order_1d[dim])
                        {
                            order_1d[dim] = order_1d[dim] - 1;
                        }
                    }

                    point_num = point_num + typeMethods.i4vec_product(dim_num, order_1d);

                    if (!more)
                    {
                        break;
                    }
                }
            }

            return point_num;
        }

        public static int sparse_grid_cc_size(int dim_num, int level_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPARSE_GRID_CC_SIZE sizes a sparse grid using Clenshaw Curtis rules.
            //
            //  Discussion:
            //
            //    The grid is defined as the sum of the product rules whose LEVEL
            //    satisfies:
            //
            //      0 <= LEVEL <= LEVEL_MAX.
            //
            //    This calculation is much faster than a previous method.  It simply
            //    computes the number of new points that are added at each level in the
            //    1D rule, and then counts the new points at a given DIM_NUM dimensional
            //    level vector as the product of the new points added in each dimension.
            //
            //    This approach will work for nested families, and may be extensible
            //    to other families, and to mixed rules.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 December 2009
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
            //    Output, int SPARSE_GRID_CC_SIZE, the number of points in the grid.
            //
        {
            int dim;
            int h;
            int j;
            int l;
            int level;
            int[] level_1d;
            bool more;
            int[] new_1d;
            int point_num;
            int t;
            int v;
            //
            //  Special case.
            //
            if (level_max < 0)
            {
                point_num = 0;
                return point_num;
            }

            if (level_max == 0)
            {
                point_num = 1;
                return point_num;
            }

            //
            //  Construct the vector that counts the new points in the 1D rule.
            //
            new_1d = new int[level_max + 1];

            new_1d[0] = 1;
            new_1d[1] = 2;

            j = 1;
            for (l = 2; l <= level_max; l++)
            {
                j = j * 2;
                new_1d[l] = j;
            }

            //
            //  Count the number of points by counting the number of new points 
            //  associated with each level vector.
            //
            level_1d = new int[dim_num];

            point_num = 0;

            for (level = 0; level <= level_max; level++)
            {
                more = false;
                h = 0;
                t = 0;

                for (;;)
                {
                    Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);

                    v = 1;
                    for (dim = 0; dim < dim_num; dim++)
                    {
                        v = v * new_1d[level_1d[dim]];
                    }

                    point_num = point_num + v;

                    if (!more)
                    {
                        break;
                    }
                }
            }

            return point_num;
        }

        public static void sparse_grid(int dim_num, int level_max, int rule, int point_num,
                ref double[] grid_weight, ref double[] grid_point)

            //***************************************************************************80 
            // 
            //  Purpose:
            //
            //    SPARSE_GRID computes a sparse grid. 
            // 
            //  Discussion: 
            //
            //    A Smolyak construction is used to create a multidimensional sparse grid.  
            // 
            //    The user specifies: 
            //    * the spatial dimension of the quadrature region, 
            //    * the level that defines the Smolyak grid. 
            //    * the 1D quadrature rule. 
            // 
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 March 2008
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
            //    Input, int LEVEL_MAX, controls the size of the final 
            //    sparse grid. 
            // 
            //    Input, int RULE, the index of the rule.
            //    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
            //    2, "F1", Fejer 1 Open Fully Nested rule.
            //    3, "F2", Fejer 2 Open Fully Nested rule.
            //    4, "GP", Gauss Patterson Open Fully Nested rule.
            //    5, "GL", Gauss Legendre Open Weakly Nested rule.
            //    6, "GH", Gauss Hermite Open Weakly Nested rule.
            //    7, "LG", Gauss Laguerre Open Non Nested rule.
            //
            //    Input, int POINT_NUM, the number of points in the grid, 
            //    as determined by LEVELS_INDEX_SIZE.
            //
            //    Output, double GRID_WEIGHT[POINT_NUM], the weights. 
            // 
            //    Output, double GRID_POINT[DIM_NUM*POINT_NUM], the points. 
            // 
        {
            if (rule == 1)
            {
                sparse_grid_cfn(dim_num, level_max, rule, point_num, ref grid_weight,
                    ref grid_point);
            }
            else if (2 <= rule && rule <= 4)
            {
                sparse_grid_ofn(dim_num, level_max, rule, point_num, ref grid_weight,
                    ref grid_point);
            }
            else if (5 <= rule && rule <= 6)
            {
                sparse_grid_own(dim_num, level_max, rule, point_num, ref grid_weight,
                    ref grid_point);
            }
            else if (7 == rule)
            {
                sparse_grid_onn(dim_num, level_max, rule, point_num, ref grid_weight,
                    ref grid_point);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("SPARSE_GRID - Fatal error!");
                Console.WriteLine("  Illegal input rule index = " + rule + "");
            }
        }

        public static void sparse_grid_cfn(int dim_num, int level_max, int rule, int point_num,
                ref double[] grid_weight, ref double[] grid_point)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPARSE_GRID_CFN computes a sparse grid based on a CFN 1D rule. 
            //
            //  Discussion:
            //
            //    The 1D quadrature rule is assumed to be Closed Fully Nested.
            //
            //    Closed Fully Nested rules include Clenshaw Curtis rules.
            //
            //    A Smolyak construction is used to create a multidimensional sparse grid.  
            // 
            //    The user specifies: 
            //    * the spatial dimension of the quadrature region, 
            //    * the level that defines the Smolyak grid. 
            //    * the quadrature rule.
            //    * the number of points. 
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
            //    Input, int LEVEL_MAX, controls the size of the final sparse grid.
            //
            //    Input, int RULE, the index of the rule.
            //    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
            //    2, "F1", Fejer 1 Open Fully Nested rule.
            //    3, "F2", Fejer 2 Open Fully Nested rule.
            //    4, "GP", Gauss Patterson Open Fully Nested rule.
            //    5, "GL", Gauss Legendre Open Weakly Nested rule.
            //    6, "GH", Gauss Hermite Open Weakly Nested rule.
            //    7, "LG", Gauss Laguerre Open Non Nested rule.
            //
            //    Input, int POINT_NUM, the number of points in the grid, as determined
            //    by SPARSE_GRID_SIZE_CFN.
            //
            //    Output, double GRID_WEIGHTS[POINT_NUM], the weights.
            //
            //    Output, double GRID_POINTS[DIM_NUM*POINT_NUM], the points.
            //
        {
            int dim;
            int[] grid_base;
            int[] grid_index;
            int order_max;
            int point;

            if (rule != 1)
            {
                Console.WriteLine("");
                Console.WriteLine("SPARSE_GRID_CFN - Fatal error!");
                Console.WriteLine("  Illegal input rule index = " + rule + "");
                return;
            }

            //
            //  Determine the index vector, relative to the full product grid,
            //  that identifies the points in the sparse grid.
            //
            grid_index = new int[dim_num * point_num];
            grid_base = new int[dim_num * point_num];

            levels_index_cfn(dim_num, level_max, point_num, ref grid_index, ref grid_base);
            //
            //  Compute the physical coordinates of the abscissas.
            //
            if (0 == level_max)
            {
                order_max = 1;
            }
            else
            {
                order_max = (int) Math.Pow(2, level_max) + 1;
            }

            for (point = 0; point < point_num; point++)
            {
                for (dim = 0; dim < dim_num; dim++)
                {
                    if (rule == 1)
                    {
                        grid_point[dim + point * dim_num] =
                            ClenshawCurtis.cc_abscissa(order_max, grid_index[dim + point * dim_num] + 1);
                    }
                }
            }

            //
            //  Gather the weights.
            //
            sparse_grid_weights_cfn(dim_num, level_max, rule, point_num, grid_index,
                ref grid_weight);
        }

        public static void sparse_grid_ofn(int dim_num, int level_max, int rule, int point_num,
                ref double[] grid_weight, ref double[] grid_point)

            //***************************************************************************80 
            // 
            //  Purpose:
            //
            //    SPARSE_GRID_OFN computes a sparse grid based on an OFN 1D rule. 
            // 
            //  Discussion: 
            //
            //    The 1D quadrature rule is assumed to be Open Fully Nested.
            //
            //    Open Fully Nested rules include Fejer 1, Fejer 2, and Gauss Patterson rules.
            //
            //    A Smolyak construction is used to create a multidimensional sparse grid.  
            // 
            //    The user specifies: 
            //    * the spatial dimension of the quadrature region, 
            //    * the level that defines the Smolyak grid. 
            //    * the 1D quadrature rule. 
            // 
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 March 2008
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
            //    Input, int LEVEL_MAX, controls the size of the final 
            //    sparse grid. 
            // 
            //    Input, int RULE, the index of the rule.
            //    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
            //    2, "F1", Fejer 1 Open Fully Nested rule.
            //    3, "F2", Fejer 2 Open Fully Nested rule.
            //    4, "GP", Gauss Patterson Open Fully Nested rule.
            //    5, "GL", Gauss Legendre Open Weakly Nested rule.
            //    6, "GH", Gauss Hermite Open Weakly Nested rule.
            //    7, "LG", Gauss Laguerre Open Non Nested rule.
            //
            //    Input, int POINT_NUM, the number of points in the grid, 
            //    as determined by LEVELS_INDEX_SIZE.
            //
            //    Output, double GRID_WEIGHT[POINT_NUM], the weights. 
            // 
            //    Output, double GRID_POINT[DIM_NUM*POINT_NUM], the points. 
            // 
        {
            int dim;
            int[] grid_base;
            int[] grid_index;
            int order_max;
            int point;

            if (rule < 2 || 4 < rule)
            {
                Console.WriteLine("");
                Console.WriteLine("SPARSE_GRID_OFN - Fatal error!");
                Console.WriteLine("  Illegal input rule index = " + rule + "");
                return;
            }

            // 
            //  Determine the index vector, relative to the full product grid, 
            //  that identifies the points in the sparse grid. 
            //
            grid_base = new int[dim_num * point_num];
            grid_index = new int[dim_num * point_num];

            levels_index_ofn(dim_num, level_max, point_num, ref grid_index, ref grid_base);
            // 
            //  Compute the physical coordinates of the abscissas. 
            //
            order_max = (int) Math.Pow(2, level_max + 1) - 1;

            for (point = 0; point < point_num; point++)
            {
                for (dim = 0; dim < dim_num; dim++)
                {
                    if (rule == 2)
                    {
                        grid_point[dim + point * dim_num] =
                            Fejer1.f1_abscissa(order_max, grid_index[dim + point * dim_num]);
                    }
                    else if (rule == 3)
                    {
                        grid_point[dim + point * dim_num] =
                            Fejer2.f2_abscissa(order_max, grid_index[dim + point * dim_num]);
                    }
                    else if (rule == 4)
                    {
                        grid_point[dim + point * dim_num] =
                            PattersonQuadrature.gp_abscissa(order_max, grid_index[dim + point * dim_num]);
                    }
                }
            }

            // 
            //  Gather the weights. 
            //
            sparse_grid_weights_ofn(dim_num, level_max, rule, point_num,
                grid_index, ref grid_weight);

        }

        public static int sparse_grid_ofn_size(int dim_num, int level_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPARSE_GRID_OFN_SIZE sizes a sparse grid using Open Fully Nested rules.
            //
            //  Discussion:
            //
            //    The grid is defined as the sum of the product rules whose LEVEL
            //    satisfies:
            //
            //      0 <= LEVEL <= LEVEL_MAX.
            //
            //    This calculation is much faster than a previous method.  It simply
            //    computes the number of new points that are added at each level in the
            //    1D rule, and then counts the new points at a given DIM_NUM dimensional
            //    level vector as the product of the new points added in each dimension.
            //
            //    This approach will work for nested families, and may be extensible
            //    to other families, and to mixed rules.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 December 2009
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
            //    Output, int SPARSE_GRID_CC_SIZE, the number of points in the grid.
            //
        {
            int dim;
            int h;
            int l;
            int level;
            int[] level_1d;
            bool more;
            int[] new_1d;
            int point_num;
            int t;
            int v;
            //
            //  Special case.
            //
            if (level_max < 0)
            {
                point_num = 0;
                return point_num;
            }

            if (level_max == 0)
            {
                point_num = 1;
                return point_num;
            }

            //
            //  Construct the vector that counts the new points in the 1D rule.
            //
            new_1d = new int[level_max + 1];

            new_1d[0] = 1;
            for (l = 1; l <= level_max; l++)
            {
                new_1d[l] = 2 * new_1d[l - 1];
            }

            //
            //  Count the number of points by counting the number of new points 
            //  associated with each level vector.
            //
            level_1d = new int[dim_num];

            point_num = 0;

            for (level = 0; level <= level_max; level++)
            {
                more = false;
                h = 0;
                t = 0;

                for (;;)
                {
                    Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);

                    v = 1;
                    for (dim = 0; dim < dim_num; dim++)
                    {
                        v = v * new_1d[level_1d[dim]];
                    }

                    point_num = point_num + v;

                    if (!more)
                    {
                        break;
                    }
                }
            }

            return point_num;
        }

        public static void sparse_grid_onn(int dim_num, int level_max, int rule, int point_num,
                ref double[] grid_weight, ref double[] grid_point)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPARSE_GRID_ONN computes a sparse grid based on a ONN 1D rule.
            //
            //  Discussion:
            //
            //    The 1D quadrature rule is assumed to be Open Non-Nested.
            //    Such rules include Gauss Laguerre rules.
            //
            //    A Smolyak construction is used to create a multidimensional sparse grid.  
            //
            //    The user specifies: 
            //    * the spatial dimension of the quadrature region, 
            //    * the level that defines the Smolyak grid. 
            //    * the quadrature rule;
            //    * the number of points in the rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 July 2008
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
            //    Input, int RULE, the index of the rule.
            //    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
            //    2, "F1", Fejer 1 Open Fully Nested rule.
            //    3, "F2", Fejer 2 Open Fully Nested rule.
            //    4, "GP", Gauss Patterson Open Fully Nested rule.
            //    5, "GL", Gauss Legendre Open Weakly Nested rule.
            //    6, "GH", Gauss Hermite Open Weakly Nested rule.
            //    7, "LG", Gauss Laguerre Open Non Nested rule.
            //
            //    Input, int POINT_NUM, the number of points in the grid, as determined
            //    by SPARSE_GRID_SIZE_ONN.
            //
            //    Output, double GRID_WEIGHT[POINT_NUM], the weights.
            //
            //    Output, double GRID_POINT[DIM_NUM*POINT_NUM], the points.
            //
        {
            double coeff;
            int dim;
            int[] grid_base2;
            int[] grid_index2;
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
                    //
                    LevelToOrder.level_to_order_open(dim_num, level_1d, ref order_1d);

                    for (dim = 0; dim < dim_num; dim++)
                    {
                        grid_base2[dim] = order_1d[dim];
                    }

                    //
                    //  The product of the 1D orders gives us the number of points in this grid.
                    //
                    order_nd = typeMethods.i4vec_product(dim_num, order_1d);
                    //
                    //  Compute the weights for this product grid.
                    //
                    grid_weight2 = Product.product_weights(dim_num, order_1d, order_nd, rule);
                    //
                    //  Now determine the coefficient of the weight.
                    //
                    coeff = typeMethods.r8_mop(level_max - level)
                            * typeMethods.r8_choose(dim_num - 1, level_max - level);
                    //
                    //  The inner (hidden) loop generates all points corresponding to given grid.
                    //  The grid indices will be between -M to +M, where 2*M + 1 = ORDER_1D(DIM).
                    //
                    grid_index2 = Multigrid.multigrid_index_onn(dim_num, order_1d, order_nd);

                    for (point = 0; point < order_nd; point++)
                    {

                        if (point_num <= point_num2)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("SPARSE_GRID_ONN - Fatal error!");
                            Console.WriteLine("  Exceeding maximum point index POINT_NUM = "
                                              + point_num + "");
                            return;
                        }

                        Legendre.QuadratureRule.lg_abscissa(dim_num, 1, grid_index2,
                            grid_base2, ref grid_point, gridIndex: + point * dim_num, gridPtIndex: + point_num2 * dim_num);

                        grid_weight[point_num2] = coeff * grid_weight2[point];

                        point_num2 = point_num2 + 1;
                    }

                    if (!more)
                    {
                        break;
                    }
                }
            }

            if (point_num2 < point_num)
            {
                Console.WriteLine("");
                Console.WriteLine("SPARSE_GRID_ONN - Fatal error!");
                Console.WriteLine("  Set fewer points than POINT_NUM = " + point_num + "");
            }

        }

        public static void sparse_grid_own(int dim_num, int level_max, int rule, int point_num,
                ref double[] grid_weight, ref double[] grid_point)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPARSE_GRID_OWN computes a sparse grid based on an OWN 1D rule.
            //
            //  Discussion:
            //
            //    The 1D quadrature rule is assumed to be Open Weakly Nested.
            //    Such rules include Gauss Hermite and Gauss Legendre rules.
            //
            //    A Smolyak construction is used to create a multidimensional sparse grid.  
            // 
            //    The user specifies: 
            //    * the spatial dimension of the quadrature region, 
            //    * the level that defines the Smolyak grid,
            //    * the rule;
            //    * the number of points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 July 2008
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
            //    Input, int RULE, the index of the rule.
            //    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
            //    2, "F1", Fejer 1 Open Fully Nested rule.
            //    3, "F2", Fejer 2 Open Fully Nested rule.
            //    4, "GP", Gauss Patterson Open Fully Nested rule.
            //    5, "GL", Gauss Legendre Open Weakly Nested rule.
            //    6, "GH", Gauss Hermite Open Weakly Nested rule.
            //    7, "LG", Gauss Laguerre Open Non Nested rule.
            //
            //    Input, int POINT_NUM, the number of points in the grid, as determined
            //    by LEVELS_INDEX_SIZE_OWN.
            //
            //    Output, double GRID_WEIGHT[POINT_NUM], the weights.
            //
            //    Output, double GRID_POINT[DIM_NUM*POINT_NUM], the points.
            //
        {
            double coeff;
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
            int level_min2;
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

            if (dim_num == 1)
            {
                level_min2 = level_min;
            }
            else
            {
                level_min2 = 0;
            }

            grid_base2 = new int[dim_num];
            level_1d = new int[dim_num];
            order_1d = new int[dim_num];

            for (level = level_min2; level <= level_max; level++)
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
                    LevelToOrder.level_to_order_open(dim_num, level_1d, ref order_1d);

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
                    grid_weight2 = Product.product_weights(dim_num, order_1d, order_nd, rule);
                    //
                    //  Now determine the coefficient of the weight.
                    //
                    coeff = typeMethods.r8_mop(level_max - level)
                            * typeMethods.r8_choose(dim_num - 1, level_max - level);
                    //
                    //  The inner (hidden) loop generates all points corresponding to given grid.
                    //  The grid indices will be between -M to +M, where 2*M + 1 = ORDER_1D(DIM).
                    //
                    grid_index2 = Multigrid.multigrid_index_own(dim_num, order_1d, order_nd);
                    //
                    //  Determine the first level of appearance of each of the points.
                    //  This allows us to flag certain points as being repeats of points
                    //  generated on a grid of lower level.  
                    //
                    //  This is SLIGHTLY tricky.
                    //
                    grid_level = LevelToOrder.index_level_own(level, level_max, dim_num, order_nd,
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

                            if (point_num <= point_num2)
                            {
                                Console.WriteLine("");
                                Console.WriteLine("SPARSE_GRID_OWN - Fatal error!");
                                Console.WriteLine("  Exceeding maximum point index POINT_NUM = "
                                                  + point_num + "");
                                return;
                            }

                            if (rule == 5)
                            {
                                GaussQuadrature.gl_abscissa(dim_num, 1, grid_index2,
                                    grid_base2, ref grid_point, gridIndex: + point * dim_num, gridPointIndex: + point_num2 * dim_num);
                            }
                            else if (rule == 6)
                            {
                                GaussHermite.gh_abscissa(dim_num, 1, grid_index2,
                                    grid_base2, ref grid_point, gridIndIndex: + point * dim_num, gridPtIndex: + point_num2 * dim_num);
                            }
                            else
                            {
                                Console.WriteLine("");
                                Console.WriteLine("SPARSE_GRID_OWN - Fatal error!");
                                Console.WriteLine("  Unrecognized rule number = " + rule + "");
                                return;
                            }

                            if (level_min <= level)
                            {
                                grid_weight[point_num2] = coeff * grid_weight2[point];
                            }

                            point_num2 = point_num2 + 1;
                        }
                        //
                        //  or an already existing point (create point temporarily, find match,
                        //  add weight to matched point's weight).
                        //
                        else
                        {
                            if (level_min <= level)
                            {
                                grid_point_temp = new double[dim_num];

                                if (rule == 5)
                                {
                                    GaussQuadrature.gl_abscissa(dim_num, 1, grid_index2,
                                        grid_base2, ref grid_point_temp, gridIndex: + point * dim_num);
                                }
                                else if (rule == 6)
                                {
                                    GaussHermite.gh_abscissa(dim_num, 1, grid_index2,
                                        grid_base2, ref grid_point_temp, gridIndIndex: + point * dim_num);
                                }
                                else
                                {
                                    Console.WriteLine("");
                                    Console.WriteLine("SPARSE_GRID_OWN - Fatal error!");
                                    Console.WriteLine("  Unrecognized rule number = " + rule + "");
                                    return;
                                }

                                for (point2 = 0; point2 < point_num2; point2++)
                                {
                                    point3 = point2;
                                    for (dim = 0; dim < dim_num; dim++)
                                    {
                                        if (grid_point[dim + point2 * dim_num] != grid_point_temp[dim])
                                        {
                                            point3 = -1;
                                            break;
                                        }
                                    }

                                    if (point3 == point2)
                                    {
                                        break;
                                    }
                                }

                                if (point3 == -1)
                                {
                                    Console.WriteLine("");
                                    Console.WriteLine("SPARSE_GRID_OWN - Fatal error!");
                                    Console.WriteLine("  Could not match point.");
                                    return;
                                }

                                grid_weight[point3] = grid_weight[point3] +
                                                      (double) (coeff) * grid_weight2[point];
                            }
                        }
                    }

                    if (!more)
                    {
                        break;
                    }
                }
            }

            if (point_num2 < point_num)
            {
                Console.WriteLine("");
                Console.WriteLine("SPARSE_GRID_OWN - Fatal error!");
                Console.WriteLine("  Set fewer points than POINT_NUM = " + point_num + "");
            }

        }

        public static void sparse_grid_weights_cfn(int dim_num, int level_max, int rule,
                int point_num, int[] grid_index, ref double[] grid_weight)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPARSE_GRID_WEIGHTS_CFN computes sparse grid weights based on a CFN 1D rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 July 2008
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
            //    Input, int RULE, the index of the rule.
            //    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
            //    2, "F1", Fejer 1 Open Fully Nested rule.
            //    3, "F2", Fejer 2 Open Fully Nested rule.
            //    4, "GP", Gauss Patterson Open Fully Nested rule.
            //    5, "GL", Gauss Legendre Open Weakly Nested rule.
            //    6, "GH", Gauss Hermite Open Weakly Nested rule.
            //    7, "LG", Gauss Laguerre Open Non Nested rule.
            //
            //    Input, int POINT_NUM, the total number of points in the grids.
            //
            //    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], a list of point indices,
            //    representing a subset of the product grid of level LEVEL_MAX,
            //    representing (exactly once) each point that will show up in a
            //    sparse grid of level LEVEL_MAX.
            //
            //    Output, double GRID_WEIGHT[POINT_NUM], the weights
            //    associated with the sparse grid points.
            //
        {
            bool all_equal;
            double coeff;
            int dim;
            int[] grid_index2;
            double[] grid_weight2;
            int h;
            int level;
            int[] level_1d;
            int level_min;
            bool more;
            int order_nd;
            int[] order_1d;
            int point;
            int point2;
            int t;

            if (level_max == 0)
            {
                for (point = 0; point < point_num; point++)
                {
                    grid_weight[point] = Math.Pow(2, dim_num);
                }

                return;
            }

            level_1d = new int[dim_num];
            order_1d = new int[dim_num];

            for (point = 0; point < point_num; point++)
            {
                grid_weight[point] = 0.0;
            }

            level_min = Math.Max(0, level_max + 1 - dim_num);

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
                    LevelToOrder.level_to_order_closed(dim_num, level_1d, ref order_1d);
                    //
                    //  The product of the 1D orders gives us the number of points in this grid.
                    //
                    order_nd = typeMethods.i4vec_product(dim_num, order_1d);
                    //
                    //  Generate the indices of the points corresponding to the grid.
                    //
                    grid_index2 = Multigrid.multigrid_index_cfn(dim_num, order_1d, order_nd);
                    //
                    //  Compute the weights for this grid.
                    //
                    grid_weight2 = Product.product_weights(dim_num, order_1d, order_nd, rule);
                    //
                    //  Adjust the grid indices to reflect LEVEL_MAX.
                    //
                    Multigrid.multigrid_scale_closed(dim_num, order_nd, level_max, level_1d,
                        ref grid_index2);
                    //
                    //  Now determine the coefficient.
                    //
                    coeff = typeMethods.r8_mop(level_max - level)
                            * typeMethods.r8_choose(dim_num - 1, level_max - level);

                    for (point2 = 0; point2 < order_nd; point2++)
                    {
                        for (point = 0; point < point_num; point++)
                        {
                            all_equal = true;
                            for (dim = 0; dim < dim_num; dim++)
                            {
                                if (grid_index2[dim + point2 * dim_num] !=
                                    grid_index[dim + point * dim_num])
                                {
                                    all_equal = false;
                                    break;
                                }
                            }

                            if (all_equal)
                            {
                                grid_weight[point] = grid_weight[point]
                                                     + coeff * grid_weight2[point2];
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

        public static void sparse_grid_weights_ofn(int dim_num, int level_max, int rule,
                int point_num, int[] grid_index, ref double[] grid_weight)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPARSE_GRID_WEIGHTS_OFN computes sparse grid weights based on a OFN 1D rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 July 2008
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
            //    Input, int RULE, the index of the rule.
            //    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
            //    2, "F1", Fejer 1 Open Fully Nested rule.
            //    3, "F2", Fejer 2 Open Fully Nested rule.
            //    4, "GP", Gauss Patterson Open Fully Nested rule.
            //    5, "GL", Gauss Legendre Open Weakly Nested rule.
            //    6, "GH", Gauss Hermite Open Weakly Nested rule.
            //    7, "LG", Gauss Laguerre Open Non Nested rule.
            //
            //    Input, int POINT_NUM, the total number of points in the grids.
            //
            //    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], a list of point indices,
            //    representing a subset of the product grid of level LEVEL_MAX,
            //    representing (exactly once) each point that will show up in a
            //    sparse grid of level LEVEL_MAX.
            //
            //    Output, double GRID_WEIGHT[POINT_NUM], the weights
            //    associated with the sparse grid points.
            //
        {
            bool all_equal;
            double coeff;
            int dim;
            int[] grid_index2;
            double[] grid_weight2;
            int h;
            int level;
            int[] level_1d;
            int level_min;
            bool more;
            int order_nd;
            int[] order_1d;
            int point;
            int point2;
            int t;

            if (level_max == 0)
            {
                for (point = 0; point < point_num; point++)
                {
                    grid_weight[point] = Math.Pow(2, dim_num);
                }

                return;
            }

            level_1d = new int[dim_num];
            order_1d = new int[dim_num];

            for (point = 0; point < point_num; point++)
            {
                grid_weight[point] = 0.0;
            }

            level_min = Math.Max(0, level_max + 1 - dim_num);

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
                    LevelToOrder.level_to_order_open(dim_num, level_1d, ref order_1d);
                    //
                    //  The product of the 1D orders gives us the number of points in this grid.
                    //
                    order_nd = typeMethods.i4vec_product(dim_num, order_1d);
                    //
                    //  Generate the indices of the points corresponding to the grid.
                    //
                    grid_index2 = Multigrid.multigrid_index_ofn(dim_num, order_1d, order_nd);
                    //
                    //  Compute the weights for this grid.
                    //
                    grid_weight2 = Product.product_weights(dim_num, order_1d, order_nd, rule);
                    //
                    //  Adjust the grid indices to reflect LEVEL_MAX.
                    //
                    Multigrid.multigrid_scale_open(dim_num, order_nd, level_max, ref level_1d,
                        ref grid_index2);
                    //
                    //  Now determine the coefficient.
                    //
                    coeff = typeMethods.r8_mop(level_max - level)
                            * typeMethods.r8_choose(dim_num - 1, level_max - level);

                    for (point2 = 0; point2 < order_nd; point2++)
                    {
                        for (point = 0; point < point_num; point++)
                        {
                            all_equal = true;
                            for (dim = 0; dim < dim_num; dim++)
                            {
                                if (grid_index2[dim + point2 * dim_num] !=
                                    grid_index[dim + point * dim_num])
                                {
                                    all_equal = false;
                                    break;
                                }
                            }

                            if (all_equal)
                            {
                                grid_weight[point] = grid_weight[point]
                                                     + coeff * grid_weight2[point2];
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

        public static int[] spgrid_open_index(int dim_num, int level_max, int point_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEVELS_OPEN_INDEX computes open grids with 0 <= LEVEL <= LEVEL_MAX.
            //
            //  Discussion:
            //
            //    The necessary dimensions of GRID_INDEX can be
            //    determined by calling SPGRID_OPEN_SIZE first.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 July 2008
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
            //    Output, int LEVELS_MAX_INDEX[DIM_NUM*POINT_NUM], a list of point indices,
            //    representing a subset of the product grid of level LEVEL_MAX,
            //    representing (exactly once) each point that will show up in a
            //    sparse grid of level LEVEL_MAX.
            //
        {
            int dim;
            int[] grid_index;
            int[] grid_index2;
            int h;
            int level;
            int[] level_1d;
            bool more;
            int[] order_1d;
            int order_nd;
            int point;
            int point_num2;
            int t;
            bool test;
            //
            //  The outer loop generates LEVELs from 0 to LEVEL_MAX.
            //
            grid_index = new int[dim_num * point_num];
            level_1d = new int[dim_num];
            order_1d = new int[dim_num];

            point_num2 = 0;

            for (level = 0; level <= level_max; level++)
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
                    LevelToOrder.level_to_order_open(dim_num, level_1d, ref order_1d);
                    //
                    //  The product of the 1D orders gives us the number of points in this grid.
                    //
                    order_nd = typeMethods.i4vec_product(dim_num, order_1d);
                    //
                    //  The inner (hidden) loop generates all points corresponding to given grid.
                    //
                    grid_index2 = Multigrid.multigrid_index1(dim_num, order_1d, order_nd);
                    //
                    //  Only keep those points which first appear on this level.
                    //  If you keep a point, it is necessary to rescale each of its components
                    //  so that we save the coordinates as they apply on the final grid.
                    //
                    for (point = 0; point < order_nd; point++)
                    {
                        test = true;
                        for (dim = 0; dim < dim_num; dim++)
                        {
                            if (grid_index2[dim + point * dim_num] % 2 == 0)
                            {
                                test = false;
                            }
                        }

                        if (test)
                        {
                            for (dim = 0; dim < dim_num; dim++)
                            {
                                grid_index[dim + point_num2 * dim_num] =
                                    (int) Math.Pow(2, level_max - level_1d[dim])
                                    * grid_index2[dim + point * dim_num];
                            }

                            point_num2 = point_num2 + 1;
                        }
                    }

                    if (!more)
                    {
                        break;
                    }
                }
            }

            return grid_index;
        }

        public static double[] spgrid_open_weights(int dim_num, int level_max, int point_num,
                int[] grid_index, int rule)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPGRID_OPEN_WEIGHTS gathers the weights.
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
            //    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], a list of point indices,
            //    representing a subset of the product grid of level LEVEL_MAX,
            //    representing (exactly once) each point that will show up in a
            //    sparse grid of level LEVEL_MAX.
            //
            //    Input, int RULE, the 1D quadrature rule being used.
            //    2, Fejer Type 2 Rule;
            //    3, Gauss-Patterson Rule,
            //    4, Newton-Cotes Open Rule,
            //    5, Newton-Cotes Open Half Rule.
            //
            //    Output, double SPGRID_OPEN_WEIGHTS[POINT_NUM], the weights
            //    associated with the sparse grid points.
            //
        {
            bool all_equal;
            int coeff;
            int dim;
            int[] grid_index2;
            double[] grid_weight;
            double[] grid_weight2;
            int h;
            int level;
            int[] level_1d;
            int level_min;
            int match;
            bool more;
            int order_nd;
            int[] order_1d;
            int point;
            int point2;
            int t;

            grid_weight = new double[point_num];

            if (level_max == 0)
            {
                for (point = 0; point < point_num; point++)
                {
                    grid_weight[point] = (int) Math.Pow(2, dim_num);
                }

                return grid_weight;
            }

            for (point = 0; point < point_num; point++)
            {
                grid_weight[point] = 0.0;
            }

            level_1d = new int[dim_num];
            order_1d = new int[dim_num];

            level_min = Math.Max(0, level_max + 1 - dim_num);

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
                    LevelToOrder.level_to_order_open(dim_num, level_1d, ref order_1d);
                    //
                    //  The product of the 1D orders gives us the number of points in this grid.
                    //
                    order_nd = typeMethods.i4vec_product(dim_num, order_1d);
                    //
                    //  Generate the indices of the points corresponding to the grid.
                    //
                    grid_index2 = Multigrid.multigrid_index1(dim_num, order_1d, order_nd);
                    //
                    //  Compute the weights for this grid.
                    //
                    grid_weight2 = Product.product_weights_open(dim_num, order_1d, order_nd, rule);
                    //
                    //  Adjust the grid indices to reflect LEVEL_MAX.
                    //
                    Multigrid.multigrid_scale_open(dim_num, order_nd, level_max, ref level_1d,
                        ref grid_index2);
                    //
                    //  Now determine the coefficient.
                    //
                    coeff = (int) Math.Pow(-1, level_max - level)
                            * typeMethods.i4_choose(dim_num - 1, level_max - level);

                    for (point2 = 0; point2 < order_nd; point2++)
                    {
                        match = -1;

                        for (point = 0; point < point_num; point++)
                        {
                            all_equal = true;
                            for (dim = 0; dim < dim_num; dim++)
                            {
                                if (grid_index2[dim + point2 * dim_num]
                                    != grid_index[dim + point * dim_num])
                                {
                                    all_equal = false;
                                    break;
                                }
                            }

                            if (all_equal)
                            {
                                grid_weight[point] = grid_weight[point]
                                                     + (double) (coeff) * grid_weight2[point2];
                                match = point;
                                break;
                            }
                        }

                        if (match == -1)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("SPGRID_OPEN_WEIGHTS - Fatal error!");
                            Console.WriteLine("  Could not match grid index.");
                            Console.WriteLine("  Point index = " + point2 + "");
                            Console.WriteLine("");
                            Console.WriteLine("  LEVEL = " + level + "");
                            Console.WriteLine("");
                            Console.WriteLine("  LEVEL_1D:");
                            string cout = "";
                            for (dim = 0; dim < dim_num; dim++)
                            {
                                cout += level_1d[dim].ToString().PadLeft(6);
                            }

                            Console.WriteLine(cout);
                            Console.WriteLine("");
                            Console.WriteLine("  ORDER_1D:");
                            cout = "";
                            for (dim = 0; dim < dim_num; dim++)
                            {
                                cout += order_1d[dim].ToString().PadLeft(6);
                            }

                            Console.WriteLine(cout);
                            Console.WriteLine("");
                            Console.WriteLine("  GRID_INDEX2");
                            cout = "";
                            for (dim = 0; dim < dim_num; dim++)
                            {
                                cout += grid_index2[dim + point2 * dim_num].ToString().PadLeft(6);
                            }

                            Console.WriteLine(cout);
                            return grid_weight;
                        }
                    }

                    if (!more)
                    {
                        break;
                    }
                }
            }

            return grid_weight;
        }

    }
}