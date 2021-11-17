using System;
using Burkardt.Composition;
using Burkardt.Types;

namespace Burkardt.ClenshawCurtisNS;

public static class ClenshawCurtisGrid
{
    public static double cc_abscissa(int order, int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_ABSCISSA returns the I-th abscissa of the Clenshaw Curtis rule.
        //
        //  Discussion:
        //
        //    Our convention is that the abscissas are numbered from left to
        //    right.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 November 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the Clenshaw Curtis rule.
        //
        //    Input, int I, the index of the desired abscissa.  1 <= I <= ORDER.
        //
        //    Output, double CC_ABSCISSA, the value of the I-th 
        //    abscissa in the Clenshaw Curtis rule of order ORDER.
        //
    {
            
        double value = 0;

        switch (order)
        {
            case < 1:
                value = -typeMethods.r8_huge();
                return value;
            case 1:
                value = 0.0;
                return value;
        }

        switch (i)
        {
            case < 1:
                value = -1.0;
                break;
            default:
            {
                if (i <= order)
                {
                    value = Math.Cos((order - i) * Math.PI
                                     / (order - 1));
                }
                else
                {
                    value = 1.0;
                }

                break;
            }
        }

        return value;
    }

    public static int[] cc_abscissa_level_1d(int level_max, int test_num, int[] test_val)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_ABSCISSA_LEVEL_1D: first level at which given abscissa is generated.
        //
        //  Discussion:
        //
        //    The Clenshaw Curtis abscissas are the cosines of angles, and hence
        //    are somewhat raggedy real values between -1 and 1.  However,
        //    it is sometimes convenient to think of them as being equally
        //    spaced integers, particularly if we are mainly concerned with the
        //    nesting property.  The question this routine can answer is, if
        //    we are going to generate a sequence of nested Clenshaw Curtis
        //    rules up to a given maximum level, on what level does a particular
        //    abscissa first show up?
        //
        //    Consider the sequence of numbers from 0 to 2**LEVEL_MAX.  Suppose
        //    that we organize them into levels.  The first two levels are
        //    somewhat arbitrary:
        //
        //      Level 0:  2**(LEVEL_MAX-1)
        //      Level 1:  0, 2**LEVEL_MAX
        //
        //    But after that, the generation rule for the next level is simply
        //    to generate a new value BETWEEN every consecutive pair of values
        //    that have already been generated ( and these new values are simply
        //    the averages of the consecutive pair.)
        //
        //    Thus, if our current set of values is 0, 8 and 16, on the next
        //    level we generate 4 and 12 to make a new set of 0, 4, 8, 12, 16.
        //    We continue this operation, level by level, until we have filled in
        //    all 2**LEVEL_MAX+1 values.
        //
        //    This routine returns, for any value I and maximum level LEVEL_MAX,
        //    the level on which the value I will first be produced.
        //
        //    For example, for LEVEL_MAX = 5, the numbers we are considering
        //    are 0 through 32, and they will be produced as follows:
        //
        //      Level
        //      0:                                 16
        //      1:  0                                                              32
        //      2:                  8                              24
        //      3:          4              12              20              28
        //      4:      2       6      10      14      18      22      26      30
        //      5:    1   3   5   7   9  11  13  15  17  19  21  23  25  27  29  31
        //
        //    Here is the list of levels for 0 through 32:
        //
        //          1 5 4 5 3 5 4 5 2 5 4 5 3 5 4 5 0 5 4 5 3 5 4 5 2 5 4 5 3 5 4 5 1
        //
        //    The purpose of this routine is, given the value 20 and the
        //    maximum level 5, to return level = 3, indicating that the value 20
        //    will first be generated on the 3rd level for a grid that ultimately
        //    reaches an order of 2**5+1 values.
        //
        //    The need for this routine arises from the necessity of understanding
        //    nested Clenshaw Curtis grids.  In particular, if we see a grid of
        //    17 points, this is the fifth in a series of nested grids, and 8 of
        //    the points are new, created specifically for the level 4 grid, while
        //    9 of the points arose earlier.  This routine can report exactly when
        //    each value was created.
        //
        //    The real need for this routine arises in multidimensional sparse grids,
        //    where we essentially have a fixed "budget" of levels we are allowed to
        //    use.  When we generate a multidimensional point, we determine its
        //    level in each single dimensional grid, add them up, and this value
        //    must be no greater than our budgeted value for the point to be included
        //    in the sparse grid.
        //
        //    Except for the behavior of the first two levels, it is true that
        //    the level of a value I is LEVEL_MAX minus the number of times I can be
        //    divided evenly by 2.  Because of a peculiarity of the definition of
        //    the grids, if this formula gives a level of 0 or 1, then the level
        //    should be replaced by 1 or 0, respectively.
        //
        //    Again, except for the first two levels, the calculation is equivalent
        //    to computing the location of the "first" nonzero bit in the representation
        //    of a number, and subtracting that from LEVEL_MAX.  This is why all the
        //    odd numbers, which have their first 1 bit in the 0-th position,
        //    are assigned a level of LEVEL_MAX.
        //
        //    This routine can also be called for values that lie outside the standard
        //    range of 0 through 2**LEVEL_MAX.  In that case, a MOD operation is
        //    applied first, to make a sensible result.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int LEVEL_MAX, determines the number of points in the grid as
        //    2**LEVEL_MAX + 1.
        //
        //    Input, int TEST_NUM, the number of points to be tested.
        //
        //    Input, int TEST_VAL[TEST_NUM], the values of the points to be tested.
        //    Normally, each value would be between 0 and 2**LEVEL_MAX+1.
        //
        //    Output, int CC_ABSCISSA_LEVEL_1D[TEST_NUM], the level at which the
        //    point would first be generated, assuming that a standard sequence of
        //    nested grids is used.
        //
    {
        int i;
        int level;
        int order;
        int t;
        int[] test_level;

        order = (int) Math.Pow(2, level_max) + 1;

        test_level = new int[test_num];

        for (i = 0; i < test_num; i++)
        {
            t = test_val[i];
            //
            //  The following MOD operation is only needed to handle cases where
            //  T is not in the expected range.
            //
            t %= order;

            switch (t)
            {
                case 0:
                    level = 0;
                    break;
                default:
                {
                    level = level_max;

                    while (t % 2 == 0)
                    {
                        t /= 2;
                        level -= 1;
                    }

                    break;
                }
            }

            level = level switch
            {
                0 => 1,
                1 => 0,
                _ => level
            };

            test_level[i] = level;
        }

        return test_level;
    }

    public static int[] cc_abscissa_level_nd(int level_max, int dim_num, int test_num,
            int[] test_val)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_ABSCISSA_LEVEL_ND: first level at which given abscissa is generated.
        //
        //  Discussion:
        //
        //    We assume an underlying product grid.  In each dimension, this product
        //    grid has order 2**LEVEL_MAX + 1.
        //
        //    We will say a sparse grid has total level LEVEL if each point in the
        //    grid has a total level of LEVEL or less.
        //
        //    The "level" of a point is determined as the sum of the levels of the
        //    point in each spatial dimension.
        //
        //    The level of a point in a single spatial dimension I is determined as
        //    the level, between 0 and LEVEL_MAX, at which the point's I'th index
        //    would have been generated.
        //
        //
        //    This description is terse and perhaps unenlightening.  Keep in mind
        //    that the product grid is the product of 1D Clenshaw Curtis grids,
        //    that the 1D Clenshaw Curtis grids are built up by levels, having
        //    orders (total number of points ) 1, 3, 5, 9, 17, 33 and so on,
        //    and that these 1D grids are nested, so that each point in a 1D grid
        //    has a first level at which it appears.
        //
        //    Our procedure for generating the points of a sparse grid, then, is
        //    to choose a value LEVEL_MAX, to generate the full product grid,
        //    but then only to keep those points on the full product grid whose
        //    LEVEL is less than or equal to LEVEL_MAX.  
        //
        //
        //    Note that this routine is really just testing out the idea of
        //    determining the level.  Our true desire is to be able to start
        //    with a value LEVEL, and determine, in a straightforward manner,
        //    all the points that are generated exactly at that level, or
        //    all the points that are generated up to and including that level.
        //
        //    This allows us to generate the new points to be added to one sparse
        //    grid to get the next, or to generate a particular sparse grid at once.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int LEVEL_MAX, controls the size of the final sparse grid.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int TEST_NUM, the number of points to be tested.
        //
        //    Input, int TEST_VAL[DIM_NUM*TEST_NUM], the indices of the points 
        //    to be tested.  Normally, each index would be between 0 and 2**LEVEL_MAX.
        //
        //    Output, int CC_ABSCISSA_LEVEL_ND[TEST_NUM], the value of LEVEL at which the
        //    point would first be generated, assuming that a standard sequence of
        //    nested grids is used.
        //
    {
        int dim;
        int i;
        int level;
        int order;
        int t;
        int[] test_level;

        order = (int) Math.Pow(2, level_max) + 1;

        test_level = new int[test_num];

        for (i = 0; i < test_num; i++)
        {
            test_level[i] = 0;
            for (dim = 0; dim < dim_num; dim++)
            {
                t = test_val[dim + i * dim_num];

                t %= order;

                switch (t)
                {
                    case 0:
                        level = 0;
                        break;
                    default:
                    {
                        level = level_max;

                        while (t % 2 == 0)
                        {
                            t /= 2;
                            level -= 1;
                        }

                        break;
                    }
                }

                level = level switch
                {
                    0 => 1,
                    1 => 0,
                    _ => level
                };

                test_level[i] += level;
            }
        }

        return test_level;
    }

    public static void cc_grid(int dim_num, int[] order_1d, int order_nd, ref double[] point, int ptIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_GRID returns a multidimensional Clenshaw-Curtis grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 November 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension of the points.
        //
        //    Input, int ORDER_1D[DIM_NUM], the order of the Clenshaw-Curtis
        //    rule in each dimension.
        //
        //    Input, int ORDER_ND, the number of points in the grid.
        //    This is equal to the product of the entries of ORDER_1D.
        //
        //    Output, double POINT[DIM_NUM*POINT_NUM], the points in
        //    the grid.  
        //
    {
        int[] a;
        int change = 0;
        int dim;
        bool done;
        int p;

        a = new int[dim_num];
        done = true;
        p = 0;

        typeMethods.VecGrayData data = new();
            
        for (;;)
        {
            typeMethods.vec_next_gray(ref data, dim_num, order_1d, ref a, ref done, ref change);

            if (done)
            {
                break;
            }

            switch (p)
            {
                case 0:
                {
                    for (dim = 0; dim < dim_num; dim++)
                    {
                        point[ptIndex + dim + p * dim_num] = cc_abscissa(order_1d[dim], 1);
                    }

                    break;
                }
                default:
                {
                    for (dim = 0; dim < dim_num; dim++)
                    {
                        point[ptIndex + dim + p * dim_num] = point[ptIndex + dim + (p - 1) * dim_num];
                    }

                    point[ptIndex + change + p * dim_num] = cc_abscissa(order_1d[change], a[change] + 1);
                    break;
                }
            }

            p += 1;
        }
    }

    public static void cc_grid_index(int dim_num, int[] order_1d, int order_nd, ref int[] indx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_GRID_INDEX returns an indexed multidimensional Clenshaw-Curtis grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension of the points.
        //
        //    Input, int ORDER_1D[DIM_NUM], the order of the Clenshaw-Curtis
        //    rule in each dimension.
        //
        //    Input, int ORDER_ND, the number of points in the grid.
        //    This is equal to the product of the entries of ORDER_1D.
        //
        //    Output, int INDX[DIM_NUM*POINT_NUM], the indices of the points in
        //    the grid.  
        //
    {
        int[] a;
        int change = 0;
        int dim;
        bool done;
        int p;

        a = new int[dim_num];
        done = true;
        p = 0;

        typeMethods.VecGrayData data = new();
            
        for (;;)
        {
            typeMethods.vec_next_gray(ref data, dim_num, order_1d, ref a, ref done, ref change);

            if (done)
            {
                break;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                indx[dim + p * dim_num] = a[dim];
            }

            p += 1;
        }
    }

    public static void cc_grids_constrained(int dim_num, double q_max, double[] alpha,
            int[] order_min, int[] order_max, int grid_num, int point_num,
            ref int[] grid_order, ref double[] grid_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_GRIDS_CONSTRAINED computes CC orders and grids satisfying a constraint.
        //
        //  Discussion:
        //
        //    The constraint on the order of the 1D Clenshaw Curtis rule in 
        //    spatial dimension I is:
        //
        //      ORDER_MIN(I) <= ORDER(I) <= ORDER_MAX(I) 
        //
        //    The constraint on the collection of orders making up a rule is:
        //
        //      Sum ( 1 <= I <= DIM_NUM ) ALPHA(I) * ORDER(I) <= Q_MAX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double Q_MAX, the maximum values of
        //    Q, the sum of the weighted orders in each spatial coordinate.
        //
        //    Input, double ALPHA[DIM_NUM], the weight factors for
        //    the orders in each spatial dimension.
        //
        //    Input, int ORDER_MIN[DIM_NUM], ORDER_MAX[DIM_NUM], the minimum
        //    and maximum values of the order of the 1D Clenshaw Curtis rule
        //    in each spatial dimension.
        //
        //    Input, int GRID_NUM, the number of Clenshaw Curtis
        //    grids in the constraint set.
        //
        //    Input, int POINT_NUM, the total number of points in the grids.
        //
        //    Output, int GRID_ORDER[DIM_NUM*GRID_NUM], contains, for each
        //    grid, the order of the Clenshaw-Curtis rule in each dimension.
        //
        //    Output, double GRID_POINT[DIM_NUM*POINT_NUM], contains
        //    a list of all the abscissas of all the rules, listed one grid at
        //    a time.  If a point occurs in several grids, it will be listed
        //    several times.
        //
    {
        int dim;
        bool more;
        int[] order_1d;
        int order_nd;

        order_1d = new int[dim_num];
        point_num = 0;
        grid_num = 0;

        more = false;

        for (;;)
        {
            typeMethods.vector_constrained_next4(dim_num, alpha, order_min, order_max,
                ref order_1d, q_max, ref more);

            if (!more)
            {
                break;
            }

            order_nd = 1;
            for (dim = 0; dim < dim_num; dim++)
            {
                order_nd *= order_1d[dim];
            }

            cc_grid(dim_num, order_1d, order_nd, ref grid_point, ptIndex: +point_num * dim_num);

            point_num += order_nd;

            for (dim = 0; dim < dim_num; dim++)
            {
                grid_order[dim + grid_num * dim_num] = order_1d[dim];
            }

            grid_num += 1;
        }
    }

    public static void cc_grids_constrained_size(int dim_num, double q_max, double[] alpha,
            int[] order_min, int[] order_max, ref int grid_num, ref int point_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_GRIDS_CONSTRAINED_SIZE counts grids for CC_GRIDS_CONSTRAINED.
        //
        //  Discussion:
        //
        //    The constraint on the order of the 1D Clenshaw Curtis rule in 
        //    spatial dimension I is:
        //
        //      ORDER_MIN(I) <= ORDER(I) <= ORDER_MAX(I) 
        //
        //    The constraint on the collection of orders making up a rule is:
        //
        //      Sum ( 1 <= I <= DIM_NUM ) ALPHA(I) * ORDER(I) <= Q_MAX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double Q_MAX, the maximum values of
        //    Q, the sum of the weighted orders in each spatial coordinate.
        //
        //    Input, double ALPHA[DIM_NUM], the weight factors for
        //    the orders in each spatial dimension.
        //
        //    Input, int ORDER_MIN[DIM_NUM], ORDER_MAX[DIM_NUM], the minimum
        //    and maximum values of the order of the 1D Clenshaw Curtis rule
        //    in each spatial dimension.
        //
        //    Output, int *GRID_NUM, the number of Clenshaw Curtis
        //    grids in the constraint set.
        //
        //    Output, int *POINT_NUM, the total number of points in the grids.
        //
    {
        int dim;
        bool more;
        int[] order_1d;
        int order_nd;

        order_1d = new int[dim_num];
        //
        //  Determine the total number of points that will be generated
        //  by "going through the motions".
        //
        point_num = 0;
        grid_num = 0;

        more = false;

        for (;;)
        {
            typeMethods.vector_constrained_next4(dim_num, alpha, order_min, order_max,
                ref order_1d, q_max, ref more);

            if (!more)
            {
                break;
            }

            order_nd = 1;
            for (dim = 0; dim < dim_num; dim++)
            {
                order_nd *= order_1d[dim];
            }

            point_num += order_nd;
            grid_num += 1;
        }
    }

    public static void cc_grids_minmax(int dim_num, int q_min, int q_max, int grid_num,
            int point_num, ref int[] grid_order, ref double[] grid_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_GRIDS_MINMAX computes CC orders and grids with Q_MIN <= Q <= Q_MAX.
        //
        //  Discussion:
        //
        //    The necessary dimensions of GRID_ORDER and GRID_POINT can be
        //    determined by calling CC_GRIDS_MINMAX first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int Q_MIN, Q_MAX, the minimum and maximum values of
        //    Q, the sum of the orders in each spatial coordinate.
        //
        //    Input, int GRID_NUM, the number of Clenshaw Curtis
        //    grids whose Q value is between Q_MIN and Q_MAX.
        //
        //    Input, int POINT_NUM, the total number of points in the grids.
        //
        //    Output, int GRID_ORDER[DIM_NUM*GRID_NUM], contains, for each
        //    grid, the order of the Clenshaw-Curtis rule in each dimension.
        //
        //    Output, double GRID_POINT[DIM_NUM*POINT_NUM], contains
        //    a list of all the abscissas of all the rules, listed one grid at
        //    a time.  If a point occurs in several grids, it will be listed
        //    several times.
        //
    {
        int dim;
        bool more;
        int[] order_1d;
        int order_nd;
        int q;
        //
        //  Outer loop generates Q's from Q_MIN to Q_MAX.
        //
        order_1d = new int[dim_num];
        point_num = 0;
        grid_num = 0;

        CompNZData data = new();
            
        for (q = q_min; q <= q_max; q++)
        {
            //
            //  Middle loop generates next partition that adds up to Q.
            //
            more = false;

            for (;;)
            {
                Comp.compnz_next(ref data, q, dim_num, ref order_1d, ref more);
                //
                //  Inner (hidden) loop generates all CC points corresponding to given grid.
                //
                order_nd = 1;
                for (dim = 0; dim < dim_num; dim++)
                {
                    order_nd *= order_1d[dim];
                }

                cc_grid(dim_num, order_1d, order_nd, ref grid_point, ptIndex: +point_num * dim_num);

                point_num += order_nd;

                for (dim = 0; dim < dim_num; dim++)
                {
                    grid_order[dim + grid_num * dim_num] = order_1d[dim];
                }

                grid_num += 1;

                if (!more)
                {
                    break;
                }
            }
        }
    }

    public static void cc_grids_minmax_size(int dim_num, int q_min, int q_max, ref int grid_num,
            ref int point_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_GRIDS_MINMAX_SIZE counts grids for CC_GRIDS_MINMAX.
        //
        //  Discussion:
        //
        //    This routine can be used to determine the necessary size to be
        //    allocated to arrays GRID_ORDER and GRID_POINT in a call to
        //    CC_GRIDS_MINMAX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int Q_MIN, Q_MAX, the minimum and maximum values of
        //    Q, the sum of the orders in each spatial coordinate.
        //
        //    Output, int *GRID_NUM, the number of Clenshaw Curtis
        //    grids whose Q value is between Q_MIN and Q_MAX.
        //
        //    Output, int *POINT_NUM, the total number of points in the grids.
        //
    {
        int dim;
        bool more;
        int[] order_1d;
        int order_nd;
        int q;
        //
        //  Determine the total number of points that will be generated
        //  by "going through the motions".
        //
        order_1d = new int[dim_num];
        point_num = 0;
        grid_num = 0;

        CompNZData data = new();
            
        for (q = q_min; q <= q_max; q++)
        {
            more = false;

            for (;;)
            {
                Comp.compnz_next(ref data, q, dim_num, ref order_1d, ref more);

                order_nd = 1;
                for (dim = 0; dim < dim_num; dim++)
                {
                    order_nd *= order_1d[dim];
                }

                point_num += order_nd;
                grid_num += 1;

                if (!more)
                {
                    break;
                }
            }
        }
    }

    public static void cc_level_to_order(int dim_num, int[] level, ref int[] order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_LEVEL_TO_ORDER converts a CC nesting level to a CC order.
        //
        //  Discussion:
        //
        //    Clenshaw Curtis grids can naturally be nested.  Except for the
        //    first case of LEVEL = 0, the relationship is
        //
        //      ORDER = 2**LEVEL + 1
        //
        //    Nesting    Order
        //    Level
        //
        //    0          1
        //    1          3 =  2 + 1
        //    2          5 =  4 + 1
        //    3          9 =  8 + 1
        //    4         17 = 16 + 1
        //    5         33 = 32 + 1
        //
        //    In this routine, we assume that a vector of levels is given,!
        //    and the corresponding orders are desired.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 November 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL[DIM_NUM], the nesting level.
        //
        //    Output, int ORDER[DIM_NUM], the order (number of points) 
        //    of the Clenshaw Curtis rule.
        //
    {
        int dim;

        for (dim = 0; dim < dim_num; dim++)
        {
            order[dim] = level[dim] switch
            {
                < 0 => -1,
                0 => 1,
                _ => (int) Math.Pow(2, level[dim]) + 1
            };
        }
    }

    public static void cc_levels_constrained(int dim_num, double q_max, double[] alpha,
            int[] level_min, int[] level_max, int grid_num, int point_num,
            ref int[] grid_level, ref double[] grid_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_LEVELS_CONSTRAINED: CC grids with constrained levels.
        //
        //  Discussion:
        //
        //    The constraint on the levels of the 1D Clenshaw Curtis rule in 
        //    spatial dimension I is:
        //
        //      LEVEL_MIN(I) <= LEVEL(I) <= LEVEL_MAX(I) 
        //
        //    The constraint on the collection of levels making up a rule is:
        //
        //      Sum ( 1 <= I <= DIM_NUM ) ALPHA(I) * LEVEL(I) <= Q_MAX.
        //
        //    The relationship of level to order is roughly 
        //
        //      ORDER = 2**LEVEL+1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 November 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double Q_MAX, the maximum values of
        //    Q, the sum of the weighted orders in each spatial coordinate.
        //
        //    Input, double ALPHA[DIM_NUM], the weight factors for
        //    the orders in each spatial dimension.
        //
        //    Input, integer LEVEL_MIN(DIM_NUM), LEVEL_MAX(DIM_NUM), the minimum
        //    and maximum values of the level of the 1D Clenshaw Curtis rule
        //    in each spatial dimension.
        //
        //    Input, int GRID_NUM, the number of Clenshaw Curtis
        //    grids in the constraint set.
        //
        //    Input, int POINT_NUM, the total number of points in the grids.
        //
        //    Output, integer GRID_LEVEL[DIM_NUM*GRID_NUM], contains, for each
        //    grid, the level of the Clenshaw-Curtis rule in each dimension.
        //
        //    Output, double GRID_POINT[DIM_NUM*POINT_NUM], contains
        //    a list of all the abscissas of all the rules, listed one grid at
        //    a time.  If a point occurs in several grids, it will be listed
        //    several times.
        //
    {
        int dim;
        int[] level_1d;
        bool more;
        int[] order_1d;
        int order_nd;

        level_1d = new int[dim_num];
        order_1d = new int[dim_num];
        point_num = 0;
        grid_num = 0;

        more = false;

        for (;;)
        {
            typeMethods.vector_constrained_next4(dim_num, alpha, level_min, level_max,
                ref level_1d, q_max, ref more);

            if (!more)
            {
                break;
            }

            cc_level_to_order(dim_num, level_1d, ref order_1d);

            order_nd = 1;
            for (dim = 0; dim < dim_num; dim++)
            {
                order_nd *= order_1d[dim];
            }

            cc_grid(dim_num, order_1d, order_nd, ref grid_point, ptIndex: +point_num * dim_num);

            point_num += order_nd;

            for (dim = 0; dim < dim_num; dim++)
            {
                grid_level[dim + grid_num * dim_num] = level_1d[dim];
            }

            grid_num += 1;
        }
    }

    public static void cc_levels_constrained_size(int dim_num, double q_max, double[] alpha,
            int[] level_min, int[] level_max, ref int grid_num, ref int point_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_LEVELS_CONSTRAINED_SIZE counts grids for CC_LEVELS_CONSTRAINED.
        //
        //  Discussion:
        //
        //    The constraint on the levels of the 1D Clenshaw Curtis rule in 
        //    spatial dimension I is:
        //
        //      LEVEL_MIN(I) <= LEVEL(I) <= LEVEL_MAX(I) 
        //
        //    The constraint on the collection of levels making up a rule is:
        //
        //      Sum ( 1 <= I <= DIM_NUM ) ALPHA(I) * LEVEL(I) <= Q_MAX.
        //
        //    The relationship of level to order is roughly 
        //
        //      ORDER = 2**LEVEL+1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 November 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double Q_MAX, the maximum values of
        //    Q, the sum of the weighted orders in each spatial coordinate.
        //
        //    Input, double ALPHA[DIM_NUM], the weight factors for
        //    the orders in each spatial dimension.
        //
        //    Input, int LEVEL_MIN[DIM_NUM], LEVEL_MAX[DIM_NUM], the minimum
        //    and maximum values of the level of the 1D Clenshaw Curtis rule
        //    in each spatial dimension.
        //
        //    Output, int *GRID_NUM, the number of Clenshaw Curtis
        //    grids in the constraint set.
        //
        //    Output, int *POINT_NUM, the total number of points in the grids.
        //
    {
        int dim;
        int[] level_1d;
        bool more;
        int[] order_1d;
        int order_nd;

        level_1d = new int[dim_num];
        order_1d = new int[dim_num];
        //
        //  Determine the total number of points that will be generated
        //  by "going through the motions".
        //
        point_num = 0;
        grid_num = 0;

        more = false;

            
        for (;;)
        {
            typeMethods.vector_constrained_next4(dim_num, alpha, level_min, level_max,
                ref level_1d, q_max, ref more);

            if (!more)
            {
                break;
            }

            cc_level_to_order(dim_num, level_1d, ref order_1d);

            order_nd = 1;
            for (dim = 0; dim < dim_num; dim++)
            {
                order_nd *= order_1d[dim];
            }

            point_num += order_nd;
            grid_num += 1;
        }
    }

    public static void cc_levels_minmax(int dim_num, int level_min, int level_max,
            int grid_num, int point_num, ref int[] grid_level, ref int[] grid_order,
            ref double[] grid_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_LEVELS_MINMAX computes CC grids with LEVEL_MIN <= LEVEL <= LEVEL_MAX.
        //
        //  Discussion:
        //
        //    The CC grids are required to have an order that is 2**LEVEL + 1.
        //
        //    The necessary dimensions of GRID_LEVEL, GRID_ORDER and GRID_POINT can be
        //    determined by calling CC_LEVELS_MINMAX first.
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
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MIN, LEVEL_MAX, the minimum and maximum values of
        //    LEVEL.
        //
        //    Input, int GRID_NUM, the number of Clenshaw Curtis
        //    grids whose LEVEL value is between LEVEL_MIN and LEVEL_MAX.
        //
        //    Input, int POINT_NUM, the total number of points in the grids.
        //
        //    Output, int GRID_LEVEL[DIM_NUM*GRID_NUM], contains, for each
        //    grid, the level of the Clenshaw-Curtis rule in each dimension.
        //
        //    Output, int GRID_ORDER[DIM_NUM*GRID_NUM], contains, for each
        //    grid, the order of the Clenshaw-Curtis rule in each dimension.
        //
        //    Output, double GRID_POINT[DIM_NUM*POINT_NUM], contains
        //    a list of all the abscissas of all the rules, listed one grid at
        //    a time.  If a point occurs in several grids, it will be listed
        //    several times.
        //
    {
        int dim;
        int h;
        int level;
        int[] level_1d;
        bool more;
        int[] order_1d;
        int order_nd;
        int t;
        //
        //  Outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
        //
        level_1d = new int[dim_num];
        order_1d = new int[dim_num];

        point_num = 0;
        grid_num = 0;

        for (level = level_min; level <= level_max; level++)
        {
            //
            //  Middle loop generates next partition that adds up to LEVEL.
            //
            more = false;
            h = 0;
            t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);
                //
                //  Inner (hidden) loop generates all CC points corresponding to given grid.
                //
                cc_level_to_order(dim_num, level_1d, ref order_1d);

                order_nd = 1;
                for (dim = 0; dim < dim_num; dim++)
                {
                    order_nd *= order_1d[dim];
                }

                cc_grid(dim_num, order_1d, order_nd, ref grid_point, ptIndex: +point_num * dim_num);

                point_num += order_nd;

                for (dim = 0; dim < dim_num; dim++)
                {
                    grid_level[dim + grid_num * dim_num] = level_1d[dim];
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    grid_order[dim + grid_num * dim_num] = order_1d[dim];
                }

                grid_num += 1;

                if (!more)
                {
                    break;
                }
            }
        }
    }

    public static void cc_levels_minmax_size(int dim_num, int level_min, int level_max,
            ref int grid_num, ref int point_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_LEVELS_MINMAX_SIZE counts grids for CC_LEVELS_MINMAX.
        //
        //  Discussion:
        //
        //    This routine can be used to determine the necessary size to be
        //    allocated to arrays GRID_LEVEL, GRID_ORDER and GRID_POINT in a call to
        //    CC_LEVELS_MINMAX.
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
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MIN, LEVEL_MAX, the minimum and maximum values of
        //    LEVEL, the sum of the levels in each spatial coordinate.
        //
        //    Output, int *GRID_NUM, the number of Clenshaw Curtis
        //    grids whose LEVEL value is between LEVEL_MIN and LEVEL_MAX.
        //
        //    Output, int *POINT_NUM, the total number of points in the grids.
        //
    {
        int dim;
        int h;
        int level;
        int[] level_1d;
        bool more;
        int[] order_1d;
        int order_nd;
        int t;
        //
        //  Determine the total number of points that will be generated
        //  by "going through the motions".
        //
        level_1d = new int[dim_num];
        order_1d = new int[dim_num];

        point_num = 0;
        grid_num = 0;

        for (level = level_min; level <= level_max; level++)
        {
            more = false;
            h = 0;
            t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);

                cc_level_to_order(dim_num, level_1d, ref order_1d);

                order_nd = 1;
                for (dim = 0; dim < dim_num; dim++)
                {
                    order_nd *= order_1d[dim];
                }

                point_num += order_nd;
                grid_num += 1;

                if (!more)
                {
                    break;
                }
            }
        }
    }

    public static double cc_weight(int order, int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_WEIGHT returns the I-th Clenshaw Curtis weight.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Input, int I, the index of the desired weight.  1 <= I <= ORDER.
        //
        //    Output, double CC_WEIGHT, the I-th weight in the
        //    Clenshaw-Curtis rule of order ORDER.
        //
    {
        double angle;
        double b;
        int j;
        double value = 0;

        switch (order)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("CC_WEIGHT - Fatal error!");
                Console.WriteLine("  ORDER < 1.");
                return 1;
            case 1:
                value = 2.0;
                return value;
        }

        value = 1.0;

        for (j = 1; j <= (order - 1) / 2; j++)
        {
            if (2 * j == order - 1)
            {
                b = 1.0;
            }
            else
            {
                b = 2.0;
            }

            angle = 2 * j * (i - 1) * Math.PI
                    / (order - 1);

            value -= b * Math.Cos(angle) / (4 * j * j - 1);
        }

        switch (i)
        {
            case 1:
                value /= order - 1;
                break;
            default:
            {
                if (i <= order - 1)
                {
                    value = 2.0 * value / (order - 1);
                }
                else if (i == order)
                {
                    value /= order - 1;
                }

                break;
            }
        }

        return value;
    }

    public static void clenshaw_curtis_compute_nd(int dim_num, int[] order_1d, ref double[] point,
            ref double[] weight )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CLENSHAW_CURTIS_COMPUTE_ND returns a multidimensional Clenshaw-Curtis rule.
        //
        //  Discussion:
        //
        //    The value ORDER_ND, used to dimension the output arrays, is
        //    simply the product of the entries in ORDER_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension of the points.
        //
        //    Input, int ORDER_1D[DIM_NUM], the order of the Clenshaw-Curtis
        //    rule in each dimension.
        //
        //    Output, double POINT[DIM_NUM*ORDER_ND], the points in
        //    the grid.
        //
        //    Output, double WEIGHT[ORDER_ND], the integration weights
        //    associated with the points.
        //
    {
        int dim;
        int i;
        int j;
        int k;
        int n1;
        int n2;
        int n3;
        int order;
        int order_nd;
        int order_old;
        int p;
        double[] w1d = null;
        double[] x1d = null;

        order_nd = typeMethods.i4vec_product(dim_num, order_1d);

        for (order = 0; order < order_nd; order++)
        {
            weight[order] = 1.0;
        }

        order = -1;

        for (dim = 0; dim < dim_num; dim++)
        {
            order_old = order;

            order = order_1d[dim];
            //
            //  For efficiency's sake, we reuse the 1D rule if we can.
            //
            if (order != order_old)
            {
                x1d = new double[order];
                w1d = new double[order];

                ClenshawCurtis.clenshaw_curtis_compute(order, ref x1d, ref w1d);
            }

            p = 0;
            n1 = typeMethods.i4vec_product(dim, order_1d);
            n2 = order_1d[dim];
            n3 = typeMethods.i4vec_product(dim_num - dim - 1, order_1d, aIndex: + dim + 1);

            for (k = 0; k < n1; k++)
            {
                for (j = 0; j < n2; j++)
                {
                    for (i = 0; i < n3; i++)
                    {
                        point[dim + p * dim_num] = x1d[j];
                        weight[p] *= w1d[j];
                        p += 1;
                    }
                }
            }
        }
    }

    public static void clenshaw_curtis_set(int order, ref double[] xtab, ref double[] weight )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CLENSHAW_CURTIS_SET sets a Clenshaw-Curtis quadrature rule.
        //
        //  Discussion:
        //
        //    The integration interval is [ -1, 1 ].
        //
        //    The integral to approximate:
        //
        //      Integral ( -1 <= X <= 1 ) F(X) dX
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 May 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Charles Clenshaw and A Curtis,
        //    A Method for Numerical Integration on an Automatic Computer,
        //    Numerische Mathematik,
        //    Volume 2, Number 1, pages 197-205.
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the rule.
        //    ORDER must be between 1 and 17, 33, 65 or 129.
        //
        //    Output, double XTAB[ORDER], the abscissas of the rule.
        //
        //    Output, double WEIGHT[ORDER], the weights of the rule.
        //    The weights are symmetric and sum to 2.
        //
    {
        switch (order)
        {
            case 1:
                xtab[1 - 1] = 0.0;
                weight[1 - 1] = 2.0;
                break;
            case 2:
                xtab[1 - 1] = -1.0;
                xtab[2 - 1] = 1.0;

                weight[1 - 1] = 1.0;
                weight[2 - 1] = 1.0;
                break;
            case 3:
                xtab[1 - 1] = -1.0;
                xtab[2 - 1] = 0.0;
                xtab[3 - 1] = 1.0;

                weight[1 - 1] = 0.33333333333333;
                weight[2 - 1] = 1.33333333333333;
                weight[3 - 1] = 0.33333333333333;
                break;
            case 4:
                xtab[1 - 1] = -1.0;
                xtab[2 - 1] = -0.5;
                xtab[3 - 1] = 0.5;
                xtab[4 - 1] = 1.0;

                weight[1 - 1] = 0.11111111111111;
                weight[2 - 1] = 0.88888888888889;
                weight[3 - 1] = 0.88888888888889;
                weight[4 - 1] = 0.11111111111111;
                break;
            case 5:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.70710678118655;
                xtab[3 - 1] = 0.00000000000000;
                xtab[4 - 1] = 0.70710678118655;
                xtab[5 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.06666666666667;
                weight[2 - 1] = 0.53333333333333;
                weight[3 - 1] = 0.80000000000000;
                weight[4 - 1] = 0.53333333333333;
                weight[5 - 1] = 0.06666666666667;
                break;
            case 6:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.80901699437495;
                xtab[3 - 1] = -0.30901699437495;
                xtab[4 - 1] = 0.30901699437495;
                xtab[5 - 1] = 0.80901699437495;
                xtab[6 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.04000000000000;
                weight[2 - 1] = 0.36074304120001;
                weight[3 - 1] = 0.59925695879999;
                weight[4 - 1] = 0.59925695879999;
                weight[5 - 1] = 0.36074304120001;
                weight[6 - 1] = 0.04000000000000;
                break;
            case 7:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.86602540378444;
                xtab[3 - 1] = -0.50000000000000;
                xtab[4 - 1] = 0.00000000000000;
                xtab[5 - 1] = 0.50000000000000;
                xtab[6 - 1] = 0.86602540378444;
                xtab[7 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.02857142857143;
                weight[2 - 1] = 0.25396825396825;
                weight[3 - 1] = 0.45714285714286;
                weight[4 - 1] = 0.52063492063492;
                weight[5 - 1] = 0.45714285714286;
                weight[6 - 1] = 0.25396825396825;
                weight[7 - 1] = 0.02857142857143;
                break;
            case 8:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.90096886790242;
                xtab[3 - 1] = -0.62348980185873;
                xtab[4 - 1] = -0.22252093395631;
                xtab[5 - 1] = 0.22252093395631;
                xtab[6 - 1] = 0.62348980185873;
                xtab[7 - 1] = 0.90096886790242;
                xtab[8 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.02040816326531;
                weight[2 - 1] = 0.19014100721821;
                weight[3 - 1] = 0.35224242371816;
                weight[4 - 1] = 0.43720840579833;
                weight[5 - 1] = 0.43720840579833;
                weight[6 - 1] = 0.35224242371816;
                weight[7 - 1] = 0.19014100721821;
                weight[8 - 1] = 0.02040816326531;
                break;
            case 9:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.92387953251129;
                xtab[3 - 1] = -0.70710678118655;
                xtab[4 - 1] = -0.38268343236509;
                xtab[5 - 1] = 0.00000000000000;
                xtab[6 - 1] = 0.38268343236509;
                xtab[7 - 1] = 0.70710678118655;
                xtab[8 - 1] = 0.92387953251129;
                xtab[9 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.01587301587302;
                weight[2 - 1] = 0.14621864921602;
                weight[3 - 1] = 0.27936507936508;
                weight[4 - 1] = 0.36171785872049;
                weight[5 - 1] = 0.39365079365079;
                weight[6 - 1] = 0.36171785872049;
                weight[7 - 1] = 0.27936507936508;
                weight[8 - 1] = 0.14621864921602;
                weight[9 - 1] = 0.01587301587302;
                break;
            case 10:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.93969262078591;
                xtab[3 - 1] = -0.76604444311898;
                xtab[4 - 1] = -0.50000000000000;
                xtab[5 - 1] = -0.17364817766693;
                xtab[6 - 1] = 0.17364817766693;
                xtab[7 - 1] = 0.50000000000000;
                xtab[8 - 1] = 0.76604444311898;
                xtab[9 - 1] = 0.93969262078591;
                xtab[10 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.01234567901235;
                weight[2 - 1] = 0.11656745657204;
                weight[3 - 1] = 0.22528432333810;
                weight[4 - 1] = 0.30194003527337;
                weight[5 - 1] = 0.34386250580414;
                weight[6 - 1] = 0.34386250580414;
                weight[7 - 1] = 0.30194003527337;
                weight[8 - 1] = 0.22528432333810;
                weight[9 - 1] = 0.11656745657204;
                weight[10 - 1] = 0.01234567901235;
                break;
            case 11:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.95105651629515;
                xtab[3 - 1] = -0.80901699437495;
                xtab[4 - 1] = -0.58778525229247;
                xtab[5 - 1] = -0.30901699437495;
                xtab[6 - 1] = 0.00000000000000;
                xtab[7 - 1] = 0.30901699437495;
                xtab[8 - 1] = 0.58778525229247;
                xtab[9 - 1] = 0.80901699437495;
                xtab[10 - 1] = 0.95105651629515;
                xtab[11 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.01010101010101;
                weight[2 - 1] = 0.09457905488370;
                weight[3 - 1] = 0.18563521442425;
                weight[4 - 1] = 0.25358833328369;
                weight[5 - 1] = 0.29921327042424;
                weight[6 - 1] = 0.31376623376623;
                weight[7 - 1] = 0.29921327042424;
                weight[8 - 1] = 0.25358833328369;
                weight[9 - 1] = 0.18563521442425;
                weight[10 - 1] = 0.09457905488370;
                weight[11 - 1] = 0.01010101010101;
                break;
            case 12:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.95949297361450;
                xtab[3 - 1] = -0.84125353283118;
                xtab[4 - 1] = -0.65486073394528;
                xtab[5 - 1] = -0.41541501300189;
                xtab[6 - 1] = -0.14231483827329;
                xtab[7 - 1] = 0.14231483827329;
                xtab[8 - 1] = 0.41541501300189;
                xtab[9 - 1] = 0.65486073394529;
                xtab[10 - 1] = 0.84125353283118;
                xtab[11 - 1] = 0.95949297361450;
                xtab[12 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.00826446280992;
                weight[2 - 1] = 0.07856015374620;
                weight[3 - 1] = 0.15504045508256;
                weight[4 - 1] = 0.21556254600087;
                weight[5 - 1] = 0.25991734106692;
                weight[6 - 1] = 0.28265504129354;
                weight[7 - 1] = 0.28265504129354;
                weight[8 - 1] = 0.25991734106692;
                weight[9 - 1] = 0.21556254600087;
                weight[10 - 1] = 0.15504045508256;
                weight[11 - 1] = 0.07856015374620;
                weight[12 - 1] = 0.00826446280992;
                break;
            case 13:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.96592582628907;
                xtab[3 - 1] = -0.86602540378444;
                xtab[4 - 1] = -0.70710678118655;
                xtab[5 - 1] = -0.50000000000000;
                xtab[6 - 1] = -0.25881904510252;
                xtab[7 - 1] = 0.00000000000000;
                xtab[8 - 1] = 0.25881904510252;
                xtab[9 - 1] = 0.50000000000000;
                xtab[10 - 1] = 0.70710678118655;
                xtab[11 - 1] = 0.86602540378444;
                xtab[12 - 1] = 0.96592582628907;
                xtab[13 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.00699300699301;
                weight[2 - 1] = 0.06605742495207;
                weight[3 - 1] = 0.13154253154253;
                weight[4 - 1] = 0.18476338476338;
                weight[5 - 1] = 0.22697302697303;
                weight[6 - 1] = 0.25267569378104;
                weight[7 - 1] = 0.26198986198986;
                weight[8 - 1] = 0.25267569378104;
                weight[9 - 1] = 0.22697302697303;
                weight[10 - 1] = 0.18476338476338;
                weight[11 - 1] = 0.13154253154253;
                weight[12 - 1] = 0.06605742495207;
                weight[13 - 1] = 0.00699300699301;
                break;
            case 14:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.97094181742605;
                xtab[3 - 1] = -0.88545602565321;
                xtab[4 - 1] = -0.74851074817110;
                xtab[5 - 1] = -0.56806474673116;
                xtab[6 - 1] = -0.35460488704254;
                xtab[7 - 1] = -0.12053668025532;
                xtab[8 - 1] = 0.12053668025532;
                xtab[9 - 1] = 0.35460488704254;
                xtab[10 - 1] = 0.56806474673116;
                xtab[11 - 1] = 0.74851074817110;
                xtab[12 - 1] = 0.88545602565321;
                xtab[13 - 1] = 0.97094181742605;
                xtab[14 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.00591715976331;
                weight[2 - 1] = 0.05646531376341;
                weight[3 - 1] = 0.11276867248986;
                weight[4 - 1] = 0.16003802611672;
                weight[5 - 1] = 0.19899241036578;
                weight[6 - 1] = 0.22590304977856;
                weight[7 - 1] = 0.23991536772235;
                weight[8 - 1] = 0.23991536772235;
                weight[9 - 1] = 0.22590304977856;
                weight[10 - 1] = 0.19899241036578;
                weight[11 - 1] = 0.16003802611672;
                weight[12 - 1] = 0.11276867248986;
                weight[13 - 1] = 0.05646531376341;
                weight[14 - 1] = 0.00591715976331;
                break;
            case 15:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.97492791218182;
                xtab[3 - 1] = -0.90096886790242;
                xtab[4 - 1] = -0.78183148246803;
                xtab[5 - 1] = -0.62348980185873;
                xtab[6 - 1] = -0.43388373911756;
                xtab[7 - 1] = -0.22252093395631;
                xtab[8 - 1] = 0.00000000000000;
                xtab[9 - 1] = 0.22252093395631;
                xtab[10 - 1] = 0.43388373911756;
                xtab[11 - 1] = 0.62348980185873;
                xtab[12 - 1] = 0.78183148246803;
                xtab[13 - 1] = 0.90096886790242;
                xtab[14 - 1] = 0.97492791218182;
                xtab[15 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.00512820512821;
                weight[2 - 1] = 0.04869938729509;
                weight[3 - 1] = 0.09782039167605;
                weight[4 - 1] = 0.13966507849560;
                weight[5 - 1] = 0.17560578900107;
                weight[6 - 1] = 0.20205146748238;
                weight[7 - 1] = 0.21888151163057;
                weight[8 - 1] = 0.22429633858205;
                weight[9 - 1] = 0.21888151163057;
                weight[10 - 1] = 0.20205146748238;
                weight[11 - 1] = 0.17560578900107;
                weight[12 - 1] = 0.13966507849560;
                weight[13 - 1] = 0.09782039167605;
                weight[14 - 1] = 0.04869938729509;
                weight[15 - 1] = 0.00512820512821;
                break;
            case 16:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.97814760073381;
                xtab[3 - 1] = -0.91354545764260;
                xtab[4 - 1] = -0.80901699437495;
                xtab[5 - 1] = -0.66913060635886;
                xtab[6 - 1] = -0.50000000000000;
                xtab[7 - 1] = -0.30901699437495;
                xtab[8 - 1] = -0.10452846326765;
                xtab[9 - 1] = 0.10452846326765;
                xtab[10 - 1] = 0.30901699437495;
                xtab[11 - 1] = 0.50000000000000;
                xtab[12 - 1] = 0.66913060635886;
                xtab[13 - 1] = 0.80901699437495;
                xtab[14 - 1] = 0.91354545764260;
                xtab[15 - 1] = 0.97814760073381;
                xtab[16 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.00444444444444;
                weight[2 - 1] = 0.04251476624753;
                weight[3 - 1] = 0.08553884025933;
                weight[4 - 1] = 0.12294010082849;
                weight[5 - 1] = 0.15573317603967;
                weight[6 - 1] = 0.18132978132978;
                weight[7 - 1] = 0.19921478132639;
                weight[8 - 1] = 0.20828410952436;
                weight[9 - 1] = 0.20828410952436;
                weight[10 - 1] = 0.19921478132639;
                weight[11 - 1] = 0.18132978132978;
                weight[12 - 1] = 0.15573317603967;
                weight[13 - 1] = 0.12294010082849;
                weight[14 - 1] = 0.08553884025933;
                weight[15 - 1] = 0.04251476624753;
                weight[16 - 1] = 0.00444444444444;
                break;
            case 17:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.98078528040323;
                xtab[3 - 1] = -0.92387953251129;
                xtab[4 - 1] = -0.83146961230255;
                xtab[5 - 1] = -0.70710678118655;
                xtab[6 - 1] = -0.55557023301960;
                xtab[7 - 1] = -0.38268343236509;
                xtab[8 - 1] = -0.19509032201613;
                xtab[9 - 1] = 0.00000000000000;
                xtab[10 - 1] = 0.19509032201613;
                xtab[11 - 1] = 0.38268343236509;
                xtab[12 - 1] = 0.55557023301960;
                xtab[13 - 1] = 0.70710678118655;
                xtab[14 - 1] = 0.83146961230255;
                xtab[15 - 1] = 0.92387953251129;
                xtab[16 - 1] = 0.98078528040323;
                xtab[17 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.00392156862745;
                weight[2 - 1] = 0.03736870283721;
                weight[3 - 1] = 0.07548233154315;
                weight[4 - 1] = 0.10890555258189;
                weight[5 - 1] = 0.13895646836823;
                weight[6 - 1] = 0.16317266428170;
                weight[7 - 1] = 0.18147378423649;
                weight[8 - 1] = 0.19251386461293;
                weight[9 - 1] = 0.19641012582189;
                weight[10 - 1] = 0.19251386461293;
                weight[11 - 1] = 0.18147378423649;
                weight[12 - 1] = 0.16317266428170;
                weight[13 - 1] = 0.13895646836823;
                weight[14 - 1] = 0.10890555258189;
                weight[15 - 1] = 0.07548233154315;
                weight[16 - 1] = 0.03736870283721;
                weight[17 - 1] = 0.00392156862745;
                break;
            case 33:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.99518472667220;
                xtab[3 - 1] = -0.98078528040323;
                xtab[4 - 1] = -0.95694033573221;
                xtab[5 - 1] = -0.92387953251129;
                xtab[6 - 1] = -0.88192126434835;
                xtab[7 - 1] = -0.83146961230255;
                xtab[8 - 1] = -0.77301045336274;
                xtab[9 - 1] = -0.70710678118655;
                xtab[10 - 1] = -0.63439328416365;
                xtab[11 - 1] = -0.55557023301960;
                xtab[12 - 1] = -0.47139673682600;
                xtab[13 - 1] = -0.38268343236509;
                xtab[14 - 1] = -0.29028467725446;
                xtab[15 - 1] = -0.19509032201613;
                xtab[16 - 1] = -0.09801714032956;
                xtab[17 - 1] = 0.00000000000000;
                xtab[18 - 1] = 0.09801714032956;
                xtab[19 - 1] = 0.19509032201613;
                xtab[20 - 1] = 0.29028467725446;
                xtab[21 - 1] = 0.38268343236509;
                xtab[22 - 1] = 0.47139673682600;
                xtab[23 - 1] = 0.55557023301960;
                xtab[24 - 1] = 0.63439328416365;
                xtab[25 - 1] = 0.70710678118655;
                xtab[26 - 1] = 0.77301045336274;
                xtab[27 - 1] = 0.83146961230255;
                xtab[28 - 1] = 0.88192126434836;
                xtab[29 - 1] = 0.92387953251129;
                xtab[30 - 1] = 0.95694033573221;
                xtab[31 - 1] = 0.98078528040323;
                xtab[32 - 1] = 0.99518472667220;
                xtab[33 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.00097751710655;
                weight[2 - 1] = 0.00939319796296;
                weight[3 - 1] = 0.01923424513268;
                weight[4 - 1] = 0.02845791667723;
                weight[5 - 1] = 0.03759434191405;
                weight[6 - 1] = 0.04626276283775;
                weight[7 - 1] = 0.05455501630398;
                weight[8 - 1] = 0.06227210954529;
                weight[9 - 1] = 0.06942757563044;
                weight[10 - 1] = 0.07588380044139;
                weight[11 - 1] = 0.08163481765494;
                weight[12 - 1] = 0.08657753844183;
                weight[13 - 1] = 0.09070611286772;
                weight[14 - 1] = 0.09394324443877;
                weight[15 - 1] = 0.09629232594549;
                weight[16 - 1] = 0.09769818820806;
                weight[17 - 1] = 0.09817857778177;
                weight[18 - 1] = 0.09769818820806;
                weight[19 - 1] = 0.09629232594549;
                weight[20 - 1] = 0.09394324443877;
                weight[21 - 1] = 0.09070611286772;
                weight[22 - 1] = 0.08657753844183;
                weight[23 - 1] = 0.08163481765494;
                weight[24 - 1] = 0.07588380044139;
                weight[25 - 1] = 0.06942757563044;
                weight[26 - 1] = 0.06227210954529;
                weight[27 - 1] = 0.05455501630398;
                weight[28 - 1] = 0.04626276283775;
                weight[29 - 1] = 0.03759434191405;
                weight[30 - 1] = 0.02845791667723;
                weight[31 - 1] = 0.01923424513268;
                weight[32 - 1] = 0.00939319796296;
                weight[33 - 1] = 0.00097751710655;
                break;
            case 65:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.99879545620517;
                xtab[3 - 1] = -0.99518472667220;
                xtab[4 - 1] = -0.98917650996478;
                xtab[5 - 1] = -0.98078528040323;
                xtab[6 - 1] = -0.97003125319454;
                xtab[7 - 1] = -0.95694033573221;
                xtab[8 - 1] = -0.94154406518302;
                xtab[9 - 1] = -0.92387953251129;
                xtab[10 - 1] = -0.90398929312344;
                xtab[11 - 1] = -0.88192126434835;
                xtab[12 - 1] = -0.85772861000027;
                xtab[13 - 1] = -0.83146961230255;
                xtab[14 - 1] = -0.80320753148064;
                xtab[15 - 1] = -0.77301045336274;
                xtab[16 - 1] = -0.74095112535496;
                xtab[17 - 1] = -0.70710678118655;
                xtab[18 - 1] = -0.67155895484702;
                xtab[19 - 1] = -0.63439328416365;
                xtab[20 - 1] = -0.59569930449243;
                xtab[21 - 1] = -0.55557023301960;
                xtab[22 - 1] = -0.51410274419322;
                xtab[23 - 1] = -0.47139673682600;
                xtab[24 - 1] = -0.42755509343028;
                xtab[25 - 1] = -0.38268343236509;
                xtab[26 - 1] = -0.33688985339222;
                xtab[27 - 1] = -0.29028467725446;
                xtab[28 - 1] = -0.24298017990326;
                xtab[29 - 1] = -0.19509032201613;
                xtab[30 - 1] = -0.14673047445536;
                xtab[31 - 1] = -0.09801714032956;
                xtab[32 - 1] = -0.04906767432742;
                xtab[33 - 1] = 0.00000000000000;
                xtab[34 - 1] = 0.04906767432742;
                xtab[35 - 1] = 0.09801714032956;
                xtab[36 - 1] = 0.14673047445536;
                xtab[37 - 1] = 0.19509032201613;
                xtab[38 - 1] = 0.24298017990326;
                xtab[39 - 1] = 0.29028467725446;
                xtab[40 - 1] = 0.33688985339222;
                xtab[41 - 1] = 0.38268343236509;
                xtab[42 - 1] = 0.42755509343028;
                xtab[43 - 1] = 0.47139673682600;
                xtab[44 - 1] = 0.51410274419322;
                xtab[45 - 1] = 0.55557023301960;
                xtab[46 - 1] = 0.59569930449243;
                xtab[47 - 1] = 0.63439328416365;
                xtab[48 - 1] = 0.67155895484702;
                xtab[49 - 1] = 0.70710678118655;
                xtab[50 - 1] = 0.74095112535496;
                xtab[51 - 1] = 0.77301045336274;
                xtab[52 - 1] = 0.80320753148064;
                xtab[53 - 1] = 0.83146961230255;
                xtab[54 - 1] = 0.85772861000027;
                xtab[55 - 1] = 0.88192126434836;
                xtab[56 - 1] = 0.90398929312344;
                xtab[57 - 1] = 0.92387953251129;
                xtab[58 - 1] = 0.94154406518302;
                xtab[59 - 1] = 0.95694033573221;
                xtab[60 - 1] = 0.97003125319454;
                xtab[61 - 1] = 0.98078528040323;
                xtab[62 - 1] = 0.98917650996478;
                xtab[63 - 1] = 0.99518472667220;
                xtab[64 - 1] = 0.99879545620517;
                xtab[65 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.00024420024420;
                weight[2 - 1] = 0.00235149067531;
                weight[3 - 1] = 0.00483146544879;
                weight[4 - 1] = 0.00719269316174;
                weight[5 - 1] = 0.00958233879528;
                weight[6 - 1] = 0.01192339471421;
                weight[7 - 1] = 0.01425206043235;
                weight[8 - 1] = 0.01653498765729;
                weight[9 - 1] = 0.01878652974180;
                weight[10 - 1] = 0.02098627442974;
                weight[11 - 1] = 0.02314069493436;
                weight[12 - 1] = 0.02523506498175;
                weight[13 - 1] = 0.02727225714147;
                weight[14 - 1] = 0.02924065319747;
                weight[15 - 1] = 0.03114129710407;
                weight[16 - 1] = 0.03296454656998;
                weight[17 - 1] = 0.03471049818093;
                weight[18 - 1] = 0.03637092028664;
                weight[19 - 1] = 0.03794545992128;
                weight[20 - 1] = 0.03942698871296;
                weight[21 - 1] = 0.04081501340036;
                weight[22 - 1] = 0.04210333111142;
                weight[23 - 1] = 0.04329151496169;
                weight[24 - 1] = 0.04437417923926;
                weight[25 - 1] = 0.04535110955166;
                weight[26 - 1] = 0.04621766751093;
                weight[27 - 1] = 0.04697395904661;
                weight[28 - 1] = 0.04761604458525;
                weight[29 - 1] = 0.04814443257251;
                weight[30 - 1] = 0.04855584485714;
                weight[31 - 1] = 0.04885125664307;
                weight[32 - 1] = 0.04902801843103;
                weight[33 - 1] = 0.04908762351494;
                weight[34 - 1] = 0.04902801843103;
                weight[35 - 1] = 0.04885125664307;
                weight[36 - 1] = 0.04855584485714;
                weight[37 - 1] = 0.04814443257251;
                weight[38 - 1] = 0.04761604458525;
                weight[39 - 1] = 0.04697395904661;
                weight[40 - 1] = 0.04621766751093;
                weight[41 - 1] = 0.04535110955166;
                weight[42 - 1] = 0.04437417923926;
                weight[43 - 1] = 0.04329151496169;
                weight[44 - 1] = 0.04210333111142;
                weight[45 - 1] = 0.04081501340036;
                weight[46 - 1] = 0.03942698871296;
                weight[47 - 1] = 0.03794545992128;
                weight[48 - 1] = 0.03637092028664;
                weight[49 - 1] = 0.03471049818093;
                weight[50 - 1] = 0.03296454656998;
                weight[51 - 1] = 0.03114129710407;
                weight[52 - 1] = 0.02924065319747;
                weight[53 - 1] = 0.02727225714147;
                weight[54 - 1] = 0.02523506498175;
                weight[55 - 1] = 0.02314069493436;
                weight[56 - 1] = 0.02098627442974;
                weight[57 - 1] = 0.01878652974180;
                weight[58 - 1] = 0.01653498765729;
                weight[59 - 1] = 0.01425206043235;
                weight[60 - 1] = 0.01192339471421;
                weight[61 - 1] = 0.00958233879528;
                weight[62 - 1] = 0.00719269316174;
                weight[63 - 1] = 0.00483146544879;
                weight[64 - 1] = 0.00235149067531;
                weight[65 - 1] = 0.00024420024420;
                break;
            case 129:
                xtab[1 - 1] = -1.00000000000000;
                xtab[2 - 1] = -0.99969881869620;
                xtab[3 - 1] = -0.99879545620517;
                xtab[4 - 1] = -0.99729045667869;
                xtab[5 - 1] = -0.99518472667220;
                xtab[6 - 1] = -0.99247953459871;
                xtab[7 - 1] = -0.98917650996478;
                xtab[8 - 1] = -0.98527764238894;
                xtab[9 - 1] = -0.98078528040323;
                xtab[10 - 1] = -0.97570213003853;
                xtab[11 - 1] = -0.97003125319454;
                xtab[12 - 1] = -0.96377606579544;
                xtab[13 - 1] = -0.95694033573221;
                xtab[14 - 1] = -0.94952818059304;
                xtab[15 - 1] = -0.94154406518302;
                xtab[16 - 1] = -0.93299279883474;
                xtab[17 - 1] = -0.92387953251129;
                xtab[18 - 1] = -0.91420975570353;
                xtab[19 - 1] = -0.90398929312344;
                xtab[20 - 1] = -0.89322430119552;
                xtab[21 - 1] = -0.88192126434835;
                xtab[22 - 1] = -0.87008699110871;
                xtab[23 - 1] = -0.85772861000027;
                xtab[24 - 1] = -0.84485356524971;
                xtab[25 - 1] = -0.83146961230255;
                xtab[26 - 1] = -0.81758481315158;
                xtab[27 - 1] = -0.80320753148064;
                xtab[28 - 1] = -0.78834642762661;
                xtab[29 - 1] = -0.77301045336274;
                xtab[30 - 1] = -0.75720884650648;
                xtab[31 - 1] = -0.74095112535496;
                xtab[32 - 1] = -0.72424708295147;
                xtab[33 - 1] = -0.70710678118655;
                xtab[34 - 1] = -0.68954054473707;
                xtab[35 - 1] = -0.67155895484702;
                xtab[36 - 1] = -0.65317284295378;
                xtab[37 - 1] = -0.63439328416365;
                xtab[38 - 1] = -0.61523159058063;
                xtab[39 - 1] = -0.59569930449243;
                xtab[40 - 1] = -0.57580819141785;
                xtab[41 - 1] = -0.55557023301960;
                xtab[42 - 1] = -0.53499761988710;
                xtab[43 - 1] = -0.51410274419322;
                xtab[44 - 1] = -0.49289819222978;
                xtab[45 - 1] = -0.47139673682600;
                xtab[46 - 1] = -0.44961132965461;
                xtab[47 - 1] = -0.42755509343028;
                xtab[48 - 1] = -0.40524131400499;
                xtab[49 - 1] = -0.38268343236509;
                xtab[50 - 1] = -0.35989503653499;
                xtab[51 - 1] = -0.33688985339222;
                xtab[52 - 1] = -0.31368174039889;
                xtab[53 - 1] = -0.29028467725446;
                xtab[54 - 1] = -0.26671275747490;
                xtab[55 - 1] = -0.24298017990326;
                xtab[56 - 1] = -0.21910124015687;
                xtab[57 - 1] = -0.19509032201613;
                xtab[58 - 1] = -0.17096188876030;
                xtab[59 - 1] = -0.14673047445536;
                xtab[60 - 1] = -0.12241067519922;
                xtab[61 - 1] = -0.09801714032956;
                xtab[62 - 1] = -0.07356456359967;
                xtab[63 - 1] = -0.04906767432742;
                xtab[64 - 1] = -0.02454122852291;
                xtab[65 - 1] = 0.00000000000000;
                xtab[66 - 1] = 0.02454122852291;
                xtab[67 - 1] = 0.04906767432742;
                xtab[68 - 1] = 0.07356456359967;
                xtab[69 - 1] = 0.09801714032956;
                xtab[70 - 1] = 0.12241067519922;
                xtab[71 - 1] = 0.14673047445536;
                xtab[72 - 1] = 0.17096188876030;
                xtab[73 - 1] = 0.19509032201613;
                xtab[74 - 1] = 0.21910124015687;
                xtab[75 - 1] = 0.24298017990326;
                xtab[76 - 1] = 0.26671275747490;
                xtab[77 - 1] = 0.29028467725446;
                xtab[78 - 1] = 0.31368174039889;
                xtab[79 - 1] = 0.33688985339222;
                xtab[80 - 1] = 0.35989503653499;
                xtab[81 - 1] = 0.38268343236509;
                xtab[82 - 1] = 0.40524131400499;
                xtab[83 - 1] = 0.42755509343028;
                xtab[84 - 1] = 0.44961132965461;
                xtab[85 - 1] = 0.47139673682600;
                xtab[86 - 1] = 0.49289819222978;
                xtab[87 - 1] = 0.51410274419322;
                xtab[88 - 1] = 0.53499761988710;
                xtab[89 - 1] = 0.55557023301960;
                xtab[90 - 1] = 0.57580819141785;
                xtab[91 - 1] = 0.59569930449243;
                xtab[92 - 1] = 0.61523159058063;
                xtab[93 - 1] = 0.63439328416365;
                xtab[94 - 1] = 0.65317284295378;
                xtab[95 - 1] = 0.67155895484702;
                xtab[96 - 1] = 0.68954054473707;
                xtab[97 - 1] = 0.70710678118655;
                xtab[98 - 1] = 0.72424708295147;
                xtab[99 - 1] = 0.74095112535496;
                xtab[100 - 1] = 0.75720884650648;
                xtab[101 - 1] = 0.77301045336274;
                xtab[102 - 1] = 0.78834642762661;
                xtab[103 - 1] = 0.80320753148064;
                xtab[104 - 1] = 0.81758481315158;
                xtab[105 - 1] = 0.83146961230255;
                xtab[106 - 1] = 0.84485356524971;
                xtab[107 - 1] = 0.85772861000027;
                xtab[108 - 1] = 0.87008699110871;
                xtab[109 - 1] = 0.88192126434836;
                xtab[110 - 1] = 0.89322430119552;
                xtab[111 - 1] = 0.90398929312344;
                xtab[112 - 1] = 0.91420975570353;
                xtab[113 - 1] = 0.92387953251129;
                xtab[114 - 1] = 0.93299279883474;
                xtab[115 - 1] = 0.94154406518302;
                xtab[116 - 1] = 0.94952818059304;
                xtab[117 - 1] = 0.95694033573221;
                xtab[118 - 1] = 0.96377606579544;
                xtab[119 - 1] = 0.97003125319454;
                xtab[120 - 1] = 0.97570213003853;
                xtab[121 - 1] = 0.98078528040323;
                xtab[122 - 1] = 0.98527764238894;
                xtab[123 - 1] = 0.98917650996478;
                xtab[124 - 1] = 0.99247953459871;
                xtab[125 - 1] = 0.99518472667220;
                xtab[126 - 1] = 0.99729045667869;
                xtab[127 - 1] = 0.99879545620517;
                xtab[128 - 1] = 0.99969881869620;
                xtab[129 - 1] = 1.00000000000000;

                weight[1 - 1] = 0.00006103888177;
                weight[2 - 1] = 0.00058807215383;
                weight[3 - 1] = 0.00120930061875;
                weight[4 - 1] = 0.00180308126695;
                weight[5 - 1] = 0.00240715327877;
                weight[6 - 1] = 0.00300345869904;
                weight[7 - 1] = 0.00360197835813;
                weight[8 - 1] = 0.00419553798719;
                weight[9 - 1] = 0.00478862143341;
                weight[10 - 1] = 0.00537724746840;
                weight[11 - 1] = 0.00596388034731;
                weight[12 - 1] = 0.00654590843862;
                weight[13 - 1] = 0.00712483332325;
                weight[14 - 1] = 0.00769875778896;
                weight[15 - 1] = 0.00826865154203;
                weight[16 - 1] = 0.00883303867470;
                weight[17 - 1] = 0.00939256583935;
                weight[18 - 1] = 0.00994602784923;
                weight[19 - 1] = 0.01049386202577;
                weight[20 - 1] = 0.01103504877427;
                weight[21 - 1] = 0.01156988348291;
                weight[22 - 1] = 0.01209748052807;
                weight[23 - 1] = 0.01261803597978;
                weight[24 - 1] = 0.01313076516694;
                weight[25 - 1] = 0.01363579321294;
                weight[26 - 1] = 0.01413241437853;
                weight[27 - 1] = 0.01462070254634;
                weight[28 - 1] = 0.01510001572479;
                weight[29 - 1] = 0.01557039073899;
                weight[30 - 1] = 0.01603123858745;
                weight[31 - 1] = 0.01648256956220;
                weight[32 - 1] = 0.01692383985846;
                weight[33 - 1] = 0.01735504125411;
                weight[34 - 1] = 0.01777566938875;
                weight[35 - 1] = 0.01818570377926;
                weight[36 - 1] = 0.01858467519567;
                weight[37 - 1] = 0.01897255587068;
                weight[38 - 1] = 0.01934890842392;
                weight[39 - 1] = 0.01971370183700;
                weight[40 - 1] = 0.02006652805198;
                weight[41 - 1] = 0.02040735612004;
                weight[42 - 1] = 0.02073580533490;
                weight[43 - 1] = 0.02105184759002;
                weight[44 - 1] = 0.02135512797426;
                weight[45 - 1] = 0.02164562356713;
                weight[46 - 1] = 0.02192300400599;
                weight[47 - 1] = 0.02218725355897;
                weight[48 - 1] = 0.02243806539723;
                weight[49 - 1] = 0.02267543270457;
                weight[50 - 1] = 0.02289907134391;
                weight[51 - 1] = 0.02310898491627;
                weight[52 - 1] = 0.02330491126131;
                weight[53 - 1] = 0.02348686571193;
                weight[54 - 1] = 0.02365460746058;
                weight[55 - 1] = 0.02380816473024;
                weight[56 - 1] = 0.02394731750477;
                weight[57 - 1] = 0.02407210792328;
                weight[58 - 1] = 0.02418233623893;
                weight[59 - 1] = 0.02427805942076;
                weight[60 - 1] = 0.02435909748928;
                weight[61 - 1] = 0.02442552306157;
                weight[62 - 1] = 0.02447717542743;
                weight[63 - 1] = 0.02451414358882;
                weight[64 - 1] = 0.02453628559651;
                weight[65 - 1] = 0.02454370750551;
                weight[66 - 1] = 0.02453628559651;
                weight[67 - 1] = 0.02451414358882;
                weight[68 - 1] = 0.02447717542743;
                weight[69 - 1] = 0.02442552306157;
                weight[70 - 1] = 0.02435909748928;
                weight[71 - 1] = 0.02427805942076;
                weight[72 - 1] = 0.02418233623893;
                weight[73 - 1] = 0.02407210792328;
                weight[74 - 1] = 0.02394731750477;
                weight[75 - 1] = 0.02380816473024;
                weight[76 - 1] = 0.02365460746058;
                weight[77 - 1] = 0.02348686571193;
                weight[78 - 1] = 0.02330491126131;
                weight[79 - 1] = 0.02310898491627;
                weight[80 - 1] = 0.02289907134391;
                weight[81 - 1] = 0.02267543270457;
                weight[82 - 1] = 0.02243806539723;
                weight[83 - 1] = 0.02218725355897;
                weight[84 - 1] = 0.02192300400599;
                weight[85 - 1] = 0.02164562356713;
                weight[86 - 1] = 0.02135512797426;
                weight[87 - 1] = 0.02105184759002;
                weight[88 - 1] = 0.02073580533490;
                weight[89 - 1] = 0.02040735612004;
                weight[90 - 1] = 0.02006652805198;
                weight[91 - 1] = 0.01971370183700;
                weight[92 - 1] = 0.01934890842392;
                weight[93 - 1] = 0.01897255587068;
                weight[94 - 1] = 0.01858467519567;
                weight[95 - 1] = 0.01818570377926;
                weight[96 - 1] = 0.01777566938875;
                weight[97 - 1] = 0.01735504125411;
                weight[98 - 1] = 0.01692383985846;
                weight[99 - 1] = 0.01648256956220;
                weight[100 - 1] = 0.01603123858745;
                weight[101 - 1] = 0.01557039073899;
                weight[102 - 1] = 0.01510001572479;
                weight[103 - 1] = 0.01462070254634;
                weight[104 - 1] = 0.01413241437853;
                weight[105 - 1] = 0.01363579321294;
                weight[106 - 1] = 0.01313076516694;
                weight[107 - 1] = 0.01261803597978;
                weight[108 - 1] = 0.01209748052807;
                weight[109 - 1] = 0.01156988348291;
                weight[110 - 1] = 0.01103504877427;
                weight[111 - 1] = 0.01049386202577;
                weight[112 - 1] = 0.00994602784923;
                weight[113 - 1] = 0.00939256583935;
                weight[114 - 1] = 0.00883303867470;
                weight[115 - 1] = 0.00826865154203;
                weight[116 - 1] = 0.00769875778896;
                weight[117 - 1] = 0.00712483332325;
                weight[118 - 1] = 0.00654590843862;
                weight[119 - 1] = 0.00596388034731;
                weight[120 - 1] = 0.00537724746840;
                weight[121 - 1] = 0.00478862143341;
                weight[122 - 1] = 0.00419553798719;
                weight[123 - 1] = 0.00360197835813;
                weight[124 - 1] = 0.00300345869904;
                weight[125 - 1] = 0.00240715327877;
                weight[126 - 1] = 0.00180308126695;
                weight[127 - 1] = 0.00120930061875;
                weight[128 - 1] = 0.00058807215383;
                weight[129 - 1] = 0.00006103888177;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("CLENSHAW_CURTIS_SET - Fatal error!");
                Console.WriteLine("  Illegal value of ORDER = " + order + "");
                Console.WriteLine("  Legal values are 1 to 17, 33, 65 or 129.");
                break;
        }
    }
}