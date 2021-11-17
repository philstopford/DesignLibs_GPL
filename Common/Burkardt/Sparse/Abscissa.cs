using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Quadrature;

namespace Burkardt.Sparse;

public static class Abscissa
{
    public static int[] abscissa_level_closed_nd(int level_max, int dim_num, int test_num,
            int[] test_val)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ABSCISSA_LEVEL_CLOSED_ND: first level at which an abscissa is generated.
        //
        //  Discussion:
        //
        //    We need this routine because the sparse grid is generated as a sum of 
        //    product grids, and many points in the sparse grid will belong to several
        //    of these product grids, and we need to do something special the very 
        //    first time we encounter such a point - namely, count it.  So this routine 
        //    determines, for any point in the full product grid, the first level 
        //    at which that point would be included.
        //
        //
        //    We assume an underlying product grid.  In each dimension, this product
        //    grid has order 2^LEVEL_MAX + 1.
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
        //    that the product grid is the product of 1D grids,
        //    that the 1D grids are built up by levels, having
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
        //    Input, int LEVEL_MAX, controls the size of the final sparse grid.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int TEST_NUM, the number of points to be tested.
        //
        //    Input, int TEST_VAL[DIM_NUM*TEST_NUM], the indices of the points 
        //    to be tested.  Normally, each index would be between 0 and 2^LEVEL_MAX.
        //
        //    Output, int ABSCISSA_LEVEL_ND[TEST_NUM], the value of LEVEL at which the
        //    point would first be generated, assuming that a standard sequence of
        //    nested grids is used.
        //
    {
        int j;
        int order;
        int[] test_level;

        test_level = new int[test_num];

        switch (level_max)
        {
            case 0:
            {
                for (j = 0; j < test_num; j++)
                {
                    test_level[j] = 0;
                }

                return test_level;
            }
        }

        order = (int)Math.Pow(2, level_max) + 1;

        for (j = 0; j < test_num; j++)
        {
            test_level[j] = ClenshawCurtis.index_to_level_closed(dim_num, test_val,
                order, level_max,  tIndex: + j * dim_num);
        }

        return test_level;
    }

    public static int[] abscissa_level_open_nd(int level_max, int dim_num, int test_num,
            int[] test_val)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ABSCISSA_LEVEL_OPEN_ND: first level at which given abscissa is generated.
        //
        //  Discussion:
        //
        //    We assume an underlying product grid.  In each dimension, this product
        //    grid has order 2**(LEVEL_MAX+1) - 1.
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
        //    that the product grid is the product of 1D grids,
        //    that the 1D grids are built up by levels, having
        //    orders (total number of points ) 1, 3, 7, 15, 31 and so on,
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
        //  Modified:
        //
        //    19 April 2007
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
        //    Input, int LEVEL_MAX, controls the size of the final sparse grid.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int TEST_NUM, the number of points to be tested.
        //
        //    Input, int TEST_VAL[DIM_NUM*TEST_NUM], the indices of the points 
        //    to be tested.  Normally, each index would be between 0 and 2**LEVEL_MAX.
        //
        //    Output, int ABSCISSA_OPEN_LEVEL_ND[TEST_NUM], the value of LEVEL at which the
        //    point would first be generated, assuming that a standard sequence of
        //    nested grids is used.
        //
    {
        int j;
        int order;
        int[] test_level;

        test_level = new int[test_num];

        switch (level_max)
        {
            case 0:
            {
                for (j = 0; j < test_num; j++)
                {
                    test_level[j] = 0;
                }

                return test_level;
            }
        }

        order = (int)Math.Pow(2, level_max) + 1;

        for (j = 0; j < test_num; j++)
        {
            test_level[j] = LevelToOrder.index_to_level_open(dim_num, test_val,
                order, level_max,  + j * dim_num);
        }

        return test_level;
    }
}