﻿using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Composition;
using Burkardt.Grid;
using Burkardt.Types;

namespace Burkardt.Sparse
{
public static class Grid_ClenshawCurtis
{
public static void sparse_grid_cc ( int dim_num, int level_max, int point_num, 
ref double[] grid_weight, ref double[] grid_point )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_CC computes a sparse grid of Clenshaw Curtis points.
//
//  Discussion:
//
//    This program computes a quadrature rule and writes it to a file.
//
//    The quadrature rule is associated with a sparse grid derived from
//    a Smolyak construction using a closed 1D quadrature rule. 
//
//    The user specifies:
//    * the spatial dimension of the quadrature region,
//    * the level that defines the Smolyak grid.
//    * the closed 1D quadrature rule (Clenshaw-Curtis or Newton-Cotes Closed).
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
//    Input, int POINT_NUM, the number of points in the grid, as determined
//    by SPARSE_GRID_CC_SIZE.
//
//    Output, double GRID_WEIGHTS[POINT_NUM], the weights.
//
//    Output, double GRID_POINTS[DIM_NUM*POINT_NUM], the points.
//
{
int dim;
int[] grid_index;
int order_max;
int point;
//
//  Determine the index vector, relative to the full product grid,
//  that identifies the points in the sparse grid.
//
grid_index = sparse_grid_cc_index ( dim_num, level_max, point_num );
//
//  Compute the physical coordinates of the abscissas.
//
if ( 0 == level_max )
{
order_max = 1;
}
else
{
order_max = (int)Math.Pow ( 2, level_max ) + 1;
}

for ( point = 0; point < point_num; point++ )
{
for ( dim = 0; dim < dim_num; dim++ )
{
grid_point[dim+point*dim_num] = 
ClenshawCurtisGrid.cc_abscissa ( order_max, grid_index[dim+point*dim_num] + 1 );
}
}
//
//  Gather the weights.
//
sparse_grid_cc_weights ( dim_num, level_max, point_num, grid_index, 
ref grid_weight );
}

public static int[] sparse_grid_cc_index ( int dim_num, int level_max, int point_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_CC_INDEX indexes the points forming a sparse grid.
//
//  Discussion:
//
//    The points forming the sparse grid are guaranteed to be a subset
//    of a certain product grid.  The product grid is formed by DIM_NUM
//    copies of a 1D rule of fixed order.  The orders of the 1D rule,
//    (called ORDER_1D) and the order of the product grid, (called ORDER)
//    are determined from the value LEVEL_MAX.
//
//    Thus, any point in the product grid can be identified by its grid index,
//    a set of DIM_NUM indices, each between 1 and ORDER_1D.
//
//    This routine creates the GRID_INDEX array, listing (uniquely) the
//    points of the sparse grid.  
//
//    An assumption has been made that the 1D rule is closed (includes
//    the interval endpoints) and nested (points that are part of a rule
//    of a given level will be part of every rule of higher level).
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
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Input, int POINT_NUM, the total number of points in the grids.
//
//    Output, int SPARSE_GRID_CC_INDEX[DIM_NUM*POINT_NUM], a list of point 
//    indices, representing a subset of the product grid of level LEVEL_MAX,
//    representing (exactly once) each point that will show up in a
//    sparse grid of level LEVEL_MAX.
//
{
int dim;
int[] grid_index;
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

grid_index = new int[dim_num*point_num];
//
//  The outer loop generates LEVELs from 0 to LEVEL_MAX.
//
point_num2 = 0;

level_1d = new int[dim_num];
order_1d = new int[dim_num];

for ( level = 0; level <= level_max; level++ )
{
//
//  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
//  that adds up to LEVEL.
//
more = false;
h = 0;
t = 0;

for ( ; ; )
{
Comp.comp_next ( level, dim_num, ref level_1d, ref more, ref h, ref t );
//
//  Transform each 1D level to a corresponding 1D order.
//
ClenshawCurtis.level_to_order_closed ( dim_num, level_1d, ref order_1d );
//
//  The product of the 1D orders gives us the number of points in this grid.
//
order_nd = typeMethods.i4vec_product ( dim_num, order_1d );
//
//  The inner (hidden) loop generates all points corresponding to given grid.
//
grid_index2 = Multigrid.multigrid_index0 ( dim_num, order_1d, order_nd );
//
//  Adjust these grid indices to reflect LEVEL_MAX.
//
Multigrid.multigrid_scale_closed ( dim_num, order_nd, level_max, level_1d, 
grid_index2 );
//
//  Determine the first level of appearance of each of the points.
//
grid_level = Abscissa.abscissa_level_closed_nd ( level_max, dim_num, order_nd, 
grid_index2 );
//
//  Only keep those points which first appear on this level.
//
for ( point = 0; point < order_nd; point++ )
{
if ( grid_level[point] == level )
{
for ( dim = 0; dim < dim_num; dim++ )
{
grid_index[dim+point_num2*dim_num] =
grid_index2[dim+point*dim_num];
}
point_num2 = point_num2 + 1;
}
}

if ( !more )
{
break;
}
}
}

return grid_index;
}

public static int sparse_grid_cc_size_old ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_CC_SIZE_OLD sizes a sparse grid of Clenshaw Curtis points.
//
//  Discussion:
//
//    This function has been replaced by a new version which is much faster.
//
//    This version is retained for historical interest.
//
//    The grid is defined as the sum of the product rules whose LEVEL
//    satisfies:
//
//      0 <= LEVEL <= LEVEL_MAX.
//
//    This routine works on an abstract set of nested grids.
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
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Output, int SPARSE_GRID_CC_SIZE, the number of points in the grid.
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
if ( level_max == 0 )
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

for ( level = 0; level <= level_max; level++ )
{
//
//  The middle loop generates the next partition that adds up to LEVEL.
//
more = false;
h = 0;
t = 0;

for ( ; ; )
{
Comp.comp_next ( level, dim_num, ref level_1d, ref more, ref h, ref t );
//
//  Transform each 1D level to a corresponding 1D order.
//
ClenshawCurtis.level_to_order_closed ( dim_num, level_1d, ref order_1d );
//
//  The product of the 1D orders gives us the number of points in this grid.
//
order_nd = typeMethods.i4vec_product ( dim_num, order_1d );
//
//  The inner (hidden) loop generates all points corresponding to given grid.
//
grid_index = Multigrid.multigrid_index0 ( dim_num, order_1d, order_nd );
//
//  Adjust these grid indices to reflect LEVEL_MAX.
//
Multigrid.multigrid_scale_closed ( dim_num, order_nd, level_max, level_1d, 
grid_index );
//
//  Determine the first level of appearance of each of the points.
//
grid_level = Abscissa.abscissa_level_closed_nd ( level_max, dim_num, order_nd, 
grid_index );
//
//  Only keep those points which first appear on this level.
//
for ( point = 0; point < order_nd; point++ )
{
if ( grid_level[point] == level )
{
point_num = point_num + 1;
}
}

if ( !more )
{
break;
}
}
}

return point_num;
}

public static void sparse_grid_cc_weights ( int dim_num, int level_max, int point_num, 
int[] grid_index, ref double[] grid_weight )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_CC_WEIGHTS gathers the weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2013
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
//    Output, double GRID_WEIGHT[POINT_NUM], the weights
//    associated with the sparse grid points.
//
{
bool all_equal;
int coeff;
int dim;
bool found;
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

if ( level_max == 0 )
{
for ( point = 0; point < point_num; point++ )
{ 
grid_weight[point] = Math.Pow ( 2.0, dim_num );
}
return;
}

level_1d = new int[dim_num];
order_1d = new int[dim_num];

for ( point = 0; point < point_num; point++ )
{ 
grid_weight[point] = 0.0;
}

level_min = Math.Max ( 0, level_max + 1 - dim_num );

for ( level = level_min; level <= level_max; level++ )
{
//
//  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
//  that adds up to LEVEL.
//
more = false;
h = 0;
t = 0;

for ( ; ; )
{
Comp.comp_next ( level, dim_num, ref level_1d, ref more, ref h, ref t );
//
//  Transform each 1D level to a corresponding 1D order.
//
ClenshawCurtis.level_to_order_closed ( dim_num, level_1d, ref order_1d );
//
//  The product of the 1D orders gives us the number of points in this grid.
//
order_nd = typeMethods.i4vec_product ( dim_num, order_1d );
//
//  Generate the indices of the points corresponding to the grid.
//
grid_index2 = Multigrid.multigrid_index0 ( dim_num, order_1d, order_nd );
//
//  Compute the weights for this grid.
//
grid_weight2 = ClenshawCurtis.product_weights_cc ( dim_num, order_1d, order_nd );
//
//  Adjust the grid indices to reflect LEVEL_MAX.
//
Multigrid.multigrid_scale_closed ( dim_num, order_nd, level_max, level_1d, 
grid_index2 );
//
//  Now determine the coefficient.
//
coeff = typeMethods.i4_mop ( level_max - level ) 
* typeMethods.i4_choose ( dim_num - 1, level_max - level );

for ( point2 = 0; point2 < order_nd; point2++ )
{
found = false;

for ( point = 0; point < point_num; point++ )
{
all_equal = true;
for ( dim = 0; dim < dim_num; dim++ )
{
if ( grid_index2[dim+point2*dim_num] !=
grid_index[dim+point*dim_num] )
{
all_equal = false;
break;
}
}
if ( all_equal )
{
grid_weight[point] = grid_weight[point] 
+ ( double ) ( coeff ) * grid_weight2[point2];
found = true;
break;
}
}
if ( !found )
{
Console.WriteLine("");
Console.WriteLine("SPARSE_GRID_CC_WEIGHTS - Fatal error!");
Console.WriteLine("  Could not find a match for a point.");
return;
}
}

if ( !more )
{
break;
}
}
}
}

public static int sparse_grid_ccs_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_CCS_SIZE sizes a sparse grid using Clenshaw Curtis Slow rules.
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
int l;
int level;
int[] level_1d;
bool more;
int[] new_1d;
int o;
int p;
int point_num;
int t;
int v;
//
//  Special case.
//
if ( level_max < 0 )
{
point_num = 0;
return point_num;
}

if ( level_max == 0 )
{
point_num = 1;
return point_num;
}
//
//  Construct the vector that counts the new points in the 1D rule.
//
new_1d = new int[level_max+1];

new_1d[0] = 1;
new_1d[1] = 2;

p = 3;
o = 3;

for ( l = 2; l <= level_max; l++ )
{
p = 2 * l + 1;
if ( o < p )
{
new_1d[l] = o - 1;
o = 2 * o - 1;
}
else
{
new_1d[l] = 0;
}
}
//
//  Count the number of points by counting the number of new points 
//  associated with each level vector.
//
level_1d = new int[dim_num];

point_num = 0;

for ( level = 0; level <= level_max; level++ )
{
more = false;
h = 0;
t = 0;

for ( ; ;)
{
Comp.comp_next ( level, dim_num, ref level_1d, ref more, ref h, ref t );

v = 1;
for ( dim = 0; dim < dim_num; dim++ )
{
v = v * new_1d[level_1d[dim]];
}

point_num = point_num + v;

if ( !more )
{
break;
}
}
}

return point_num;
}

public static int sparse_grid_cfn_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_CFN_SIZE sizes a sparse grid using Closed Fully Nested rules.
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
//    Output, int SPARSE_GRID_CFN_SIZE, the number of points in the grid.
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
if ( level_max < 0 )
{
point_num = 0;
return point_num;
}

if ( level_max == 0 )
{
point_num = 1;
return point_num;
}
//
//  Construct the vector that counts the new points in the 1D rule.
//
new_1d = new int[level_max+1];

new_1d[0] = 1;
new_1d[1] = 2;

j = 1;
for ( l = 2; l <= level_max; l++ )
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

for ( level = 0; level <= level_max; level++ )
{
more = false;
h = 0;
t = 0;

for ( ; ;)
{
Comp.comp_next ( level, dim_num, ref level_1d, ref more, ref h, ref t );

v = 1;
for ( dim = 0; dim < dim_num; dim++ )
{
v = v * new_1d[level_1d[dim]];
}

point_num = point_num + v;

if ( !more )
{
break;
}
}
}

return point_num;
}
}
}