# include "sandia_rules.hpp"
# include "sgmg.hpp"

# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>

namespace webbur 
{
//****************************************************************************80

void product_mixed_growth_weight ( int dim_num, int order_1d[], int order_nd, 
  int rule[], int np[], double p[], 
  void ( *gw_compute_weights[] ) ( int order, int np, double p[], double w[] ),
  double weight_nd[] )

//****************************************************************************80
//
//  Purpose:
//
//    PRODUCT_MIXED_GROWTH_WEIGHT computes the weights of a mixed product rule.
//
//  Discussion:
//
//    This routine computes the weights for a quadrature rule which is
//    a product of 1D rules of varying order and kind.
//
//    The user must preallocate space for the output array WEIGHT_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 June 2010
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
//    Input, int ORDER_1D[DIM_NUM], the order of the 1D rules.
//
//    Input, int ORDER_ND, the order of the product rule.
//
//    Input, int RULE[DIM_NUM], the rule in each dimension.
//     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
//     2, "F2",  Fejer Type 2, Open Fully Nested.
//     3, "GP",  Gauss Patterson, Open Fully Nested.
//     4, "GL",  Gauss Legendre, Open Weakly Nested.
//     5, "GH",  Gauss Hermite, Open Weakly Nested.
//     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
//     7, "LG",  Gauss Laguerre, Open Non Nested.
//     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
//     9, "GJ",  Gauss Jacobi, Open Non Nested.
//    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
//    11, "UO",  User supplied Open, presumably Non Nested.
//    12, "UC",  User supplied Closed, presumably Non Nested.
//
//    Input, int NP[DIM_NUM], the number of parameters used by each rule.
//
//    Input, double P[sum(NP[*])], the parameters needed by each rule.
//
//    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int np, double p[], double w[] ),
//    an array of pointers to functions which return the 1D quadrature weights 
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
//    Output, double WEIGHT_ND[ORDER_ND], the product rule weights.
//
{
  int dim;
  int i;
  int p_index;
  double *weight_1d;

  for ( i = 0; i < order_nd; i++ )
  {
    weight_nd[i] = 1.0;
  }

  p_index = 0;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    weight_1d = new double[order_1d[dim]];

    if ( rule[dim] == 1 )
    {
      webbur::clenshaw_curtis_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 2 )
    {
      webbur::fejer2_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 3 )
    {
      webbur::patterson_lookup_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 4 )
    {
      webbur::legendre_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 5 )
    {
      webbur::hermite_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 6 )
    {
      webbur::gen_hermite_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 7 )
    {
      webbur::laguerre_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 8 )
    {
      webbur::gen_laguerre_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 9 )
    {
      webbur::jacobi_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 10 )
    {
      webbur::hermite_genz_keister_lookup_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 11 )
    {
      gw_compute_weights[dim] ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 12 )
    {
      gw_compute_weights[dim] ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else
    {
      std::cerr << "\n";
      std::cerr << "PRODUCT_MIXED_GROWTH_WEIGHT - Fatal error!\n";
      std::cerr << "  Unexpected value of RULE[" << dim << "] = " 
           << rule[dim] << ".\n";
      std::exit ( 1 );
    }

    p_index = p_index + np[dim];

    webbur::r8vec_direct_product2 ( dim, order_1d[dim], weight_1d, 
      dim_num, order_nd, weight_nd );

    delete [] weight_1d;
  }
  return;
}
//****************************************************************************80

void sgmg_index ( int dim_num, int level_max, int rule[], 
  int point_num, int point_total_num, int sparse_unique_index[], int growth[],
  int sparse_order[], int sparse_index[] )

//****************************************************************************80
//
//  Purpose:
//
//    SGMG_INDEX indexes a sparse grid of mixed 1D rules.
//
//  Discussion:
//
//    For each "unique" point in the sparse grid, we return its INDEX and ORDER.
//
//    That is, for the I-th unique point P, we determine the product grid which
//    first generated this point, and  and we return in SPARSE_ORDER the orders 
//    of the 1D rules in that grid, and  and in SPARSE_INDEX the component 
//    indexes in those rules that generated this specific point.
//
//    For instance, say P was first generated by a rule which was a 3D product
//    of a 9th order CC rule and  and a 15th order GL rule, and  and that to 
//    generate P, we used the 7-th point of the CC rule and  and the 3rh point 
//    of the GL rule.  Then the SPARSE_ORDER information would be (9,15) and
//    the SPARSE_INDEX information would be (7,3).  This, combined with the 
//    information in RULE, is enough to regenerate the value of P.
//
//    The user must preallocate space for the output arrays SPARSE_ORDER and
//    SPARSE_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 June 2010
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
//    Input, int RULE[DIM_NUM], the rule in each dimension.
//     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
//     2, "F2",  Fejer Type 2, Open Fully Nested.
//     3, "GP",  Gauss Patterson, Open Fully Nested.
//     4, "GL",  Gauss Legendre, Open Weakly Nested.
//     5, "GH",  Gauss Hermite, Open Weakly Nested.
//     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
//     7, "LG",  Gauss Laguerre, Open Non Nested.
//     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
//     9, "GJ",  Gauss Jacobi, Open Non Nested.
//    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
//    11, "UO",  User supplied Open, presumably Non Nested.
//    12, "UC",  User supplied Closed, presumably Non Nested.
//
//    Input, int POINT_NUM, the number of unique points 
//    in the grid. 
//
//    Input, int POINT_TOTAL_NUM, the total number of points in the grid.
//
//    Input, int SPARSE_UNIQUE_INDEX[POINT_TOTAL_NUM], associates each
//    point in the grid with its unique representative.
//
//    Input, int GROWTH[DIM_NUM], the growth rule in each dimension. 
//    0, "DF", default growth associated with this quadrature rule;
//    1, "SL", slow linear, L+1;
//    2  "SO", slow linear odd, O=1+2((L+1)/2)
//    3, "ML", moderate linear, 2L+1;
//    4, "SE", slow exponential;
//    5, "ME", moderate exponential;
//    6, "FE", full exponential.
//
//    Output, int SPARSE_ORDER[DIM_NUM*POINT_NUM], lists, 
//    for each point, the order of the 1D rules used in the grid that 
//    generated it.
//
//    Output, int SPARSE_INDEX[DIM_NUM*POINT_NUM)] lists, for 
//    each point, its index in each of the 1D rules in the grid that generated 
//    it.  The indices are 1-based.
//
{
  int dim;
  int h;
  int level;
  int *level_1d;
  int level_min;
  bool more_grids;
  bool more_points;
  int *order_1d;
  int point;
  int point_count;
  int *point_index;
  int point_unique;
  int t;
//
//  Special cases.
//
  if ( level_max < 0 )
  {
    return;
  }

  if ( level_max == 0 )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_order[dim+0*dim_num] = 1;
      sparse_index[dim+0*dim_num] = 1;
    }
    return;
  }

  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_order[dim+point*dim_num] = -1;
      sparse_index[dim+point*dim_num] = -1;
    }
  }

  point_count = 0;
//
//  The outer loop generates values of LEVEL.
//
  level_1d = new int[dim_num];
  order_1d = new int[dim_num];
  point_index = new int[dim_num];

  level_min = webbur::i4_max ( 0, level_max + 1 - dim_num );

  for ( level = level_min; level <= level_max; level++ )
  {
//
//  The middle loop generates a GRID, 
//  based on the next partition that adds up to LEVEL.
//
    more_grids = false;
    h = 0;
    t = 0;

    for ( ; ; )
    {
      webbur::comp_next ( level, dim_num, level_1d, &more_grids, &h, &t );

      webbur::level_growth_to_order ( dim_num, level_1d, rule, growth, order_1d );
//
//  The inner loop generates a POINT of the GRID of the LEVEL.
//
      more_points = false;

      for ( ; ; )
      {
        webbur::vec_colex_next3 ( dim_num, order_1d, point_index, &more_points );

        if ( !more_points )
        {
          break;
        }
        point_unique = sparse_unique_index[point_count];
        for ( dim = 0; dim < dim_num; dim++ )
        {
          sparse_order[dim+point_unique*dim_num] = order_1d[dim];
        }
        for ( dim = 0; dim < dim_num; dim++ )
        {
          sparse_index[dim+point_unique*dim_num] = point_index[dim];
        }
        point_count = point_count + 1;
      }

      if ( !more_grids )
      {
        break;
      }
    }
  }

  delete [] level_1d;
  delete [] order_1d;
  delete [] point_index;

  return;
}
//****************************************************************************80

void sgmg_point ( int dim_num, int level_max, int rule[], 
  int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double x[] ),
  int point_num, int sparse_order[], int sparse_index[], int growth[],
  double sparse_point[] )

//****************************************************************************80
//
//  Purpose:
//
//    SGMG_POINT computes the points of a sparse grid rule.
//
//  Discussion:
//
//    The sparse grid is the logical sum of low degree product rules.
//
//    Each product rule is the product of 1D factor rules.
//
//    The user specifies:
//    * the spatial dimension of the quadrature region,
//    * the level that defines the Smolyak grid.
//    * the quadrature rules.
//    * the number of points.
//
//    The user must preallocate space for the output array SPARSE_POINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 March 2010
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
//    Input, int RULE[DIM_NUM], the rule in each dimension.
//     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
//     2, "F2",  Fejer Type 2, Open Fully Nested.
//     3, "GP",  Gauss Patterson, Open Fully Nested.
//     4, "GL",  Gauss Legendre, Open Weakly Nested.
//     5, "GH",  Gauss Hermite, Open Weakly Nested.
//     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
//     7, "LG",  Gauss Laguerre, Open Non Nested.
//     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
//     9, "GJ",  Gauss Jacobi, Open Non Nested.
//    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
//    11, "UO",  User supplied Open, presumably Non Nested.
//    12, "UC",  User supplied Closed, presumably Non Nested.
//
//    Input, int NP[DIM_NUM], the number of parameters used by each rule.
//
//    Input, double P[sum(NP[*])], the parameters needed by each rule.
//
//    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int np, double p[], double x[] ),
//    an array of pointers to functions which return the 1D quadrature points 
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
//    Input, int POINT_NUM, the number of points in the grid,
//    as determined by SGMG_SIZE.
//
//    Input, int SPARSE_ORDER[DIM_NUM*POINT_NUM], lists, for each point,
//    the order of the 1D rules used in the grid that generated it.
//
//    Input, int SPARSE_INDEX[DIM_NUM*POINT_NUM], lists, for each point,
//    its index in each of the 1D rules in the grid that generated it.
//    The indices are 1-based.
//
//    Input, int GROWTH[DIM_NUM], the growth rule in each dimension. 
//    0, "DF", default growth associated with this quadrature rule;
//    1, "SL", slow linear, L+1;
//    2  "SO", slow linear odd, O=1+2((L+1)/2)
//    3, "ML", moderate linear, 2L+1;
//    4, "SE", slow exponential;
//    5, "ME", moderate exponential;
//    6, "FE", full exponential.
//
//    Output, double SPARSE_POINT[DIM_NUM*POINT_NUM], the points.
//
{
  int dim;
  int level;
  int order;
  int p_index;
  int point;
  double *points;

  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_point[dim+point*dim_num] = webbur::r8_huge ( );
    }
  }
//
//  Compute the point coordinates.
//
  p_index = 0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    for ( level = 0; level <= level_max; level++ )
    {
      webbur::level_growth_to_order ( 1, &level, rule+dim, growth+dim, &order );

      points = new double[order];

      if ( rule[dim] == 1 )
      {
        webbur::clenshaw_curtis_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 2 )
      {
        webbur::fejer2_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 3 )
      {
        webbur::patterson_lookup_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 4 )
      {
        webbur::legendre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 5 )
      {
        webbur::hermite_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 6 )
      {
        webbur::gen_hermite_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 7 )
      {
        webbur::laguerre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 8 )
      {
        webbur::gen_laguerre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 9 )
      {
        webbur::jacobi_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 10 )
      {
        webbur::hermite_genz_keister_lookup_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 11 )
      {
        gw_compute_points[dim] (
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 12 )
      {
        gw_compute_points[dim] (
          order, np[dim], p+p_index, points );
      }
      else
      {
        std::cerr << "\n";
        std::cerr << "SGMG_POINT - Fatal error!\n";
        std::cerr << "  Unexpected value of RULE[" << dim << "] = " 
             << rule[dim] << ".\n";
        std::exit ( 1 );
      }

      for ( point = 0; point < point_num; point++ )
      {
        if ( sparse_order[dim+point*dim_num] == order )
        {
          sparse_point[dim+point*dim_num] = 
            points[sparse_index[dim+point*dim_num]-1];
        }
      }
      delete [] points;
    }
    p_index = p_index + np[dim];
  }

  return;
}
//****************************************************************************80

int sgmg_size ( int dim_num, int level_max, int rule[], 
  int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double x[] ),
  double tol, int growth[] )

//****************************************************************************80
//
//  Purpose:
//
//    SGMG_SIZE sizes a sparse grid, discounting duplicates.
//
//  Discussion:
//
//    The sparse grid is the logical sum of product grids with total LEVEL
//    between LEVEL_MIN and LEVEL_MAX.
//
//    Depending on the 1D rules involved, there may be many duplicate points
//    in the sparse grid.
//
//    This routine counts the unique points in the sparse grid.  It does this
//    in a straightforward way, by actually generating all the points, and
//    comparing them, with a tolerance for equality.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 July 2010
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
//    Input, int RULE[DIM_NUM], the rule in each dimension.
//     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
//     2, "F2",  Fejer Type 2, Open Fully Nested.
//     3, "GP",  Gauss Patterson, Open Fully Nested.
//     4, "GL",  Gauss Legendre, Open Weakly Nested.
//     5, "GH",  Gauss Hermite, Open Weakly Nested.
//     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
//     7, "LG",  Gauss Laguerre, Open Non Nested.
//     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
//     9, "GJ",  Gauss Jacobi, Open Non Nested.
//    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
//    11, "UO",  User supplied Open, presumably Non Nested.
//    12, "UC",  User supplied Closed, presumably Non Nested.
//
//    Input, int NP[DIM_NUM], the number of parameters used by each rule.
//
//    Input, double P[sum(NP[*])], the parameters needed by each rule.
//
//    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int np, double p[], double x[] ),
//    an array of pointers to functions which return the 1D quadrature points 
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
//    Input, double TOL, a tolerance for point equality.
//
//    Input, int GROWTH[DIM_NUM], the growth rule in each dimension. 
//    0, "DF", default growth associated with this quadrature rule;
//    1, "SL", slow linear, L+1;
//    2  "SO", slow linear odd, O=1+2((L+1)/2)
//    3, "ML", moderate linear, 2L+1;
//    4, "SE", slow exponential;
//    5, "ME", moderate exponential;
//    6, "FE", full exponential.
//
//    Output, int SGMG_SIZE, the number of unique points.
//
{
  int dim;
  int h;
  int level;
  int *level_1d;
  int level_min;
  bool more_grids;
  bool more_points;
  int order;
  int *order_1d;
  int p_index;
  int point;
  int *point_index;
  int point_num;
  int point_total_num;
  int point_total_num2;
  double *points;
  int seed;
  int *sparse_total_index;
  int *sparse_total_order;
  double *sparse_total_point;
  int t;
//
//  Special cases.
//
  if ( level_max < 0 )
  {
    point_num = -1;
    return point_num;
  }

  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  Get total number of points, including duplicates.
//
  point_total_num = webbur::sgmg_size_total ( dim_num, 
    level_max, rule, growth );
//
//  Generate SPARSE_TOTAL_ORDER and SPARSE_TOTAL_INDEX arrays 
//  for the TOTAL set of points.
//
  sparse_total_order = new int[dim_num*point_total_num];
  sparse_total_index = new int[dim_num*point_total_num];

  point_total_num2 = 0;
//
//  The outer loop generates values of LEVEL.
//
  level_1d = new int[dim_num];
  order_1d = new int[dim_num];
  point_index = new int[dim_num];

  level_min = webbur::i4_max ( 0, level_max + 1 - dim_num );

  for ( level = level_min; level <= level_max; level++ )
  {
//
//  The middle loop generates a GRID, 
//  based on the next partition that adds up to LEVEL.
//
    more_grids = false;
    h = 0;
    t = 0;

    for ( ; ; )
    {
      webbur::comp_next ( level, dim_num, level_1d, &more_grids, &h, &t );

      webbur::level_growth_to_order ( dim_num, level_1d, rule, growth, order_1d );
//
//  The inner loop generates a POINT of the GRID of the LEVEL.
//
      more_points = false;

      for ( ; ; )
      {
        webbur::vec_colex_next3 ( dim_num, order_1d, point_index, &more_points );

        if ( !more_points )
        {
          break;
        }
        for ( dim = 0; dim < dim_num; dim++ )
        {
          sparse_total_order[dim+point_total_num2*dim_num] = order_1d[dim];
        }
        for ( dim = 0; dim < dim_num; dim++ )
        {
          sparse_total_index[dim+point_total_num2*dim_num] = point_index[dim];
        }
        point_total_num2 = point_total_num2 + 1;
      }

      if ( !more_grids )
      {
        break;
      }
    }
  }
  delete [] level_1d;
  delete [] order_1d;
  delete [] point_index;
//
//  Now compute the coordinates of the TOTAL set of points.
//
  sparse_total_point = new double[dim_num*point_total_num];

  for ( point = 0; point < point_total_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_total_point[dim+point*dim_num] = webbur::r8_huge ( );
    }
  }
//
//  Compute the point coordinates.
//
  p_index = 0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    for ( level = 0; level <= level_max; level++ )
    {
      webbur::level_growth_to_order ( 1, &level, rule+dim, growth+dim, &order );

      points = new double[order];

      if ( rule[dim] == 1 )
      {
        webbur::clenshaw_curtis_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 2 )
      {
        webbur::fejer2_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 3 )
      {
        webbur::patterson_lookup_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 4 )
      {
        webbur::legendre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 5 )
      {
        webbur::hermite_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 6 )
      {
        webbur::gen_hermite_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 7 )
      {
        webbur::laguerre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 8 )
      {
        webbur::gen_laguerre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 9 )
      {
        webbur::jacobi_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 10 )
      {
        webbur::hermite_genz_keister_lookup_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 11 )
      {
        gw_compute_points[dim] (
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 12 )
      {
        gw_compute_points[dim] (
          order, np[dim], p+p_index, points );
      }
      else
      {
        std::cerr << "\n";
        std::cerr << "SGMG_SIZE - Fatal error!\n";
        std::cerr << "  Unexpected value of RULE[" << dim << "] = " 
             << rule[dim] << ".\n";
        std::exit ( 1 );
      }

      for ( point = 0; point < point_total_num; point++ )
      {
        if ( sparse_total_order[dim+point*dim_num] == order )
        {
          sparse_total_point[dim+point*dim_num] = 
            points[sparse_total_index[dim+point*dim_num]-1];
        }
      }
      delete [] points;
    }
    p_index = p_index + np[dim];
  }
//
//  Count the tolerably unique columns. 
//
  seed = 123456789;

  point_num = webbur::point_radial_tol_unique_count ( dim_num, point_total_num, 
    sparse_total_point, tol, &seed );

  delete [] sparse_total_index;
  delete [] sparse_total_order;
  delete [] sparse_total_point;

  return point_num;
}
//****************************************************************************80

int sgmg_size_total ( int dim_num, int level_max, int rule[], int growth[] )

//****************************************************************************80
//
//  Purpose:
//
//    SGMG_SIZE_TOTAL sizes a sparse grid, counting duplicates.
//
//  Discussion:
//
//    The sparse grid is the logical sum of product grids with total LEVEL
//    between LEVEL_MIN and LEVEL_MAX.
//
//    In some cases, the same point may occur in different product grids
//    used to form the sparse grid.
//
//    This routine counts the total number of points used to construct the 
//    sparse grid; if the same point occurs several times, each occurrence 
//    is added to the sum.
//
//    This computation is useful in order to be able to allocate enough
//    space for the full set of points, before they are compressed by removing
//    duplicates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 June 2010
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
//    Input, int RULE[DIM_NUM], the rule in each dimension.
//     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
//     2, "F2",  Fejer Type 2, Open Fully Nested.
//     3, "GP",  Gauss Patterson, Open Fully Nested.
//     4, "GL",  Gauss Legendre, Open Weakly Nested.
//     5, "GH",  Gauss Hermite, Open Weakly Nested.
//     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
//     7, "LG",  Gauss Laguerre, Open Non Nested.
//     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
//     9, "GJ",  Gauss Jacobi, Open Non Nested.
//    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
//    11, "UO",  User supplied Open, presumably Non Nested.
//    12, "UC",  User supplied Closed, presumably Non Nested.
//
//    Input, int GROWTH[DIM_NUM], the growth rule in each dimension. 
//    0, "DF", default growth associated with this quadrature rule;
//    1, "SL", slow linear, L+1;
//    2  "SO", slow linear odd, O=1+2((L+1)/2)
//    3, "ML", moderate linear, 2L+1;
//    4, "SE", slow exponential;
//    5, "ME", moderate exponential;
//    6, "FE", full exponential.
//
//    Output, int SGMG_SIZE_TOTAL, the number of points
//    including repetitions.
//
{
  int h;
  int level;
  int *level_1d;
  int level_min;
  bool more_grids;
  int *order_1d;
  int point_total_num;
  int t;
//
//  Special case.
//
  if ( level_max == 0 )
  {
    point_total_num = 1;
    return point_total_num;
  }

  point_total_num = 0;

  level_1d = new int[dim_num];
  order_1d = new int[dim_num];
//
//  The outer loop generates values of LEVEL.
//
  level_min = webbur::i4_max ( 0, level_max + 1 - dim_num );

  for ( level = level_min; level <= level_max; level++ )
  {
//
//  The middle loop generates a GRID, 
//  based on the next partition that adds up to LEVEL.
//
    more_grids = false;
    h = 0;
    t = 0;

    for ( ; ; )
    {
      webbur::comp_next ( level, dim_num, level_1d, &more_grids, &h, &t );

      webbur::level_growth_to_order ( dim_num, level_1d, rule, growth, order_1d );

      point_total_num = point_total_num + webbur::i4vec_product ( dim_num, 
        order_1d );

      if ( !more_grids )
      {
        break;
      }
    }
  }
  delete [] level_1d;
  delete [] order_1d;

  return point_total_num;
}
//****************************************************************************80

void sgmg_unique_index ( int dim_num, int level_max, 
  int rule[], int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double x[] ),
  double tol, int point_num, int point_total_num, int growth[], 
  int sparse_unique_index[] )

//****************************************************************************80
//
//  Purpose:
//
//    SGMG_UNIQUE_INDEX maps nonunique to unique points.
//
//  Discussion:
//
//    The sparse grid usually contains many points that occur in more
//    than one product grid.
//
//    When generating the point locations, it is easy to realize that a point
//    has already been generated.
//
//    But when it's time to compute the weights of the sparse grids, it is
//    necessary to handle situations in which weights corresponding to 
//    the same point generated in multiple grids must be collected together.
//
//    This routine generates ALL the points, including their multiplicities,
//    and figures out a mapping from them to the collapsed set of unique points.
//
//    This mapping can then be used during the weight calculation so that
//    a contribution to the weight gets to the right place.
//
//    The user must preallocate space for the output array SPARSE_UNIQUE_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2010
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
//    Input, int RULE[DIM_NUM], the rule in each dimension.
//     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
//     2, "F2",  Fejer Type 2, Open Fully Nested.
//     3, "GP",  Gauss Patterson, Open Fully Nested.
//     4, "GL",  Gauss Legendre, Open Weakly Nested.
//     5, "GH",  Gauss Hermite, Open Weakly Nested.
//     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
//     7, "LG",  Gauss Laguerre, Open Non Nested.
//     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
//     9, "GJ",  Gauss Jacobi, Open Non Nested.
//    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
//    11, "UO",  User supplied Open, presumably Non Nested.
//    12, "UC",  User supplied Closed, presumably Non Nested.
//
//    Input, int NP[DIM_NUM], the number of parameters used by each rule.
//
//    Input, double P[sum(NP[*])], the parameters needed by each rule.
//
//    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int np, double p[], double x[] ),
//    an array of pointers to functions which return the 1D quadrature points 
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
//    Input, double TOL, a tolerance for point equality.
//
//    Input, int POINT_NUM, the number of unique points 
//    in the grid. 
//
//    Input, int POINT_TOTAL_NUM, the total number of points 
//    in the grid. 
//
//    Input, int GROWTH[DIM_NUM], the growth rule in each dimension. 
//    0, "DF", default growth associated with this quadrature rule;
//    1, "SL", slow linear, L+1;
//    2  "SO", slow linear odd, O=1+2((L+1)/2)
//    3, "ML", moderate linear, 2L+1;
//    4, "SE", slow exponential;
//    5, "ME", moderate exponential;
//    6, "FE", full exponential.
//
//    Output, int SPARSE UNIQUE_INDEX[POINT_TOTAL_NUM], lists, 
//    for each (nonunique) point, the corresponding index of the same point in 
//    the unique listing.
//
{
  int dim;
  int h;
  int level;
  int *level_1d;
  int level_min;
  bool more_grids;
  bool more_points;
  int order;
  int *order_1d;
  int p_index;
  int point;
  int *point_index;
  int point_total_num2;
  double *points;
  int rep;
  int seed;
  int *sparse_total_index;
  int *sparse_total_order;
  double *sparse_total_point;
  int t;
  int *undx;
//
//  Special cases.
//
  if ( level_max < 0 )
  {
    return;
  }

  if ( level_max == 0 )
  {
    sparse_unique_index[0] = 0;
    return;
  }
//
//  Generate SPARSE_TOTAL_ORDER and SPARSE_TOTAL_INDEX arrays 
//  for the TOTAL set of points.
//
  sparse_total_order = new int[dim_num*point_total_num];
  sparse_total_index = new int[dim_num*point_total_num];

  point_total_num2 = 0;
//
//  The outer loop generates values of LEVEL.
//
  level_1d = new int[dim_num];
  order_1d = new int[dim_num];
  point_index = new int[dim_num];

  level_min = webbur::i4_max ( 0, level_max + 1 - dim_num );

  for ( level = level_min; level <= level_max; level++ )
  {
//
//  The middle loop generates a GRID, 
//  based on the next partition that adds up to LEVEL.
//
    more_grids = false;
    h = 0;
    t = 0;

    for ( ; ; )
    {
      webbur::comp_next ( level, dim_num, level_1d, &more_grids, &h, &t );

      webbur::level_growth_to_order ( dim_num, level_1d, rule, growth, order_1d );
//
//  The inner loop generates a POINT of the GRID of the LEVEL.
//
      more_points = false;

      for ( ; ; )
      {
        webbur::vec_colex_next3 ( dim_num, order_1d, point_index, &more_points );

        if ( !more_points )
        {
          break;
        }
        for ( dim = 0; dim < dim_num; dim++ )
        {
          sparse_total_order[dim+point_total_num2*dim_num] = order_1d[dim];
        }
        for ( dim = 0; dim < dim_num; dim++ )
        {
          sparse_total_index[dim+point_total_num2*dim_num] = point_index[dim];
        }
        point_total_num2 = point_total_num2 + 1;
      }

      if ( !more_grids )
      {
        break;
      }
    }
  }
  delete [] level_1d;
  delete [] order_1d;
  delete [] point_index;
//
//  Now compute the coordinates of the TOTAL set of points.
//
  sparse_total_point = new double[dim_num*point_total_num];

  for ( point = 0; point < point_total_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_total_point[dim+point*dim_num] = webbur::r8_huge ( );
    }
  }
//
//  Compute the point coordinates.
//
  p_index = 0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    for ( level = 0; level <= level_max; level++ )
    {
      webbur::level_growth_to_order ( 1, &level, rule+dim, growth+dim, &order );

      points = new double[order];

      if ( rule[dim] == 1 )
      {
        webbur::clenshaw_curtis_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 2 )
      {
        webbur::fejer2_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 3 )
      {
        webbur::patterson_lookup_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 4 )
      {
        webbur::legendre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 5 )
      {
        webbur::hermite_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 6 )
      {
        webbur::gen_hermite_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 7 )
      {
        webbur::laguerre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 8 )
      {
        webbur::gen_laguerre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 9 )
      {
        webbur::jacobi_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 10 )
      {
        webbur::hermite_genz_keister_lookup_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 11 )
      {
        gw_compute_points[dim] (
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 12 )
      {
        gw_compute_points[dim] (
          order, np[dim], p+p_index, points );
      }
      else
      {
        std::cerr << "\n";
        std::cerr << "SGMG_UNIQUE_INDEX - Fatal error!\n";
        std::cerr << "  Unexpected value of RULE[" << dim << "] = "
             << rule[dim] << ".\n";
        std::exit ( 1 );
      }

      for ( point = 0; point < point_total_num; point++ )
      {
        if ( sparse_total_order[dim+point*dim_num] == order )
        {
          sparse_total_point[dim+point*dim_num] = 
            points[sparse_total_index[dim+point*dim_num]-1];
        }
      }
      delete [] points;
    }
    p_index = p_index + np[dim];
  }
//
//  Merge points that are too close.
//
  seed = 123456789;

  undx = new int[point_num];

  webbur::point_radial_tol_unique_index ( dim_num, point_total_num, 
    sparse_total_point, tol, &seed, undx, sparse_unique_index );

  for ( point = 0; point < point_total_num; point++ )
  {
    rep = undx[sparse_unique_index[point]];
    if ( point != rep )
    {
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_total_point[dim+point*dim_num] = sparse_total_point[dim+rep*dim_num];
      }
    }
  }
//
//  Construct an index that indicates the "rank" of the unique points.
//
  webbur::point_unique_index ( dim_num, point_total_num, sparse_total_point,
    point_num, undx, sparse_unique_index );

  delete [] undx;

  delete [] sparse_total_index;
  delete [] sparse_total_order;
  delete [] sparse_total_point;

  return;
}
//****************************************************************************80

void sgmg_weight ( int dim_num, int level_max, int rule[], 
  int np[], double p[], 
  void ( *gw_compute_weights[] ) ( int order, int np, double p[], double w[] ),
  int point_num, int point_total_num, int sparse_unique_index[], int growth[],
  double sparse_weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    SGMG_WEIGHT: sparse grid weights for a mix of 1D rules.
//
//  Discussion:
//
//    The user must preallocate space for the output array SPARSE_WEIGHT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 June 2010
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
//    Input, int RULE[DIM_NUM], the rule in each dimension.
//     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
//     2, "F2",  Fejer Type 2, Open Fully Nested.
//     3, "GP",  Gauss Patterson, Open Fully Nested.
//     4, "GL",  Gauss Legendre, Open Weakly Nested.
//     5, "GH",  Gauss Hermite, Open Weakly Nested.
//     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
//     7, "LG",  Gauss Laguerre, Open Non Nested.
//     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
//     9, "GJ",  Gauss Jacobi, Open Non Nested.
//    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
//    11, "UO",  User supplied Open, presumably Non Nested.
//    12, "UC",  User supplied Closed, presumably Non Nested.
//
//    Input, int NP[DIM_NUM], the number of parameters used by each rule.
//
//    Input, double P[sum(NP[*])], the parameters needed by each rule.
//
//    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int np, double p[], double w[] ),
//    an array of pointers to functions which return the 1D quadrature weights 
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
//    Input, int POINT_NUM, the number of unique points 
//    in the grid. 
//
//    Input, int POINT_TOTAL_NUM, the total number of points 
//    in the grid. 
//
//    Input, int SPARSE UNIQUE_INDEX[POINT_TOTAL_NUM], lists, 
//    for each (nonunique) point, the corresponding index of the same point in 
//    the unique listing.
//
//    Input, int GROWTH[DIM_NUM], the growth rule in each dimension. 
//    0, "DF", default growth associated with this quadrature rule;
//    1, "SL", slow linear, L+1;
//    2  "SO", slow linear odd, O=1+2((L+1)/2)
//    3, "ML", moderate linear, 2L+1;
//    4, "SE", slow exponential;
//    5, "ME", moderate exponential;
//    6, "FE", full exponential.
//
//    Output, double SPARSE_WEIGHT[POINT_NUM], the weights
//    associated with the sparse grid points.
//
{
  double coeff;
  double *grid_weight;
  int h;
  int level;
  int *level_1d;
  int level_min;
  bool more_grids;
  int order;
  int *order_1d;
  int order_nd;
  int point;
  int point_total;
  int point_unique;
  int t;

  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

  for ( point = 0; point < point_num; point++ )
  {
    sparse_weight[point] = 0.0;
  }

  point_total = 0;

  level_min = webbur::i4_max ( 0, level_max + 1 - dim_num );

  for ( level = level_min; level <= level_max; level++ )
  {
//
//  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
//  that adds up to LEVEL.
//
    more_grids = false;
    h = 0;
    t = 0;

    for ( ; ; )
    {
      webbur::comp_next ( level, dim_num, level_1d, &more_grids, &h, &t );
//
//  Transform each 1D level to a corresponding 1D order.
//
      webbur::level_growth_to_order ( dim_num, level_1d, rule, growth, order_1d );
//
//  The product of the 1D orders gives us the number of points in this grid.
//
      order_nd = webbur::i4vec_product ( dim_num, order_1d );
//
//  Compute the weights for this grid.
//
//  The correct transfer of data from the product grid to the sparse grid
//  depends on the fact that the product rule weights are stored under colex
//  order of the points, and this is the same ordering implicitly used in
//  generating the SPARSE_UNIQUE_INDEX array.
//
      grid_weight = new double[order_nd];

      webbur::product_mixed_growth_weight ( dim_num, order_1d, order_nd, rule, 
        np, p, gw_compute_weights, grid_weight );
//
//  Compute Smolyak's binomial coefficient for this grid.
//
      coeff = webbur::r8_mop ( level_max - level ) 
        * webbur::r8_choose ( dim_num - 1, level_max - level );
//
//  Add these weights to the rule.
//
      for ( order = 0; order < order_nd; order++ )
      {
        point_unique = sparse_unique_index[point_total];

        point_total = point_total + 1;

        sparse_weight[point_unique] = sparse_weight[point_unique] 
          + coeff * grid_weight[order];
      }

      delete [] grid_weight;

      if ( !more_grids )
      {
        break;
      }
    }
  }

  delete [] level_1d;
  delete [] order_1d;

  return;
}
//****************************************************************************80

void sgmg_write ( int dim_num, int rule[], int np[],
  double p[], int point_num, double sparse_weight[], double sparse_point[],
  std::string file_name )

//****************************************************************************80
//
//  Purpose:
//
//    SGMG_WRITE writes a sparse grid rule to five files.
//
//  Discussion:
//
//    The files are:
//    * the "N" file stores the NP values, as a DIM_NUM x 1 list.
//    * the "P" file stores the P values, as a sum(NP[*]) x 1 list.
//    * the "X" file stores the abscissas as a DIM_NUM x POINT_NUM list;
//    * the "W" file stores the weights as a POINT_NUM list;
//    * the "R" file stores the region, as a DIM_NUM x 2 list.
//
//    The entries in the "R" file are the two corners of the DIM_NUM dimensional
//    rectangle that constitutes the integration region.  Coordinates that
//    should be infinite are set to 1.0E+30.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 June 2010
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
//    Input, int RULE[DIM_NUM], the rule in each dimension.
//     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
//     2, "F2",  Fejer Type 2, Open Fully Nested.
//     3, "GP",  Gauss Patterson, Open Fully Nested.
//     4, "GL",  Gauss Legendre, Open Weakly Nested.
//     5, "GH",  Gauss Hermite, Open Weakly Nested.
//     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
//     7, "LG",  Gauss Laguerre, Open Non Nested.
//     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
//     9, "GJ",  Gauss Jacobi, Open Non Nested.
//    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
//    11, "UO",  User supplied Open, presumably Non Nested.
//    12, "UC",  User supplied Closed, presumably Non Nested.
//
//    Input, int NP[DIM_NUM], the number of parameters used by each rule.
//
//    Input, double P[sum(NP[*])], the parameters needed by each rule.
//
//    Input, int POINT_NUM, the number of unique points 
//    in the grid. 
//
//    Input, double SPARSE_WEIGHT[POINT_NUM], the weights.
//
//    Input, double SPARSE_POINT[DIM_NUM*POINT_NUM], the points.
//
//    Input, string FILE_NAME, the main part of the file name.
//
{
  int dim;
  std::string file_name_n;
  std::string file_name_p;
  std::string file_name_r;
  std::string file_name_w;
  std::string file_name_x;
  int np_sum;
  int point;
  double *sparse_region;
  double t1;
  double t2;

  sparse_region = new double[dim_num*2];

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( rule[dim] == 1 )
    {
      sparse_region[dim+0*dim_num] = -1.0;
      sparse_region[dim+1*dim_num] = +1.0;
    }
    else if ( rule[dim] == 2 )
    {
      sparse_region[dim+0*dim_num] = -1.0;
      sparse_region[dim+1*dim_num] = +1.0;
    }
    else if ( rule[dim] == 3 )
    {
      sparse_region[dim+0*dim_num] = -1.0;
      sparse_region[dim+1*dim_num] = +1.0;
    }
    else if ( rule[dim] == 4 )
    {
      sparse_region[dim+0*dim_num] = -1.0;
      sparse_region[dim+1*dim_num] = +1.0;
    }
    else if ( rule[dim] == 5 )
    {
      sparse_region[dim+0*dim_num] = - webbur::r8_huge ( );
      sparse_region[dim+1*dim_num] = + webbur::r8_huge ( );
    }
    else if ( rule[dim] == 6 )
    {
      sparse_region[dim+0*dim_num] = - webbur::r8_huge ( );
      sparse_region[dim+1*dim_num] = + webbur::r8_huge ( );
    }
    else if ( rule[dim] == 7 )
    {
      sparse_region[dim+0*dim_num] = 0.0;
      sparse_region[dim+1*dim_num] = webbur::r8_huge ( );
    }
    else if ( rule[dim] == 8 )
    {
      sparse_region[dim+0*dim_num] = 0.0;
      sparse_region[dim+1*dim_num] = webbur::r8_huge ( );
    }
    else if ( rule[dim] == 9 )
    {
      sparse_region[dim+0*dim_num] = -1.0;
      sparse_region[dim+1*dim_num] = +1.0;
    }
    else if ( rule[dim] == 10 )
    {
      sparse_region[dim+0*dim_num] = - webbur::r8_huge ( );
      sparse_region[dim+1*dim_num] = + webbur::r8_huge ( );
    }
//
//  Best guess as to region extent for rules of type 11 or 12.
//
    else if ( rule[dim] == 11 )
    {
      t1 =   webbur::r8_huge ( );
      t2 = - webbur::r8_huge ( );
      for ( point = 0; point < point_num; point++ )
      {
        t1 = webbur::r8_min ( t1, sparse_point[dim+point*dim_num] );
        t2 = webbur::r8_max ( t2, sparse_point[dim+point*dim_num] );
      }
      sparse_region[dim+0*dim_num] = t1;
      sparse_region[dim+1*dim_num] = t2;
    }
    else if ( rule[dim] == 12 )
    {
      t1 =   webbur::r8_huge ( );
      t2 = - webbur::r8_huge ( );
      for ( point = 0; point < point_num; point++ )
      {
        t1 = webbur::r8_min ( t1, sparse_point[dim+point*dim_num] );
        t2 = webbur::r8_max ( t2, sparse_point[dim+point*dim_num] );
      }
      sparse_region[dim+0*dim_num] = t1;
      sparse_region[dim+1*dim_num] = t2;
    }
    else 
    {
      std::cerr << "\n";
      std::cerr << "SGMG_WRITE - Fatal error!\n";
      std::cerr << "  Unexpected value of RULE[" << dim << "] = " 
           << rule[dim] << ".\n";
      std::exit ( 1 );
    }
  }
  std::cout << "\n";
  std::cout << "SGMG_WRITE:\n";

  file_name_n = file_name + "_n.txt";
  webbur::i4mat_write ( file_name_n, dim_num, 1, np );
  std::cout << "  Wrote the N file = \"" << file_name_n << "\".\n";

  np_sum = i4vec_sum ( dim_num, np );
  file_name_p = file_name + "_p.txt";
  webbur::r8mat_write ( file_name_p, 1, np_sum, p );
  std::cout << "  Wrote the P file = \"" << file_name_p << "\".\n";

  file_name_r = file_name + "_r.txt";
  webbur::r8mat_write ( file_name_r, dim_num, 2, sparse_region );
  std::cout << "  Wrote the R file = \"" << file_name_r << "\".\n";

  file_name_w = file_name + "_w.txt";
  webbur::r8mat_write ( file_name_w, 1, point_num, sparse_weight );
  std::cout << "  Wrote the W file = \"" << file_name_w << "\".\n";

  file_name_x = file_name + "_x.txt";
  webbur::r8mat_write ( file_name_x, dim_num, point_num, sparse_point );
  std::cout << "  Wrote the X file = \"" << file_name_x << "\".\n";

  delete [] sparse_region;

  return;
}
//
//  This final curly bracket closes the namespace.
//
}
