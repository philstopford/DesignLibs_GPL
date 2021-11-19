using System;
using System.Linq;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace Burkardt.Sparse;

public static class SGMGAniso
{
    public static double[] sgmga_aniso_balance(double alpha_max, int dim_num,
            double[] level_weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_ANISO_BALANCE "balances" an anisotropic weight vector.
        //
        //  Discussion:
        //
        //    The entries in LEVEL_WEIGHT are essentially arbitrary nonnegative numbers.
        //
        //    The ratio between two entries indicates their relative importance.
        //    For instance,
        //
        //      LEVEL_WEIGHT(1) / LEVEL_WEIGHT(2) = 10
        //
        //    means that variable 2 is 10 times more important than variable 1.
        //    Here, being 10 times more important means that we will generate 10 levels
        //    of sparse grid in direction 2 as we generate 1 level in direction 1.
        //
        //    Under this interpretation, a ratio of 10 already indicates an extreme 
        //    imbalanace in the variables, since 10 sparse grid levels in the second
        //    variable corresponds roughly to approximating x^1 only, and 
        //    all of y^1 through y^10.  A ratio higher than this seems unreasonable.
        //
        //    Therefore, this function tries to take a somewhat arbitrary level weight
        //    vector, and produce a "balanced" level weight vector with the properties
        //    that the mininum entry is 1 (representing the item of most importance)
        //    and the maximum entry is ALPHA_MAX.  A reasonable value of ALPHA_MAX
        //    might be 10 or even 5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 February 2010
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, double ALPHA_MAX, the maximum legal value of 
        //    LEVEL_WEIGHT, after all entries have been divided by the minimum 
        //    nonzero entry.  1 <= ALPHA_MAX.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.  
        //    The values must be positive.  
        //
        //    Output, double SGMGA_ANISO_BALANCE[DIM_NUM], the balanced 
        //    anisotropic weights.  The smallest nonzero entry is 1.0 and 
        //    no entry is greater than ALPHA_MAX.
        //
    {
        int dim;

        switch (alpha_max)
        {
            case < 1.0:
                Console.WriteLine("");
                Console.WriteLine("SGMGA_ANISO_BALANCE - Fatal error!");
                Console.WriteLine("  ALPHA_MAX < 1.0");
                return null;
        }

        //
        //  Find the smallest nonzero entry.
        //
        double level_weight_min = typeMethods.r8_huge();
        int nonzero_num = 0;

        for (dim = 0; dim < dim_num; dim++)
        {
            switch (level_weight[dim])
            {
                case > 0.0:
                {
                    if (level_weight[dim] < level_weight_min)
                    {
                        level_weight_min = level_weight[dim];
                        nonzero_num += 1;
                    }

                    break;
                }
            }
        }

        switch (nonzero_num)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("SGMGA_ANISO_BALANCE - Fatal error!");
                Console.WriteLine("  Could not find a positive entry in LEVEL_WEIGHT.");
                return null;
        }

        //
        //  Rescale so the smallest nonzero entry is 1.
        //
        double[] level_weight2 = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            level_weight2[dim] = level_weight[dim] / level_weight_min;
        }

        //
        //  Set the maximum entry to no more than ALPHA_MAX.
        //
        for (dim = 0; dim < dim_num; dim++)
        {
            level_weight2[dim] = Math.Min(alpha_max, level_weight2[dim]);
        }

        return level_weight2;
    }

    public static void sgmga_aniso_normalize(int option, int dim_num, ref double[] level_weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_ANISO_NORMALIZE normalizes the SGMGA anisotropic weight vector.
        //
        //  Discussion:
        //
        //    It is convenient for the user to initialize the anisotropic weight
        //    vector with any set of positive values.  These values are to be used
        //    as coefficients of the 1D levels, to evaluate an expression which 
        //    determines which 1D levels will be included in a given rule.
        //
        //    This means that a relatively LARGE coefficient forces the corresponding 
        //    level to be relatively SMALL.  This is perhaps the opposite of what
        //    a user might expect.  If a user wishes to use an importance vector,
        //    so that a relatively large importance should correspond to more levels,
        //    and hence more points, in that dimension, then the function
        //    SGMGA_IMPORTANCE_TO_ANISO should be called first!
        //
        //    Since the weights only represent the relative importance of the
        //    components, they may be multiplied by any (positive) scale factor.
        //    Nonetheless, it may be convenient to choose a particular normalization
        //    for the weights.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 November 2009
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int OPTION, the normalization option.
        //    0, no scaling is applied.
        //    1, the weights are scaled so that the minimum nonzero entry is 1.
        //    2, the weights are scaled so that they sum to DIM_NUM.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input/output, double LEVEL_WEIGHT[DIM_NUM], the anisotropic
        //    weights.  The input values must be strictly positive.  
        //    On output, these have been normalized.
        //
    {
        int dim;
        switch (option)
        {
            //
            //  Option 0, no normalization.
            //
            case 0:
                break;
            //
            //  Option 1, the minimum nonzero entry is 1.
            //
            case 1:
            {
                double level_weight_min = typeMethods.r8_huge();
                int found = 0;
                for (dim = 0; dim < dim_num; dim++)
                {
                    switch (level_weight[dim])
                    {
                        case > 0.0:
                        {
                            if (level_weight[dim] < level_weight_min)
                            {
                                level_weight_min = level_weight[dim];
                                found += 1;
                            }

                            break;
                        }
                    }
                }

                switch (found)
                {
                    case 0:
                        Console.WriteLine("");
                        Console.WriteLine("SGMGA_ANISO_NORMALIZE - Fatal error!");
                        Console.WriteLine("  Could not find a positive entry in LEVEL_WEIGHT.");
                        return;
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    level_weight[dim] /= level_weight_min;
                }

                break;
            }
            //
            //  Option 2, rescale so sum of weights is DIM_NUM.
            //
            case 2:
            {
                double level_weight_sum = typeMethods.r8vec_sum(dim_num, level_weight);

                switch (level_weight_sum)
                {
                    case <= 0.0:
                        Console.WriteLine("");
                        Console.WriteLine("SGMGA_ANISO_NORMALIZE - Fatal error!");
                        Console.WriteLine("  Sum of level weights is not positive.");
                        return;
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    level_weight[dim] = dim_num * level_weight[dim]
                                        / level_weight_sum;
                }

                break;
            }
        }
    }

    public static void sgmga_importance_to_aniso(int dim_num, double[] importance,
            ref double[] level_weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_IMPORTANCE_TO_ANISO: importance vector to anisotropic weight vector.
        //
        //  Discussion:
        //
        //    To specify the anisotropy of a multidimensional problem, the user is
        //    allowed to specify an "importance vector".  This vector can contain
        //    any set of positive values.  These values represent the relative
        //    importance of each dimension.  These values, with a suitable normalization,
        //    will be used to evaluate a constraint of the following form:
        //
        //      QMIN < Level(1) / Importance(1) + Level(2) / Importance(2) + ...
        //             Level(N) / Importance(N) <= QMAX
        //
        //    and a set of levels that satisfies this constraint will then be included
        //    in a given anistotropic sparse grid rule.  Thus, increasing the
        //    importance value of a particular dimension allows larger level values
        //    in that dimension to satisfy the constraint.
        //
        //    The program actually works with coefficients LEVEL_WEIGHT that are
        //    the inverse of the importance vector entries, with a suitable
        //    normalization.  This function is supplied to convert between the
        //    more natural "importance vector" and the internally useful 
        //    "level_weight" vector.
        //
        //    This function converts the importance vector to an unnormalized 
        //    anisotropy weight vector.
        //
        //    Note that some (but not all) of the IMPORTANCE vector entries may be zero.
        //    This indicates that the corresponding dimension is of "zero" or
        //    rather "minimal" importance.  In such a case, only a one-point quadrature
        //    rule will be applied for that dimension, no matter what sparse grid
        //    level is requested for the overall problem.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 November 2009
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double IMPORTANCE[DIM_NUM], the importance vector.
        //    All entries must be nonnegative, and at least one must be positive.
        //
        //    Output, double LEVEL_WEIGHT[DIM_NUM], the anisotropic
        //    weights.
        //
    {
        int dim;

        for (dim = 0; dim < dim_num; dim++)
        {
            switch (importance[dim])
            {
                case < 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("SGMGA_IMPORTANCE_TO_ANISO - Fatal error!");
                    Console.WriteLine("  Some IMPORTANCE entries are not positive.");
                    return;
            }
        }

        int found = 0;

        for (dim = 0; dim < dim_num; dim++)
        {
            switch (importance[dim])
            {
                case > 0.0:
                    level_weight[dim] = 1.0 / importance[dim];
                    found += 1;
                    break;
                default:
                    level_weight[dim] = 0.0;
                    break;
            }
        }

        switch (found)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("SGMGA_IMPORTANCE_TO_ANISO - Fatal error!");
                Console.WriteLine("  No importance entry is positive.");
                return;
        }
    }

    public static void sgmga_index(int dim_num, double[] level_weight, int level_max,
            int[] rule, int point_num, int point_total_num, int[] sparse_unique_index,
            int[] growth, ref int[] sparse_order, ref int[] sparse_index)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_INDEX indexes an SGMGA grid.
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
        //    26 April 2011
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
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
        //    Input, int GROWTH[DIM_NUM], the desired growth in each dimension.
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
        int point;
        SGMGAData data = new();
        switch (level_max)
        {
            //
            //  Special cases.
            //
            case < 0:
                return;
            case 0:
            {
                point = 0;
                for (dim = 0; dim < dim_num; dim++)
                {
                    sparse_order[dim + point * dim_num] = 1;
                    sparse_index[dim + point * dim_num] = 1;
                }

                return;
            }
        }

        //
        //  Initialize the INDEX and ORDER arrays to -1 to help catch errors.
        //
        for (point = 0; point < point_num; point++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                sparse_order[dim + point * dim_num] = -1;
                sparse_index[dim + point * dim_num] = -1;
            }
        }

        int point_count = 0;

        int[] level_1d = new int[dim_num];
        int[] level_1d_max = new int[dim_num];
        int[] order_1d = new int[dim_num];
        int[] point_index = new int[dim_num];
        //
        //  Initialization for SGMGA_VCN_ORDERED.
        //
        double level_weight_min_pos = typeMethods.r8vec_min_pos(dim_num, level_weight);
        double q_min = level_max * level_weight_min_pos
                       - typeMethods.r8vec_sum(dim_num, level_weight);
        double q_max = level_max * level_weight_min_pos;
        for (dim = 0; dim < dim_num; dim++)
        {
            switch (level_weight[dim])
            {
                case > 0.0:
                {
                    level_1d_max[dim] = (int) Math.Floor(q_max / level_weight[dim]) + 1;
                    if (q_max <= (level_1d_max[dim] - 1) * level_weight[dim])
                    {
                        level_1d_max[dim] -= 1;
                    }

                    break;
                }
                default:
                    level_1d_max[dim] = 0;
                    break;
            }
        }

        bool more_grids = false;
        //
        //  Seek all vectors LEVEL_1D which satisfy the constraint.
        //
        //    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
        //      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
        //      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
        //
        for (;;)
        {
            sgmga_vcn_ordered(ref data, dim_num, level_weight, level_1d_max,
                ref level_1d, q_min, q_max, ref more_grids);

            if (!more_grids)
            {
                break;
            }

            //
            //  Compute the combinatorial coefficient.
            //
            double coef = sgmga_vcn_coef(dim_num, level_weight, level_1d, q_max);

            switch (coef)
            {
                case 0.0:
                    continue;
            }

            //
            //  Transform each 1D level to a corresponding 1D order.
            //
            LevelToOrder.level_growth_to_order(dim_num, level_1d, rule, growth, ref order_1d);
            //
            //  The inner loop generates a POINT of the GRID of the LEVEL.
            //
            bool more_points = false;

            for (;;)
            {
                typeMethods.vec_colex_next3(dim_num, order_1d, ref point_index, ref more_points);

                if (!more_points)
                {
                    break;
                }

                int point_unique = sparse_unique_index[point_count];
                for (dim = 0; dim < dim_num; dim++)
                {
                    sparse_order[dim + point_unique * dim_num] = order_1d[dim];
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    sparse_index[dim + point_unique * dim_num] = point_index[dim];
                }

                point_count += 1;
            }
        }
    }

    public static void sgmga_point(int dim_num, double[] level_weight, int level_max,
            int[] rule, int[] np, double[] p,
            Func<int, int, double[], double[], double[]>[] gw_compute_points,
            int point_num, int[] sparse_order, int[] sparse_index, int[] growth,
            ref double[] sparse_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_POINT computes the points of an SGMGA rule.
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
        //    26 April 2011
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
        //
        //    Input, int LEVEL_MAX, controls the size of the final sparse grid.
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
        //    as determined by SGMGA_SIZE.
        //
        //    Input, int SPARSE_ORDER[DIM_NUM*POINT_NUM], lists, for each point,
        //    the order of the 1D rules used in the grid that generated it.
        //
        //    Input, int SPARSE_INDEX[DIM_NUM*POINT_NUM], lists, for each point,
        //    its index in each of the 1D rules in the grid that generated it.
        //    The indices are 1-based.
        //
        //    Input, int GROWTH[DIM_NUM], the desired growth in each dimension.
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
        int[] orders = null;
        int point;

        for (point = 0; point < point_num; point++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                sparse_point[dim + point * dim_num] = -typeMethods.r8_huge();
            }
        }

        //
        //  Compute the point coordinates.
        //
        int[] level_1d_max = new int[dim_num];
        int[] levels = new int[dim_num];
        double level_weight_min_pos = typeMethods.r8vec_min_pos(dim_num, level_weight);
        double q_max = level_max * level_weight_min_pos;

        int p_index = 0;

        for (dim = 0; dim < dim_num; dim++)
        {
            switch (level_weight[dim])
            {
                case > 0.0:
                {
                    level_1d_max[dim] = (int) Math.Floor(q_max / level_weight[dim]) + 1;
                    if (q_max <= (level_1d_max[dim] - 1) * level_weight[dim])
                    {
                        level_1d_max[dim] -= 1;
                    }

                    break;
                }
                default:
                    level_1d_max[dim] = 0;
                    break;
            }

            int level;
            for (level = 0; level <= level_1d_max[dim]; level++)
            {
                LevelToOrder.level_growth_to_order(1, levels, rule.Skip(+dim).ToArray(),
                    growth.Skip(+dim).ToArray(), ref orders);

                int order = orders.Length;

                double[] points = new double[order];

                switch (rule[dim])
                {
                    case 1:
                        ClenshawCurtis.clenshaw_curtis_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 2:
                        Fejer2.fejer2_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 3:
                        PattersonQuadrature.patterson_lookup_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 4:
                        Legendre.QuadratureRule.legendre_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 5:
                        HermiteQuadrature.hermite_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 6:
                        HermiteQuadrature.gen_hermite_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 7:
                        Laguerre.QuadratureRule.laguerre_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 8:
                        Laguerre.QuadratureRule.gen_laguerre_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 9:
                        JacobiQuadrature.jacobi_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 10:
                        HermiteQuadrature.hermite_genz_keister_lookup_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 11:
                    case 12:
                        gw_compute_points[dim](
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    default:
                        Console.WriteLine("");
                        Console.WriteLine("SGMGA_POINT - Fatal error!");
                        Console.WriteLine("  Unexpected value of RULE[" + dim + "] = "
                                          + rule[dim] + ".");
                        return;
                }

                for (point = 0; point < point_num; point++)
                {
                    if (sparse_order[dim + point * dim_num] == order)
                    {
                        sparse_point[dim + point * dim_num] =
                            points[sparse_index[dim + point * dim_num] - 1];
                    }
                }
            }

            p_index += np[dim];
        }

        //
        //  Check to see if we missed any points.
        //
        for (point = 0; point < point_num; point++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                if (!(Math.Abs(sparse_point[dim + point * dim_num] - -typeMethods.r8_huge()) <= double.Epsilon))
                {
                    continue;
                }

                Console.WriteLine("");
                Console.WriteLine("SGMGA_POINT - Fatal error!");
                Console.WriteLine("  At least one point component was not assigned.");
                Console.WriteLine("  POINT = " + point + "");
                Console.WriteLine("  DIM = " + dim + "");
                Console.WriteLine("  SPARSE_ORDER(DIM,POINT) = "
                                  + sparse_order[dim + point * dim_num] + "");
                Console.WriteLine("  LEVEL_WEIGHT(DIM) = " + level_weight[dim] + "");
                return;
            }
        }
    }

    public static void sgmga_product_weight(int dim_num, int[] order_1d, int order_nd,
            int[] rule, int[] np, double[] p,
            Func<int, int, double[], double[], double[]>[] gw_compute_weights,
            ref double[] weight_nd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_PRODUCT_WEIGHT computes the weights of a mixed product rule.
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
        //    09 June 2010
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
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
        typeMethods.r8vecDPData data = new();

        for (i = 0; i < order_nd; i++)
        {
            weight_nd[i] = 1.0;
        }

        int p_index = 0;

        for (dim = 0; dim < dim_num; dim++)
        {
            double[] weight_1d = new double[order_1d[dim]];

            switch (rule[dim])
            {
                case 1:
                    ClenshawCurtis.clenshaw_curtis_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(p_index).ToArray(), weight_1d);
                    break;
                case 2:
                    Fejer2.fejer2_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(p_index).ToArray(), weight_1d);
                    break;
                case 3:
                    PattersonQuadrature.patterson_lookup_weights_np(
                        order_1d[dim], np[dim], p.Skip(p_index).ToArray(), weight_1d);
                    break;
                case 4:
                    Legendre.QuadratureRule.legendre_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(p_index).ToArray(), weight_1d);
                    break;
                case 5:
                    HermiteQuadrature.hermite_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(p_index).ToArray(), weight_1d);
                    break;
                case 6:
                    HermiteQuadrature.gen_hermite_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(p_index).ToArray(), weight_1d);
                    break;
                case 7:
                    Laguerre.QuadratureRule.laguerre_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(p_index).ToArray(), weight_1d);
                    break;
                case 8:
                    Laguerre.QuadratureRule.gen_laguerre_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(p_index).ToArray(), weight_1d);
                    break;
                case 9:
                    JacobiQuadrature.jacobi_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(p_index).ToArray(), weight_1d);
                    break;
                case 10:
                    HermiteQuadrature.hermite_genz_keister_lookup_weights_np(
                        order_1d[dim], np[dim], p.Skip(p_index).ToArray(), weight_1d);
                    break;
                case 11:
                case 12:
                    gw_compute_weights[dim](
                        order_1d[dim], np[dim], p.Skip(p_index).ToArray(), weight_1d);
                    break;
                default:
                    Console.WriteLine("");
                    Console.WriteLine("SGMGA_PRODUCT_WEIGHT - Fatal error!");
                    Console.WriteLine("  Unexpected value of RULE[" + dim + "] = "
                                      + rule[dim] + ".");
                    return;
            }

            p_index += np[dim];

            typeMethods.r8vec_direct_product2(ref data, dim, order_1d[dim], weight_1d,
                dim_num, order_nd, ref weight_nd);

        }
    }

    public static int sgmga_size(int dim_num, double[] level_weight, int level_max, int[] rule,
            int[] np, double[] p,
            Func<int, int, double[], double[], double[]>[] gw_compute_points,
            double tol, int[] growth)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_SIZE sizes an SGMGA grid, discounting duplicates.
        //
        //  Discussion:
        //
        //    The sparse grid is the logical sum of product grids that satisfy
        //    a particular constraint.
        //
        //    Depending on the 1D rules involved, there may be many duplicate points
        //    in the sparse grid.
        //
        //    This function counts the unique points in the sparse grid.  It does this
        //    in a straightforward way, by actually generating all the points, and
        //    comparing them, with a tolerance for equality.
        //
        //    This function has been modified to automatically omit points for which
        //    the "combinatorial coefficient" is zero, since such points would have
        //    a weight of zero in the grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 April 2011
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
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
        //    Input, int GROWTH[DIM_NUM], the desired growth in each dimension.
        //    0, "DF", default growth associated with this quadrature rule;
        //    1, "SL", slow linear, L+1;
        //    2  "SO", slow linear odd, O=1+2((L+1)/2)
        //    3, "ML", moderate linear, 2L+1;
        //    4, "SE", slow exponential;
        //    5, "ME", moderate exponential;
        //    6, "FE", full exponential.
        //
        //    Output, int SGMGA_SIZE, the number of unique points.
        //
    {
        int dim;
        int[] orders = null;
        int point;
        int point_num;
        SGMGAData data = new();
        switch (level_max)
        {
            //
            //  Special cases.
            //
            case < 0:
                point_num = -1;
                return point_num;
            case 0:
                point_num = 1;
                return point_num;
        }

        //
        //  Get total number of points, including duplicates.
        //
        int point_total_num = sgmga_size_total(dim_num, level_weight, level_max,
            rule, growth);
        //
        //  Generate SPARSE_TOTAL_ORDER and SPARSE_TOTAL_INDEX arrays 
        //  for the TOTAL set of points.
        //
        int[] sparse_total_order = new int[dim_num * point_total_num];
        int[] sparse_total_index = new int[dim_num * point_total_num];

        int point_total_num2 = 0;

        int[] level_1d = new int[dim_num];
        int[] levels = new int[dim_num];
        int[] level_1d_max = new int[dim_num];
        int[] order_1d = new int[dim_num];
        int[] point_index = new int[dim_num];
        //
        //  Initialization for SGMGA_VCN_ORDERED.
        //
        double level_weight_min_pos = typeMethods.r8vec_min_pos(dim_num, level_weight);
        double q_min = level_max * level_weight_min_pos
                       - typeMethods.r8vec_sum(dim_num, level_weight);
        double q_max = level_max * level_weight_min_pos;
        for (dim = 0; dim < dim_num; dim++)
        {
            switch (level_weight[dim])
            {
                case > 0.0:
                {
                    level_1d_max[dim] = (int) Math.Floor(q_max / level_weight[dim]) + 1;
                    if (q_max <= (level_1d_max[dim] - 1) * level_weight[dim])
                    {
                        level_1d_max[dim] -= 1;
                    }

                    break;
                }
                default:
                    level_1d_max[dim] = 0;
                    break;
            }
        }

        bool more_grids = false;
        //
        //  Seek all vectors LEVEL_1D which satisfy the constraint.
        //
        //    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
        //      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
        //      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
        //
        for (;;)
        {
            sgmga_vcn_ordered(ref data, dim_num, level_weight, level_1d_max,
                ref level_1d, q_min, q_max, ref more_grids);

            if (!more_grids)
            {
                break;
            }

            //
            //  Compute the combinatorial coefficient.
            //
            double coef = sgmga_vcn_coef(dim_num, level_weight, level_1d, q_max);

            switch (coef)
            {
                case 0.0:
                    continue;
            }

            //
            //  Transform each 1D level to a corresponding 1D order.
            //
            LevelToOrder.level_growth_to_order(dim_num, level_1d, rule, growth, ref order_1d);
            //
            //  The inner loop generates a POINT of the GRID of the LEVEL.
            //
            bool more_points = false;

            for (;;)
            {
                typeMethods.vec_colex_next3(dim_num, order_1d, ref point_index, ref more_points);

                if (!more_points)
                {
                    break;
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    sparse_total_order[dim + point_total_num2 * dim_num] = order_1d[dim];
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    sparse_total_index[dim + point_total_num2 * dim_num] = point_index[dim];
                }

                point_total_num2 += 1;
            }
        }

        //
        //  Now compute the coordinates of the TOTAL set of points.
        //
        double[] sparse_total_point = new double[dim_num * point_total_num];

        for (point = 0; point < point_total_num; point++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                sparse_total_point[dim + point * dim_num] = typeMethods.r8_huge();
            }
        }

        //
        //  Compute the point coordinates.
        //
        level_1d_max = new int[dim_num];
        level_weight_min_pos = typeMethods.r8vec_min_pos(dim_num, level_weight);
        q_max = level_max * level_weight_min_pos;

        int p_index = 0;
        for (dim = 0; dim < dim_num; dim++)
        {
            switch (level_weight[dim])
            {
                case > 0.0:
                {
                    level_1d_max[dim] = (int) Math.Floor(q_max / level_weight[dim]) + 1;
                    if (q_max <= (level_1d_max[dim] - 1) * level_weight[dim])
                    {
                        level_1d_max[dim] -= 1;
                    }

                    break;
                }
                default:
                    level_1d_max[dim] = 0;
                    break;
            }

            int level;
            for (level = 0; level <= level_1d_max[dim]; level++)
            {
                LevelToOrder.level_growth_to_order(1, levels, rule.Skip(+dim).ToArray(),
                    growth.Skip(+dim).ToArray(), ref orders);

                int order = orders.Length;

                double[] points = new double[order];

                switch (rule[dim])
                {
                    case 1:
                        ClenshawCurtis.clenshaw_curtis_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 2:
                        Fejer2.fejer2_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 3:
                        PattersonQuadrature.patterson_lookup_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 4:
                        Legendre.QuadratureRule.legendre_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 5:
                        HermiteQuadrature.hermite_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 6:
                        HermiteQuadrature.gen_hermite_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 7:
                        Laguerre.QuadratureRule.laguerre_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 8:
                        Laguerre.QuadratureRule.gen_laguerre_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 9:
                        JacobiQuadrature.jacobi_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 10:
                        HermiteQuadrature.hermite_genz_keister_lookup_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 11:
                    case 12:
                        gw_compute_points[dim](
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    default:
                        Console.WriteLine("");
                        Console.WriteLine("SGMGA_SIZE - Fatal error!");
                        Console.WriteLine("  Unexpected value of RULE[" + dim + "] = "
                                          + rule[dim] + ".");
                        return 1;
                }

                for (point = 0; point < point_total_num; point++)
                {
                    if (sparse_total_order[dim + point * dim_num] == order)
                    {
                        sparse_total_point[dim + point * dim_num] =
                            points[sparse_total_index[dim + point * dim_num] - 1];
                    }
                }
            }

            p_index += np[dim];
        }

        //
        //  Count the tolerably unique columns. 
        //
        int seed = 123456789;

        point_num = typeMethods.point_radial_tol_unique_count(dim_num, point_total_num,
            sparse_total_point, tol, ref seed);

        return point_num;
    }

    public static int sgmga_size_total(int dim_num, double[] level_weight, int level_max,
            int[] rule, int[] growth)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_SIZE_TOTAL sizes an SGMGA grid, counting duplicates.
        //
        //  Discussion:
        //
        //    This routine returns the total point count for an SGMGA
        //    ( Sparse Grid of Mixed type with Growth rule and Anisotropic weights).
        //
        //    The sparse grid is the logical sum of product grids.
        //
        //    The sparse grid has an associated integer index LEVEL_MAX, whose lowest 
        //    value is 0.  LEVEL_MAX = 0 indicates the sparse grid made up of one product 
        //    grid, which in turn is the product of 1D factor grids of the lowest level.
        //    This usually means the sparse grid with LEVEL_MAX equal to 0 is a
        //    one point grid.
        //
        //    We can assign a level to each factor grid, and hence a LEVEL vector
        //    to the corresponding product grid, and a weighted index
        //    LEVEL_GRID (which will in general be a real number):
        //
        //      LEVEL_GRID = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * LEVEL(I)
        //
        //    The product grid will participate in the formation of the sparse grid
        //    if it satisfies the following weighted constraint:
        //
        //      LEVEL_MAX - DIM_NUM < LEVEL_GRID <= LEVEL_MAX
        //
        //    This routine determines the total number of abscissas in all the 
        //    product rules used to form the SGMGA associated with the index LEVEL_MAX.
        //    The count disregards duplication.  If the same multidimensional abcsissa
        //    occurs in two different product rules that are part of the SGMGA, then
        //    that single abcissa is counted twice. 
        //
        //    This computation is useful in cases where the entire set of abscissas
        //    is going to be generated, preparatory to compression to finding, indexing
        //    and merging the duplicate abcissass.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 April 2011
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
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
        //    Input, int GROWTH[DIM_NUM], the desired growth in each dimension.
        //    0, "DF", default growth associated with this quadrature rule;
        //    1, "SL", slow linear, L+1;
        //    2  "SO", slow linear odd, O=1+2((L+1)/2)
        //    3, "ML", moderate linear, 2L+1;
        //    4, "SE", slow exponential;
        //    5, "ME", moderate exponential;
        //    6, "FE", full exponential.
        //
        //    Output, int SGMGA_SIZE_TOTAL, the number of points
        //    including repetitions.
        //
    {
        int dim;
        int point_total_num;
        SGMGAData data = new();
        switch (level_max)
        {
            //
            //  Special case.
            //
            case 0:
                point_total_num = 1;
                return point_total_num;
        }

        point_total_num = 0;

        int[] level_1d = new int[dim_num];
        int[] level_1d_max = new int[dim_num];
        int[] order_1d = new int[dim_num];
        //
        //  Initialization for SGMGA_VCN_ORDERED.
        //
        double level_weight_min_pos = typeMethods.r8vec_min_pos(dim_num, level_weight);
        double q_min = level_max * level_weight_min_pos
                       - typeMethods.r8vec_sum(dim_num, level_weight);
        double q_max = level_max * level_weight_min_pos;
        for (dim = 0; dim < dim_num; dim++)
        {
            switch (level_weight[dim])
            {
                case > 0.0:
                {
                    level_1d_max[dim] = (int) Math.Floor(q_max / level_weight[dim]) + 1;
                    if (q_max <= (level_1d_max[dim] - 1) * level_weight[dim])
                    {
                        level_1d_max[dim] -= 1;
                    }

                    break;
                }
                default:
                    level_1d_max[dim] = 0;
                    break;
            }
        }

        bool more_grids = false;
        //
        //  Seek all vectors LEVEL_1D which satisfy the constraint.
        //
        //    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
        //      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
        //      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
        //
        for (;;)
        {
            sgmga_vcn_ordered(ref data, dim_num, level_weight, level_1d_max,
                ref level_1d, q_min, q_max, ref more_grids);

            if (!more_grids)
            {
                break;
            }

            //
            //  Compute the combinatorlal coefficient.
            //
            double coef = sgmga_vcn_coef(dim_num, level_weight, level_1d, q_max);

            switch (coef)
            {
                case 0.0:
                    continue;
                default:
                    //
                    //  Transform each 1D level to a corresponding 1D order.
                    //
                    LevelToOrder.level_growth_to_order(dim_num, level_1d, rule, growth, ref order_1d);

                    point_total_num += typeMethods.i4vec_product(dim_num,
                        order_1d);
                    break;
            }
        }

        return point_total_num;
    }

    public static void sgmga_unique_index(int dim_num, double[] level_weight, int level_max,
            int[] rule, int[] np, double[] p,
            Func<int, int, double[], double[], double[]>[] gw_compute_points,
            double tol, int point_num, int point_total_num, int[] growth,
            ref int[] sparse_unique_index)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_UNIQUE_INDEX maps nonunique to unique points.
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
        //    26 April 2011
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the maximum value of LEVEL.
        //
        //    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
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
        //    Input, int POINT_NUM, the number of unique points in the grid. 
        //
        //    Input, int POINT_TOTAL_NUM, the total number of points 
        //    in the grid. 
        //
        //    Input, int GROWTH[DIM_NUM], the desired growth in each dimension.
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
        int[] orders = null;
        int point;
        SGMGAData data = new();
        switch (level_max)
        {
            //
            //  Special cases.
            //
            case < 0:
                return;
            case 0:
                sparse_unique_index[0] = 0;
                return;
        }

        //
        //  Generate SPARSE_TOTAL_ORDER and SPARSE_TOTAL_INDEX arrays 
        //  for the TOTAL set of points.
        //
        int[] sparse_total_order = new int[dim_num * point_total_num];
        int[] sparse_total_index = new int[dim_num * point_total_num];

        int[] level_1d = new int[dim_num];
        int[] level_1d_max = new int[dim_num];
        int[] order_1d = new int[dim_num];
        int[] point_index = new int[dim_num];

        int point_total_num2 = 0;
        //
        //  Initialization for SGMGA_VCN_ORDERED.
        //
        double level_weight_min_pos = typeMethods.r8vec_min_pos(dim_num, level_weight);
        double q_min = level_max * level_weight_min_pos
                       - typeMethods.r8vec_sum(dim_num, level_weight);
        double q_max = level_max * level_weight_min_pos;
        for (dim = 0; dim < dim_num; dim++)
        {
            switch (level_weight[dim])
            {
                case > 0.0:
                {
                    level_1d_max[dim] = (int) Math.Floor(q_max / level_weight[dim]) + 1;
                    if (q_max <= (level_1d_max[dim] - 1) * level_weight[dim])
                    {
                        level_1d_max[dim] -= 1;
                    }

                    break;
                }
                default:
                    level_1d_max[dim] = 0;
                    break;
            }
        }

        bool more_grids = false;
        //
        //  Seek all vectors LEVEL_1D which satisfy the constraint:
        //
        //    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
        //      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
        //      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
        //
        for (;;)
        {
            sgmga_vcn_ordered(ref data, dim_num, level_weight, level_1d_max,
                ref level_1d, q_min, q_max, ref more_grids);

            if (!more_grids)
            {
                break;
            }

            //
            //  Compute the combinatorial coefficient.
            //
            double coef = sgmga_vcn_coef(dim_num, level_weight, level_1d, q_max);

            switch (coef)
            {
                case 0.0:
                    continue;
            }

            //
            //  Transform each 1D level to a corresponding 1D order.
            //
            LevelToOrder.level_growth_to_order(dim_num, level_1d, rule, growth, ref order_1d);
            //
            //  The inner loop generates a POINT of the GRID of the LEVEL.
            //
            bool more_points = false;

            for (;;)
            {
                typeMethods.vec_colex_next3(dim_num, order_1d, ref point_index, ref more_points);

                if (!more_points)
                {
                    break;
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    sparse_total_order[dim + point_total_num2 * dim_num] = order_1d[dim];
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    sparse_total_index[dim + point_total_num2 * dim_num] = point_index[dim];
                }

                point_total_num2 += 1;
            }
        }

        //
        //  Now compute the coordinates of the TOTAL set of points.
        //
        double[] sparse_total_point = new double[dim_num * point_total_num];

        for (point = 0; point < point_total_num; point++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                sparse_total_point[dim + point * dim_num] = typeMethods.r8_huge();
            }
        }

        //
        //  Compute the point coordinates.
        //
        level_1d_max = new int[dim_num];
        int[] levels = new int[dim_num];
        level_weight_min_pos = typeMethods.r8vec_min_pos(dim_num, level_weight);
        q_max = level_max * level_weight_min_pos;

        int p_index = 0;
        for (dim = 0; dim < dim_num; dim++)
        {
            switch (level_weight[dim])
            {
                case > 0.0:
                {
                    level_1d_max[dim] = (int) Math.Floor(q_max / level_weight[dim]) + 1;
                    if (q_max <= (level_1d_max[dim] - 1) * level_weight[dim])
                    {
                        level_1d_max[dim] -= 1;
                    }

                    break;
                }
                default:
                    level_1d_max[dim] = 0;
                    break;
            }

            int level;
            for (level = 0; level <= level_1d_max[dim]; level++)
            {
                LevelToOrder.level_growth_to_order(1, levels, rule.Skip(+dim).ToArray(),
                    growth.Skip(+dim).ToArray(), ref orders);

                int order = orders.Length;

                double[] points = new double[order];

                switch (rule[dim])
                {
                    case 1:
                        ClenshawCurtis.clenshaw_curtis_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 2:
                        Fejer2.fejer2_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 3:
                        PattersonQuadrature.patterson_lookup_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 4:
                        Legendre.QuadratureRule.legendre_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 5:
                        HermiteQuadrature.hermite_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 6:
                        HermiteQuadrature.gen_hermite_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 7:
                        Laguerre.QuadratureRule.laguerre_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 8:
                        Laguerre.QuadratureRule.gen_laguerre_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 9:
                        JacobiQuadrature.jacobi_compute_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 10:
                        HermiteQuadrature.hermite_genz_keister_lookup_points_np(
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    case 11:
                    case 12:
                        gw_compute_points[dim](
                            order, np[dim], p.Skip(p_index).ToArray(), points);
                        break;
                    default:
                        Console.WriteLine("");
                        Console.WriteLine("SGMGA_UNIQUE_INDEX - Fatal error!");
                        Console.WriteLine("  Unexpected value of RULE[" + dim + "] = "
                                          + rule[dim] + ".");
                        return;
                }

                for (point = 0; point < point_total_num; point++)
                {
                    if (sparse_total_order[dim + point * dim_num] == order)
                    {
                        sparse_total_point[dim + point * dim_num] =
                            points[sparse_total_index[dim + point * dim_num] - 1];
                    }
                }
            }

            p_index += np[dim];
        }

        //
        //  Merge points that are too close.
        //
        int seed = 123456789;

        int[] undx = new int[point_num];

        typeMethods.point_radial_tol_unique_index(dim_num, point_total_num,
            sparse_total_point, tol, ref seed, ref undx, ref sparse_unique_index);

        for (point = 0; point < point_total_num; point++)
        {
            int rep = undx[sparse_unique_index[point]];
            if (point == rep)
            {
                continue;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                sparse_total_point[dim + point * dim_num] = sparse_total_point[dim + rep * dim_num];
            }
        }

        //
        //  Construct an index that indicates the "rank" of the unique points.
        //
        typeMethods.point_unique_index(dim_num, point_total_num, sparse_total_point,
            point_num, ref undx, ref sparse_unique_index);

    }

    public class SGMGAData
    {
        public int dir;
        public int n2;
        public int nstart;
        public int[] xmax;
        public int xmin;
        public double q_max2;
        public double q_min2;
    }

    public static void sgmga_vcn(ref SGMGAData data, int n, double[] w, ref int[] x, double q_min, double q_max,
            ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_VCN returns the next constrained vector.
        //
        //  Discussion:
        //
        //    This function is intended to replace the "naive" version, now called
        //    SGMGA_VCN_NAIVE, which is too slow for high dimensional problems.
        //
        //    For nonnegative vectors X of dimension N, and nonnegative
        //    weights W, we define:
        //
        //      Q = sum ( 1 <= I <= N ) W(I) * X(I)
        //
        //    and seek X satisfying the constraint:
        //
        //      Q_MIN < Q <= Q_MAX
        //
        //    This routine returns, one at a time exactly those X which satisfy
        //    the constraint.  No attempt is made to return the X values in 
        //    any particular order as far as Q goes.  
        //
        //  Example:
        // 
        //        W               4.0 3.0 5.0       
        //      MIN     16.0       0   0   0
        //      ---     ----      -----------
        //        1     20.0       5   0   0
        //        2     19.0       4   1   0
        //        3     18.0       3   2   0
        //        4     17.0       2   3   0
        //        5     20.0       2   4   0
        //        6     19.0       1   5   0
        //        7     18.0       0   6   0
        //        8     17.0       3   0   1
        //        9     20.0       3   1   1
        //       10     19.0       2   2   1
        //       11     18.0       1   3   1
        //       12     17.0       0   4   1
        //       13     20.0       0   5   1
        //       14     18.0       2   0   2
        //       15     17.0       1   1   2
        //       16     20.0       1   2   2
        //       17     19.0       0   3   2
        //       18     19.0       1   0   3
        //       19     18.0       0   1   3
        //       20     20.0       0   0   4
        //      ---     ----      ----------
        //      MAX     20.0       6   7   5         
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 May 2010
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vector.
        //
        //    Input, double W[N], the weights, which should be nonnegative.
        //    At least one weight must be positive.
        //
        //    Input/output, int X[N].  On first call, with 
        //    MORE = FALSE, the input value of X is not important.  On subsequent calls, 
        //    the input value of X should be the output value from the previous call.
        //    On output, (with MORE = TRUE), the value of X will be the "next"
        //    vector in the reverse lexicographical list of vectors that satisfy
        //    the condition.  However, on output with MORE = FALSE, the vector
        //    X is meaningless, because there are no more vectors in the list.
        //
        //    Input, double Q_MIN, Q_MAX, the lower and upper limits on the sum.
        //
        //    Input/output, bool *MORE.  On input, if the user has set MORE
        //    FALSE, the user is requesting the initiation of a new sequence
        //    of values.  If MORE is TRUE, then the user is requesting "more"
        //    values in the current sequence.  On output, if MORE is TRUE,
        //    then another value was found and returned in X, but if MORE is
        //    FALSE, then there are no more values in the sequence, and X is
        //    NOT the next value.
        //
    {
        int i;
        switch (more)
        {
            //
            //  Initialization for first call.
            //
            //  Allocate XMAX to remember the currently maximum possible value for each X.
            //
            //  Locate NSTART, the index of the first nonzero weight.
            //  The algorithm is easier to program if the last index we look at
            //  has a nonzero weight, so that it can always make up the remainder.
            //
            case false:
            {
                data.xmax = new int[n];

                data.nstart = -1;

                for (i = 0; i < n; i++)
                {
                    if (!(0.0 < w[i]))
                    {
                        continue;
                    }

                    data.nstart = i;
                    break;
                }

                switch (data.nstart)
                {
                    //
                    //  Theoretically, we could even handle the case where all weights are zero.
                    //  That case is ruled out elsewhere in this software, so I will not try
                    //  to deal with it here for now.
                    //
                    case -1:
                        Console.WriteLine("");
                        Console.WriteLine(" SGMGA_VCN - Fatal error!");
                        Console.WriteLine("  No weight is positive.");
                        return;
                }

                //
                //  Initialize X to zero, even the indices we ignore.
                //
                for (i = 0; i < n; i++)
                {
                    x[i] = 0;
                }

                //
                //  N2 points to our current index of interest.
                //
                data.n2 = n;
                data.dir = -1;

                more = true;
                break;
            }
        }

        //
        //  Look for the next solution vector X.
        //
        for (;;)
        {
            //
            //  If no more, the search is terminated.
            //
            if (!more)
            {
                break;
            }
            //
            //  DIR = -1, decrement N2, and, if possible, set X[N2] to XMIN.
            //  DIR =  0, hold N2 at current value, and see if we can increment X[N2].
            //

            if (data.dir is -1 or 0)
            {
                switch (data.dir)
                {
                    case -1:
                        data.n2 -= 1;
                        break;
                }

                switch (w[data.n2])
                {
                    case 0.0:
                        data.xmin = 0;
                        data.xmax[data.n2] = 0;
                        break;
                    default:
                    {
                        double t;
                        if (data.nstart < data.n2)
                        {
                            data.xmin = 0;
                            t = q_max;
                            for (i = data.n2 + 1; i < n; i++)
                            {
                                t -= w[i] * x[i];
                            }

                            data.xmax[data.n2] = (int) Math.Floor(t / w[data.n2]);
                        }
                        else if (data.n2 == data.nstart && data.dir == -1)
                        {
                            t = q_min;
                            for (i = data.n2 + 1; i < n; i++)
                            {
                                t -= w[i] * x[i];
                            }

                            data.xmin = data.xmin switch
                            {
                                < 0 => 0,
                                _ => (int) Math.Ceiling(t / w[data.n2])
                            };

                            t = 0.0;
                            for (i = 0; i < data.n2; i++)
                            {
                                t += w[i] * x[i];
                            }

                            t += w[data.n2] * data.xmin;
                            for (i = data.n2 + 1; i < n; i++)
                            {
                                t += w[i] * x[i];
                            }

                            if (t <= q_min)
                            {
                                data.xmin += 1;
                            }

                            x[data.n2] = data.xmin;
                            t = q_max;
                            for (i = data.n2 + 1; i < n; i++)
                            {
                                t -= w[i] * x[i];
                            }

                            data.xmax[data.n2] = (int) Math.Floor(t / w[data.n2]);
                        }

                        break;
                    }
                }

                if (data.xmax[data.n2] < data.xmin)
                {
                    data.dir = +1;
                }
                else
                {
                    if (data.n2 == data.nstart)
                    {
                        if (data.dir == -1)
                        {
                            data.dir = 0;
                            break;
                        }

                        if (data.dir != 0)
                        {
                            continue;
                        }

                        x[data.n2] += 1;
                        if (x[data.n2] <= data.xmax[data.n2])
                        {
                            break;
                        }

                        data.dir = +1;
                    }
                    else
                    {
                        x[data.n2] = data.xmin;
                    }
                }
            }
            else
            {
                switch (data.dir)
                {
                    //
                    //  DIR = + 1:
                    //  Try moving backwards to find an index N2 whose X we can increment.
                    //
                    case +1:
                    {
                        for (;;)
                        {
                            if (data.n2 == n - 1)
                            {
                                data.dir = 0;
                                more = false;
                                break;
                            }

                            data.n2 += 1;

                            if (!(0.0 < w[data.n2]))
                            {
                                continue;
                            }

                            if (x[data.n2] >= data.xmax[data.n2])
                            {
                                continue;
                            }

                            x[data.n2] += 1;
                            data.dir = -1;
                            break;
                        }

                        break;
                    }
                }
            }
        }
    }

    public static double sgmga_vcn_coef(int dim_num, double[] level_weight, int[] x,
            double q_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_VCN_COEF returns the "next" constrained vector's coefficient.
        //
        //  Discussion:
        //
        //    The related code "SGMGA_VCN_COEF_NAIVE" represents a "naive" approach to
        //    this calculation.  This code carries out the same calculation, but tries
        //    to avoid the potential explosion in work that is exponential in the
        //    spatial dimension for the naive approach.
        //
        //    We are considering nonnegative integer vectors X of dimension DIM_NUM 
        //    for which the functional
        //
        //      Q(X) = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)
        //
        //   satisfies the "Q" constraint:
        //
        //      Q_MIN < Q(X) <= Q_MAX
        //
        //    where LEVEL_WEIGHT is a vector of (essentially) positive weights.
        //    Some, but not all of the entries of LEVEL_WEIGHT might be zero;
        //    in that case, the corresponding values of X never vary, and do not
        //    play a part in the following computation.
        //
        //    Supposing we have a suitable vector X, we now wish to count the 
        //    number of distinct vectors Y which also satisfy the Q constraint
        //    as well as the following "binary" constraint:
        //
        //      Y(I) = X(I) + B(I)
        //
        //    where every entry of B is 0 or 1.
        //
        //    Clearly, there are 2^DIM_NUM vectors Y which satisfy the binary
        //    constraint, and a naive calculation would simply generate each 
        //    possible Y, evaluate Q(Y), and if Y satisfies the Q constraint,
        //    add it to the count.
        //
        //    But if we are considering even moderately large values of DIM_NUM, 
        //    say 20 <= DIM_NUM, then the mere task of generating all possible 
        //    Y vectors is burdensome.  If there are in fact likely to be only 
        //    a few satisfactory Y vectors, (which depends on the values of 
        //    Q_MIN, Q_MAX, and LEVEL_WEIGHT, of course) then it may be possible to
        //    find and count them more rapidly.
        //
        //    This function attempts a more rapid computation by carrying out the
        //    search in a particular order, and realizing that, in certain cases,
        //    if a particular value Y* does not satisfy the Q constraint, then
        //    a consecutive sequence of Y's following Y* also cannot satisfy the
        //    constraint, and hence the search can jump over them.
        //
        //  Example:
        //
        //    DIM_NUM = 3
        //    LEVEL_WEIGHT    3.0  2.0  1.0
        //    Q_MAX    6.0
        //
        //    U = unsigned count 
        //    S =   signed count returned as COEF
        //                 
        //    #   U  S   X1 X2 X3
        //
        //    1   8  0    0  0  0
        //    2   7  1    0  0  1
        //    3   6  0    0  0  2
        //    4   5 -1    0  0  3
        //    5   3 -1    0  0  4
        //    6   2  0    0  0  5
        //    7   1  1    0  0  6
        //    8   6  0    0  1  0
        //    9   5 -1    0  1  1
        //   10   3 -1    0  1  2
        //   11   2  0    0  1  3
        //   12   1  1    0  1  4
        //   13   3 -1    0  2  0
        //   14   2  0    0  2  1
        //   15   1  1    0  2  2
        //   16   1  1    0  3  0
        //   17   5 -1    1  0  0
        //   18   3 -1    1  0  1
        //   19   2  0    1  0  2
        //   20   1  1    1  0  3
        //   21   2  0    1  1  0
        //   22   1  1    1  1  1
        //   23   1  1    2  0  0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 May 2010
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double LEVEL_WEIGHT[DIM_NUM], the weights.
        //
        //    Input, int X[DIM_NUM], satisfies the Q constraint.
        //
        //    Input, double Q_MAX, the Q constraint maximum.
        //
        //    Output, double SGMGA_VCN_COEF, the combinatorial coefficient.
        //
    {
        int i;

        int c = 0;
        int[] b = new int[dim_num];

        for (i = 0; i < dim_num; i++)
        {
            b[i] = 0;
        }

        for (;;)
        {
            //
            //  Generate the next binary perturbation.
            //
            i = -1;

            int j;
            while (i < dim_num - 1)
            {
                i += 1;
                //
                //  If LEVEL_WEIGHT(I) == 0, B(I) is fixed at 0.  Next I.
                //
                if (level_weight[i] == 0.0)
                {
                }
                //
                //  If B(I) is 1, set it to 0.  Next I.
                //
                else if (b[i] == 1)
                {
                    b[i] = 0;
                }
                //
                //  B(I) is 0.  Convert it to 1.
                //
                else
                {
                    b[i] = 1;

                    for (;;)
                    {
                        // 
                        //  Does X + B satisfy the Q_MAX constraint?
                        //
                        double q = 0.0;
                        for (j = 0; j < dim_num; j++)
                        {
                            q += level_weight[j] * (x[j] + b[j]);
                        }

                        if (q <= q_max)
                        {
                            break;
                        }

                        //
                        //  If Q(X+B) now exceeds QMAX, B is rejected.  But we can also skip
                        //  all perturbations which agree with B through the I-th position.
                        //  To skip to the next "interesting" candidate, we essentially carry
                        //  out binary addition between B and a vector B' which has a single 1
                        //  in the I-th position.
                        //
                        b[i] = 0;

                        while (i < dim_num - 1)
                        {
                            i += 1;

                            if (level_weight[i] == 0.0)
                            {
                            }
                            else if (b[i] == 1)
                            {
                                b[i] = 0;
                            }
                            else
                            {
                                b[i] = 1;
                                break;
                            }
                        }
                    }

                    break;
                }
            }

            int b_sum = 0;
            for (j = 0; j < dim_num; j++)
            {
                b_sum += b[j];
            }

            //
            //  X+B is another solution to be counted.
            //
            c = c + 1 - 2 * (b_sum % 2);
            //
            //  We're done if we've got back to 0.
            //
            if (b_sum == 0)
            {
                break;
            }
        }

        double coef = c;

        return coef;
    }

    public static double sgmga_vcn_coef_naive(int dim_num, double[] level_weight, int[] x_max,
            int[] x, double q_min, double q_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_VCN_COEF_NAIVE returns the "next" constrained vector's coefficient.
        //
        //  Discussion:
        //
        //    This function uses a naive approach to the computation, resulting in
        //    a set of 2^DIM_NUM tasks.  Hence it is not suitable for cases where
        //    DIM_NUM is moderately large.  The function SGMGA_VCN_COEF carries out
        //    a more complicated but more efficient algorithm for the same computation.
        //
        //    We are given a vector X of dimension DIM_NUM which satisfies:
        //
        //      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).
        //
        //    and the following constraint:
        //
        //      sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I) <= Q_MAX
        //
        //    This routine computes the appropriate coefficient for X in the
        //    anisotropic sparse grid scheme.
        //
        //    The coefficient is calculated as follows:
        //
        //      Let B be a binary vector of length DIM_NUM, and let ||B|| represent
        //      the sum of the entries of B.
        //
        //      Coef = sum ( all B such that X+B satisfies constraints ) (-1)^||B||
        //
        //    Since X+0 satisfies the constraint, there is always at least one 
        //    summand.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 May 2010
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the number of components in the vector.
        //
        //    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
        //
        //    Input, int X_MAX[DIM_NUM], the maximum values allowed in each component.
        //
        //    Input, int X[DIM_NUM], a point which satisifies the constraints.
        //
        //    Input, double Q_MIN, Q_MAX, the lower and upper limits on the sum.
        //
        //    Output, double SGMGA_VCN_COEF_NAIVE, the combinatorial coefficient.
        //
    {
        int i;

        int[] b = new int[dim_num];

        for (i = 0; i < dim_num; i++)
        {
            b[i] = 0;
        }

        double coef = 1.0;

        for (;;)
        {
            //
            //  Generate the next binary perturbation.
            //
            typeMethods.binary_vector_next(dim_num, ref b);
            int b_sum = typeMethods.i4vec_sum(dim_num, b);
            //
            //  We're done if we've got back to 0.
            //
            if (b_sum == 0)
            {
                break;
            }

            //
            //  Does it satisfy the XMAX constraint?
            //  (THIS CHECK IS SURPRISINGLY NECESSARY, PARTLY BECAUSE OF ZERO WEIGHT).
            //
            bool too_big = false;
            for (i = 0; i < dim_num; i++)
            {
                if (x_max[i] >= x[i] + b[i])
                {
                    continue;
                }

                too_big = true;
                break;
            }

            switch (too_big)
            {
                case true:
                    continue;
            }

            //
            //  Does it satisfy the Q_MAX constraint?
            //
            double q = 0.0;
            for (i = 0; i < dim_num; i++)
            {
                q += level_weight[i] * (x[i] + b[i]);
            }

            if (q <= q_max)
            {
                coef += typeMethods.r8_mop(b_sum);
            }
        }

        return coef;
    }

    public static void sgmga_vcn_naive(int dim_num, double[] level_weight, int[] x_max, int[] x,
            double q_min, double q_max, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_VCN_NAIVE returns the next constrained vector.
        //
        //  Discussion:
        //
        //    This function uses a naive algorithm, which quickly becomes unsuitable
        //    for higher dimensions.  The function SGMGA_VCN is an attempt at a more
        //    efficient calculation of the same quantities.
        //
        //    We consider vectors X of dimension DIM_NUM satisfying:
        //
        //      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).
        //
        //    and define
        //
        //      Q = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)
        //
        //    and seek X satisfying the constraint:
        //
        //      Q_MIN < Q <= Q_MAX
        //
        //    For sparse grid applications, we compute
        //
        //      LEVEL_WEIGHT_MIN_POS = minimum positive entry in LEVEL_WEIGHT
        //
        //    and assume there is an underlying LEVEL used to index the sets of 
        //    constrained vectors, and that 
        //
        //      Q_MAX = LEVEL * LEVEL_WEIGHT_MIN_POS
        //      Q_MIN = LEVEL - LEVEL_WEIGHT_MIN_POS * sum ( LEVEL_WEIGHT(:) )
        //      X_MAX(I) = LEVEL * LEVEL_WEIGHT_MIN_POS / LEVEL_WEIGHT(I)
        //
        //    This routine returns, one at a time exactly those X which satisfy
        //    the constraint.  No attempt is made to return the X values in 
        //    any particular order as far as Q goes.  
        //
        //  Example:
        //
        //    LEVEL_WEIGHT:          1.000000        1.000000
        //
        //    Q_MIN:        0.000000
        //    Q_MAX:        2.000000
        //    X_MAX:                         2         2
        //
        //         1        1.000000         1         0
        //         2        2.000000         2         0
        //         3        1.000000         0         1
        //         4        2.000000         1         1
        //         5        2.000000         0         2
        //
        //    LEVEL_WEIGHT:          1.000000        2.000000
        //
        //    Q_MIN:       -1.000000
        //    Q_MAX:        2.000000
        //    X_MAX:                         2         1
        //
        //         1        0.000000         0         0
        //         2        1.000000         1         0
        //         3        2.000000         2         0
        //         4        2.000000         0         1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 October 2009
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the number of components in the vector.
        //
        //    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
        //
        //    Input, int X_MAX[DIM_NUM], the maximum values allowed in each component.
        //
        //    Input/output, int X[DIM_NUM].  On first call (with MORE = FALSE),
        //    the input value of X is not important.  On subsequent calls, the
        //    input value of X should be the output value from the previous call.
        //    On output, (with MORE = TRUE), the value of X will be the "next"
        //    vector in the reverse lexicographical list of vectors that satisfy
        //    the condition.  However, on output with MORE = FALSE, the vector
        //    X is meaningless, because there are no more vectors in the list.
        //
        //    Input, double Q_MIN, Q_MAX, the lower and upper
        //    limits on the sum.
        //
        //    Input/output, bool *MORE.  On input, if the user has set MORE
        //    FALSE, the user is requesting the initiation of a new sequence
        //    of values.  If MORE is TRUE, then the user is requesting "more"
        //    values in the current sequence.  On output, if MORE is TRUE,
        //    then another value was found and returned in X, but if MORE is
        //    FALSE, then there are no more values in the sequence, and X is
        //    NOT the next value.
        //
    {
        int i;
        double q;

        switch (more)
        {
            case false:
            {
                more = true;
                for (i = 0; i < dim_num; i++)
                {
                    x[i] = 0;
                }

                q = 0.0;
                for (i = 0; i < dim_num; i++)
                {
                    q += level_weight[i] * x[i];
                }

                if (q_min < q && q <= q_max)
                {
                    return;
                }

                break;
            }
        }

        for (;;)
        {
            int j = 0;

            for (;;)
            {
                if (x[j] < x_max[j])
                {
                    break;
                }

                if (dim_num - 1 <= j)
                {
                    more = false;
                    return;
                }

                j += 1;
            }

            x[j] += 1;
            for (i = 0; i < j; i++)
            {
                x[i] = 0;
            }

            q = 0.0;
            for (i = 0; i < dim_num; i++)
            {
                q += level_weight[i] * x[i];
            }

            if (q_min < q && q <= q_max)
            {
                break;
            }
        }
    }


    public static void sgmga_vcn_ordered(ref SGMGAData data, int dim_num, double[] level_weight, int[] x_max,
            ref int[] x, double q_min, double q_max, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_VCN_ORDERED returns the "next" constrained vector, with ordering.
        //
        //  Discussion:
        //
        //    We consider vectors X of dimension DIM_NUM satisfying:
        //
        //      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).
        //
        //    and define
        //
        //      Q = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)
        //
        //    and seek X's satisfying the constraint:
        //
        //      Q_MIN < Q <= Q_MAX
        //
        //    For sparse grid applications, we compute
        //
        //      LEVEL_WEIGHT_MIN_POS = minimum positive entry in LEVEL_WEIGHT
        //
        //    and assume there is an underlying LEVEL used to index the sets of 
        //    constrained vectors, and that 
        //
        //      Q_MAX = LEVEL * LEVEL_WEIGHT_MIN_POS
        //      Q_MIN = LEVEL - LEVEL_WEIGHT_MIN_POS * sum ( LEVEL_WEIGHT(:) )
        //      X_MAX(I) = LEVEL * LEVEL_WEIGHT_MIN_POS / LEVEL_WEIGHT(I)
        //
        //    This function returns, one at a time exactly those X which satisfy
        //    the constraint.
        //
        //    A weak ordering is imposed on the solution vectors.  This function 
        //    subdivides the range Q_MIN through Q_MAX into subintervals of width 1, so 
        //    that the X vectors returned are roughly sorted (or at least binned) 
        //    by Q value.
        //
        //  Example:
        //
        //    If the weights are also integral, then the X vectors are in fact SORTED 
        //    by Q value:
        //
        //    LEVEL_WEIGHT:          1.000000        1.000000
        //    Q_MIN:        0.000000
        //    Q_MAX:        2.000000
        //    X_MAX:                         2         2
        //
        //         1        1.000000         1         0
        //         2        1.000000         0         1
        //         3        2.000000         2         0
        //         4        2.000000         1         1
        //         5        2.000000         0         2
        //
        //    When the weights are not integral, then the X values are only BINNED
        //    by Q value, that is, we first get all X's with Q values between Q_MIN
        //    and Q_MIN+1, then Q_MIN+1 to Q_MIN+2 and so on, as demonstrated here:
        //
        //    LEVEL_WEIGHT:             1.5               1
        //    Q_MIN:  0.5
        //    Q_MAX:  3
        //    X_MAX:                           2         3
        //
        //           1             1.5         1         0
        //           2               1         0         1
        //           3             2.5         1         1
        //           4               2         0         2
        //           5               3         2         0
        //           6               3         0         3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 May 2010
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the number of components in the vector.
        //
        //    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
        //
        //    Input, int X_MAX[DIM_NUM], the maximum values allowed in each component.
        //
        //    Input/output, int X[DIM_NUM].  On first call (with MORE = FALSE),
        //    the input value of X is not important.  On subsequent calls, the
        //    input value of X should be the output value from the previous call.
        //    On output, (with MORE = TRUE), the value of X will be the "next"
        //    vector in the reverse lexicographical list of vectors that satisfy
        //    the condition.  However, on output with MORE = FALSE, the vector
        //    X is meaningless, because there are no more vectors in the list.
        //
        //    Input, double Q_MIN, Q_MAX, the lower and upper
        //    limits on the sum.
        //
        //    Input/output, bool *MORE.  On input, if the user has set MORE
        //    FALSE, the user is requesting the initiation of a new sequence
        //    of values.  If MORE is TRUE, then the user is requesting "more"
        //    values in the current sequence.  On output, if MORE is TRUE,
        //    then another value was found and returned in X, but if MORE is
        //    FALSE, then there are no more values in the sequence, and X is
        //    NOT the next value.
        //
    {
        switch (more)
        {
            //
            //  On first call, initialize the subrange.
            //
            case false:
                data.q_min2 = q_min;
                data.q_max2 = Math.Min(q_min + 1.0, q_max);
                break;
        }

        //
        //  Call a lower level function to search the subrange.
        //
        for (;;)
        {
            sgmga_vcn(ref data, dim_num, level_weight, ref x, data.q_min2, data.q_max2, ref more);
            switch (more)
            {
                //
                //  If another solution was found, return it.
                //
                case true:
                    return;
            }

            //
            //  If the current subrange is exhausted, try to move to the next one.
            //
            if (data.q_max2 < q_max)
            {
                data.q_min2 = data.q_max2;
                data.q_max2 = Math.Min(data.q_max2 + 1.0, q_max);
            }
            //
            //  If there are no more subranges, we're done.
            //
            else
            {
                break;
            }
        }
    }

    public static void sgmga_vcn_ordered_naive(ref SGMGAData data, int dim_num, double[] level_weight, int[] x_max,
            int[] x, double q_min, double q_max, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_VCN_ORDERED_NAIVE returns the "next" constrained vector, with ordering.
        //
        //  Discussion:
        //
        //    We consider vectors X of dimension DIM_NUM satisfying:
        //
        //      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).
        //
        //    and define
        //
        //      Q = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)
        //
        //    and seek X's satisfying the constraint:
        //
        //      Q_MIN < Q <= Q_MAX
        //
        //    For sparse grid applications, we compute
        //
        //      LEVEL_WEIGHT_MIN_POS = minimum positive entry in LEVEL_WEIGHT
        //
        //    and assume there is an underlying LEVEL used to index the sets of 
        //    constrained vectors, and that 
        //
        //      Q_MAX = LEVEL * LEVEL_WEIGHT_MIN_POS
        //      Q_MIN = LEVEL - LEVEL_WEIGHT_MIN_POS * sum ( LEVEL_WEIGHT(:) )
        //      X_MAX(I) = LEVEL * LEVEL_WEIGHT_MIN_POS / LEVEL_WEIGHT(I)
        //
        //    This function returns, one at a time exactly those X which satisfy
        //    the constraint.
        //
        //    A weak ordering is imposed on the solution vectors.  This function 
        //    subdivides the range Q_MIN through Q_MAX into subintervals of width 1, so 
        //    that the X vectors returned are roughly sorted (or at least binned) 
        //    by Q value.
        //
        //  Example:
        //
        //    If the weights are also integral, then the X vectors are in fact SORTED 
        //    by Q value:
        //
        //    LEVEL_WEIGHT:          1.000000        1.000000
        //    Q_MIN:        0.000000
        //    Q_MAX:        2.000000
        //    X_MAX:                         2         2
        //
        //         1        1.000000         1         0
        //         2        1.000000         0         1
        //         3        2.000000         2         0
        //         4        2.000000         1         1
        //         5        2.000000         0         2
        //
        //    When the weights are not integral, then the X values are only BINNED
        //    by Q value, that is, we first get all X's with Q values between Q_MIN
        //    and Q_MIN+1, then Q_MIN+1 to Q_MIN+2 and so on, as demonstrated here:
        //
        //    LEVEL_WEIGHT:             1.5               1
        //    Q_MIN:  0.5
        //    Q_MAX:  3
        //    X_MAX:                           2         3
        //
        //           1             1.5         1         0
        //           2               1         0         1
        //           3             2.5         1         1
        //           4               2         0         2
        //           5               3         2         0
        //           6               3         0         3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2009
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the number of components in the vector.
        //
        //    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
        //
        //    Input, int X_MAX[DIM_NUM], the maximum values allowed in each component.
        //
        //    Input/output, int X[DIM_NUM].  On first call (with MORE = FALSE),
        //    the input value of X is not important.  On subsequent calls, the
        //    input value of X should be the output value from the previous call.
        //    On output, (with MORE = TRUE), the value of X will be the "next"
        //    vector in the reverse lexicographical list of vectors that satisfy
        //    the condition.  However, on output with MORE = FALSE, the vector
        //    X is meaningless, because there are no more vectors in the list.
        //
        //    Input, double Q_MIN, Q_MAX, the lower and upper
        //    limits on the sum.
        //
        //    Input/output, bool *MORE.  On input, if the user has set MORE
        //    FALSE, the user is requesting the initiation of a new sequence
        //    of values.  If MORE is TRUE, then the user is requesting "more"
        //    values in the current sequence.  On output, if MORE is TRUE,
        //    then another value was found and returned in X, but if MORE is
        //    FALSE, then there are no more values in the sequence, and X is
        //    NOT the next value.
        //
    {
        switch (more)
        {
            //
            //  On first call, initialize the subrange.
            //
            case false:
                data.q_min2 = q_min;
                data.q_max2 = Math.Min(q_min + 1.0, q_max);
                break;
        }

        //
        //  Call a lower level function to search the subrange.
        //
        for (;;)
        {
            sgmga_vcn_naive(dim_num, level_weight, x_max, x, data.q_min2, data.q_max2,
                ref more);
            switch (more)
            {
                //
                //  If another solution was found, return it.
                //
                case true:
                    return;
            }

            //
            //  If the current subrange is exhausted, try to move to the next one.
            //
            if (data.q_max2 < q_max)
            {
                data.q_min2 = data.q_max2;
                data.q_max2 = Math.Min(data.q_max2 + 1.0, q_max);
            }
            //
            //  If there are no more subranges, we're done.
            //
            else
            {
                break;
            }
        }
    }

    public static void sgmga_weight(int dim_num, double[] level_weight, int level_max,
            int[] rule, int[] np, double[] p,
            Func<int, int, double[], double[], double[]>[] gw_compute_weights,
            int point_num, int point_total_num, int[] sparse_unique_index,
            int[] growth, double[] sparse_weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_WEIGHT computes weights for an SGMGA grid.
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
        //    26 April 2011
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
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
        //    Input, int GROWTH[DIM_NUM], the desired growth in each dimension.
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
        int dim;
        int point;
        SGMGAData data = new();

        for (point = 0; point < point_num; point++)
        {
            sparse_weight[point] = 0.0;
        }

        int point_total = 0;

        int[] level_1d = new int[dim_num];
        int[] order_1d = new int[dim_num];
        int[] level_1d_max = new int[dim_num];
        //
        //  Initialization for SGMGA_VCN_ORDERED.
        //
        double level_weight_min_pos = typeMethods.r8vec_min_pos(dim_num, level_weight);
        double q_min = level_max * level_weight_min_pos
                       - typeMethods.r8vec_sum(dim_num, level_weight);
        double q_max = level_max * level_weight_min_pos;
        for (dim = 0; dim < dim_num; dim++)
        {
            switch (level_weight[dim])
            {
                case > 0.0:
                {
                    level_1d_max[dim] = (int) Math.Floor(q_max / level_weight[dim]) + 1;
                    if (q_max <= (level_1d_max[dim] - 1) * level_weight[dim])
                    {
                        level_1d_max[dim] -= 1;
                    }

                    break;
                }
                default:
                    level_1d_max[dim] = 0;
                    break;
            }
        }

        bool more_grids = false;
        //
        //  Seek all vectors LEVEL_1D which satisfy the constraint:
        //
        //    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
        //      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
        //      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
        //
        for (;;)
        {
            sgmga_vcn_ordered(ref data, dim_num, level_weight, level_1d_max,
                ref level_1d, q_min, q_max, ref more_grids);

            if (!more_grids)
            {
                break;
            }

            //
            //  Compute the combinatorial coefficient.
            //
            double coef = sgmga_vcn_coef(dim_num, level_weight, level_1d, q_max);

            switch (coef)
            {
                case 0.0:
                    continue;
            }

            //
            //  Transform each 1D level to a corresponding 1D order.
            //
            LevelToOrder.level_growth_to_order(dim_num, level_1d, rule, growth, ref order_1d);
            //
            //  The product of the 1D orders gives us the number of points in this grid.
            //
            int order_nd = typeMethods.i4vec_product(dim_num, order_1d);
            //
            //  Compute the weights for this grid.
            //
            //  The correct transfer of data from the product grid to the sparse grid
            //  depends on the fact that the product rule weights are stored under colex
            //  order of the points, and this is the same ordering implicitly used in
            //  generating the SPARSE_UNIQUE_INDEX array.
            //
            double[] grid_weight = new double[order_nd];

            sgmga_product_weight(dim_num, order_1d, order_nd, rule,
                np, p, gw_compute_weights, ref grid_weight);
            //
            //  Add these weights to the rule.
            //
            int order;
            for (order = 0; order < order_nd; order++)
            {
                int point_unique = sparse_unique_index[point_total];

                point_total += 1;

                sparse_weight[point_unique] += coef * grid_weight[order];
            }

        }
    }

    public static void sgmga_write(int dim_num, double[] level_weight, int[] rule, int[] np,
            double[] p, int point_num, double[] sparse_weight, double[] sparse_point,
            string file_name)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_WRITE writes an SGMGA rule to six files.
        //
        //  Discussion:
        //
        //    The files are:
        //    * the "A" file stores the anisotropic weights, as a DIM_NUM x 1 list.
        //    * the "N" file stores the NP values, as a DIM_NUM x 1 list.
        //    * the "P" file stores the P values, as a sum(NP[*]) x 1 list.
        //    * the "R" file stores the region, as a DIM_NUM x 2 list;
        //    * the "W" file stores the weights as a POINT_NUM list;
        //    * the "X" file stores the abscissas as a DIM_NUM x POINT_NUM list.
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
        //    09 June 2010
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
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
        //    Differential Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2411-2442.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
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

        double[] sparse_region = new double[dim_num * 2];

        for (dim = 0; dim < dim_num; dim++)
        {
            int point;
            double t1;
            double t2;
            switch (rule[dim])
            {
                case 1:
                    sparse_region[dim + 0 * dim_num] = -1.0;
                    sparse_region[dim + 1 * dim_num] = +1.0;
                    break;
                case 2:
                    sparse_region[dim + 0 * dim_num] = -1.0;
                    sparse_region[dim + 1 * dim_num] = +1.0;
                    break;
                case 3:
                    sparse_region[dim + 0 * dim_num] = -1.0;
                    sparse_region[dim + 1 * dim_num] = +1.0;
                    break;
                case 4:
                    sparse_region[dim + 0 * dim_num] = -1.0;
                    sparse_region[dim + 1 * dim_num] = +1.0;
                    break;
                case 5:
                    sparse_region[dim + 0 * dim_num] = -typeMethods.r8_huge();
                    sparse_region[dim + 1 * dim_num] = +typeMethods.r8_huge();
                    break;
                case 6:
                    sparse_region[dim + 0 * dim_num] = -typeMethods.r8_huge();
                    sparse_region[dim + 1 * dim_num] = +typeMethods.r8_huge();
                    break;
                case 7:
                    sparse_region[dim + 0 * dim_num] = 0.0;
                    sparse_region[dim + 1 * dim_num] = typeMethods.r8_huge();
                    break;
                case 8:
                    sparse_region[dim + 0 * dim_num] = 0.0;
                    sparse_region[dim + 1 * dim_num] = typeMethods.r8_huge();
                    break;
                case 9:
                    sparse_region[dim + 0 * dim_num] = -1.0;
                    sparse_region[dim + 1 * dim_num] = +1.0;
                    break;
                case 10:
                    sparse_region[dim + 0 * dim_num] = -typeMethods.r8_huge();
                    sparse_region[dim + 1 * dim_num] = +typeMethods.r8_huge();
                    break;
                //
                //  Best guess as to region extent for rules of type 11 or 12.
                //
                case 11:
                {
                    t1 = typeMethods.r8_huge();
                    t2 = -typeMethods.r8_huge();
                    for (point = 0; point < point_num; point++)
                    {
                        t1 = Math.Min(t1, sparse_point[dim + point * dim_num]);
                        t2 = Math.Max(t2, sparse_point[dim + point * dim_num]);
                    }

                    sparse_region[dim + 0 * dim_num] = t1;
                    sparse_region[dim + 1 * dim_num] = t2;
                    break;
                }
                case 12:
                {
                    t1 = typeMethods.r8_huge();
                    t2 = -typeMethods.r8_huge();
                    for (point = 0; point < point_num; point++)
                    {
                        t1 = Math.Min(t1, sparse_point[dim + point * dim_num]);
                        t2 = Math.Max(t2, sparse_point[dim + point * dim_num]);
                    }

                    sparse_region[dim + 0 * dim_num] = t1;
                    sparse_region[dim + 1 * dim_num] = t2;
                    break;
                }
                default:
                    Console.WriteLine("");
                    Console.WriteLine("SGMGA_WRITE - Fatal error!");
                    Console.WriteLine("  Unexpected value of RULE[" + dim + "] = "
                                      + rule[dim] + ".");
                    return;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("SGMGA_WRITE:");

        string file_name_a = file_name + "_a.txt";
        typeMethods.r8mat_write(file_name_a, 1, dim_num, level_weight);
        Console.WriteLine("  Wrote the A file = \"" + file_name_a + "\".");

        string file_name_n = file_name + "_n.txt";
        typeMethods.i4mat_write(file_name_n, 1, dim_num, np);
        Console.WriteLine("  Wrote the N file = \"" + file_name_n + "\".");

        int np_sum = typeMethods.i4vec_sum(dim_num, np);
        string file_name_p = file_name + "_p.txt";
        typeMethods.r8mat_write(file_name_p, 1, np_sum, p);
        Console.WriteLine("  Wrote the P file = \"" + file_name_p + "\".");

        string file_name_r = file_name + "_r.txt";
        typeMethods.r8mat_write(file_name_r, dim_num, 2, sparse_region);
        Console.WriteLine("  Wrote the R file = \"" + file_name_r + "\".");

        string file_name_w = file_name + "_w.txt";
        typeMethods.r8mat_write(file_name_w, 1, point_num, sparse_weight);
        Console.WriteLine("  Wrote the W file = \"" + file_name_w + "\".");

        string file_name_x = file_name + "_x.txt";
        typeMethods.r8mat_write(file_name_x, dim_num, point_num, sparse_point);
        Console.WriteLine("  Wrote the X file = \"" + file_name_x + "\".");

    }
}