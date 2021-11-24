using System;
using System.Globalization;
using Burkardt.Composition;
using Burkardt.Quadrature;
using Burkardt.Sparse;

namespace SandiaSparseTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SANDIA_SPARSE_TEST.
        //
        //  Discussion:
        //
        //    SANDIA_SPARSE_TEST tests the SANDIA_SPARSE library.
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
    {
        int dim_num;
        int level_max;
        int rule;
        const int rule_max = 7;

        Console.WriteLine("");
        Console.WriteLine("SANDIA_SPARSE_TEST");
        Console.WriteLine("  Test the SANDIA_SPARSE library.");
        switch (true)
        {
            //
            //  Test LEVELS_INDEX_SIZE for one example each of CFN, OFN, OWN and ONN rules.
            //
            case true:
                rule = 1;

                int dim_min = 1;
                int dim_max = 1;
                int level_max_min = 0;
                int level_max_max = 10;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                dim_min = 1;
                dim_max = 6;
                level_max_min = 0;
                level_max_max = 6;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                dim_min = 6;
                dim_max = 10;
                level_max_min = 0;
                level_max_max = 5;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                dim_min = 100;
                dim_max = 100;
                level_max_min = 0;
                level_max_max = 2;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                rule = 2;

                dim_min = 1;
                dim_max = 1;
                level_max_min = 0;
                level_max_max = 10;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                dim_min = 1;
                dim_max = 6;
                level_max_min = 0;
                level_max_max = 6;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                dim_min = 6;
                dim_max = 10;
                level_max_min = 0;
                level_max_max = 5;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                dim_min = 100;
                dim_max = 100;
                level_max_min = 0;
                level_max_max = 2;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                rule = 5;

                dim_min = 1;
                dim_max = 1;
                level_max_min = 0;
                level_max_max = 10;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                dim_min = 1;
                dim_max = 6;
                level_max_min = 0;
                level_max_max = 6;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                dim_min = 6;
                dim_max = 10;
                level_max_min = 0;
                level_max_max = 5;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                dim_min = 100;
                dim_max = 100;
                level_max_min = 0;
                level_max_max = 2;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                rule = 7;

                dim_min = 1;
                dim_max = 1;
                level_max_min = 0;
                level_max_max = 10;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                dim_min = 1;
                dim_max = 6;
                level_max_min = 0;
                level_max_max = 6;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                dim_min = 6;
                dim_max = 10;
                level_max_min = 0;
                level_max_max = 5;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);

                dim_min = 100;
                dim_max = 100;
                level_max_min = 0;
                level_max_max = 2;

                levels_index_size_test(rule, dim_min, dim_max, level_max_min,
                    level_max_max);
                break;
        }

        switch (true)
        {
            //
            //  Test LEVELS_INDEX for one example each of CFN, OFN, OWN and ONN rules.
            //
            case true:
                rule = 1;

                dim_num = 2;
                level_max = 1;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 2;
                level_max = 3;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 3;
                level_max = 0;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 3;
                level_max = 2;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 6;
                level_max = 2;
                levels_index_test(rule, dim_num, level_max);

                rule = 2;

                dim_num = 2;
                level_max = 1;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 2;
                level_max = 3;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 3;
                level_max = 0;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 3;
                level_max = 2;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 6;
                level_max = 2;
                levels_index_test(rule, dim_num, level_max);

                rule = 5;

                dim_num = 2;
                level_max = 1;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 2;
                level_max = 3;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 3;
                level_max = 0;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 3;
                level_max = 2;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 6;
                level_max = 2;
                levels_index_test(rule, dim_num, level_max);

                rule = 7;

                dim_num = 2;
                level_max = 1;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 2;
                level_max = 3;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 3;
                level_max = 0;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 3;
                level_max = 2;
                levels_index_test(rule, dim_num, level_max);

                dim_num = 6;
                level_max = 2;
                levels_index_test(rule, dim_num, level_max);
                break;
        }

        switch (true)
        {
            //
            //  Test SPARSE_GRID by having it compute a few sparse grids based on each rule.
            //
            case true:
            {
                for (rule = 1; rule <= rule_max; rule++)
                {
                    dim_num = 2;
                    level_max = 1;
                    sparse_grid_compute_test(rule, dim_num, level_max);

                    dim_num = 2;
                    level_max = 2;
                    sparse_grid_compute_test(rule, dim_num, level_max);

                    dim_num = 3;
                    level_max = 1;
                    sparse_grid_compute_test(rule, dim_num, level_max);
                }

                break;
            }
        }

        switch (true)
        {
            //
            //  Test SPARSE_GRID by having it compute a few sparse grids based on each rule,
            //  and comparing the sum of the quadrature weights to the expected sum.
            //
            case true:
            {
                for (rule = 1; rule <= rule_max; rule++)
                {
                    dim_num = 2;
                    level_max = 4;
                    sparse_grid_weight_test(rule, dim_num, level_max);

                    dim_num = 3;
                    level_max = 0;
                    sparse_grid_weight_test(rule, dim_num, level_max);

                    dim_num = 3;
                    level_max = 1;
                    sparse_grid_weight_test(rule, dim_num, level_max);

                    dim_num = 3;
                    level_max = 6;
                    sparse_grid_weight_test(rule, dim_num, level_max);

                    dim_num = 10;
                    level_max = 3;
                    sparse_grid_weight_test(rule, dim_num, level_max);
                }

                break;
            }
        }

        switch (true)
        {
            //
            //  Test SPARSE_GRID by having it compute a few sparse grids based on each rule,
            //  and comparing estimated and exact monomial integrals.
            //
            case true:
            {
                for (rule = 1; rule <= rule_max; rule++)
                {
                    dim_num = 2;
                    level_max = 0;
                    int degree_max = 3;
                    sparse_grid_monomial_test(rule, dim_num, level_max, degree_max);

                    dim_num = 2;
                    level_max = 1;
                    degree_max = 5;
                    sparse_grid_monomial_test(rule, dim_num, level_max, degree_max);

                    dim_num = 2;
                    level_max = 2;
                    degree_max = 7;
                    sparse_grid_monomial_test(rule, dim_num, level_max, degree_max);

                    dim_num = 2;
                    level_max = 3;
                    degree_max = 9;
                    sparse_grid_monomial_test(rule, dim_num, level_max, degree_max);

                    dim_num = 3;
                    level_max = 0;
                    degree_max = 2;
                    sparse_grid_monomial_test(rule, dim_num, level_max, degree_max);

                    dim_num = 3;
                    level_max = 1;
                    degree_max = 4;
                    sparse_grid_monomial_test(rule, dim_num, level_max, degree_max);

                    dim_num = 3;
                    level_max = 2;
                    degree_max = 6;
                    sparse_grid_monomial_test(rule, dim_num, level_max, degree_max);
                }

                break;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("SANDIA_SPARSE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void levels_index_size_test(int rule, int dim_min, int dim_max,
            int level_max_min, int level_max_max)

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    LEVELS_INDEX_SIZE_TEST tests LEVELS_INDEX_SIZE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
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
        //    Input, int DIM_MIN, the minimum spatial dimension.
        //
        //    Input, int DIM_MAX, the maximum spatial dimension.
        //
        //    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX.
        //
        //    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX.
        //
    {
        int dim_num;
        int level_max;

        Console.WriteLine("");
        Console.WriteLine("LEVELS_INDEX_SIZE_TEST");
        Console.WriteLine("  LEVELS_INDEX_SIZE returns the number of distinct");
        Console.WriteLine("  points in a sparse grid derived from a 1D rule.");
        Console.WriteLine("");
        Console.WriteLine("  We are looking at rules like rule " + rule + "");
        Console.WriteLine("");
        Console.WriteLine("  Each sparse grid is of spatial dimension DIM,");
        Console.WriteLine("  and is made up of product grids such that");
        Console.WriteLine("  LEVEL_MIN <= LEVEL <= LEVEL_MAX.");

        Console.WriteLine("");
        string cout = "   DIM: ";
        for (dim_num = dim_min; dim_num <= dim_max; dim_num++)
        {
            cout += "  " + dim_num.ToString().PadLeft(8);
        }

        Console.WriteLine(cout);
        Console.WriteLine("   LEVEL_MAX");
        Console.WriteLine("   ---------");

        for (level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            cout = "    " + level_max.ToString().PadLeft(4);
            for (dim_num = dim_min; dim_num <= dim_max; dim_num++)
            {
                cout += "  " + Grid.levels_index_size(dim_num, level_max, rule).ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
        }

    }

    private static void levels_index_test(int rule, int dim_num, int level_max)

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    LEVELS_INDEX_TEST tests LEVELS_INDEX.
        //
        //  Discussion:
        //
        //    The routine computes the indices of the unique points used in a sparse 
        //    multidimensional grid whose size is controlled by a parameter LEVEL_MAX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
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
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the level.
        //
    {
        int point;

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        Console.WriteLine("");
        Console.WriteLine("LEVELS_INDEX_TEST");
        Console.WriteLine("  LEVELS_INDEX returns all grid indexes");
        Console.WriteLine("  whose level value satisfies");
        Console.WriteLine("    LEVEL_MIN <= LEVEL <= LEVEL_MAX.");
        Console.WriteLine("  Here, LEVEL is the sum of the levels of the 1D rules,");
        Console.WriteLine("  and the order of the rule is 2^LEVEL + 1.");
        Console.WriteLine("");
        Console.WriteLine("  We are looking at rules like rule " + rule + "");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        Console.WriteLine("  LEVEL_MIN =                 " + level_min + "");
        Console.WriteLine("  LEVEL_MAX =                 " + level_max + "");

        int point_num = Grid.levels_index_size(dim_num, level_max, rule);

        Console.WriteLine("  Unique points in the grid = " + point_num + "");
        //
        //  Compute the orders and points.
        //
        int[] grid_base = new int[dim_num * point_num];
        int[] grid_index = new int[dim_num * point_num];

        Grid.levels_index(dim_num, level_max, rule, point_num, ref grid_index,
            ref grid_base);
        //
        //  Now we're done.  Print the merged grid data.
        //
        Console.WriteLine("");
        Console.WriteLine(" Point     Grid indices:");
        Console.WriteLine("           Grid bases:");
        Console.WriteLine("");
        for (point = 0; point < point_num; point++)
        {
            string cout = "  " + point.ToString().PadLeft(4) + "  ";
            int dim;
            for (dim = 0; dim < dim_num; dim++)
            {
                cout += grid_index[dim + point * dim_num].ToString().PadLeft(6);
            }

            Console.WriteLine(cout);
            cout = "        ";
            for (dim = 0; dim < dim_num; dim++)
            {
                cout += grid_base[dim + point * dim_num].ToString().PadLeft(6);
            }

            Console.WriteLine(cout);
        }
    }

    private static void sparse_grid_compute_test(int rule, int dim_num, int level_max)

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    SPARSE_GRID_COMPUTE_TEST computes and prints a sparse grid rule.
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
        //  Parameters:
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
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the level.
        //
    {
        int point;

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        Console.WriteLine("");
        Console.WriteLine("SPARSE_GRID_COMPUTE_TEST:");
        Console.WriteLine("  SPARSE_GRID can make a sparse grid.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        Console.WriteLine("  LEVEL_MIN =                 " + level_min + "");
        Console.WriteLine("  LEVEL_MAX =                 " + level_max + "");
        Console.WriteLine("  1D quadrature index RULE =  " + rule + "");
        //
        //  Determine the number of points.
        //
        int point_num = Grid.levels_index_size(dim_num, level_max, rule);

        Console.WriteLine("  Unique points in the grid = " + point_num + "");
        //
        //  Allocate space for the weights and points.
        //
        double[] grid_weight = new double[point_num];
        double[] grid_point = new double[dim_num * point_num];
        //
        //  Compute the weights and points.
        //
        Grid.sparse_grid(dim_num, level_max, rule, point_num, ref grid_weight,
            ref grid_point);
        //
        //  Print them out.
        //
        Console.WriteLine("");
        Console.WriteLine("  Grid weights:");
        Console.WriteLine("");
        for (point = 0; point < point_num; point++)
        {
            Console.WriteLine("  " + point.ToString().PadLeft(4)
                                   + "  " + grid_weight[point].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Grid points:");
        Console.WriteLine("");
        for (point = 0; point < point_num; point++)
        {
            string cout = "  " + point + "  ";
            int dim;
            for (dim = 0; dim < dim_num; dim++)
            {
                cout += grid_point[dim + point * dim_num].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }
    }

    private static void sparse_grid_weight_test(int rule, int dim_num, int level_max)

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    SPARSE_GRID_WEIGHT_TEST checks the sum of the quadrature weights.
        //
        //  Discussion:
        //
        //    This routine gets the sparse grid indices and determines the 
        //    corresponding sparse grid abscissas.
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
        //  Parameters:
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
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the level.
        //
    {
        int point;
        double weight_sum_exact = 0;

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        weight_sum_exact = rule switch
        {
            >= 1 and <= 5 => Math.Pow(2.0, dim_num),
            6 => Math.Sqrt(Math.Pow(Math.PI, dim_num)),
            7 => 1.0,
            _ => weight_sum_exact
        };

        Console.WriteLine("");
        Console.WriteLine("SPARSE_GRID_WEIGHT_TEST:");
        Console.WriteLine("  Compute the weights of a sparse grid.");
        Console.WriteLine("");
        Console.WriteLine("  As a simple test, sum these weights.");
        Console.WriteLine("  They should sum to exactly " + weight_sum_exact + "");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        Console.WriteLine("  LEVEL_MIN =                 " + level_min + "");
        Console.WriteLine("  LEVEL_MAX =                 " + level_max + "");
        Console.WriteLine("  1D quadrature index RULE =  " + rule + "");
        //
        //  Determine the number of points.
        //
        int point_num = Grid.levels_index_size(dim_num, level_max, rule);

        Console.WriteLine("  Unique points in the grid = " + point_num + "");
        //
        //  Allocate space for the weights and points.
        //
        double[] grid_weight = new double[point_num];
        double[] grid_point = new double[dim_num * point_num];
        //
        //  Compute the weights and points.
        //
        Grid.sparse_grid(dim_num, level_max, rule, point_num, ref grid_weight,
            ref grid_point);
        //
        //  Sum the weights.
        //
        double weight_sum = 0.0;
        for (point = 0; point < point_num; point++)
        {
            weight_sum += grid_weight[point];
        }

        double weight_sum_error = Math.Abs(weight_sum - weight_sum_exact);

        Console.WriteLine("");
        Console.WriteLine("    Weight sum  Expected sum    Difference");
        Console.WriteLine("");
        Console.WriteLine("  " + weight_sum.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + weight_sum_exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + weight_sum_error.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void sparse_grid_monomial_test(int rule, int dim_num, int level_max,
            int degree_max)

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    SPARSE_GRID_MONOMIAL_TEST tests monomial exactness of the sparse grid rules.
        //
        //  Discussion:
        //
        //    This test is going to check EVERY monomial of total degree DEGREE_MAX
        //    or less.  Even for a moderately high dimension of DIM_NUM = 10, you
        //    do NOT want to use a large value of DEGREE_MAX, since there are
        //
        //      1         monomials of total degree 0,
        //      DIM_NUM   monomials of total degree 1,
        //      DIM_NUM^2 monomials of total degree 2,
        //      DIM_NUM^3 monomials of total degree 3, and so on.
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
        //  Parameters:
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
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the level.
        //
        //    Input, int DEGREE_MAX, the maximum monomial total 
        //    degree to check.
        //
    {
        int degree;
        int h = 0;
        int t = 0;

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        Console.WriteLine("");
        Console.WriteLine("SPARSE_GRID_MONOMIAL_TEST");
        Console.WriteLine("  Check the exactness of a sparse grid quadrature rule,");
        Console.WriteLine("  applied to all monomials of orders 0 to DEGREE_MAX.");
        Console.WriteLine("");
        Console.WriteLine("  For cases where the dimension is greater than 1,");
        Console.WriteLine("  many sparse grid of this level have accuracy through");
        Console.WriteLine("  monomials of total degree   " + 2 * level_max + 1 + "");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        Console.WriteLine("  LEVEL_MIN =                 " + level_min + "");
        Console.WriteLine("  LEVEL_MAX =                 " + level_max + "");
        Console.WriteLine("  1D quadrature index RULE =  " + rule + "");
        Console.WriteLine("  Check up to DEGREE_MAX =    " + degree_max + "");
        //
        //  Determine the number of points in the rule.
        //
        int point_num = Grid.levels_index_size(dim_num, level_max, rule);

        Console.WriteLine("  Unique points in the grid = " + point_num + "");
        //
        //  Allocate space for the weights and points.
        //
        double[] grid_weight = new double[point_num];
        double[] grid_point = new double[dim_num * point_num];
        //
        //  Compute the weights and points.
        //
        Grid.sparse_grid(dim_num, level_max, rule, point_num, ref grid_weight, ref grid_point);
        //
        //  Compare exact and estimated values of the integrals of various monomials.
        //
        int[] expon = new int[dim_num];

        Console.WriteLine("");
        Console.WriteLine("      Error      Total   Monomial");
        Console.WriteLine("                 Degree  Exponents");
        Console.WriteLine("");

        for (degree = 0; degree <= degree_max; degree++)
        {
            bool more = false;

            for (;;)
            {
                Comp.comp_next(degree, dim_num, ref expon, ref more, ref h, ref t);

                double quad_error = MonomialQuadrature.monomial_quadrature(dim_num, expon, point_num,
                    grid_weight, grid_point, rule);

                string cout = "  " + quad_error.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + degree.ToString().PadLeft(2)
                                   + "  ";

                int dim;
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += expon[dim].ToString().PadLeft(2);
                }

                Console.WriteLine(cout);

                if (!more)
                {
                    break;
                }
            }

            switch (dim_num)
            {
                case > 1:
                    Console.WriteLine("");
                    break;
            }
        }
    }
}