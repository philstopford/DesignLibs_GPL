using System;
using System.Globalization;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace ClenshawCurtisGridTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80*
        //
        //  Purpose:
        //
        //    MAIN is the main program for CLENSHAW_CURTIS_GRID_TEST.
        //
        //  Discussion:
        //
        //    CLENSHAW_CURTIS_TEST tests the CLENSHAW_CURTIS routines.
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
    {
        Console.WriteLine("");
        Console.WriteLine("CLENSHAW_CURTIS_GRID_TEST");
        Console.WriteLine("  Test the routines in the CLENSHAW_CURTIS_GRID library.");

        test005();
        test01();
        test015();
        test02();
        test025();
        test03();
        test035();
        test04();
        test045();
        test05();
        test06();
        test07();
        test08();
        test09();
        test10();
        test11();
        test12();
        test13();

        Console.WriteLine("");
        Console.WriteLine("CLENSHAW_CURTIS_GRID_TEST");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }

    private static void test005()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST005 calls CC_GRID for the 1D problem.
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
    {
        const int DIM_NUM = 1;

        int i;
        int[] order_1d = new int[DIM_NUM];

        Console.WriteLine("");
        Console.WriteLine("TEST005:");
        Console.WriteLine("  CC_GRID returns a grid of Clenshaw-Curtis points.");
        Console.WriteLine("  Here, we simply call for grids in the 1D case");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension of grid = " + DIM_NUM + "");

        for (i = 1; i <= 10; i++)
        {
            Console.WriteLine("");

            order_1d[0] = i;
            int order_nd = order_1d[0];

            double[] grid_point = new double[order_nd];

            ClenshawCurtisGrid.cc_grid(DIM_NUM, order_1d, order_nd, ref grid_point);

            int j;
            for (j = 0; j < order_nd; j++)
            {
                string cout = "  " + (j + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8);
                int dim;
                for (dim = 0; dim < DIM_NUM; dim++)
                {
                    cout += "  " + grid_point[dim + j * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12);
                }

                Console.WriteLine(cout);
            }
        }
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 simply calls CC_GRID once.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 3;

        int dim;
        int[] order_1d =  {
                3, 4, 2
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  CC_GRID returns a grid of Clenshaw-Curtis points.");
        Console.WriteLine("  Here, we simply call for a specific grid.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension of grid = " + DIM_NUM + "");

        int order_nd = 1;
        for (dim = 0; dim < DIM_NUM; dim++)
        {
            order_nd *= order_1d[dim];
        }

        double[] grid_point = new double[DIM_NUM * order_nd];

        Console.WriteLine("");
        Console.WriteLine("  Total number of points in the grid = " + order_nd + "");
        Console.WriteLine("");

        ClenshawCurtisGrid.cc_grid(DIM_NUM, order_1d, order_nd, ref grid_point);

        int j = 1;
        int q = 0;
        for (dim = 0; dim < DIM_NUM; dim++)
        {
            q += order_1d[dim];
        }

        Console.WriteLine("");
        Console.WriteLine("         I         Q          Grid orders:");
        Console.WriteLine("");

        string cout = "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                           + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(8);
        for (dim = 0; dim < DIM_NUM; dim++)
        {
            cout += "  " + order_1d[dim].ToString(CultureInfo.InvariantCulture).PadLeft(8);
        }

        Console.WriteLine(cout);

        Console.WriteLine("");
        Console.WriteLine("  Grid points:");
        Console.WriteLine("");

        for (j = 0; j < order_nd; j++)
        {
            cout = "  " + (j + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8);
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                cout += "  " + grid_point[dim + j * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test015()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST015 calls CC_GRID_INDEX.
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
    {
        const int DIM_NUM = 3;

        int dim;
        int[] order_1d =  {
                3, 4, 2
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST015:");
        Console.WriteLine("  CC_GRID_INDEX returns an indexed grid of Clenshaw-Curtis points.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension of grid = " + DIM_NUM + "");

        int order_nd = 1;
        for (dim = 0; dim < DIM_NUM; dim++)
        {
            order_nd *= order_1d[dim];
        }

        int[] grid_index = new int[DIM_NUM * order_nd];

        Console.WriteLine("");
        Console.WriteLine("  Total number of points in the grid = " + order_nd + "");
        Console.WriteLine("");

        ClenshawCurtisGrid.cc_grid_index(DIM_NUM, order_1d, order_nd, ref grid_index);

        int j = 1;
        int q = 0;
        for (dim = 0; dim < DIM_NUM; dim++)
        {
            q += order_1d[dim];
        }

        Console.WriteLine("");
        Console.WriteLine("         I         Q          Grid orders:");
        Console.WriteLine("");

        string cout = "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                           + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  ";
        for (dim = 0; dim < DIM_NUM; dim++)
        {
            cout += "  " + order_1d[dim].ToString(CultureInfo.InvariantCulture).PadLeft(4);
        }

        Console.WriteLine(cout);

        Console.WriteLine("");
        Console.WriteLine("  Grid indexed points:");
        Console.WriteLine("");

        for (j = 0; j < order_nd; j++)
        {
            cout = "  " + (j + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8) + "            ";
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                cout += "  " + grid_index[dim + j * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 calls CC_GRIDS_MINMAX to collect all points on 2D grids for Q = 3 to 5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 2;

        int dim;
        int grid_num = 0;
        int j;
        int point_num = 0;
        const int q_max = 5;
        const int q_min = 3;

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  CC_GRIDS_MINMAX returns all Clenshaw Curtis grids");
        Console.WriteLine("  whose Q value satisfies Q_MIN <= Q <= Q_MAX.");
        Console.WriteLine("  Here, Q is the sum of the orders of the 1D rules, and");
        Console.WriteLine("  Q_MIN = " + q_min + "");
        Console.WriteLine("  Q_MAX = " + q_max + "");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension of grids = " + DIM_NUM + "");

        ClenshawCurtisGrid.cc_grids_minmax_size(DIM_NUM, q_min, q_max, ref grid_num, ref point_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of grids = " + grid_num + "");
        Console.WriteLine("  Number of points in the grids = " + point_num + "");
        //
        //  Allocate the space.
        //
        int[] grid_order = new int[DIM_NUM * grid_num];
        double[] grid_point = new double[DIM_NUM * point_num];
        //
        //  Compute the orders and points.
        //
        ClenshawCurtisGrid.cc_grids_minmax(DIM_NUM, q_min, q_max, grid_num, point_num,
            ref grid_order, ref grid_point);
        //
        //  Now we're done.  Print the merged grid data.
        //
        Console.WriteLine("");
        Console.WriteLine("         I         Q          Grid orders:");
        Console.WriteLine("");
        for (j = 0; j < grid_num; j++)
        {
            int q = 0;
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                q += grid_order[dim + j * DIM_NUM];
            }

            string cout = "  " + (j + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(8);
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                cout += "  " + grid_order[dim + j * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(8);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Grid points:");
        Console.WriteLine("");

        for (j = 0; j < point_num; j++)
        {
            string cout = "  " + (j + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8);
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                cout += "  " + grid_point[dim + j * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test025()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST025 calls CC_LEVELS_MINMAX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 November 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 2;
        const int TEST_NUM = 3;

        int grid_num = 0;
        int[] level_max_test =  {
                2, 3, 3
            }
            ;
        int[] level_min_test =  {
                2, 0, 3
            }
            ;
        int point_num = 0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST025:");
        Console.WriteLine("  CC_LEVELS_MINMAX returns all Clenshaw Curtis grids");
        Console.WriteLine("  whose LEVEL value satisfies");
        Console.WriteLine("    LEVEL_MIN <= LEVEL <= LEVEL_MAX.");
        Console.WriteLine("  Here, LEVEL is the sum of the levels of the 1D rules,");
        Console.WriteLine("  and the order of the rule is 2**LEVEL - 1.");

        for (test = 0; test < TEST_NUM; test++)
        {
            int level_min = level_min_test[test];
            int level_max = level_max_test[test];

            Console.WriteLine("");
            Console.WriteLine("  LEVEL_MIN = " + level_min + "");
            Console.WriteLine("  LEVEL_MAX = " + level_max + "");
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension of grids = " + DIM_NUM + "");

            ClenshawCurtisGrid.cc_levels_minmax_size(DIM_NUM, level_min, level_max, ref grid_num,
                ref point_num);

            Console.WriteLine("");
            Console.WriteLine("  Number of grids = " + grid_num + "");
            Console.WriteLine("  Number of points in the grids = " + point_num + "");
            //
            //  Allocate the space.
            //
            int[] grid_level = new int[DIM_NUM * grid_num];
            int[] grid_order = new int[DIM_NUM * grid_num];
            double[] grid_point = new double[DIM_NUM * point_num];
            //
            //  Compute the orders and points.
            //
            ClenshawCurtisGrid.cc_levels_minmax(DIM_NUM, level_min, level_max, grid_num, point_num,
                ref grid_level, ref grid_order, ref grid_point);
            //
            //  Now we're done.  Print the merged grid data.
            //
            Console.WriteLine("");
            Console.WriteLine("      Grid     Level           Grid Levels         Grid orders:");
            Console.WriteLine("      ----     -----          ------------        ------------");
            Console.WriteLine("");
            int dim;
            int j;
            for (j = 0; j < grid_num; j++)
            {
                int level = 0;
                for (dim = 0; dim < DIM_NUM; dim++)
                {
                    level += grid_level[dim + j * DIM_NUM];
                }

                string cout = "  " + (j + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + level.ToString(CultureInfo.InvariantCulture).PadLeft(8);
                for (dim = 0; dim < DIM_NUM; dim++)
                {
                    cout += "  " + grid_level[dim + j * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(8);
                }

                for (dim = 0; dim < DIM_NUM; dim++)
                {
                    cout += "  " + grid_order[dim + j * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(8);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  Grid points:");
            Console.WriteLine("");

            for (j = 0; j < point_num; j++)
            {
                string cout = "  " + (j + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8);
                for (dim = 0; dim < DIM_NUM; dim++)
                {
                    cout += "  " + grid_point[dim + j * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12);
                }

                Console.WriteLine(cout);
            }
        }
    }

    private static void test03()

        //****************************************************************************80*
        //
        //  Purpose:
        //
        //    TEST03 calls CC_GRIDS_CONSTRAINED to collect constrained grids.
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
    {
        const int DIM_NUM = 2;

        double[] alpha =  {
                2.0, 3.0
            }
            ;
        int dim;
        int grid_num = 0;
        int j;
        int[] order_max = new int[DIM_NUM];
        int[] order_min = new int[DIM_NUM];
        int point_num = 0;
        const double q_max = 13.0;

        Console.WriteLine("");
        Console.WriteLine("TEST03:");
        Console.WriteLine("  CC_GRIDS_CONSTRAINED returns all Clenshaw Curtis grids");
        Console.WriteLine("  satisfying a set of constraints.");
        Console.WriteLine("");
        Console.WriteLine("  ORDER(I), the order of the 1D rule in dimension I,");
        Console.WriteLine("  is constrained by ");
        Console.WriteLine("");
        Console.WriteLine("    ORDER_MIN(I) <= ORDER(I) <= ORDER_MAX(I)");
        Console.WriteLine("");
        Console.WriteLine("  We also define the total weighted order Q");
        Console.WriteLine("");
        Console.WriteLine("    Q = ALPHA(1) * ORDER(1) + ... + ALPHA(N) * ORDER(N)");
        Console.WriteLine("");
        Console.WriteLine("  and further constrain our grids to satisfy");
        Console.WriteLine("");
        Console.WriteLine("    Q <= Q_MAX = " + q_max + "");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension of grids = " + DIM_NUM + "");

        for (dim = 0; dim < DIM_NUM; dim++)
        {
            order_min[dim] = 1;
        }

        for (dim = 0; dim < DIM_NUM; dim++)
        {
            order_max[dim] = 5;
        }

        Console.WriteLine("");
        Console.WriteLine(" Dimension Order_min Order_max     Alpha");
        Console.WriteLine("");

        for (dim = 0; dim < DIM_NUM; dim++)
        {
            Console.WriteLine("  " + (dim + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + order_min[dim].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + order_max[dim].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + alpha[dim].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        ClenshawCurtisGrid.cc_grids_constrained_size(DIM_NUM, q_max, alpha,
            order_min, order_max, ref grid_num, ref point_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of grids = " + grid_num + "");
        Console.WriteLine("  Number of points in the grids = " + point_num + "");
        //
        //  Allocate the space.
        //
        int[] grid_order = new int[DIM_NUM * grid_num];
        double[] grid_point = new double[DIM_NUM * point_num];

        ClenshawCurtisGrid.cc_grids_constrained(DIM_NUM, q_max, alpha,
            order_min, order_max, grid_num, point_num, ref grid_order, ref grid_point);
        //
        //  Now we're done.  Print the merged grid data.
        //
        Console.WriteLine("");
        Console.WriteLine("         I               Q          Grid orders:");
        Console.WriteLine("");
        for (j = 0; j < grid_num; j++)
        {
            double q = 0.0;
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                q += alpha[dim] * grid_order[dim + j * DIM_NUM];
            }

            string cout = "  " + (j + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                cout += "  " + grid_order[dim + j * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(8);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Grid points:");
        Console.WriteLine("");
        for (j = 0; j < point_num; j++)
        {
            string cout = "  " + (j + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8);
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                cout += "  " + grid_point[dim + j * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test035()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST035 calls CC_LEVELS_CONSTRAINED to collect constrained grids.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 November 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 2;

        double[] alpha =  {
                2.0, 3.0
            }
            ;
        int dim;
        int grid_num = 0;
        int j;
        int[] level_max = new int[DIM_NUM];
        int[] level_min = new int[DIM_NUM];
        int point_num = 0;
        const double q_max = 13.0;

        Console.WriteLine("");
        Console.WriteLine("TEST035:");
        Console.WriteLine("  CC_LEVELS_CONSTRAINED returns all Clenshaw Curtis grids");
        Console.WriteLine("  satisfying a set of constraints.");
        Console.WriteLine("");
        Console.WriteLine("  The constraint on the levels of the 1D Clenshaw Curtis");
        Console.WriteLine("  rule in spatial dimension I is:");
        Console.WriteLine("");
        Console.WriteLine("    LEVEL_MIN(I) <= LEVEL(I) <= LEVEL_MAX(I) ");
        Console.WriteLine("");
        Console.WriteLine("  The constraint on the levels making up a rule is:");
        Console.WriteLine("");
        Console.WriteLine("    Sum ( 1 <= I <= DIM_NUM ) ALPHA(I) * LEVEL(I) <= Q_MAX.");
        Console.WriteLine("");
        Console.WriteLine("    where Q_MAX = " + q_max + "");
        Console.WriteLine("");
        Console.WriteLine("  The relationship of level to order is roughly ");
        Console.WriteLine("");
        Console.WriteLine("    ORDER = 2^LEVEL+1.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension of grids = " + DIM_NUM + "");

        for (dim = 0; dim < DIM_NUM; dim++)
        {
            level_min[dim] = 1;
        }

        for (dim = 0; dim < DIM_NUM; dim++)
        {
            level_max[dim] = 5;
        }

        Console.WriteLine("");
        Console.WriteLine(" Dimension Level_min Level_max     Alpha");
        Console.WriteLine("");

        for (dim = 0; dim < DIM_NUM; dim++)
        {
            Console.WriteLine("  " + (dim + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + level_min[dim].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + level_max[dim].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + alpha[dim].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        ClenshawCurtisGrid.cc_levels_constrained_size(DIM_NUM, q_max, alpha,
            level_min, level_max, ref grid_num, ref point_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of grids = " + grid_num + "");
        Console.WriteLine("  Number of points in the grids = " + point_num + "");
        //
        //  Allocate the space.
        //
        int[] grid_level = new int[DIM_NUM * grid_num];
        double[] grid_point = new double[DIM_NUM * point_num];

        ClenshawCurtisGrid.cc_levels_constrained(DIM_NUM, q_max, alpha,
            level_min, level_max, grid_num, point_num, ref grid_level, ref grid_point);
        //
        //  Now we're done.  Print the merged grid data.
        //
        Console.WriteLine("");
        Console.WriteLine("         I               Q          Grid levels:");
        Console.WriteLine("");
        for (j = 0; j < grid_num; j++)
        {
            double q = 0.0;
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                q += alpha[dim] * grid_level[dim + j * DIM_NUM];
            }

            string cout = "  " + (j + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                cout += "  " + grid_level[dim + j * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(8);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Grid points:");
        Console.WriteLine("");
        for (j = 0; j < point_num; j++)
        {
            string cout = "  " + (j + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8);
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                cout += "  " + grid_point[dim + j * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests CLENSHAW_CURTIS_COMPUTE
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int ORDER_MAX = 16;

        int order;
        double[] w = new double[ORDER_MAX];
        double[] x = new double[ORDER_MAX];

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  CLENSHAW_CURTIS_COMPUTE computes");
        Console.WriteLine("  a Clenshaw-Curtis quadrature rule over [-1,1]");
        Console.WriteLine("  of given order.");

        Console.WriteLine("");
        Console.WriteLine("    Order  W             X");
        Console.WriteLine("");

        for (order = 1; order <= 10; order++)
        {
            ClenshawCurtis.clenshaw_curtis_compute(order, ref x, ref w);

            Console.WriteLine("");
            Console.WriteLine("  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

            int i;
            for (i = 0; i < order; i++)
            {
                Console.WriteLine("          "
                                  + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    private static void test045()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST045 tests CC_ABSCISSA and CC_WEIGHT.
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
    {
        int i;
        const int order = 10;

        Console.WriteLine("");
        Console.WriteLine("TEST045");
        Console.WriteLine("  To compute a single Clenshaw Curtis weight or abscissa,");
        Console.WriteLine("  CC_ABSCISSA computes one abscissa,");
        Console.WriteLine("  CC_WEIGHT computes one weight.");
        Console.WriteLine("");
        Console.WriteLine("  We use these routines wastefully,");
        Console.WriteLine("  to compute the order 10 rule one value at a time.");
        Console.WriteLine("");
        Console.WriteLine("     Order       W               X");
        Console.WriteLine("");
        Console.WriteLine("  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

        for (i = 1; i <= order; i++)
        {
            double x = ClenshawCurtisGrid.cc_abscissa(order, i);
            double w = ClenshawCurtisGrid.cc_weight(order, i);

            Console.WriteLine("          "
                              + "  " + w.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests CLENSHAW_CURTIS_COMPUTE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 May 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int ORDER_MAX = 16;

        int j;
        int order;
        double[] result = new double[3];
        double[] weight = new double[ORDER_MAX];
        double[] xtab = new double[ORDER_MAX];

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  CLENSHAW_CURTIS_COMPUTE computes a Clenshaw-Curtis rule;");
        Console.WriteLine("");
        Console.WriteLine("  The integration interval is [-1,1].");
        Console.WriteLine("  Quadrature order will vary.");
        Console.WriteLine("  Integrand will vary.");
        Console.WriteLine("");
        Console.WriteLine("      Order     F1              F2              F3");
        Console.WriteLine("");

        for (order = 1; order <= ORDER_MAX; order++)
        {
            ClenshawCurtis.clenshaw_curtis_compute(order, ref xtab, ref weight);

            result[0] = 0.0;
            int i;
            for (i = 0; i < order; i++)
            {
                result[0] += weight[i] * f1(xtab[i]);
            }

            result[1] = 0.0;
            for (i = 0; i < order; i++)
            {
                result[1] += weight[i] * f2(xtab[i]);
            }

            result[2] = 0.0;
            for (i = 0; i < order; i++)
            {
                result[2] += weight[i] * f3(xtab[i]);
            }

            string cout = "  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(6);
            for (j = 0; j < 3; j++)
            {
                cout += "  " + result[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");

        result[0] = 46.0 * Math.Sinh(1.0) / 25.0 - 2.0 * Math.Sin(1.0);
        result[1] = 1.5822329637296729331;
        result[2] = (Math.Sqrt(2.0) + 3.0 * Math.Sqrt(6.0)) / 6.0;

        string cout2 = "  " + "Exact ";
        for (j = 0; j < 3; j++)
        {
            cout2 += "  " + result[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout2);
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests CLENSHAW_CURTIS_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 May 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int ORDER_MAX = 16;

        int j;
        int order;
        double[] result = new double[3];
        double[] weight = new double[ORDER_MAX];
        double[] xtab = new double[ORDER_MAX];

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  CLENSHAW_CURTIS_SET sets up a Clenshaw-Curtis rule;");
        Console.WriteLine("");
        Console.WriteLine("  The integration interval is [-1,1].");
        Console.WriteLine("  Quadrature order will vary.");
        Console.WriteLine("  Integrand will vary.");
        Console.WriteLine("");
        Console.WriteLine("  Order     F1              F2              F3");
        Console.WriteLine("");

        for (order = 1; order <= ORDER_MAX; order++)
        {
            ClenshawCurtisGrid.clenshaw_curtis_set(order, ref xtab, ref weight);

            result[0] = 0.0;
            int i;
            for (i = 0; i < order; i++)
            {
                result[0] += weight[i] * f1(xtab[i]);
            }

            result[1] = 0.0;
            for (i = 0; i < order; i++)
            {
                result[1] += weight[i] * f2(xtab[i]);
            }

            result[2] = 0.0;
            for (i = 0; i < order; i++)
            {
                result[2] += weight[i] * f3(xtab[i]);
            }

            string cout = "  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(6);
            for (j = 0; j < 3; j++)
            {
                cout += "  " + result[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);

        }

        Console.WriteLine("");

        result[0] = 46.0 * Math.Sinh(1.0) / 25.0 - 2.0 * Math.Sin(1.0);
        result[1] = 1.5822329637296729331;
        result[2] = (Math.Sqrt(2.0) + 3.0 * Math.Sqrt(6.0)) / 6.0;

        string cout2 = "  " + "Exact ";
        for (j = 0; j < 3; j++)
        {
            cout2 += "  " + result[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout2);
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests DTABLE_WRITE0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 1;
        const int ORDER = 9;

        const string r_file = "cc_r_d1_o9.txt";
        double[] r = new double[2];
        const string w_file = "cc_w_d1_o9.txt";
        double[] w = new double[ORDER];
        const string x_file = "cc_x_d1_o9.txt";
        double[] x = new double[ORDER];

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  DTABLE_WRITE0 writes a Clenshaw-Curtis");
        Console.WriteLine("  quadrature rule to a file.");

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension = " + DIM_NUM + "");
        Console.WriteLine("  Computing the rule of order = " + ORDER + "");

        ClenshawCurtis.clenshaw_curtis_compute(ORDER, ref x, ref w);

        Console.WriteLine("");
        Console.WriteLine("  Write abscissas to file \"" + x_file + "\".");

        typeMethods.dtable_write0(x_file, DIM_NUM, ORDER, x);

        Console.WriteLine("  Write weights to file \"" + w_file + "\".");

        typeMethods.dtable_write0(w_file, 1, ORDER, w);

        Console.WriteLine("  Write region to file \"" + r_file + "\".");

        r[0 + 0 * DIM_NUM] = -1.0;
        r[0 + 1 * DIM_NUM] = +1.0;

        typeMethods.dtable_write0(r_file, DIM_NUM, 2, r);
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests CLENSHAW_CURTIS_COMPUTE_ND.
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
    {
        int dim;
        int order;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  CLENSHAW_CURTIS_COMPUTE_ND computes");
        Console.WriteLine("  a multidimensional Clenshaw-Curtis quadrature rule");
        Console.WriteLine("  over the hypercube [-1,1]^ND of given");
        Console.WriteLine("  (possibly different) orders in each dimension.");

        const int dim_num = 2;

        int[] order_1d = new int[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            order_1d[dim] = 5;
        }

        int order_nd = typeMethods.i4vec_product(dim_num, order_1d);

        double[] point = new double[dim_num * order_nd];
        double[] weight = new double[order_nd];

        Console.WriteLine("");
        Console.WriteLine("  In this example, we use the SAME ORDER");
        Console.WriteLine("  in all dimensions.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        string cout = "  1D orders = ";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + order_1d[dim];
        }

        Console.WriteLine(cout);
        Console.WriteLine("  Number of points = " + order_nd + "");

        ClenshawCurtisGrid.clenshaw_curtis_compute_nd(dim_num, order_1d, ref point, ref weight);

        Console.WriteLine("");
        Console.WriteLine("      Weight           X(1)           X(2)");
        Console.WriteLine("  --------------  --------------  --------------");
        Console.WriteLine("");

        for (order = 0; order < order_nd; order++)
        {
            string cout2 = "  " + weight[order].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (dim = 0; dim < dim_num; dim++)
            {
                cout2 += "  " + point[dim + order * dim_num].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout2);
        }

        double weight_sum = typeMethods.r8vec_sum(order_nd, weight);

        Console.WriteLine("");
        Console.WriteLine("  " + weight_sum.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests CLENSHAW_CURTIS_COMPUTE_ND.
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
    {
        int dim;
        int order;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  CLENSHAW_CURTIS_COMPUTE_ND computes");
        Console.WriteLine("  a multidimensional Clenshaw-Curtis quadrature rule");
        Console.WriteLine("  over the hypercube [-1,1]^ND of given");
        Console.WriteLine("  (possibly different) orders in each dimension.");

        const int dim_num = 3;

        int[] order_1d = new int[dim_num];

        order_1d[0] = 2;
        order_1d[1] = 4;
        order_1d[2] = 3;

        int order_nd = typeMethods.i4vec_product(dim_num, order_1d);

        double[] point = new double[dim_num * order_nd];
        double[] weight = new double[order_nd];

        Console.WriteLine("");
        Console.WriteLine("  In this example, we use DIFFERENT ORDERS");
        Console.WriteLine("  in each dimension.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        string cout = "  1D orders = ";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + order_1d[dim];
        }

        Console.WriteLine(cout);
        Console.WriteLine("  Number of points = " + order_nd + "");

        ClenshawCurtisGrid.clenshaw_curtis_compute_nd(dim_num, order_1d, ref point, ref weight);

        Console.WriteLine("");
        Console.WriteLine("      Weight           X(1)           X(2)           X(3)");
        Console.WriteLine("  --------------  --------------  --------------  --------------");
        Console.WriteLine("");

        for (order = 0; order < order_nd; order++)
        {
            string cout2 = "  " + weight[order].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (dim = 0; dim < dim_num; dim++)
            {
                cout2 += "  " + point[dim + order * dim_num].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout2);
        }

        double weight_sum = typeMethods.r8vec_sum(order_nd, weight);

        Console.WriteLine("");
        Console.WriteLine("  " + weight_sum.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 uses DTABLE_WRITE0 to write out a multidimensional CC rule.
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
    {
        int dim;
        const string r_file = "cc_r_d4_o81.txt";
        const string w_file = "cc_w_d4_o81.txt";
        const string x_file = "cc_x_d4_o81.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  Use DTABLE_WRITE0 to write out a multidimensional rule.");
        Console.WriteLine("  CLENSHAW_CURTIS_COMPUTE_ND computes");
        Console.WriteLine("  a multidimensional Clenshaw-Curtis quadrature rule");
        Console.WriteLine("  over the hypercube [-1,1]^ND of given");
        Console.WriteLine("  (possibly different) orders in each dimension.");

        const int dim_num = 4;

        int[] order_1d = new int[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            order_1d[dim] = 3;
        }

        int order_nd = typeMethods.i4vec_product(dim_num, order_1d);

        double[] point = new double[dim_num * order_nd];
        double[] weight = new double[order_nd];
        double[] r = new double[dim_num * 2];

        Console.WriteLine("");
        Console.WriteLine("  In this example, we use the SAME ORDER");
        Console.WriteLine("  in all dimensions.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        string cout = "  1D orders = ";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + order_1d[dim];
        }

        Console.WriteLine(cout);
        Console.WriteLine("  Number of points = " + order_nd + "");

        ClenshawCurtisGrid.clenshaw_curtis_compute_nd(dim_num, order_1d, ref point, ref weight);

        for (dim = 0; dim < dim_num; dim++)
        {
            r[dim + 0 * dim_num] = -1.0;
            r[dim + 1 * dim_num] = +1.0;
        }

        Console.WriteLine("");
        Console.WriteLine("  Write abscissas to file \"" + x_file + "\".");

        typeMethods.dtable_write0(x_file, dim_num, order_nd, point);

        Console.WriteLine("  Write weights to file \"" + w_file + "\".");

        typeMethods.dtable_write0(w_file, 1, order_nd, weight);

        Console.WriteLine("  Write region to file \"" + r_file + "\".");

        typeMethods.dtable_write0(r_file, dim_num, 2, r);
    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests CC_ABSCISSA_LEVEL_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  CC_ABSCISSA_LEVEL_1D reports the level on which");
        Console.WriteLine("  a Clenshaw Curtis abscissa of given index will first");
        Console.WriteLine("  be generated, assuming a series of grids that grow");
        Console.WriteLine("  in order as 2^LEVEL+1.");

        const int base_ = 5;
        int order = (int)Math.Pow(2, base_) + 1;
        int test_num = (int)Math.Pow(2, base_) + 1;

        Console.WriteLine("");
        Console.WriteLine("  Base B = " + base_ + "");
        Console.WriteLine("  ORDER 2^B+1 = " + order + "");

        int[] test_val = new int[test_num];

        for (i = 0; i < test_num; i++)
        {
            test_val[i] = i;
        }

        int[] test_level = ClenshawCurtisGrid.cc_abscissa_level_1d(base_, test_num, test_val);

        Console.WriteLine("");
        Console.WriteLine("         I  Level(I)");
        Console.WriteLine("");

        for (i = 0; i < test_num; i++)
        {
            Console.WriteLine("  " + test_val[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + test_level[i].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }
    }

    private static void test12()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tests CC_ABSCISSA_LEVEL_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST12");
        Console.WriteLine("  CC_ABSCISSA_LEVEL_1D can also be called for values");
        Console.WriteLine("  outside the standard range of 0 through 2^LEVEL_MAX.");
        Console.WriteLine("  In that case, a MOD operation is applied first,");
        Console.WriteLine("  to make a sensible result.");

        const int base_ = 5;
        int order = (int)Math.Pow(2, base_) + 1;
        int seed = 123456789;
        const int test_num = 20;

        Console.WriteLine("");
        Console.WriteLine("  Base B = " + base_ + "");
        Console.WriteLine("  ORDER = 2^B+1 = " + order + "");

        int[] test_val = UniformRNG.i4vec_uniform(test_num, -20, 100, ref seed);

        int[] test_level = ClenshawCurtisGrid.cc_abscissa_level_1d(base_, test_num, test_val);

        Console.WriteLine("");
        Console.WriteLine("         I  Mod(I,O)  Level(I)");
        Console.WriteLine("");

        for (i = 0; i < test_num; i++)
        {
            Console.WriteLine("  " + test_val[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + (test_val[i] % order).ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + test_level[i].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }
    }

    private static void test13()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST13 tests CC_ABSCISSA_LEVEL_ND.
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
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST13");
        Console.WriteLine("  CC_ABSCISSA_LEVEL_ND reports the level on which");
        Console.WriteLine("  a Clenshaw Curtis abscissa of given index will first");
        Console.WriteLine("  be generated, assuming a series of grids that grow");
        Console.WriteLine("  in order as 2^LEVEL+1.");
        Console.WriteLine("");
        Console.WriteLine("  This routine is applied for multidimensional cases.");

        const int base_ = 3;
        int order = (int)Math.Pow(2, base_) + 1;
        const int dim_num = 2;
        int test_num = order * order;

        Console.WriteLine("");
        Console.WriteLine("  Base B = " + base_ + "");
        Console.WriteLine("  ORDER 2^B+1 = " + order + "");
        Console.WriteLine("  DIM_NUM = " + dim_num + "");

        int[] test_val = new int[dim_num * test_num];

        int k = 0;
        for (i = 0; i < order; i++)
        {
            int j;
            for (j = 0; j < order; j++)
            {
                test_val[0 + k * dim_num] = i;
                test_val[1 + k * dim_num] = j;
                k += 1;
            }
        }

        int[] test_level = ClenshawCurtisGrid.cc_abscissa_level_nd(base_, dim_num, test_num, test_val);

        Console.WriteLine("");
        Console.WriteLine("         I         J  Level(I,J)");
        Console.WriteLine("");

        for (k = 0; k < test_num; k++)
        {
            Console.WriteLine("  " + test_val[0 + k * dim_num].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + test_val[1 + k * dim_num].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + test_level[k].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }
    }

    private static double f1(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F1 evaluates F1(X) = 23 * Math.Cosh ( x ) / 25 - Math.Cos ( x ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Charles Clenshaw, Alan Curtis,
        //    A Method for Numerical Integration on an Automatic Computer,
        //    Numerische Mathematik,
        //    Volume 2, Number 1, December 1960, pages 197-205.
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double F1, the value of the function.
        //
    {
        double value = 23.0 * Math.Cosh(x) / 25.0 - Math.Cos(x);

        return value;
    }

    private static double f2(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F2 evaluates F2(X) = 1 / ( x^4 + x^2 + 0.9 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Charles Clenshaw, Alan Curtis,
        //    A Method for Numerical Integration on an Automatic Computer,
        //    Numerische Mathematik,
        //    Volume 2, Number 1, December 1960, pages 197-205.
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double F2, the value of the function.
        //
    {
        double value = 1.0 / (Math.Pow(x, 4) + Math.Pow(x, 2) + 0.9);

        return value;
    }

    private static double f3(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F3 evaluates F3(X) = Math.Sqrt ( abs ( x + 1/2 ) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Charles Clenshaw, Alan Curtis,
        //    A Method for Numerical Integration on an Automatic Computer,
        //    Numerische Mathematik,
        //    Volume 2, Number 1, December 1960, pages 197-205.
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double F3, the value of the function.
        //
    {
        double value = Math.Sqrt(Math.Abs(x + 0.5));

        return value;
    }
}