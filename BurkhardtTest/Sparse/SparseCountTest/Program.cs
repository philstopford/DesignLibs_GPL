using System;
using Burkardt.Sparse;

namespace SparseCountTest;

internal static class Program
{
    private static void Main()
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPARSE_COUNT_TEST.
//
//  Discussion:
//
//    SPARSE_COUNT_TEST tests the SPARSE_COUNT library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2012
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
        Console.WriteLine("");
        Console.WriteLine("SPARSE_COUNT_TEST");
        Console.WriteLine("  Test the SPARSE_COUNT library.");
//
//  CC_SE
//
        int dim_min = 1;
        int dim_max = 5;
        int level_max_min = 0;
        int level_max_max = 10;
        test01(dim_min, dim_max, level_max_min, level_max_max);

        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 10;
        test01(dim_min, dim_max, level_max_min, level_max_max);
//
//  CFN_E
//
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 10;
        test02(dim_min, dim_max, level_max_min, level_max_max);

        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 10;
        test02(dim_min, dim_max, level_max_min, level_max_max);
//
//  F2_SE
//
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 10;
        test03(dim_min, dim_max, level_max_min, level_max_max);

        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 10;
        test03(dim_min, dim_max, level_max_min, level_max_max);
//
//  GP_SE
//
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 10;
        test04(dim_min, dim_max, level_max_min, level_max_max);

        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 10;
        test04(dim_min, dim_max, level_max_min, level_max_max);
//
//  OFN_E
//
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 10;
        test05(dim_min, dim_max, level_max_min, level_max_max);

        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 10;
        test05(dim_min, dim_max, level_max_min, level_max_max);
//
//  ONN_E
//
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 10;
        test06(dim_min, dim_max, level_max_min, level_max_max);

        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 10;
        test06(dim_min, dim_max, level_max_min, level_max_max);
//
//  ONN_L
//
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 10;
        test07(dim_min, dim_max, level_max_min, level_max_max);

        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 10;
        test07(dim_min, dim_max, level_max_min, level_max_max);
//
//  OWN_E
//
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 10;
        test08(dim_min, dim_max, level_max_min, level_max_max);

        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 10;
        test08(dim_min, dim_max, level_max_min, level_max_max);
//
//  OWN_L2
//
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 10;
        test09(dim_min, dim_max, level_max_min, level_max_max);

        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 10;
        test09(dim_min, dim_max, level_max_min, level_max_max);
//
//  OWN_LS
//
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 10;
        test10(dim_min, dim_max, level_max_min, level_max_max);

        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 10;
        test10(dim_min, dim_max, level_max_min, level_max_max);

        Console.WriteLine("");
        Console.WriteLine("SPARSE_COUNT_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int dim_min, int dim_max, int level_max_min, int level_max_max)
//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests CC_SE_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_MIN, the minimum spatial dimension to consider.
//
//    Input, int DIM_MAX, the maximum spatial dimension to consider.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
//
    {
        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  CC_SE_SIZE returns the number of ");
        Console.WriteLine("  distinct points in a sparse grid made from ");
        Console.WriteLine("  * CC_SE, Clenshaw Curtis Slow Exponential Growth family.");
        Console.WriteLine("");
        string cout = "   DIM: ";
        for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
        {
            cout += "  " + dim_num.ToString().PadLeft(10);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("   LEVEL_MAX");
        Console.WriteLine("");

        for (int level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            cout = "    " + level_max.ToString().PadLeft(4);
            for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
            {
                int point_num = SparseCount.cc_se_size(dim_num, level_max);
                cout += "  " + point_num.ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test02(int dim_min, int dim_max, int level_max_min, int level_max_max)
//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests CFN_E_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_MIN, the minimum spatial dimension to consider.
//
//    Input, int DIM_MAX, the maximum spatial dimension to consider.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
//
    {
        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  CFN_E_SIZE returns the number of ");
        Console.WriteLine("  distinct points in a CFN_E sparse grid made from ");
        Console.WriteLine("  any closed fully nested family of 1D quadrature");
        Console.WriteLine("  rules with exponential growth, including:");
        Console.WriteLine("  * CC_E, the Clenshaw Curtis Exponential Growth family;");
        Console.WriteLine("  * NCC_E, the Newton Cotes Closed Exponential Growth family.");
        Console.WriteLine("");
        string cout = "   DIM: ";
        for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
        {
            cout += "  " + dim_num.ToString().PadLeft(10);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("   LEVEL_MAX");
        Console.WriteLine("");

        for (int level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            cout = "    " + level_max.ToString().PadLeft(4);
            for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
            {
                int point_num = SparseCount.cfn_e_size(dim_num, level_max);
                cout += "  " + point_num.ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test03(int dim_min, int dim_max, int level_max_min, int level_max_max)
//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests F2_SE_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_MIN, the minimum spatial dimension to consider.
//
//    Input, int DIM_MAX, the maximum spatial dimension to consider.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
//
    {
        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  F2_SE_SIZE returns the number of ");
        Console.WriteLine("  distinct points in a sparse grid made from ");
        Console.WriteLine("  * F2_SE, Fejer Type 2 Slow Exponential Growth family.");
        Console.WriteLine("");
        string cout = "   DIM: ";
        for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
        {
            cout += "  " + dim_num.ToString().PadLeft(10);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("   LEVEL_MAX");
        Console.WriteLine("");

        for (int level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            cout = "    " + level_max.ToString().PadLeft(4);
            for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
            {
                int point_num = SparseCount.f2_se_size(dim_num, level_max);
                cout += "  " + point_num.ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test04(int dim_min, int dim_max, int level_max_min, int level_max_max)
//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests GP_SE_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_MIN, the minimum spatial dimension to consider.
//
//    Input, int DIM_MAX, the maximum spatial dimension to consider.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
//
    {
        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  GP_SE_SIZE returns the number of ");
        Console.WriteLine("  distinct points in a sparse grid made from ");
        Console.WriteLine("  * GP_SE, Gauss Patterson Slow Exponential Growth family.");
        Console.WriteLine("");
        string cout = "   DIM: ";
        for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
        {
            cout += "  " + dim_num.ToString().PadLeft(10);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("   LEVEL_MAX");
        Console.WriteLine("");

        for (int level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            cout = "    " + level_max.ToString().PadLeft(4);
            for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
            {
                int point_num = SparseCount.gp_se_size(dim_num, level_max);
                cout += "  " + point_num.ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test05(int dim_min, int dim_max, int level_max_min, int level_max_max)
//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests OFN_E_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_MIN, the minimum spatial dimension to consider.
//
//    Input, int DIM_MAX, the maximum spatial dimension to consider.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
//
    {
        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  OFN_E_SIZE returns the number of ");
        Console.WriteLine("  distinct points in an OFN_E sparse grid made from ");
        Console.WriteLine("  product grids formed from open fully nested ");
        Console.WriteLine("  quadrature rules with Exponential Growth, including");
        Console.WriteLine("  * F2_E, the Fejer Type 2 Exponential Growth Family;");
        Console.WriteLine("  * GP_E, the Gauss Patterson Exponential Growth Family;");
        Console.WriteLine("  * NCO_E, the Newton Cotes Open Exponential Growth Family.");
        Console.WriteLine("");
        string cout = "   DIM: ";
        for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
        {
            cout += "  " + dim_num.ToString().PadLeft(10);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("   LEVEL_MAX");
        Console.WriteLine("");

        for (int level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            cout = "    " + level_max.ToString().PadLeft(4);
            for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
            {
                int point_num = SparseCount.ofn_e_size(dim_num, level_max);
                cout += "  " + point_num.ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test06(int dim_min, int dim_max, int level_max_min, int level_max_max)
//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests ONN_E_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_MIN, the minimum spatial dimension to consider.
//
//    Input, int DIM_MAX, the maximum spatial dimension to consider.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
//
    {
        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  ONN_E_SIZE returns the number of ");
        Console.WriteLine("  distinct points in an ONN_E sparse grid made from ");
        Console.WriteLine("  product grids formed from open non-nested ");
        Console.WriteLine("  quadrature rules with exponential growth, including:");
        Console.WriteLine("  * LG_E, the Gauss Laguerre Exponential Growth Family;");
        Console.WriteLine("  * GJ_E, the Gauss Jacobi Exponential Growth Family;");
        Console.WriteLine("  * GLG_E, the Generalized Gauss Laguerre Exponential Growth Family");
        Console.WriteLine("  * GW_E, any Golub Welsch Exponential Growth Family;");
        Console.WriteLine("");
        string cout = "   DIM: ";
        for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
        {
            cout += "  " + dim_num.ToString().PadLeft(10);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("   LEVEL_MAX");
        Console.WriteLine("");

        for (int level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            cout = "    " + level_max.ToString().PadLeft(4);
            for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
            {
                int point_num = SparseCount.onn_e_size(dim_num, level_max);
                cout += "  " + point_num.ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test07(int dim_min, int dim_max, int level_max_min, int level_max_max)
//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests ONN_L_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_MIN, the minimum spatial dimension to consider.
//
//    Input, int DIM_MAX, the maximum spatial dimension to consider.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
//
    {
        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  ONN_L_SIZE returns the number of ");
        Console.WriteLine("  distinct points in an ONN_L sparse grid made from ");
        Console.WriteLine("  product grids formed from open non-nested ");
        Console.WriteLine("  quadrature rules with linear growth, including:");
        Console.WriteLine("  * LG_L, the Gauss Laguerre Linear Growth Family;");
        Console.WriteLine("  * GJ_L, the Gauss Jacobi Linear Growth Family;");
        Console.WriteLine("  * GLG_L, the Generalized Gauss Laguerre Linear Growth Family;");
        Console.WriteLine("  * GW_L, any Golub Welsch Linear Growth Family;");
        Console.WriteLine("");
        string cout = "   DIM: ";
        for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
        {
            cout += "  " + dim_num.ToString().PadLeft(10);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("   LEVEL_MAX");
        Console.WriteLine("");

        for (int level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            cout = "    " + level_max.ToString().PadLeft(4);
            for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
            {
                int point_num = SparseCount.onn_l_size(dim_num, level_max);
                cout += "  " + point_num.ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test08(int dim_min, int dim_max, int level_max_min, int level_max_max)
//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests OWN_E_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_MIN, the minimum spatial dimension to consider.
//
//    Input, int DIM_MAX, the maximum spatial dimension to consider.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
//
    {
        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  OWN_E_SIZE returns the number of ");
        Console.WriteLine("  distinct points in an OWN_E sparse grid made from ");
        Console.WriteLine("  product grids formed from open weakly nested ");
        Console.WriteLine("  quadrature rules with exponential growth, including:");
        Console.WriteLine("  * GGH_E, the Generalized Gauss-Hermite Exponential Growth Family;");
        Console.WriteLine("  * GH_E, the Gauss-Hermite Exponential Growth Family;");
        Console.WriteLine("  * LG_E, the Gauss-Legendre Exponential Growth Family;");
        Console.WriteLine("");
        string cout = "   DIM: ";
        for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
        {
            cout += "  " + dim_num.ToString().PadLeft(10);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("   LEVEL_MAX");
        Console.WriteLine("");

        for (int level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            cout = "    " + level_max.ToString().PadLeft(4);
            for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
            {
                int point_num = SparseCount.own_e_size(dim_num, level_max);
                cout += "  " + point_num.ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test09(int dim_min, int dim_max, int level_max_min, int level_max_max)
//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests OWN_L2_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_MIN, the minimum spatial dimension to consider.
//
//    Input, int DIM_MAX, the maximum spatial dimension to consider.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
//
    {
        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  OWN_L2_SIZE returns the number of ");
        Console.WriteLine("  distinct points in an OWN_L2 sparse grid made from ");
        Console.WriteLine("  product grids formed from open weakly nested ");
        Console.WriteLine("  quadrature rules with linear growth, including:");
        Console.WriteLine("  * GGH_L2, the Generalized Gauss-Hermite Linear Growth Family;");
        Console.WriteLine("  * GH_L2, the Gauss-Hermite Linear Growth Family;");
        Console.WriteLine("  * LG_L2, the Gauss-Legendre Linear Growth Family;");
        Console.WriteLine("");
        string cout = "   DIM: ";
        for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
        {
            cout += "  " + dim_num.ToString().PadLeft(10);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("   LEVEL_MAX");
        Console.WriteLine("");

        for (int level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            cout =  "    " + level_max.ToString().PadLeft(4);
            for (int dim_num = dim_min; dim_num <= dim_max; dim_num++)
            {
                int point_num = SparseCount.own_l2_size(dim_num, level_max);
                cout += "  " + point_num.ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test10(int dim_min, int dim_max, int level_max_min, int level_max_max)
//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests OWN_O_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_MIN, the minimum spatial dimension to consider.
//
//    Input, int DIM_MAX, the maximum spatial dimension to consider.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
//
    {
        int dim_num;
        int level_max;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  OWN_O_SIZE returns the number of ");
        Console.WriteLine("  distinct points in an OWN_O sparse grid made from ");
        Console.WriteLine("  product grids formed from open weakly nested ");
        Console.WriteLine("  quadrature rules with odd growth, including:");
        Console.WriteLine("  * GGH_L, the Generalized Gauss-Hermite Odd Growth Family;");
        Console.WriteLine("  * GH_L, the Gauss-Hermite Odd Growth Family;");
        Console.WriteLine("  * LG_L, the Gauss-Legendre Odd Growth Family;");
        Console.WriteLine("");
        string cout = "   DIM: ";
        for (dim_num = dim_min; dim_num <= dim_max; dim_num++)
        {
            cout += "  " + dim_num.ToString().PadLeft(10);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("   LEVEL_MAX");
        Console.WriteLine("");

        for (level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            cout = "    " + level_max.ToString().PadLeft(4);
            for (dim_num = dim_min; dim_num <= dim_max; dim_num++)
            {
                int point_num = SparseCount.own_o_size(dim_num, level_max);
                cout += "  " + point_num.ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }
    }
}