using System;
using Burkardt.Sparse;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SGMGAnisoNormalizeTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SGMGA_ANISO_NORMALIZE_TEST.
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
    {
        Console.WriteLine("");
        Console.WriteLine("SGMGA_ANISO_NORMALIZE_TEST:");
        Console.WriteLine("  Test the SGMGA_ANISO_NORMALIZE and");
        Console.WriteLine("  SGMGA_IMPORTANCE_TO_ANISO functions.");

        sgmga_aniso_balance_tests();

        sgmga_aniso_normalize_tests();

        sgmga_importance_to_aniso_tests();

        Console.WriteLine("");
        Console.WriteLine("SGMGA_ANISO_NORMALIZE_TEST:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }

    private static void sgmga_aniso_balance_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_ANISO_BALANCE_TESTS call SGMGA_ANISO_BALANCE_TEST.
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
    {
        double alpha_max;
        int dim;
        int dim_num;
        double[] level_weight;
        int seed;
        int test;
        int test_num;

        Console.WriteLine("");
        Console.WriteLine("SGMGA_ANISO_BALANCE_TESTS");
        Console.WriteLine("  Call SGMGA_ANISO_BALANCE_TEST with various arguments.");

        test_num = 5;
        dim_num = 5;
        level_weight = new double[dim_num];

        alpha_max = 10.0;
        seed = 123456789;
        for (test = 1; test <= test_num; test++)
        {
            UniformRNG.r8vec_uniform_01(dim_num, ref seed, ref level_weight);
            for (dim = 0; dim < dim_num; dim++)
            {
                level_weight[dim] = 10.0 * level_weight[dim];
            }

            sgmga_aniso_balance_test(alpha_max, dim_num, level_weight);
        }

        alpha_max = 5.0;
        seed = 123456789;
        for (test = 1; test <= test_num; test++)
        {
            UniformRNG.r8vec_uniform_01(dim_num, ref seed, ref level_weight);
            for (dim = 0; dim < dim_num; dim++)
            {
                level_weight[dim] = 10.0 * level_weight[dim];
            }

            sgmga_aniso_balance_test(alpha_max, dim_num, level_weight);
        }
    }

    private static void sgmga_aniso_balance_test(double alpha_max, int dim_num,
            double[] level_weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_ANISO_BALANCE_TEST calls SGMGA_ANISO_BALANCE.
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
    {
        int dim;
        double[] level_weight2;

        Console.WriteLine("");
        Console.WriteLine("SGMGA_ANISO_BALANCE_TEST");
        Console.WriteLine("  ALPHA_MAX = " + alpha_max + "");
        Console.WriteLine("  Input weight sum: "
                          + typeMethods.r8vec_sum(dim_num, level_weight) + "");
        string cout = "";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_weight[dim].ToString(CultureInfo.InvariantCulture).PadLeft(12);
        }

        Console.WriteLine(cout);

        level_weight2 = SGMGAniso.sgmga_aniso_balance(alpha_max, dim_num, level_weight);

        Console.WriteLine("  Output weight sum: "
                          + typeMethods.r8vec_sum(dim_num, level_weight2) + "");
        cout = "";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_weight2[dim].ToString(CultureInfo.InvariantCulture).PadLeft(12);
        }

        Console.WriteLine(cout);

    }

    private static void sgmga_aniso_normalize_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_ANISO_NORMALIZE_TESTS call SGMGA_ANISO_NORMALIZE_TEST.
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
    {
        int dim_num;
        double[] level_weight;

        Console.WriteLine("");
        Console.WriteLine("SGMGA_ANISO_NORMALIZE_TESTS");
        Console.WriteLine("  Call SGMGA_ANISO_NORMALIZE_TEST with various arguments.");

        dim_num = 2;
        level_weight = new double[dim_num];
        level_weight[0] = 1.0;
        level_weight[1] = 1.0;
        sgmga_aniso_normalize_test(dim_num, level_weight);

        dim_num = 2;
        level_weight = new double[dim_num];
        level_weight[0] = 10.0;
        level_weight[1] = 10.0;
        sgmga_aniso_normalize_test(dim_num, level_weight);

        dim_num = 2;
        level_weight = new double[dim_num];
        level_weight[0] = 10.0;
        level_weight[1] = 2.0;
        sgmga_aniso_normalize_test(dim_num, level_weight);

        dim_num = 2;
        level_weight = new double[dim_num];
        level_weight[0] = 1.0;
        level_weight[1] = 2.0;
        sgmga_aniso_normalize_test(dim_num, level_weight);

        dim_num = 3;
        level_weight = new double[dim_num];
        level_weight[0] = 1.0;
        level_weight[1] = 2.0;
        level_weight[2] = 3.0;
        sgmga_aniso_normalize_test(dim_num, level_weight);
        //
        //  Try a case in which one variable has 0 weight.
        //
        dim_num = 3;
        level_weight = new double[dim_num];
        level_weight[0] = 2.0;
        level_weight[1] = 0.0;
        level_weight[2] = 1.5;
        sgmga_aniso_normalize_test(dim_num, level_weight);

        dim_num = 4;
        level_weight = new double[dim_num];
        level_weight[0] = 1.0;
        level_weight[1] = 2.0;
        level_weight[2] = 3.0;
        level_weight[3] = 4.0;
        sgmga_aniso_normalize_test(dim_num, level_weight);

    }

    private static void sgmga_aniso_normalize_test(int dim_num, double[] level_weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_ANISO_NORMALIZE_TEST calls SGMGA_ANISO_NORMALIZE.
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
    {
        int dim;
        int option;

        Console.WriteLine("");
        Console.WriteLine("SGMGA_ANISO_NORMALIZE_TEST");
        Console.WriteLine("  Input weight sum: "
                          + typeMethods.r8vec_sum(dim_num, level_weight) + "");
        string cout = "";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_weight[dim].ToString(CultureInfo.InvariantCulture).PadLeft(12);
        }

        Console.WriteLine(cout);

        for (option = 0; option <= 2; option++)
        {
            SGMGAniso.sgmga_aniso_normalize(option, dim_num, ref level_weight);

            Console.WriteLine("  For OPTION = " + option
                                                + "  Normalized weight sum: "
                                                + typeMethods.r8vec_sum(dim_num, level_weight) + "");
            cout = "";
            for (dim = 0; dim < dim_num; dim++)
            {
                cout += "  " + level_weight[dim].ToString(CultureInfo.InvariantCulture).PadLeft(12);
            }

            Console.WriteLine(cout);
        }
    }

    private static void sgmga_importance_to_aniso_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_IMPORTANCE_TO_ANISO_TESTS call SGMGA_IMPORTANCE_TO_ANISO_TEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num;
        double[] importance;
        double[] level_weight;

        Console.WriteLine("");
        Console.WriteLine("SGMGA_IMPORTANCE_TO_ANISO_TESTS");
        Console.WriteLine("  Call SGMGA_IMPORTANCE_TO_ANISO_TEST with various arguments.");

        dim_num = 2;
        importance = new double[dim_num];
        level_weight = new double[dim_num];
        importance[0] = 1.0;
        importance[1] = 1.0;
        sgmga_importance_to_aniso_test(dim_num, importance, level_weight);

        dim_num = 2;
        importance = new double[dim_num];
        level_weight = new double[dim_num];
        importance[0] = 10.0;
        importance[1] = 10.0;
        sgmga_importance_to_aniso_test(dim_num, importance, level_weight);

        dim_num = 2;
        importance = new double[dim_num];
        level_weight = new double[dim_num];
        importance[0] = 10.0;
        importance[1] = 2.0;
        sgmga_importance_to_aniso_test(dim_num, importance, level_weight);

        dim_num = 2;
        importance = new double[dim_num];
        level_weight = new double[dim_num];
        importance[0] = 1.0;
        importance[1] = 2.0;
        sgmga_importance_to_aniso_test(dim_num, importance, level_weight);

        dim_num = 3;
        importance = new double[dim_num];
        level_weight = new double[dim_num];
        importance[0] = 1.0;
        importance[1] = 2.0;
        importance[2] = 3.0;
        sgmga_importance_to_aniso_test(dim_num, importance, level_weight);
        //
        //  Try a case in which one variable has 0 importance.
        //
        dim_num = 3;
        importance = new double[dim_num];
        level_weight = new double[dim_num];
        importance[0] = 2.0;
        importance[1] = 0.0;
        importance[2] = 1.5;
        sgmga_importance_to_aniso_test(dim_num, importance, level_weight);

        dim_num = 4;
        importance = new double[dim_num];
        level_weight = new double[dim_num];
        importance[0] = 1.0;
        importance[1] = 2.0;
        importance[2] = 3.0;
        importance[3] = 4.0;
        sgmga_importance_to_aniso_test(dim_num, importance, level_weight);
    }

    private static void sgmga_importance_to_aniso_test(int dim_num, double[] importance,
            double[] level_weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_IMPORTANCE_TO_ANISO_TEST calls SGMGA_IMPORTANCE_TO_ANISO.
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
    {
        int dim;

        Console.WriteLine("");
        Console.WriteLine("SGMGA_IMPORTANCE_TO_ANISO_TEST");
        Console.WriteLine("  Importances:");
        string cout = "";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + importance[dim].ToString(CultureInfo.InvariantCulture).PadLeft(12);
        }

        Console.WriteLine(cout);

        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);

        Console.WriteLine("  Anisotropic coefficients:");
        cout = "";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_weight[dim].ToString(CultureInfo.InvariantCulture).PadLeft(12);
        }

        Console.WriteLine(cout);
    }
}