using System;
using Burkardt.Sparse;
using Burkardt.Types;

namespace SGMGAnisoVCNCoefTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SGMGA_VCN_COEF_TEST.
        //
        //  Discussion:
        //
        //    SGMGA_VCN_COEF_TEST tests the SGMGA_VCN_COEF function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SGMGA_VCN_COEF_TEST");
        Console.WriteLine("  Test the SGMGA_VCN_COEF function.");
        //
        //  Compute examples of the combinatorial coefficent.
        //
        sgmga_vcn_coef_tests();
        //
        //  That's all.
        //
        Console.WriteLine("");
        Console.WriteLine("SGMGA_VCN_COEF_TEST");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }

    private static void sgmga_vcn_coef_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_VCN_COEF_TESTS calls SGMGA_VCN_COEF_TEST.
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
        int dim_num;
        double[] importance;
        int level_max_max;
        int level_max_min;
        double[] level_weight;

        Console.WriteLine("");
        Console.WriteLine("SGMGA_VCN_COEF_TESTS");
        Console.WriteLine("  calls SGMGA_VCN_COEF_TEST.");

        dim_num = 2;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = 1.0;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 4;
        sgmga_vcn_coef_test(dim_num, importance, level_weight, level_max_min,
            level_max_max);

        dim_num = 2;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = dim + 1;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 4;
        sgmga_vcn_coef_test(dim_num, importance, level_weight, level_max_min,
            level_max_max);

        dim_num = 3;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = 1.0;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 4;
        sgmga_vcn_coef_test(dim_num, importance, level_weight, level_max_min,
            level_max_max);

        dim_num = 3;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = dim + 1;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 4;
        sgmga_vcn_coef_test(dim_num, importance, level_weight, level_max_min,
            level_max_max);

        dim_num = 4;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = dim + 1;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 3;
        sgmga_vcn_coef_test(dim_num, importance, level_weight, level_max_min,
            level_max_max);
        //
        //  Try a case with a dimension of "0 importance".
        //
        dim_num = 3;
        importance = new double[dim_num];
        importance[0] = 1.0;
        importance[1] = 0.0;
        importance[2] = 1.0;
        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 3;
        sgmga_vcn_coef_test(dim_num, importance, level_weight, level_max_min,
            level_max_max);

    }

    private static void sgmga_vcn_coef_test(int dim_num, double[] importance,
            double[] level_weight, int level_max_min, int level_max_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_VCN_COEF_TEST tests SGMGA_VCN_COEF.
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
    {
        double coef1;
        double coef1_sum;
        double coef2;
        double coef2_sum;
        int dim;
        int i;
        int[] level_1d;
        int[] level_1d_max;
        int[] level_1d_min;
        int level_max;
        double level_weight_min_pos;
        bool more_grids;
        double q;
        double q_max;
        double q_min;
        SGMGAniso.SGMGAData data = new();
        string cout = "";

        level_1d = new int[dim_num];
        level_1d_max = new int[dim_num];
        level_1d_min = new int[dim_num];

        Console.WriteLine("");
        Console.WriteLine("SGMGA_VCN_COEF_TEST");
        Console.WriteLine("  For anisotropic problems, a \"combinatorial coefficent\"");
        Console.WriteLine("  must be computed for each component product grid.");
        Console.WriteLine("  SGMGA_VCN_COEF_NAIVE does this in a simple, inefficient way.");
        Console.WriteLine("  SGMGA_VCN_COEF tries to be more efficient.");
        Console.WriteLine("  Here, we simply compare COEF1 and COEF2, the same");
        Console.WriteLine("  coefficient computed by the naive and efficient ways.");
        Console.WriteLine("");
        Console.WriteLine("  IMPORTANCE:");
        cout = "";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + importance[dim].ToString().PadLeft(14);
        }

        Console.WriteLine(cout);
        Console.WriteLine("  LEVEL_WEIGHT:");
        cout = "";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_weight[dim].ToString().PadLeft(14);
        }

        Console.WriteLine(cout);

        for (level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            i = 0;
            coef1_sum = 0.0;
            coef2_sum = 0.0;
            //
            //  Initialization.
            //
            level_weight_min_pos = typeMethods.r8vec_min_pos(dim_num, level_weight);
            q_min = level_max * level_weight_min_pos
                    - typeMethods.r8vec_sum(dim_num, level_weight);
            q_max = level_max * level_weight_min_pos;
            for (dim = 0; dim < dim_num; dim++)
            {
                level_1d_min[dim] = 0;
            }

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

            more_grids = false;


            Console.WriteLine("");
            Console.WriteLine("     I               Q       Coef1       Coef2   X");
            cout = "   MIN" + "  " + q_min.ToString().PadLeft(14)
                   + "                        ";
            for (dim = 0; dim < dim_num; dim++)
            {
                cout += "  " + level_1d_min[dim].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
            //
            //  Seek all vectors LEVEL_1D which satisfy the constraint:
            //
            //    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
            //      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
            //      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
            //
            for (;;)
            {
                SGMGAniso.sgmga_vcn_ordered_naive(ref data, dim_num, level_weight, level_1d_max,
                    level_1d, q_min, q_max, ref more_grids);

                if (!more_grids)
                {
                    break;
                }

                //
                //  Compute the combinatorial coefficient.
                //
                coef1 = SGMGAniso.sgmga_vcn_coef_naive(dim_num, level_weight, level_1d_max,
                    level_1d, q_min, q_max);

                coef2 = SGMGAniso.sgmga_vcn_coef(dim_num, level_weight, level_1d, q_max);

                i += 1;

                q = 0.0;
                for (dim = 0; dim < dim_num; dim++)
                {
                    q += level_weight[dim] * level_1d[dim];
                }

                coef1_sum += coef1;
                coef2_sum += coef2;

                cout = "  " + i.ToString().PadLeft(4)
                            + "  " + q.ToString().PadLeft(14)
                            + "  " + coef1.ToString().PadLeft(10)
                            + "  " + coef2.ToString().PadLeft(10);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + level_1d[dim].ToString().PadLeft(2);
                }

                Console.WriteLine(cout);
            }

            cout = "   MAX" + "  " + q_max.ToString().PadLeft(14)
                   + "                        ";
            for (dim = 0; dim < dim_num; dim++)
            {
                cout += "  " + level_1d_max[dim].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
            Console.WriteLine("   SUM                "
                              + "  " + coef1_sum.ToString().PadLeft(10)
                              + "  " + coef2_sum.ToString().PadLeft(10) + "");
        }
    }
}