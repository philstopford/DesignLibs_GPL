using System;
using System.Globalization;
using Burkardt.Sparse;
using Burkardt.Types;

namespace SGMGAnisoVCNTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SGMGA_VCN_TEST.
        //
        //  Discussion:
        //
        //    SGMGA_VCN_TEST tests SGMGA_VCN.
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
    {
        Console.WriteLine("");
        Console.WriteLine("SGMGA_VCN_TEST");
        Console.WriteLine("  Test the SGMGA_VCN and SGMGA_VCN_ORDERED functions.");

        sgmga_vcn_tests();

        sgmga_vcn_timing_tests();

        sgmga_vcn_ordered_tests();

        Console.WriteLine("");
        Console.WriteLine("SGMGA_VCN_TEST");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }

    private static void sgmga_vcn_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_VCN_TESTS calls SGMGA_VCN_TEST.
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
    {
        int dim;
        int dim_num;
        int[] dim_num_array =
        {
            2, 2, 2, 2, 2,
            3, 3, 3, 3, 3,
            4, 4
        };
        double[] importance;
        int level_max;
        int[] level_max_array =
        {
            0, 1, 2, 3, 4,
            0, 1, 2, 3, 4,
            2, 3
        };
        double[] level_weight;
        double q_max;
        double q_min;
        int test;
        const int test_num = 12;

        Console.WriteLine("");
        Console.WriteLine("SGMGA_VCN_TESTS");
        Console.WriteLine("  calls SGMGA_VCN_TEST.");
        //
        //  Isotropic examples.
        //
        for (test = 0; test < test_num; test++)
        {
            dim_num = dim_num_array[test];
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = 1.0;
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = level_max_array[test];
            q_min = level_max - typeMethods.r8vec_sum(dim_num, level_weight);
            q_max = level_max;

            sgmga_vcn_test(dim_num, importance, level_weight, q_min, q_max);

        }

        //
        //  Zero weight example.
        //
        dim_num = 3;
        importance = new double[dim_num];
        importance[0] = 1.0;
        importance[1] = 0.0;
        importance[2] = 1.0;
        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max = 2;
        q_min = level_max - typeMethods.r8vec_sum(dim_num, level_weight);
        q_max = level_max;

        sgmga_vcn_test(dim_num, importance, level_weight, q_min, q_max);

        //
        //  Anisotropic examples.
        //
        for (test = 0; test < test_num; test++)
        {
            dim_num = dim_num_array[test];
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = dim + 1;
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = level_max_array[test];
            q_min = level_max - typeMethods.r8vec_sum(dim_num, level_weight);
            q_max = level_max;

            sgmga_vcn_test(dim_num, importance, level_weight, q_min, q_max);

        }

    }

    private static void sgmga_vcn_test(int dim_num, double[] importance, double[] level_weight,
            double q_min, double q_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_VCN_TEST tests SGMGA_VCN.
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
    {
        int dim;
        double q;
        SGMGAniso.SGMGAData data = new();
        string cout;

        int[] level_1d = new int[dim_num];
        int[] level_1d_max = new int[dim_num];
        int[] level_1d_min = new int[dim_num];

        Console.WriteLine("");
        Console.WriteLine("SGMGA_VCN_TEST");
        Console.WriteLine("  Consider vectors 0 <= LEVEL_1D(1:N) <= LEVEL_1D_MAX(1:N),");
        Console.WriteLine("  Set Q = sum ( LEVEL_WEIGHT(1:N) * LEVEL_1D(1:N) )");
        Console.WriteLine("  Accept vectors for which Q_MIN < Q <= Q_MAX");
        Console.WriteLine("  No particular order is imposed on the LEVEL_1D values.");
        Console.WriteLine("  SGMGA_VCN_NAIVE uses a naive approach;");
        Console.WriteLine("  SGMGA_VCN tries to be more efficient.");
        Console.WriteLine("  Here, we just compare the results.");

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

        for (dim = 0; dim < dim_num; dim++)
        {
            level_1d_min[dim] = 0;
        }

        bool more_grids = false;

        Console.WriteLine("");
        Console.WriteLine("  IMPORTANCE:");
        cout = "";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + importance[dim].ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);
        Console.WriteLine("  LEVEL_WEIGHT:");
        cout = "";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_weight[dim].ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("  SGMGA_VCN_NAIVE");
        Console.WriteLine("     I               Q   X");
        cout = "   MIN" + "  " + q_min.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_1d_min[dim].ToString().PadLeft(2);
        }

        Console.WriteLine(cout);

        int i = 0;

        for (;;)
        {
            SGMGAniso.sgmga_vcn_naive(dim_num, level_weight, level_1d_max, level_1d,
                q_min, q_max, ref more_grids);

            if (!more_grids)
            {
                break;
            }

            q = 0.0;
            for (dim = 0; dim < dim_num; dim++)
            {
                q += level_weight[dim] * level_1d[dim];
            }

            i += 1;
            cout = "  " + i.ToString().PadLeft(4)
                        + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (dim = 0; dim < dim_num; dim++)
            {
                cout += "  " + level_1d[dim].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
        }

        cout = "   MAX" + "  " + q_max.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_1d_max[dim].ToString().PadLeft(2);
        }

        Console.WriteLine(cout);

        Console.WriteLine("");
        Console.WriteLine("  SGMGA_VCN");
        Console.WriteLine("     I               Q   X");
        cout = "   MIN" + "  " + q_min.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_1d_min[dim].ToString().PadLeft(2);
        }

        Console.WriteLine(cout);

        i = 0;

        for (;;)
        {
            SGMGAniso.sgmga_vcn(ref data, dim_num, level_weight, ref level_1d, q_min, q_max,
                ref more_grids);

            if (!more_grids)
            {
                break;
            }

            q = 0.0;
            for (dim = 0; dim < dim_num; dim++)
            {
                q += level_weight[dim] * level_1d[dim];
            }

            i += 1;
            cout = "  " + i.ToString().PadLeft(4)
                        + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (dim = 0; dim < dim_num; dim++)
            {
                cout += "  " + level_1d[dim].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
        }

        cout = "   MAX" + "  " + q_max.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_1d_max[dim].ToString().PadLeft(2);
        }

        Console.WriteLine(cout);
    }

    private static void sgmga_vcn_timing_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_VCN_TIMING_TESTS calls SGMGA_VCN_TIMING_TEST.
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
    {
        int dim;
        double[] importance;
        int level_max;
        double[] level_weight;
        double q_max;
        double q_min;
        int test;

        Console.WriteLine("");
        Console.WriteLine("SGMGA_VCN_TIMING_TESTS");
        Console.WriteLine("  calls SGMGA_VCN_TIMING_TEST.");
        //
        //  Isotropic examples.
        //
        int dim_num = 2;

        for (test = 0; test < 2; test++)
        {
            dim_num *= 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = 1.0;
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 2;
            q_min = level_max - typeMethods.r8vec_sum(dim_num, level_weight);
            q_max = level_max;

            sgmga_vcn_timing_test(dim_num, importance, level_weight, q_min, q_max);

        }

        //
        //  Anisotropic examples.
        //
        dim_num = 2;

        for (test = 0; test < 2; test++)
        {
            dim_num *= 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = dim + 1;
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 2;
            q_min = level_max - typeMethods.r8vec_sum(dim_num, level_weight);
            q_max = level_max;

            sgmga_vcn_timing_test(dim_num, importance, level_weight, q_min, q_max);

        }

    }

    private static void sgmga_vcn_timing_test(int dim_num, double[] importance,
            double[] level_weight, double q_min, double q_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_VCN_TIMING_TEST times SGMGA_VCN and SGMGA_VCN_NAIVE.
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
    {
        int dim;
        SGMGAniso.SGMGAData data = new();
        string cout;

        int[] level_1d = new int[dim_num];
        int[] level_1d_max = new int[dim_num];
        int[] level_1d_min = new int[dim_num];

        Console.WriteLine("");
        Console.WriteLine("SGMGA_VCN_TIMING_TEST");
        Console.WriteLine("  Consider vectors 0 <= LEVEL_1D(1:N) <= LEVEL_1D_MAX(1:N),");
        Console.WriteLine("  Set Q = sum ( LEVEL_WEIGHT(1:N) * LEVEL_1D(1:N) )");
        Console.WriteLine("  Accept vectors for which Q_MIN < Q <= Q_MAX");
        Console.WriteLine("  No particular order is imposed on the LEVEL_1D values.");
        Console.WriteLine("  SGMGA_VCN_NAIVE uses a naive approach;");
        Console.WriteLine("  SGMGA_VCN tries to be more efficient.");
        Console.WriteLine("  Here, we compare the timings.");

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

        for (dim = 0; dim < dim_num; dim++)
        {
            level_1d_min[dim] = 0;
        }

        bool more_grids = false;

        Console.WriteLine("");
        Console.WriteLine("  IMPORTANCE:");
        cout = "";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + importance[dim].ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);
        Console.WriteLine("  LEVEL_WEIGHT:");
        cout = "";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_weight[dim].ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("  SGMGA_VCN_NAIVE");
        Console.WriteLine("     I               Q   X");
        cout = "   MIN" + "  " + q_min.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_1d_min[dim].ToString().PadLeft(2);
        }

        Console.WriteLine(cout);

        DateTime t1 = DateTime.Now;
        for (;;)
        {
            SGMGAniso.sgmga_vcn_naive(dim_num, level_weight, level_1d_max, level_1d,
                q_min, q_max, ref more_grids);

            if (!more_grids)
            {
                break;
            }
        }

        DateTime t2 = DateTime.Now;
        cout = "   MAX" + "  " + q_max.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_1d_max[dim].ToString().PadLeft(2);
        }

        Console.WriteLine(cout);
        Console.WriteLine("  TIME" + "  "
                                   + (t2 - t1).TotalSeconds.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        Console.WriteLine("");
        Console.WriteLine("  SGMGA_VCN");
        Console.WriteLine("     I               Q   X");
        cout = "   MIN" + "  " + q_min.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_1d_min[dim].ToString().PadLeft(2);
        }

        Console.WriteLine(cout);

        t1 = DateTime.Now;
        for (;;)
        {
            SGMGAniso.sgmga_vcn(ref data, dim_num, level_weight, ref level_1d, q_min, q_max,
                ref more_grids);

            if (!more_grids)
            {
                break;
            }
        }

        t2 = DateTime.Now;
        cout = "   MAX" + "  " + q_max.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_1d_max[dim].ToString().PadLeft(2);
        }

        Console.WriteLine(cout);
        Console.WriteLine("  TIME" + "  "
                                   + (t2 - t1).TotalSeconds.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void sgmga_vcn_ordered_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_VCN_ORDERED_TESTS calls SGMGA_VCN_ORDERED_TEST.
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
    {
        int dim;
        int dim_num;
        int[] dim_num_array =
        {
            2, 2, 2, 2, 2,
            3, 3, 3, 3, 3,
            4, 4
        };
        double[] importance;
        int level_max;
        int[] level_max_array =
        {
            0, 1, 2, 3, 4,
            0, 1, 2, 3, 4,
            2, 3
        };
        double[] level_weight;
        double q_max;
        double q_min;
        int test;
        const int test_num = 12;

        Console.WriteLine("");
        Console.WriteLine("SGMGA_VCN_ORDERED_TESTS");
        Console.WriteLine("  calls SGMGA_VCN_ORDERED_TEST.");
        //
        //  Isotropic examples.
        //
        for (test = 0; test < test_num; test++)
        {
            dim_num = dim_num_array[test];
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = 1.0;
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = level_max_array[test];
            q_min = level_max - typeMethods.r8vec_sum(dim_num, level_weight);
            q_max = level_max;

            sgmga_vcn_ordered_test(dim_num, importance, level_weight, q_min, q_max);

        }

        //
        //  Anisotropic examples.
        //
        for (test = 0; test < test_num; test++)
        {
            dim_num = dim_num_array[test];
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = dim + 1;
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = level_max_array[test];
            q_min = level_max - typeMethods.r8vec_sum(dim_num, level_weight);
            q_max = level_max;

            sgmga_vcn_ordered_test(dim_num, importance, level_weight, q_min, q_max);

        }

    }

    private static void sgmga_vcn_ordered_test(int dim_num, double[] importance,
            double[] level_weight, double q_min, double q_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_VCN_ORDERED_TEST tests SGMGA_VCN_ORDERED and SGMGA_VCN_ORDERED_NAIVE.
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
    {
        int dim;
        double q;
        SGMGAniso.SGMGAData data = new();
        string cout;

        int[] level_1d = new int[dim_num];
        int[] level_1d_max = new int[dim_num];
        int[] level_1d_min = new int[dim_num];

        Console.WriteLine("");
        Console.WriteLine("SGMGA_VCN_ORDERED_TEST");

        Console.WriteLine("  Consider vectors 0 <= LEVEL_1D(1:N) <= LEVEL_1D_MAX(1:N),");
        Console.WriteLine("  Set Q = sum ( LEVEL_WEIGHT(1:N) * LEVEL_1D(1:N) )");
        Console.WriteLine("  Accept only vectors for which Q_MIN < Q <= Q_MAX");
        Console.WriteLine("  The solutions are weakly ordered by the value of Q.");
        Console.WriteLine("  SGMGA_VCN_ORDERED_NAIVE calls SGMGA_VCN_NAIVE;");
        Console.WriteLine("  SGMGA_VCN_ORDERED calls SGMGA_VCN.");

        for (dim = 0; dim < dim_num; dim++)
        {
            level_1d_max[dim] = level_weight[dim] switch
            {
                > 0.0 => (int) Math.Floor(q_max / level_weight[dim]) + 1,
                _ => 0
            };
        }

        for (dim = 0; dim < dim_num; dim++)
        {
            level_1d_min[dim] = 0;
        }

        bool more_grids = false;

        Console.WriteLine("");
        Console.WriteLine("  IMPORTANCE:");
        cout = "";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + importance[dim].ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);
        Console.WriteLine("  LEVEL_WEIGHT:");
        cout = "";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_weight[dim].ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("  SGMGA_VCN_ORDERED_NAIVE:");
        Console.WriteLine("     I               Q   X");
        cout = "   MIN" + "  " + q_min.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_1d_min[dim].ToString().PadLeft(2);
        }

        Console.WriteLine(cout);

        int i = 0;

        for (;;)
        {
            SGMGAniso.sgmga_vcn_ordered_naive(ref data, dim_num, level_weight, level_1d_max, level_1d,
                q_min, q_max, ref more_grids);

            if (!more_grids)
            {
                break;
            }

            q = 0.0;
            for (dim = 0; dim < dim_num; dim++)
            {
                q += level_weight[dim] * level_1d[dim];
            }

            i += 1;
            cout = "  " + i.ToString().PadLeft(4)
                        + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (dim = 0; dim < dim_num; dim++)
            {
                cout += "  " + level_1d[dim].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
        }

        cout = "   MAX" + "  " + q_max.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_1d_max[dim].ToString().PadLeft(2);
        }

        Console.WriteLine(cout);

        Console.WriteLine("");
        Console.WriteLine("  SGMGA_VCN_ORDERED:");
        Console.WriteLine("     I               Q   X");
        cout = "   MIN" + "  " + q_min.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_1d_min[dim].ToString().PadLeft(2);
        }

        Console.WriteLine(cout);

        i = 0;

        for (;;)
        {
            SGMGAniso.sgmga_vcn_ordered(ref data, dim_num, level_weight, level_1d_max, ref level_1d,
                q_min, q_max, ref more_grids);

            if (!more_grids)
            {
                break;
            }

            q = 0.0;
            for (dim = 0; dim < dim_num; dim++)
            {
                q += level_weight[dim] * level_1d[dim];
            }

            i += 1;
            cout = "  " + i.ToString().PadLeft(4)
                        + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14);
            for (dim = 0; dim < dim_num; dim++)
            {
                cout += "  " + level_1d[dim].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
        }

        cout = "   MAX" + "  " + q_max.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_1d_max[dim].ToString().PadLeft(2);
        }

        Console.WriteLine(cout);
    }
}