using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest;

internal partial class Program
{
    private static void knapsack_01_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KNAPSACK_01_TEST tests KNAPSACK_01.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        int i;
        double mass = 0;
        double mass_limit = 26.0;
        int n = N;
        double[] p =  {
                24.0, 13.0, 23.0, 15.0, 16.0
            }
            ;
        double profit = 0;
        double[] w =  {
                12.0, 7.0, 11.0, 8.0, 9.0
            }
            ;
        double[] x = new double[N];

        Console.WriteLine("");
        Console.WriteLine("KNAPSACK_01_TEST");
        Console.WriteLine("  KNAPSACK_01 solves the 0/1 knapsack problem.");

        Console.WriteLine("");
        Console.WriteLine("  Object, Profit, Mass, Profit Density");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(6)
                                   + "  " + p[i].ToString().PadLeft(7)
                                   + "  " + w[i].ToString().PadLeft(7)
                                   + "  " + (p[i] / w[i]).ToString().PadLeft(7) + "");
        }

        Ranking.knapsack_reorder(n, ref p, ref w);

        Console.WriteLine("");
        Console.WriteLine("  After reordering by Profit Density:");
        Console.WriteLine("");
        Console.WriteLine("  Object, Profit, Mass, Profit Density");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(6)
                                   + "  " + p[i].ToString().PadLeft(7)
                                   + "  " + w[i].ToString().PadLeft(7)
                                   + "  " + (p[i] / w[i]).ToString().PadLeft(7) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Total mass restriction is " + mass_limit + "");

        Ranking.knapsack_01(n, mass_limit, ref p, ref w, ref x, ref mass, ref profit);

        Console.WriteLine("");
        Console.WriteLine("  Object, Density, Choice, Profit, Mass");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(6)
                                   + "  " + (p[i] / w[i]).ToString().PadLeft(7)
                                   + "  " + x[i].ToString().PadLeft(7)
                                   + "  " + (x[i] * p[i]).ToString().PadLeft(7)
                                   + "  " + (x[i] * w[i]).ToString().PadLeft(7) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Total:            " + profit
                                                 + "  " + mass + "");
    }

    private static void knapsack_rational_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KNAPSACK_RATIONAL_TEST tests KNAPSACK_RATIONAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        int i;
        double mass = 0;
        double mass_limit = 26.0;
        int n = N;
        double[] p =  {
                24.0, 13.0, 23.0, 15.0, 16.0
            }
            ;
        double profit = 0;
        double[] w =  {
                12.0, 7.0, 11.0, 8.0, 9.0
            }
            ;
        double[] x = new double[N];

        Console.WriteLine("");
        Console.WriteLine("KNAPSACK_RATIONAL_TEST");
        Console.WriteLine("  KNAPSACK_RATIONAL solves the rational knapsack problem.");

        Console.WriteLine("");
        Console.WriteLine("  Object, Profit, Mass, Profit Density");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + (i + 1).ToString().PadLeft(4)
                                   + "  " + p[i].ToString().PadLeft(7)
                                   + "  " + w[i].ToString().PadLeft(7)
                                   + "  " + (p[i] / w[i]).ToString().PadLeft(7) + "");
        }

        Ranking.knapsack_reorder(n, ref p, ref w);

        Console.WriteLine("");
        Console.WriteLine("  After reordering by Profit Density:");
        Console.WriteLine("");
        Console.WriteLine("  Object, Profit, Mass, Profit Density");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + (i + 1).ToString().PadLeft(4)
                                   + "  " + p[i].ToString().PadLeft(7)
                                   + "  " + w[i].ToString().PadLeft(7)
                                   + "  " + (p[i] / w[i]).ToString().PadLeft(7) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Total mass restriction is " + mass_limit + "");

        Ranking.knapsack_rational(n, mass_limit, p, w, ref x, ref mass, ref profit);

        Console.WriteLine("");
        Console.WriteLine("  Object, Density, Choice, Profit, Mass");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + (i + 1).ToString().PadLeft(4)
                                   + "  " + (p[i] / w[i]).ToString().PadLeft(7)
                                   + "  " + (x[i] * p[i]).ToString().PadLeft(7)
                                   + "  " + (x[i] * w[i]).ToString().PadLeft(7) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Total:            " + profit
                                                 + "  " + mass + "");
    }

    private static void knapsack_reorder_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KNAPSACK_REORDER_TEST tests KNAPSACK_REORDER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n = 5;
        double[] p =  {
                24.0, 13.0, 23.0, 15.0, 16.0
            }
            ;
        double[] w =  {
                12.0, 7.0, 11.0, 8.0, 9.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("KNAPSACK_REORDER_TEST");
        Console.WriteLine("  KNAPSACK_REORDER reorders the knapsack data.");

        Console.WriteLine("");
        Console.WriteLine("  Object, Profit, Mass, Profit Density");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(6)
                                   + "  " + p[i].ToString().PadLeft(7)
                                   + "  " + w[i].ToString().PadLeft(7)
                                   + "  " + (p[i] / w[i]).ToString().PadLeft(7) + "");
        }

        Ranking.knapsack_reorder(n, ref p, ref w);

        Console.WriteLine("");
        Console.WriteLine("  After reordering by Profit Density:");
        Console.WriteLine("");
        Console.WriteLine("  Object, Profit, Mass, Profit Density");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(6)
                                   + "  " + p[i].ToString().PadLeft(7)
                                   + "  " + w[i].ToString().PadLeft(7)
                                   + "  " + (p[i] / w[i]).ToString().PadLeft(7) + "");
        }
    }
}