using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest;

internal partial class Program
{
    private static void tableau_check_test()

        //****************************************************************************80
        //
        // TABLEAU_CHECK_TEST tests TABLEAU_CHECK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        bool check;
        int n = 0;
        int[] t = new int[1];
        int[] t1 =  {
            }
            ;
        int[] t2 =  {
                1, 2,
                2, 4,
                3, 7,
                4, 9
            }
            ;
        int[] t3 =  {
                1, 2,
                3, 4,
                5, 5,
                3, 3
            }
            ;
        int[] t4 =  {
                1, 2,
                3, 4,
                4, 5,
                5, 3
            }
            ;
        int[] t5 =  {
                1, 3,
                3, 4,
                6, 7,
                7, 8
            }
            ;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TABLEAU_CHECK TEST");
        Console.WriteLine("  TABLEAU_CHECK checks a 2xN tableau.");
        Console.WriteLine("");
        Console.WriteLine("  Check?");
        Console.WriteLine("");

        for (test = 1; test <= 5; test++)
        {
            switch (test)
            {
                case 1:
                    n = 0;
                    t = typeMethods.i4mat_copy_new(2, n, t1);
                    break;
                case 2:
                    n = 4;
                    t = typeMethods.i4mat_copy_new(2, n, t2);
                    break;
                case 3:
                    n = 4;
                    t = typeMethods.i4mat_copy_new(2, n, t3);
                    break;
                case 4:
                    n = 4;
                    t = typeMethods.i4mat_copy_new(2, n, t4);
                    break;
                case 5:
                    n = 4;
                    t = typeMethods.i4mat_copy_new(2, n, t5);
                    break;
            }

            Console.WriteLine("");
            check = Ranking.tableau_check(n, t);
            Console.WriteLine("      Check = " + check + "");
            typeMethods.i4mat_print(2, n, t, "  Tableau:");
        }
    }

    private static void tableau_enum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TABLEAU_ENUM_TEST tests TABLEAU_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        int tableau_num;

        Console.WriteLine("");
        Console.WriteLine("TABLEAU_ENUM_TEST");
        Console.WriteLine("  TABLEAU_ENUM enumerates tableaus on N nodes.");

        for (n = 0; n <= 10; n++)
        {
            tableau_num = Ranking.tableau_enum(n);
            Console.WriteLine("  " + n.ToString().PadLeft(2)
                                   + "  " + tableau_num.ToString().PadLeft(6) + "");
        }
    }

    private static void tableau_to_bal_seq_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TABLEAU_TO_BAL_SEQ_TEST tests TABLEAU_TO_BAL_SEQ.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 4;
        int[] t;
        int[] tab =  {
                1, 3, 2, 4, 5, 7, 6, 8
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TABLEAU_TO_BAL_SEQ_TEST");
        Console.WriteLine("  TABLEAU_TO_BAL_SEQ converts a tableau");
        Console.WriteLine("  to a balanced sequence.");

        typeMethods.i4mat_print(2, n, tab, "  Tableau");

        t = Ranking.tableau_to_bal_seq(n, tab);

        typeMethods.i4vec_transpose_print(2 * n, t, "  Balanced sequence:");
    }
}