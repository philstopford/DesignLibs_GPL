﻿using System;
using Burkardt.Composition;
using Burkardt.Types;

namespace LegendreProductPolynomialTest;

public static class compTest
{
    public static void comp_enum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMP_ENUM_TEST tests COMP_ENUM;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("COMP_ENUM_TEST");
        Console.WriteLine("  COMP_ENUM counts compositions;");
        Console.WriteLine("");
        string cout = "";
        for (n = 0; n <= 10; n++)
        {
            int k;
            for (k = 1; k <= 10; k++)
            {
                int num = Comp.comp_enum(n, k);
                cout += "  " + num.ToString().PadLeft(6);
            }

            Console.WriteLine(cout);
        }

    }

    public static void comp_next_grlex_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMP_NEXT_GRLEX_TEST tests COMP_NEXT_GRLEX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int kc = 3;
        int rank;
        int[] xc = new int[3];

        Console.WriteLine("");
        Console.WriteLine("COMP_NEXT_GRLEX_TEST");
        Console.WriteLine("  A COMP is a composition of an integer N into K parts.");
        Console.WriteLine("  Each part is nonnegative.  The order matters.");
        Console.WriteLine("  COMP_NEXT_GRLEX determines the next COMP in");
        Console.WriteLine("  graded lexicographic (grlex) order.");

        Console.WriteLine("");
        Console.WriteLine("  Rank:     NC       COMP");
        Console.WriteLine("  ----:     --   ------------");

        for (rank = 1; rank <= 71; rank++)
        {
            int j;
            switch (rank)
            {
                case 1:
                {
                    for (j = 0; j < kc; j++)
                    {
                        xc[j] = 0;
                    }

                    break;
                }
                default:
                    Comp.comp_next_grlex(kc, ref xc);
                    break;
            }

            int nc = typeMethods.i4vec_sum(kc, xc);

            string cout = "   " + rank.ToString().PadLeft(3) + ": ";
            cout += "    " + nc.ToString().PadLeft(2) + " = ";
            for (j = 0; j < kc - 1; j++)
            {
                cout += xc[j].ToString().PadLeft(2) + " + ";
            }

            Console.WriteLine(cout + xc[kc - 1].ToString().PadLeft(2) + "");
            //
            //  When XC(1) == NC, we have completed the compositions associated with
            //  a particular integer, and are about to advance to the next integer.
            //
            if (xc[0] == nc)
            {
                Console.WriteLine("  ----:     --   ------------");
            }
        }
    }

    public static void comp_random_grlex_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMP_RANDOM_GRLEX_TEST tests COMP_RANDOM_GRLEX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int rank = 0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("COMP_RANDOM_GRLEX_TEST");
        Console.WriteLine("  A COMP is a composition of an integer N into K parts.");
        Console.WriteLine("  Each part is nonnegative.  The order matters.");
        Console.WriteLine("  COMP_RANDOM_GRLEX selects a random COMP in");
        Console.WriteLine("  graded lexicographic (grlex) order between indices RANK1 and RANK2.");
        Console.WriteLine("");

        const int kc = 3;
        const int rank1 = 20;
        const int rank2 = 60;
        int seed = 123456789;

        for (test = 1; test <= 5; test++)
        {
            int[] xc = Comp.comp_random_grlex(kc, rank1, rank2, ref seed, ref rank);
            int nc = typeMethods.i4vec_sum(kc, xc);

            string cout = "   " + rank.ToString().PadLeft(3) + ": ";
            cout += "    " + nc.ToString().PadLeft(2) + " = ";
            int j;
            for (j = 0; j < kc - 1; j++)
            {
                cout += xc[j].ToString().PadLeft(2) + " + ";
            }

            Console.WriteLine(cout + xc[kc - 1].ToString().PadLeft(2) + "");
        }

    }

    public static void comp_rank_grlex_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMP_RANK_GRLEX_TEST tests COMP_RANK_GRLEX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int rank3 = 0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("COMP_RANK_GRLEX_TEST");
        Console.WriteLine("  A COMP is a composition of an integer N into K parts.");
        Console.WriteLine("  Each part is nonnegative.  The order matters.");
        Console.WriteLine("  COMP_RANK_GRLEX determines the rank of a COMP");
        Console.WriteLine("  from its parts.");
        Console.WriteLine("");
        Console.WriteLine("        Actual  Inferred");
        Console.WriteLine("  Test    Rank      Rank");
        Console.WriteLine("");

        const int kc = 3;
        const int rank1 = 20;
        const int rank2 = 60;
        int seed = 123456789;

        for (test = 1; test <= 5; test++)
        {
            int[] xc = Comp.comp_random_grlex(kc, rank1, rank2, ref seed, ref rank3);
            int rank4 = Comp.comp_rank_grlex(kc, xc);

            Console.WriteLine("  " + test.ToString().PadLeft(4) +
                              "  " + rank3.ToString().PadLeft(6) +
                              "  " + rank4.ToString().PadLeft(8) + "");
        }
    }

    public static void comp_unrank_grlex_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMP_UNRANK_GRLEX_TEST tests COMP_UNRANK_GRLEX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int kc = 3;
        int rank1;

        Console.WriteLine("");
        Console.WriteLine("COMP_UNRANK_GRLEX_TEST");
        Console.WriteLine("  A COMP is a composition of an integer N into K parts.");
        Console.WriteLine("  Each part is nonnegative.  The order matters.");
        Console.WriteLine("  COMP_UNRANK_GRLEX determines the parts");
        Console.WriteLine("  of a COMP from its rank.");

        Console.WriteLine("");
        Console.WriteLine("  Rank: ->  NC       COMP");
        Console.WriteLine("  ----:     --   ------------");

        for (rank1 = 1; rank1 <= 71; rank1++)
        {
            int[] xc = Comp.comp_unrank_grlex(kc, rank1);
            int nc = typeMethods.i4vec_sum(kc, xc);

            string cout = "   " + rank1.ToString().PadLeft(3) + ": ";
            cout += "    " + nc.ToString().PadLeft(2) + " = ";
            int j;
            for (j = 0; j < kc - 1; j++)
            {
                cout += xc[j].ToString().PadLeft(2) + " + ";
            }

            Console.WriteLine(cout + xc[kc - 1].ToString().PadLeft(2) + "");
            //
            //  When XC(1) == NC, we have completed the compositions associated with
            //  a particular integer, and are about to advance to the next integer.
            //
            if (xc[0] == nc)
            {
                Console.WriteLine("  ----:     --   ------------");
            }
        }
    }

}