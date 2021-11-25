﻿using System;
using Burkardt.Composition;
using Burkardt.SubsetNS;
using Burkardt.Types;

namespace SubsetTestNS;

public static class CompTest
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
            
        for (n = 0; n <= 10; n++)
        {
            string cout = "";
            int k;
            for (k = 1; k <= 10; k++)
            {
                int num = Comp.comp_enum(n, k);
                cout += "  " + num.ToString().PadLeft(6);
            }

            Console.WriteLine(cout);
        }
    }

    public static void comp_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMP_NEXT_TEST tests COMP_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int K = 3;

        int[] a = new int[K];
        int h = 0;
        const int n = 6;
        int t = 0;

        Console.WriteLine("");
        Console.WriteLine("COMP_NEXT_TEST");
        Console.WriteLine("  COMP_NEXT produces compositions.");
        Console.WriteLine("");
        Console.WriteLine("  Seeking all compositions of N = " + n + "");
        Console.WriteLine("  using " + K + " parts.");
        Console.WriteLine("");

        bool more = false;
        int index = 0;

        for (;;)
        {
            Comp.comp_next(n, K, ref a, ref more, ref h, ref t);

            index += 1;
            string cout = "  ";
            cout += "  " + index.ToString().PadLeft(4) + "  ";
            int i;
            for (i = 0; i < K; i++)
            {
                cout += a[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);

            if (!more)
            {
                break;
            }
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

    public static void comp_random_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMP_RANDOM_TEST tests COMP_RANDOM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int K = 5;

        int[] a = new int[K];
        int j;
        const int n = 10;

        Console.WriteLine("");
        Console.WriteLine("COMP_RANDOM_TEST");
        Console.WriteLine("  COMP_RANDOM produces compositions at random.");
        Console.WriteLine("");
        Console.WriteLine("  Seeking random compositions of N = " + n + "");
        Console.WriteLine("  using " + K + " parts.");
        Console.WriteLine("");

        int seed = 123456789;

        for (j = 1; j <= 5; j++)
        {
            Comp.comp_random(n, K, ref seed, ref a);

            string cout = "  ";
            int i;
            for (i = 0; i < K; i++)
            {
                cout += a[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);
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

    public static void comp_to_ksub_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMP_TO_KSUB_TEST tests COMP_TO_KSUB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] ac = new int[5];
        int[] as_ = new int [4];
        int i;
        int ks = 0;
        int ns = 0;

        Console.WriteLine("");
        Console.WriteLine("COMP_TO_KSUB_TEST");
        Console.WriteLine("  COMP_TO_KSUB returns the K subset corresponding to a composition.");

        int nc = 10;
        int kc = 5;
        int seed = 123456789;

        for (i = 1; i <= 5; i++)
        {
            Console.WriteLine("");

            Comp.comp_random(nc, kc, ref seed, ref ac);
            string cout = "  COMP:";
            int j;
            for (j = 0; j < kc; j++)
            {
                cout += ac[j].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);

            Comp.comp_to_ksub(nc, kc, ac, ref ns, ref ks,  ref as_ );
            cout = "  KSUB:";
            for (j = 0; j < ks; j++)
            {
                cout += as_[j].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);

            Ksub.ksub_to_comp(ns, ks,  as_, ref nc, ref kc, ref ac );
            cout = "  COMP:";
            for (j = 0; j < kc; j++)
            {
                cout += ac[j].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
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

    public static void compnz_to_ksub_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMPNZ_TO_KSUB_TEST tests COMPNZ_TO_KSUB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] ac = new int[5];
        int[] as_ = new int[4];
        int i;
        int ks = 0;
        int ns = 0;

        Console.WriteLine("");
        Console.WriteLine("COMPNZ_TO_KSUB_TEST");
        Console.WriteLine("  COMPNZ_TO_KSUB returns the K subset corresponding");
        Console.WriteLine("  to a nonzero composition.");

        int nc = 10;
        int kc = 5;
        int seed = 123456789;

        for (i = 1; i <= 5; i++)
        {
            Console.WriteLine("");

            Comp.compnz_random(nc, kc, ref seed, ref ac);
            string cout = "  COMPNZ:";
            int j;
            for (j = 0; j < kc; j++)
            {
                cout += ac[j].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);

            Comp.compnz_to_ksub(nc, kc, ac, ref ns, ref ks,  ref as_ );
            cout = "  KSUB:  ";
            for (j = 0; j < ks; j++)
            {
                cout += as_[j].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);

            Ksub.ksub_to_compnz(ns, ks,  as_, ref nc, ref kc, ref ac );
            cout = "  COMPNZ:";
            for (j = 0; j < kc; j++)
            {
                cout += ac[j].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    public static void compnz_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMPNZ_NEXT_TEST tests COMPNZ_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a = new int[3];

        const int n = 6;
        const int k = 3;
        bool more = false;

        CompNZData data = new();

        Console.WriteLine("");
        Console.WriteLine("COMPNZ_NEXT_TEST");
        Console.WriteLine("  COMPNZ_NEXT produces compositions using nonzero parts.");
        Console.WriteLine("");
        Console.WriteLine("  Seeking all compositions of N = " + n + "");
        Console.WriteLine("  using " + k + " nonzero parts.");
        Console.WriteLine("");

        for (;;)
        {
            Comp.compnz_next(ref data, n, k, ref a, ref more);

            string cout = "  ";
            int i;
            for (i = 0; i < k; i++)
            {
                cout += a[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);

            if (!more)
            {
                break;
            }

        }
    }

    public static void compnz_random_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMPNZ_RANDOM_TEST tests COMPNZ_RANDOM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int K = 5;

        int[] a = new int[K];
        int j;
        const int n = 10;

        Console.WriteLine("");
        Console.WriteLine("COMPNZ_RANDOM_TEST");
        Console.WriteLine("  COMPNZ_RANDOM produces compositions at random");
        Console.WriteLine("  with only nonzero parts.");
        Console.WriteLine("");
        Console.WriteLine("  Seeking random compositions of N = " + n + "");
        Console.WriteLine("  using " + K + " nonzero parts.");
        Console.WriteLine("");

        int seed = 123456789;

        for (j = 1; j <= 5; j++)
        {
            Comp.compnz_random(n, K, ref seed, ref a);

            string cout = "  ";
            int i;
            for (i = 0; i < K; i++)
            {
                cout += a[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);
        }
    }

}