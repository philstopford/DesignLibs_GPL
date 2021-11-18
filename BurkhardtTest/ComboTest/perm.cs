using System;
using Burkardt;
using Burkardt.Types;

namespace ComboTest;

internal partial class Program
{
    private static void perm_check_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_CHECK_TEST tests PERM_CHECK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        bool check;
        int n = 0;
        int[] s1 =  {
                5, 1, 8, 3, 4
            }
            ;
        int[] s2 =  {
                5, 1, 4, 3, 4
            }
            ;
        int[] s3 =  {
                5, 1, 2, 3, 4
            }
            ;
        int[] s = new int[1];
        int test;

        Console.WriteLine("");
        Console.WriteLine("PERM_CHECK TEST");
        Console.WriteLine("  PERM_CHECK checks a permutation.");

        for (test = 1; test <= 3; test++)
        {
            switch (test)
            {
                case 1:
                    n = 5;
                    s = typeMethods.i4vec_copy_new(n, s1);
                    break;
                case 2:
                    n = 5;
                    s = typeMethods.i4vec_copy_new(n, s2);
                    break;
                case 3:
                    n = 5;
                    s = typeMethods.i4vec_copy_new(n, s3);
                    break;
            }

            check = Permutation.perm_check(n, s);
            typeMethods.i4vec_transpose_print(n, s, "  Permutation:");
            Console.WriteLine("  Check = " + check + "");
        }
    }

    private static void perm_enum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_ENUM_TEST tests PERM_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        int perm_num;

        Console.WriteLine("");
        Console.WriteLine("PERM_ENUM_TEST");
        Console.WriteLine("  PERM_ENUM enumerates permutations of N items.");

        for (n = 0; n <= 10; n++)
        {
            perm_num = Permutation.perm_enum(n);
            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + perm_num.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }
    }

    private static void perm_to_cycle_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_TO_CYCLE_TEST tests PERM_TO_CYCLE.
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
        int i;
        int j;
        int jlo;
        int[] index;
        int n = 7;
        int ncycle = 0;
        int[] p =  {
                4, 5, 1, 2, 3, 6, 7
            }
            ;
        int[] t;

        Console.WriteLine("");
        Console.WriteLine("PERM_TO_CYCLE_TEST");
        Console.WriteLine("  PERM_TO_CYCLE converts a permutation from");
        Console.WriteLine("  array to cycle form.");

        Permutation.perm_print(n, p, "  Permutation:");
        //
        //  Convert the permutation to cycle form.
        //
        t = new int[n];
        index = new int[n];

        Permutation.perm_to_cycle(n, p, ref ncycle, ref t, ref index);

        Console.WriteLine("");
        Console.WriteLine("  Corresponding cycle form:");
        Console.WriteLine("  Number of cycles is " + ncycle + "");
        Console.WriteLine("");
        jlo = 0;
        for (i = 1; i <= ncycle; i++)
        {
            for (j = jlo + 1; j <= jlo + index[i - 1]; j++)
            {
            }

            Console.WriteLine("");
            jlo += index[i - 1];
        }
    }

    private static void perm_inv_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_INV_TEST tests PERM_INV.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 4;
        int[] p =  {
                3, 1, 2, 4
            }
            ;
        int[] q;
        int[] r;

        Console.WriteLine("");
        Console.WriteLine("PERM_INV_TEST");
        Console.WriteLine("  PERM_INV computes an inverse permutation,");

        Permutation.perm_print(n, p, "  The permutation P:");
        //
        //  Invert.
        //
        q = Permutation.perm_inv(n, p);

        Permutation.perm_print(n, q, "  The inverse permutation Q:");
        //
        //  Multiply.
        //
        r = Permutation.perm_mul(n, p, q);

        Permutation.perm_print(n, r, "  The product R = P * Q:");
    }

    private static void perm_lex_rank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_LEX_RANK_TEST tests PERM_LEX_RANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 4;
        int[] pi =  {
                3, 1, 2, 4
            }
            ;
        int rank;

        Console.WriteLine("");
        Console.WriteLine("PERM_LEX_RANK_TEST");
        Console.WriteLine("  PERM_LEX_RANK ranks");
        Console.WriteLine("  permutations of the integers,");
        Console.WriteLine("  using the lexicographic ordering:");

        Permutation.perm_print(n, pi, "  Element to be ranked:");

        rank = Permutation.perm_lex_rank(n, pi);

        Console.WriteLine("");
        Console.WriteLine("  Rank is computed as " + rank + "");
    }

    private static void perm_lex_successor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_LEX_SUCCESSOR_TEST tests PERM_LEX_SUCCESSOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n = 4;
        int[] pi;
        int rank;
        int rank_old;

        Console.WriteLine("");
        Console.WriteLine("PERM_LEX_SUCCESSOR_TEST");
        Console.WriteLine("  PERM_LEX_SUCCESSOR lists");
        Console.WriteLine("  permutations of the integers,");
        Console.WriteLine("  using the lexicographic ordering:");

        pi = new int[n];

        rank = -1;

        for (;;)
        {
            rank_old = rank;

            Permutation.perm_lex_successor(n, ref pi, ref rank);

            if (rank <= rank_old)
            {
                break;
            }

            string cout = "  " + rank.ToString(CultureInfo.InvariantCulture).PadLeft(4);
            for (i = 0; i < n; i++)
            {
                cout += "  " + pi[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void perm_lex_unrank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_LEX_UNRANK_TEST tests PERM_LEX_UNRANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 4;
        int[] pi;
        int rank;

        Console.WriteLine("");
        Console.WriteLine("PERM_LEX_UNRANK_TEST");
        Console.WriteLine("  PERM_LEX_UNRANK unranks");
        Console.WriteLine("  permutations of the integers,");
        Console.WriteLine("  using the lexicographic ordering:");

        rank = 12;

        pi = Permutation.perm_lex_unrank(rank, n);

        Permutation.perm_print(n, pi, "  The element of rank 12:");
    }

    private static void perm_mul_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_MUL_TEST tests PERM_MUL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 4;
        int[] p =  {
                3, 1, 2, 4
            }
            ;
        int[] q =  {
                2, 3, 1, 4
            }
            ;
        int[] r;

        Console.WriteLine("");
        Console.WriteLine("PERM_MUL_TEST");
        Console.WriteLine("  PERM_MUL multiplies two permutations.");

        Permutation.perm_print(n, p, "  The permutation P:");

        Permutation.perm_print(n, q, "  The permutation Q:");

        r = Permutation.perm_mul(n, p, q);

        Permutation.perm_print(n, r, "  The product R = P * Q:");
    }

    private static void perm_parity_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_PARITY_TEST tests PERM_PARITY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        int[] p;
        int parity;
        int seed;
        int test;

        Console.WriteLine("");
        Console.WriteLine("PERM_PARITY_TEST");
        Console.WriteLine("  PERM_PARITY computes the parity of a permutation.");

        n = 5;
        seed = 123456789;

        for (test = 1; test <= 5; test++)
        {
            p = Permutation.perm_random(n, ref seed);
            Permutation.perm_print(n, p, "  The permutation P:");
            parity = Permutation.perm_parity(n, p);
            Console.WriteLine("");
            Console.WriteLine("  The parity is " + parity + "");
        }
    }

    private static void perm_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_PRINT_TEST tests PERM_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 7;
        int[] p =  {
                7, 2, 4, 1, 5, 3, 6
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("PERM_PRINT_TEST");
        Console.WriteLine("  PERM_PRINT prints a permutation of (1,...,N).");

        Permutation.perm_print(n, p, "  The 1-based permutation:");
    }

    private static void perm_random_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_RANDOM_TEST tests PERM_RANDOM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        int[] p;
        int seed;
        int test;

        Console.WriteLine("");
        Console.WriteLine("PERM_RANDOM_TEST:");
        Console.WriteLine("  PERM_RANDOM randomly selects a permutation of 1, ..., N.");

        n = 5;
        seed = 123456789;

        for (test = 1; test <= 5; test++)
        {
            p = Permutation.perm_random(n, ref seed);
            typeMethods.i4vec_transpose_print(n, p, "");
        }
    }

    private static void perm_tj_rank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_TJ_RANK_TEST tests PERM_TJ_RANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 4;
        int[] pi =  {
                4, 3, 2, 1
            }
            ;
        int rank;

        Console.WriteLine("");
        Console.WriteLine("PERM_TJ_RANK_TEST");
        Console.WriteLine("  PERM_TJ_RANK ranks permutations of the integers");
        Console.WriteLine("  using the Trotter-Johnson ordering.");

        Permutation.perm_print(n, pi, "  The element to be ranked:");

        rank = Permutation.perm_tj_rank(n, pi);

        Console.WriteLine("");
        Console.WriteLine("  Rank is computed as " + rank + "");
    }

    private static void perm_tj_successor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_TJ_SUCCESSOR_TEST tests PERM_TJ_SUCCESSOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n = 4;
        int[] pi;
        int rank;
        int rank_old;

        Console.WriteLine("");
        Console.WriteLine("PERM_TJ_SUCCESOR_TEST");
        Console.WriteLine("  PERM_TJ_SUCCESSOR lists");
        Console.WriteLine("  permutations of the integers");
        Console.WriteLine("  using the Trotter-Johnson ordering.");

        pi = new int[n];

        rank = -1;

        for (;;)
        {
            rank_old = rank;

            Permutation.perm_tj_successor(n, ref pi, ref rank);

            if (rank <= rank_old)
            {
                break;
            }

            string cout = "  " + rank.ToString(CultureInfo.InvariantCulture).PadLeft(4);
            for (i = 0; i < n; i++)
            {
                cout += "  " + pi[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void perm_tj_unrank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_TJ_UNRANK_TEST tests PERM_TJ_UNRANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        int[] pi;
        int rank;

        Console.WriteLine("");
        Console.WriteLine("PERM_TJ_UNRANK_TEST");
        Console.WriteLine("  PERM_TJ_UNRANK unranks");
        Console.WriteLine("  permutations of the integers");
        Console.WriteLine("  using the Trotter-Johnson ordering.");

        rank = 12;
        n = 4;

        pi = Permutation.perm_tj_unrank(rank, n);

        Permutation.perm_print(n, pi, "  The element of rank 12:");
    }
}