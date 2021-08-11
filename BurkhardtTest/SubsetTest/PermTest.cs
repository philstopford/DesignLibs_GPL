using System;
using Burkardt;
using Burkardt.SubsetNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SubsetTest
{
    public static class PermTest
    {
        public static void multiperm_enum_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MULTIPERM_ENUM_TEST tests MULTIPERM_ENUM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 5;

            int[] counts = new int[N];
            int i;
            int k;
            int number;
            int seed = 123456789;
            int test;
            int test_num = 5;

            Console.WriteLine("");
            Console.WriteLine("MULTIPERM_ENUM_TEST");
            Console.WriteLine("  MULTIPERM_ENUM enumerates multipermutations.");
            Console.WriteLine("");
            Console.WriteLine("  N is the number of objects to be permuted.");
            Console.WriteLine("  K is the number of distinct types of objects.");
            Console.WriteLine("  COUNTS is the number of objects of each type.");
            Console.WriteLine("  NUMBER is the number of multipermutations.");
            Console.WriteLine("");
            Console.WriteLine("  Number       N       K       Counts(1:K)");
            Console.WriteLine("");

            for (test = 1; test <= test_num; test++)
            {
                k = UniformRNG.i4_uniform_ab(1, N, ref seed);

                Comp.compnz_random(N, k, ref seed, ref counts);

                number = Permutation.multiperm_enum(N, k, counts);

                string cout = "  " + number.ToString().PadLeft(6)
                                   + "  " + N.ToString().PadLeft(6)
                                   + "  " + k.ToString().PadLeft(6);
                for (i = 0; i < k; i++)
                {
                    cout += "  " + counts[i].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

        public static void multiperm_next_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MULTIPERM_NEXT_TEST tests MULTIPERM_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 March 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 6;

            int[] a = { 1, 2, 2, 3, 3, 3 };
            int i;
            bool more;
            int tally;

            Console.WriteLine("");
            Console.WriteLine("MULTIPERM_NEXT_TEST");
            Console.WriteLine("  MULTIPERM_NEXT computes multipermutations in");
            Console.WriteLine("  lexical order.");
            Console.WriteLine("");

            tally = 0;
            more = true;

            while (more)
            {
                tally = tally + 1;

                string cout = "  " + tally.ToString().PadLeft(4);
                for (i = 0; i < N; i++)
                {
                    cout += "  " + a[i].ToString().PadLeft(2);
                }

                Console.WriteLine(cout);

                Permutation.multiperm_next(N, ref a, ref more);
            }

            return;

        }

        public static void perm_ascend_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_ASCEND_TEST tests PERM_ASCEND.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 9;

            int length = 0;
            int[] p = { 1, 2, 8, 5, 6, 7, 4, 3, 0 };
            int[] subseq = new int[N];

            Console.WriteLine("");
            Console.WriteLine("PERM_ASCEND_TEST");
            Console.WriteLine("  PERM_ASCEND determines the length of the longest");
            Console.WriteLine("  increasing subsequence in a permutation.");

            Permutation.perm0_print(N, p, "  The permutation:");

            Permutation.perm_ascend(N, p, ref length, ref subseq);

            Console.WriteLine("");
            Console.WriteLine("  The length of the longest increasing subsequence is " +
                              length + "");

            typeMethods.i4vec1_print(length, subseq, "  A longest increasing subsequence:");


        }

        public static void perm_fixed_enum_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_FIXED_ENUM_TEST tests PERM_FIXED_ENUM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int fnm;
            int m;
            int n = 10;

            Console.WriteLine("");
            Console.WriteLine("PERM_FIXED_ENUM_TEST");
            Console.WriteLine("  PERM_FIXED_ENUM enumerates the permutations");
            Console.WriteLine("  of N objects that leave M unchanged.");
            Console.WriteLine("");
            Console.WriteLine("  For this test, N = " + n + "");
            Console.WriteLine("");
            Console.WriteLine("    M    F(N,M)");
            Console.WriteLine("");

            for (m = 0; m <= n; m++)
            {
                fnm = Permutation.perm_fixed_enum(n, m);

                Console.WriteLine("  "
                                  + m.ToString().PadLeft(3) + "  "
                                  + fnm.ToString().PadLeft(8) + "");
            }
        }

        public static void perm0_print_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_PRINT_TEST tests PERM0_PRINT.
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
            int[] p = { 6, 1, 3, 0, 4, 2, 5 };

            Console.WriteLine("");
            Console.WriteLine("PERM0_PRINT_TEST");
            Console.WriteLine("  PERM0_PRINT prints a permutation of (0,...,N-1).");

            Permutation.perm0_print(n, p, "  The 0-based permutation:");
        }

        public static void perm0_break_count_test()

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    PERM0_BREAK_COUNT_TEST tests PERM0_BREAK_COUNT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int break_count;
            int n = 6;
            int[] p = { 3, 4, 1, 0, 5, 2 };

            Console.WriteLine("");
            Console.WriteLine("PERM0_BREAK_COUNT_TEST");
            Console.WriteLine("  PERM0_BREAK_COUNT counts the breaks in a permutation.");

            Permutation.perm0_print(n, p, "  The permutation:");

            break_count = Permutation.perm0_break_count(n, p);

            Console.WriteLine("");
            Console.WriteLine("  The number of breaks is " + break_count + "");
        }

        public static void perm0_check_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_CHECK_TEST tests PERM0_CHECK.
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
            int n = 5;
            int[] p1 = { 5, 2, 3, 4, 1 };
            int[] p2 = { 4, 1, 3, 0, 2 };
            int[] p3 = { 0, 2, 1, 3, 2 };

            Console.WriteLine("");
            Console.WriteLine("PERM0_CHECK_TEST");
            Console.WriteLine("  PERM0_CHECK checks a permutation of (0,...,N-1).");
            Console.WriteLine("");

            Permutation.perm0_print(n, p1, "  Permutation 1:");
            Permutation.perm0_check(n, p1);

            Permutation.perm0_print(n, p2, "  Permutation 2:");
            Permutation.perm0_check(n, p2);

            Permutation.perm0_print(n, p3, "  Permutation 3:");
            Permutation.perm0_check(n, p3);
        }

        public static void perm0_cycle_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_CYCLE_TEST tests PERM0_CYCLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 9;

            int iopt = 1;
            int isgn = 0;
            int ncycle = 0;
            int[] p = { 1, 2, 8, 5, 6, 7, 4, 3, 0 };

            Console.WriteLine("");
            Console.WriteLine("PERM0_CYCLE_TEST");
            Console.WriteLine("  PERM0_CYCLE analyzes a permutation of (0,...,N-1).");

            Permutation.perm0_print(N, p, "  The permutation:");

            Permutation.perm0_cycle(N, p, ref isgn, ref ncycle, iopt);

            Console.WriteLine("");
            Console.WriteLine("  NCYCLE = " + ncycle + "");
            Console.WriteLine("  ISGN =   " + isgn + "");

            Permutation.perm0_print(N, p, "  The permutation in cycle form:");
        }

        public static void perm0_distance_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_DISTANCE_TEST tests PERM0_DISTANCE
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            int k11;
            int k12;
            int k13;
            int k21;
            int k23;
            int[] p1 = new int[N];
            int[] p2 = new int[N];
            int[] p3 = new int[N];
            int seed;

            Console.WriteLine("");
            Console.WriteLine("PERM0_DISTANCE_TEST");
            Console.WriteLine("  PERM0_DISTANCE computes the Ulam metric distance");
            Console.WriteLine("  between two permutations of (0,...,N-1).");

            seed = 123456789;

            Permutation.perm0_random2(N, ref seed, ref p1);
            Permutation.perm0_random2(N, ref seed, ref p2);
            Permutation.perm0_random2(N, ref seed, ref p3);

            Permutation.perm0_print(N, p1, "  Permutation P1");
            Permutation.perm0_print(N, p2, "  Permutation P2");
            Permutation.perm0_print(N, p3, "  Permutation P3");

            k11 = Permutation.perm0_distance(N, p1, p1);
            k12 = Permutation.perm0_distance(N, p1, p2);
            k21 = Permutation.perm0_distance(N, p2, p1);
            k13 = Permutation.perm0_distance(N, p1, p3);
            k23 = Permutation.perm0_distance(N, p2, p3);

            Console.WriteLine("");
            Console.WriteLine("  K(P1,P1) should be 0.");
            Console.WriteLine("");
            Console.WriteLine("  K(P1,P1) = " + k11 + "");
            Console.WriteLine("");
            Console.WriteLine("  K(P1,P2) should equal K(P2,P1).");
            Console.WriteLine("");
            Console.WriteLine("  K(P1,P2) = " + k12 + "");
            Console.WriteLine("  K(P2,P1) = " + k21 + "");
            Console.WriteLine("");
            Console.WriteLine("  K(P1,P2) + K(P2,P3) >= K(P1,P3).");
            Console.WriteLine("");
            Console.WriteLine("  K(P1,P3) = " + k13 + "");
            Console.WriteLine("  K(P1,P2) = " + k12 + "");
            Console.WriteLine("  K(P2,P3) = " + k23 + "");
            Console.WriteLine("  K(P1,P2) + K(P2,P3) = " + k12 + k23 + "");
        }

        public static void perm0_free_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_FREE_TEST tests PERM0_FREE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int[] ifree = new int[5];
            int[] ipart = new int[5];
            int n = 5;
            int nfree;
            int npart;
            int[] p = { 4, 1, 2, 3, 0 };

            Console.WriteLine("");
            Console.WriteLine("PERM0_FREE_TEST");
            Console.WriteLine("  PERM0_FREE returns the unused values in a partial permutation");
            Console.WriteLine("  of (0,...,N-1).");

            for (npart = 0; npart <= n; npart++)
            {
                for (i = 0; i < npart; i++)
                {
                    ipart[i] = p[i];
                }

                nfree = n - npart;
                Permutation.perm0_free(npart, ipart, nfree, ref ifree);
                typeMethods.i4vec_transpose_print(npart, ipart, "  Partial permutation:");
                typeMethods.i4vec_transpose_print(nfree, ifree, "  Values not yet used:");
            }
        }

        public static void perm0_inverse_test()

            //****************************************************************************
            //
            //  Purpose:
            //
            //    PERM0_INVERSE_TEST tests PERM0_INVERSE;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n = 7;
            int[] p1 = { 3, 2, 4, 0, 6, 5, 1 };
            int[] p2;

            Console.WriteLine("");
            Console.WriteLine("PERM0_INVERSE_TEST");
            Console.WriteLine("  PERM0_INVERSE inverts a permutation of (0,...,N-1);");
            Console.WriteLine("");

            Permutation.perm0_print(n, p1, "  The original permutation:");

            p2 = Permutation.perm0_inverse(n, p1);

            Permutation.perm0_print(n, p2, "  The inverted permutation:");
        }

        public static void perm0_inverse2_test()

            //****************************************************************************
            //
            //  Purpose:
            //
            //    PERM0_INVERSE2_TEST tests PERM0_INVERSE2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 7;

            int[] p = { 3, 2, 4, 0, 6, 5, 1 };

            Console.WriteLine("");
            Console.WriteLine("PERM0_INVERSE2_TEST");
            Console.WriteLine("  PERM0_INVERSE2 inverts a permutation of (0,...,N-1).");

            Permutation.perm0_print(N, p, "  The original permutation:");

            Permutation.perm0_inverse2(N, p);

            Permutation.perm0_print(N, p, "  The inverted permutation:");

        }

        public static void perm0_inverse3_new_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_INVERSE3_NEW_TEST tests PERM0_INVERSE3_NEW.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 May 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 7;

            int[] perm = { 3, 2, 4, 0, 6, 5, 1 };
            int[] perm_inv;

            Console.WriteLine("");
            Console.WriteLine("PERM0_INVERSE3_NEW_TEST");
            Console.WriteLine("  PERM0_INVERSE3_NEW inverts a permutation of (0,...,N-1).");

            Permutation.perm0_print(N, perm, "  The original permutation:");

            perm_inv = Permutation.perm0_inverse3_new(N, perm);

            Permutation.perm0_print(N, perm_inv, "  The inverted permutation:");
        }

        public static void perm0_lex_next_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_LEX_NEXT_TEST tests PERM0_LEX_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 4;

            bool more;
            int[] p = new int[N];

            Console.WriteLine("");
            Console.WriteLine("PERM0_LEX_NEXT_TEST");
            Console.WriteLine("  PERM0_LEX_NEXT generates permutations in order.");
            Console.WriteLine("");

            more = false;

            for (;;)
            {
                Permutation.perm0_lex_next(N, p, ref more);

                if (!more)
                {
                    break;
                }

                Permutation.perm0_print(N, p, " ");

            }

        }

        public static void perm0_mul_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_MUL_TEST tests PERM0_MUL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 5;

            int[] p1 = new int[N];
            int[] p2 = new int[N];
            int[] p3 = new int[N];
            int seed;

            Console.WriteLine("");
            Console.WriteLine("PERM0_MUL_TEST");
            Console.WriteLine("  PERM0_MUL multiplies two permutations of (0,...,N-1).");
            Console.WriteLine("");

            seed = 123456789;

            Permutation.perm0_random(N, ref seed, ref p1);
            Permutation.perm0_random(N, ref seed, ref p2);

            Permutation.perm0_print(N, p1, "  Permutation P1:");

            Permutation.perm0_print(N, p2, "  Permutation P2:");

            Permutation.perm0_mul(N, p1, p2, ref p3);

            Permutation.perm0_print(N, p3, "  Product permutation: P3");

        }

        public static void perm0_next_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_NEXT_TEST tests PERM0_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            bool even = false;
            bool more = false;
            int n = 4;
            int[] p = new int[4];

            Console.WriteLine("");
            Console.WriteLine("PERM0_NEXT_TEST");
            Console.WriteLine("  PERM0_NEXT generates permutations of (0,...,N-1).");
            Console.WriteLine("");

            more = false;

            for (;;)
            {
                Permutation.perm0_next(n, p, ref more, ref even);

                Permutation.perm0_print(n, p, "");

                if (!more)
                {
                    break;
                }
            }
        }

        public static void perm0_next2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_NEXT2_TEST tests PERM0_NEXT2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 4;

            bool done = false;
            int[] p = new int[N];
            Permutation.PermNext2Data data = new Permutation.PermNext2Data();

            Console.WriteLine("");
            Console.WriteLine("PERM0_NEXT2_TEST");
            Console.WriteLine("  PERM0_NEXT2 generates permutations of (0,...,N-1).");
            Console.WriteLine("");

            done = true;

            for (;;)
            {
                Permutation.perm0_next2(ref data, N, p, ref done);

                if (done)
                {
                    break;
                }

                Permutation.perm0_print(N, p, "");

            }
        }

        public static void perm0_next3_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_NEXT3_TEST tests PERM0_NEXT3.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 November 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            bool more;
            int n;
            int[] p;
            int rank;

            Console.WriteLine("");
            Console.WriteLine("PERM0_NEXT3_TEST");
            Console.WriteLine("  PERM0_NEXT3 generates permutations of (0,...,N-1).");
            Console.WriteLine("");

            n = 4;
            p = new int[n];
            more = false;
            rank = 0;

            for (;;)
            {
                Permutation.perm0_next3(n, p, ref more, ref rank);

                if (!more)
                {
                    break;
                }

                Permutation.perm0_print(n, p, "");
            }
        }

        public static void perm0_random_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_RANDOM_TEST tests PERM0_RANDOM;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 4;

            int i;
            int[] p = new int[N];
            int seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("PERM0_RANDOM_TEST");
            Console.WriteLine("  PERM0_RANDOM produces a random permutation of (0,...,N-1);");
            Console.WriteLine("  For this test, N = " + N + "");
            Console.WriteLine("");

            for (i = 1; i <= 5; i++)
            {
                Permutation.perm0_random(N, ref seed, ref p);
                Permutation.perm0_print(N, p, "");
            }
        }

        public static void perm0_random2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_RANDOM2_TEST tests PERM0_RANDOM2;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 4;

            int i;
            int[] p = new int[N];
            int seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("PERM0_RANDOM2_TEST");
            Console.WriteLine("  PERM0_RANDOM2 produces a random permutation of (0,...,N-1);");
            Console.WriteLine("  For this test, N = " + N + "");
            Console.WriteLine("");

            for (i = 1; i <= 5; i++)
            {
                Permutation.perm0_random2(N, ref seed, ref p);
                Permutation.perm0_print(N, p, "");
            }
        }

        public static void perm0_rank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_RANK_TEST tests PERM0_RANK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] invers = new int[4];
            int n = 4;
            int[] p = { 0, 3, 1, 2 };
            int rank;

            Console.WriteLine("");
            Console.WriteLine("PERM0_RANK_TEST");
            Console.WriteLine("  PERM0_RANK ranks a permutation of (0,...,N-1).");

            Permutation.perm0_print(n, p, "  The permutation:");

            rank = Permutation.perm0_rank(n, p, ref invers);

            Console.WriteLine("");
            Console.WriteLine("  The rank is " + rank + "");
        }

        public static void perm0_sign_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_SIGN_TEST tests PERM0_SIGN.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 4;

            int i;
            bool more;
            int[] p = new int[N];
            int rank;
            int p_sign;

            Console.WriteLine("");
            Console.WriteLine("PERM0_SIGN_TEST");
            Console.WriteLine("  PERM0_SIGN computes the sign of a permutation of (0,...,N-1).");
            Console.WriteLine("");
            Console.WriteLine("  RANK  SIGN  Permutation");
            Console.WriteLine("");

            more = false;
            rank = 0;

            for (;;)
            {
                Permutation.perm0_lex_next(N, p, ref more);

                p_sign = Permutation.perm0_sign(N, p);
                if (!more)
                {
                    break;
                }

                string cout = rank.ToString().PadLeft(4) + "  "
                                                         + p_sign.ToString().PadLeft(4) + "  ";

                for (i = 0; i < N; i++)
                {
                    cout += p[i].ToString().PadLeft(4) + "  ";
                }

                Console.WriteLine(cout);

                rank = rank + 1;
            }
        }

        public static void perm0_to_equiv_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_TO_EQUIV_TEST tests PERM0_TO_EQUIV.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] a = new int[9];
            int[] jarray = new int[9];
            int n = 9;
            int npart = 0;
            int[] p = { 1, 2, 8, 5, 6, 7, 4, 3, 0 };

            Console.WriteLine("");
            Console.WriteLine("PERM0_TO_EQUIV_TEST");
            Console.WriteLine("  PERM0_TO_EQUIV returns the set partition");
            Console.WriteLine("  or equivalence classes determined by a");
            Console.WriteLine("  permutation of (0,...,N-1).");

            Permutation.perm0_print(n, p, "  The input permutation:");

            Permutation.perm0_to_equiv(n, p, ref npart, ref jarray, ref a);

            Equiv.equiv_print(n, a, "  The partition:");

        }

        public static void perm0_to_inversion_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_TO_INVERSION_TEST tests PERM0_TO_INVERSION.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 5;

            int[] ins = new int [N];
            int[] perm = { 2, 4, 0, 3, 1 };
            int[] perm2 = new int [N];

            Console.WriteLine("");
            Console.WriteLine("PERM0_TO_INVERSION_TEST");
            Console.WriteLine("  PERM0_TO_INVERSION: permutation (0,...,N-1) to inversion.");
            Console.WriteLine("");

            typeMethods.i4vec1_print(N, perm, "  The permutation:");

            Permutation.perm0_to_inversion(N, perm, ins);

            typeMethods.i4vec1_print(N, ins, "  The inversion sequence:");

            Permutation.inversion_to_perm0(N, ins, ref perm2);

            typeMethods.i4vec1_print(N, perm2, "  The recovered permutation:");
        }

        public static void perm0_to_ytb_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_TO_YTB_TEST tests PERM0_TO_YTB.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 7;

            int[] a = new int[N];
            int[] lambda = new int[N];
            int[] p = { 7, 2, 4, 1, 5, 3, 6 };

            Console.WriteLine("");
            Console.WriteLine("PERM0_TO_YTB_TEST");
            Console.WriteLine("  PERM0_TO_YTB converts a permutation of (0,...,N-1) to a");
            Console.WriteLine("  Young tableau.");

            Permutation.perm0_print(N, p, "  The permutation:");

            Permutation.perm0_to_ytb(N, p, ref lambda, ref a);

            YoungTableau.ytb_print(N, a, "  The Young tableau:");

        }

        public static void perm0_unrank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_UNRANK_TEST tests PERM0_UNRANK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n = 4;
            int[] p = new int[4];
            int rank = 6;

            Console.WriteLine("");
            Console.WriteLine("PERM0_UNRANK_TEST");
            Console.WriteLine("  PERM0_UNRANK, given a rank, computes the");
            Console.WriteLine("  corresponding permutation of (0,...,N-1).");
            Console.WriteLine("");
            Console.WriteLine("  The requested rank is " + rank + "");

            Permutation.perm0_unrank(n, rank, p);

            Permutation.perm0_print(n, p, "  The permutation:");
        }

        public static void perm1_canon_to_cycle_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM1_CANON_TO_CYCLE_TEST tests PERM1_CANON_TO_CYCLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n = 6;
            int[] p1 = { 4, 5, 2, 1, 6, 3 };
            int[] p2 = new int[6];

            Console.WriteLine("");
            Console.WriteLine("PERM1_CANON_TO_CYCLE_TEST");
            Console.WriteLine("  PERM1_CANON_TO_CYCLE converts a permutation of (1,...,N) from");
            Console.WriteLine("  canonical to cycle form.");

            Permutation.perm1_print(n, p1, "  The permutation in canonical form:");

            Permutation.perm1_canon_to_cycle(n, p1, ref p2);

            Permutation.perm1_print(n, p2, "  The permutation in cycle form:");
        }

        public static void perm1_check_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM1_CHECK_TEST tests PERM1_CHECK.
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
            int n = 5;
            int[] p1 = { 5, 2, 3, 4, 1 };
            int[] p2 = { 4, 1, 3, 0, 2 };
            int[] p3 = { 0, 2, 1, 3, 2 };

            Console.WriteLine("");
            Console.WriteLine("PERM1_CHECK_TEST");
            Console.WriteLine("  PERM1_CHECK checks a permutation of (1,...,N).");
            Console.WriteLine("");

            Permutation.perm1_print(n, p1, "  Permutation 1:");
            Permutation.perm1_check(n, p1);

            Permutation.perm1_print(n, p2, "  Permutation 2:");
            Permutation.perm1_check(n, p2);

            Permutation.perm1_print(n, p3, "  Permutation 3:");
            Permutation.perm1_check(n, p3);
        }

        public static void perm1_cycle_to_canon_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM1_CYCLE_TO_CANON_TEST tests PERM1_CYCLE_TO_CANON.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n = 6;
            int[] p1 = { -6, 3, 1, -5, 4, -2 };
            int[] p2 = new int[6];

            Console.WriteLine("");
            Console.WriteLine("PERM1_CYCLE_TO_CANON_TEST");
            Console.WriteLine("  PERM1_CYCLE_TO_CANON converts a permutation of (1,...,N) from");
            Console.WriteLine("  cycle to canonical form.");

            Permutation.perm1_print(n, p1, "  The permutation in cycle form:");

            Permutation.perm1_cycle_to_canon(n, p1, ref p2);

            Permutation.perm1_print(n, p2, "  The permutation in canonical form:");
        }

        public static void perm1_cycle_to_index_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM1_CYCLE_TO_INDEX_TEST tests PERM1_CYCLE_TO_INDEX.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n = 9;
            int[] p1 = { 2, 3, 9, 6, 7, 8, 5, 4, 1 };
            int[] p2 = new int[9];
            int[] p3 = new int[9];

            Console.WriteLine("");
            Console.WriteLine("PERM1_CYCLE_TO_INDEX_TEST");
            Console.WriteLine("  PERM1_CYCLE_TO_INDEX converts a permutation of (1,...,N) from");
            Console.WriteLine("  cycle to standard index form.");

            Permutation.perm1_print(n, p1, "  The standard index form permutation:");

            Permutation.perm1_index_to_cycle(n, p1, p2);

            Permutation.perm1_print(n, p2, "  The permutation in cycle form:");

            Permutation.perm1_cycle_to_index(n, p2, p3);

            Permutation.perm1_print(n, p3, "  The standard index form permutation:");
        }

        public static void perm1_index_to_cycle_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM1_INDEX_TO_CYCLE_TEST tests PERM1_INDEX_TO_CYCLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n = 9;
            int[] p1 = { 2, 3, 9, 6, 7, 8, 5, 4, 1 };
            int[] p2 = new int[9];
            int[] p3 = new int[9];

            Console.WriteLine("");
            Console.WriteLine("PERM1_INDEX_TO_CYCLE_TEST");
            Console.WriteLine("  PERM1_INDEX_TO_CYCLE converts a permutation of (1,...,N) from");
            Console.WriteLine("  standard index to cycle form.");

            Permutation.perm1_print(n, p1, "  The standard index form permutation:");

            Permutation.perm1_index_to_cycle(n, p1, p2);

            Permutation.perm1_print(n, p2, "  The permutation in cycle form:");

            Permutation.perm1_cycle_to_index(n, p2, p3);

            Permutation.perm1_print(n, p3, "  The standard index form permutation:");
        }

        public static void perm1_print_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM1_PRINT_TEST tests PERM1_PRINT.
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
            int[] p = { 7, 2, 4, 1, 5, 3, 6 };

            Console.WriteLine("");
            Console.WriteLine("PERM1_PRINT_TEST");
            Console.WriteLine("  PERM1_PRINT prints a permutation of (1,...,N).");

            Permutation.perm1_print(n, p, "  The 1-based permutation:");
        }
    }
}