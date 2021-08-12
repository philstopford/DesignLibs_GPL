using System;
using Burkardt;
using Burkardt.SubsetNS;
using Burkardt.Types;

namespace SubsetTestNS
{
    public static class KsubTest
    {
        public static void ksub_next_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_NEXT_TEST tests KSUB_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] a = new int[3];
            int i;
            int k;
            int m;
            int m2;
            bool more;
            int n;
            int rank;

            Console.WriteLine("");
            Console.WriteLine("KSUB_NEXT_TEST");
            Console.WriteLine("  KSUB_NEXT generates all K subsets of an N set");
            Console.WriteLine("  in lexicographic order.");
            Console.WriteLine("");

            rank = 0;

            n = 5;
            k = 3;
            for (i = 0; i < k; i++)
            {
                a[i] = 0;
            }

            more = false;
            m = 0;
            m2 = 0;

            for (;;)
            {
                Ksub.ksub_next(n, k, ref a, ref more, ref m, ref m2);

                rank = rank + 1;

                string cout = rank.ToString().PadLeft(4) + "    ";
                for (i = 0; i < k; i++)
                {
                    cout += a[i].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);

                if (!more)
                {
                    break;
                }

            }
        }

        public static void ksub_next2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_NEXT2_TEST tests KSUB_NEXT2.
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
            int K = 3;

            int[] a = new int[K];
            int i;
            int i_in;
            int i_out;
            bool more;
            int n = 5;
            int rank;

            Console.WriteLine("");
            Console.WriteLine("KSUB_NEXT2_TEST");
            Console.WriteLine("  KSUB_NEXT2 generates the next K subset of an");
            Console.WriteLine("  N set by the revolving door method.");
            Console.WriteLine("");
            Console.WriteLine("Rank  Subset  Added  Removed");
            Console.WriteLine("");
            //
            //  KSUB_NEXT2 does not have a good way of stopping.  
            //  We will save the starting subset, and stop when the
            //  new subset is the same as the starting one.
            //
            i_in = 0;
            i_out = 0;
            rank = 0;

            typeMethods.i4vec_indicator1(K, ref a);

            for (;;)
            {
                rank = rank + 1;
                string cout = rank.ToString().PadLeft(2) + "  ";
                for (i = 0; i < K; i++)
                {
                    cout += a[i].ToString().PadLeft(2) + "  ";
                }

                cout += "   ";
                cout += i_in.ToString().PadLeft(2) + "  ";
                Console.WriteLine(cout + i_out.ToString().PadLeft(2) + "");

                Ksub.ksub_next2(n, K, ref a, ref i_in, ref i_out);

                more = false;

                for (i = 1; i <= K; i++)
                {
                    if (a[i - 1] != i)
                    {
                        more = true;
                    }
                }

                if (!more)
                {
                    break;
                }

            }
        }

        public static void ksub_next3_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_NEXT3_TEST tests KSUB_NEXT3.
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
            int K = 3;

            int[] a = new int[K];
            int i;
            int i_in = 0;
            int i_out = 0;
            bool more;
            int n = 5;
            int rank;

            Console.WriteLine("");
            Console.WriteLine("KSUB_NEXT3_TEST");
            Console.WriteLine("  KSUB_NEXT3 generates all K subsets of an N set");
            Console.WriteLine("  using the revolving door method.");
            Console.WriteLine("");
            Console.WriteLine("Rank    Subset  Added Removed");
            Console.WriteLine("");

            rank = 0;
            more = false;

            for (;;)
            {
                Ksub.ksub_next3(n, K, ref a, ref more, ref i_in, ref i_out);

                rank = rank + 1;
                string cout = rank.ToString().PadLeft(4) + "  ";
                for (i = 0; i < K; i++)
                {
                    cout += a[i].ToString().PadLeft(2) + "  ";
                }

                cout += "   ";
                cout += i_in.ToString().PadLeft(2) + "  ";
                Console.WriteLine(cout + i_out.ToString().PadLeft(2) + "");

                if (!more)
                {
                    break;
                }

            }
        }

        public static void ksub_next4_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_NEXT4_TEST tests KSUB_NEXT4.
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
            int K = 3;

            int[] a = new int[K];
            bool done;
            int i;
            int n = 5;
            int rank;

            Console.WriteLine("");
            Console.WriteLine("KSUB_NEXT4_TEST");
            Console.WriteLine("  KSUB_NEXT4 generates K subsets of an N set.");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("  K = " + K + "");
            Console.WriteLine("");
            Console.WriteLine("Rank    Subset");
            Console.WriteLine("");

            done = true;
            rank = 0;

            for (;;)
            {
                Ksub.ksub_next4(n, K, ref a, ref done);

                if (done)
                {
                    break;
                }

                rank = rank + 1;
                string cout = rank.ToString().PadLeft(4) + "  ";
                cout += "  ";
                for (i = 0; i < K; i++)
                {
                    cout += a[i].ToString().PadLeft(4) + "  ";
                }

                Console.WriteLine(cout);

            }
        }

        public static void ksub_random_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_RANDOM_TEST tests KSUB_RANDOM.
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
            int K = 3;

            int[] a = new int[K];
            int i;
            int j;
            int n = 5;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("KSUB_RANDOM_TEST");
            Console.WriteLine("  KSUB_RANDOM generates a random K subset of an N set.");
            Console.WriteLine("  Set size is N =    " + n + "");
            Console.WriteLine("  Subset size is K = " + K + "");
            Console.WriteLine("");

            seed = 123456789;

            for (i = 1; i <= 5; i++)
            {
                string cout = "";
                Ksub.ksub_random(n, K, ref seed, ref a);
                for (j = 0; j < K; j++)
                {
                    cout += "  " + a[j].ToString().PadLeft(3);
                }

                Console.WriteLine(cout);
            }
        }

        public static void ksub_random2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_RANDOM2_TEST tests KSUB_RANDOM2.
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
            int K = 3;

            int[] a = new int[K];
            int i;
            int j;
            int n = 5;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("KSUB_RANDOM2_TEST");
            Console.WriteLine("  KSUB_RANDOM2 generates a random K subset of an N set.");
            Console.WriteLine("  Set size is N =    " + n + "");
            Console.WriteLine("  Subset size is K = " + K + "");
            Console.WriteLine("");

            seed = 123456789;

            for (i = 1; i <= 5; i++)
            {
                string cout = "";
                Ksub.ksub_random2(n, K, ref seed, ref a);
                for (j = 0; j < K; j++)
                {
                    cout += "  " + a[j].ToString().PadLeft(3);
                }

                Console.WriteLine(cout);
            }
        }

        public static void ksub_random3_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_RANDOM3_TEST tests KSUB_RANDOM3.
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
            int K = 3;
            int N = 5;

            int[] a = new int[N];
            int i;
            int j;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("KSUB_RANDOM3_TEST");
            Console.WriteLine("  KSUB_RANDOM3 generates a random K subset of an N set.");
            Console.WriteLine("  Set size is N =    " + N + "");
            Console.WriteLine("  Subset size is K = " + K + "");
            Console.WriteLine("");

            seed = 123456789;

            for (i = 1; i <= 10; i++)
            {
                string cout = "";
                Ksub.ksub_random3(N, K, ref seed, ref a);
                for (j = 0; j < N; j++)
                {
                    cout += "  " + a[j].ToString().PadLeft(3);
                }

                Console.WriteLine(cout);
            }
        }

        public static void ksub_random4_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_RANDOM4_TEST tests KSUB_RANDOM4.
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
            int K = 3;
            int N = 5;

            int[] a = new int[N];
            int i;
            int j;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("KSUB_RANDOM4_TEST");
            Console.WriteLine("  KSUB_RANDOM4 generates a random K subset of an N set.");
            Console.WriteLine("  Set size is N =    " + N + "");
            Console.WriteLine("  Subset size is K = " + K + "");
            Console.WriteLine("");

            seed = 123456789;

            for (i = 1; i <= 10; i++)
            {
                string cout = "";
                Ksub.ksub_random4(N, K, ref seed, ref a);
                for (j = 0; j < K; j++)
                {
                    cout += "  " + a[j].ToString().PadLeft(3);
                }

                Console.WriteLine(cout);
            }
        }

        public static void ksub_random5_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_RANDOM5_TEST tests KSUB_RANDOM5.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 June 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] a;
            int i;
            int j;
            int k = 5;
            int n = 52;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("KSUB_RANDOM5_TEST");
            Console.WriteLine("  KSUB_RANDOM5 generates a random K subset of an N set.");
            Console.WriteLine("  Set size is N =    " + n + "");
            Console.WriteLine("  Subset size is K = " + k + "");
            Console.WriteLine("");

            seed = 123456789;

            for (i = 1; i <= 5; i++)
            {
                string cout = "";
                a = Ksub.ksub_random5(n, k, ref seed);
                for (j = 0; j < k; j++)
                {
                    cout += "  " + a[j].ToString().PadLeft(3);
                }

                Console.WriteLine(cout);
            }
        }

        public static void ksub_rank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_RANK_TEST tests KSUB_RANK.
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
            int K = 3;

            int[] a = { 1, 3, 5 };
            int i;
            int rank = 0;

            Console.WriteLine("");
            Console.WriteLine("KSUB_RANK_TEST");
            Console.WriteLine("  KSUB_RANK: determine the rank of a K subset of an N set.");
            Console.WriteLine("");
            Console.WriteLine("  For N = " + N + "");
            Console.WriteLine("  and K = " + K + "");
            Console.WriteLine("  the subset is:");
            Console.WriteLine("");

            string cout = "";
            for (i = 0; i < K; i++)
            {
                cout += a[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);

            Ksub.ksub_rank(K, a, ref rank);

            Console.WriteLine("");
            Console.WriteLine("  The computed rank is " + rank + "");
        }

        public static void ksub_to_comp_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_TO_COMP_TEST tests KSUB_TO_COMP.
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
            int[] as_ = new int[4];
            int i = 0;
            int j = 0;
            int kc = 0;
            int ks = 0;
            int nc = 0;
            int ns = 0;
            int seed = 0;

            Console.WriteLine("");
            Console.WriteLine("KSUB_TO_COMP_TEST");
            Console.WriteLine("  KSUB_TO_COMP returns the composition corresponding to a K subset.");

            nc = 10;
            kc = 5;
            seed = 123456789;

            for (i = 1; i <= 5; i++)
            {
                Console.WriteLine("");

                Comp.comp_random(nc, kc, ref seed, ref ac);
                string cout = "  COMP:";
                for (j = 0; j < kc; j++)
                {
                    cout += ac[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);

                Comp.comp_to_ksub(nc, kc, ac, ref ns, ref ks, ref as_);
                cout = "  KSUB:";
                for (j = 0; j < ks; j++)
                {
                    cout += as_[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);

                Ksub.ksub_to_comp(ns, ks, as_, ref nc, ref kc, ref ac);
                cout = "  COMP:";
                for (j = 0; j < kc; j++)
                {
                    cout += ac[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

        public static void ksub_to_compnz_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_TO_COMPNZ_TEST tests KSUB_TO_COMPNZ.
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
            int i = 0;
            int j = 0;
            int kc = 0;
            int ks = 0;
            int nc = 0;
            int ns = 0;
            int seed = 0;

            Console.WriteLine("");
            Console.WriteLine("KSUB_TO_COMPNZ_TEST");
            Console.WriteLine("  KSUB_TO_COMPNZ returns the nonzero composition ");
            Console.WriteLine("  corresponding to a K subset.");

            nc = 10;
            kc = 5;
            seed = 123456789;

            for (i = 1; i <= 5; i++)
            {
                Console.WriteLine("");

                Comp.compnz_random(nc, kc, ref seed, ref ac);
                string cout = "  COMPNZ:";
                for (j = 0; j < kc; j++)
                {
                    cout += ac[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);

                Comp.compnz_to_ksub(nc, kc, ac, ref ns, ref ks, ref as_);
                cout = "  KSUB:  ";
                for (j = 0; j < ks; j++)
                {
                    cout += as_[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);

                Ksub.ksub_to_compnz(ns, ks, as_, ref nc, ref kc, ref ac);
                cout = "  COMPNZ:";
                for (j = 0; j < kc; j++)
                {
                    cout += ac[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

        public static void ksub_unrank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KSUB_UNRANK_TEST tests KSUB_UNRANK.
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
            int K = 3;

            int i;
            int[] a = new int[K];
            int n = 5;
            int rank;

            rank = 8;

            Console.WriteLine("");
            Console.WriteLine("KSUB_UNRANK_TEST");
            Console.WriteLine("  KSUB_UNRANK: find the K-subset of an N set");
            Console.WriteLine("  of a given rank.");
            Console.WriteLine("");
            Console.WriteLine("  For N = " + n + "");
            Console.WriteLine("  and K = " + K + "");
            Console.WriteLine("  and the desired rank is " + rank + "");

            Ksub.ksub_unrank(K, rank, ref a);

            Console.WriteLine("");
            Console.WriteLine("  The subset of the given rank is:");
            Console.WriteLine("");

            string cout = "";
            for (i = 0; i < K; i++)
            {
                cout += a[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }

    }
}