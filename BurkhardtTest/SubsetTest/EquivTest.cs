using System;
using Burkardt.SubsetNS;

namespace SubsetTestNS
{
    public static class EquivTest
    {
        public static void equiv_next_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EQUIV_NEXT_TEST tests EQUIV_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 4;

            int[] a = new int[N];
            int i;
            int[] jarray = new int[N];
            bool more = false;
            int npart = 0;
            int rank = 0;

            Console.WriteLine("");
            Console.WriteLine("EQUIV_NEXT_TEST");
            Console.WriteLine("  EQUIV_NEXT generates all partitions of a set.");
            Console.WriteLine("");
            Console.WriteLine("  Rank//element:");
            Console.WriteLine("");
            string cout = "      ";
            for (i = 1; i <= N; i++)
            {
                cout += i.ToString().PadLeft(2) + "  ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("");

            rank = 0;
            more = false;

            for (;;)
            {
                Equiv.equiv_next(N, ref npart, ref jarray, ref a, ref more);

                rank = rank + 1;

                cout = "  "
                       + rank.ToString().PadLeft(2) + "  ";
                for (i = 0; i < N; i++)
                {
                    cout += a[i].ToString().PadLeft(2) + "  ";
                }

                Console.WriteLine(cout);

                if (!more)
                {
                    break;
                }

            }
        }

        public static void equiv_next2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EQUIV_NEXT2_TEST tests EQUIV_NEXT2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 4;

            int[] a = new int[N];
            bool done;
            int i;
            int rank;

            Console.WriteLine("");
            Console.WriteLine("EQUIV_NEXT2_TEST");
            Console.WriteLine("  EQUIV_NEXT2 generates all partitions of a set.");
            Console.WriteLine("");
            Console.WriteLine("  Rank//element:");
            Console.WriteLine("");
            string cout = "      ";
            for (i = 1; i <= N; i++)
            {
                cout += i.ToString().PadLeft(2) + "  ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("");

            rank = 0;
            done = true;

            for (;;)
            {
                Equiv.equiv_next2(ref done, ref a, N);

                if (done)
                {
                    break;
                }

                rank = rank + 1;

                cout = "  "
                       + rank.ToString().PadLeft(2) + "  ";
                for (i = 0; i < N; i++)
                {
                    cout += a[i].ToString().PadLeft(2) + "  ";
                }

                Console.WriteLine(cout);

            }

        }

        public static void equiv_print_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EQUIV_PRINT_TEST tests EQUIV_PRINT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] a = new int[4];
            int i;
            int n = 4;
            int npart = 0;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("EQUIV_PRINT_TEST");
            Console.WriteLine("  EQUIV_PRINT prints a set partition.");

            seed = 123456789;

            for (i = 1; i <= 5; i++)
            {
                Equiv.equiv_random(n, ref seed, ref npart, ref a);

                Equiv.equiv_print(n, a, "  The partition:");
            }
        }

        public static void equiv_print2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EQUIV_PRINT2_TEST tests EQUIV_PRINT2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] a = new int[4];
            int i;
            int n = 4;
            int npart = 0;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("EQUIV_PRINT2_TEST");
            Console.WriteLine("  EQUIV_PRINT2 prints a set partition.");

            seed = 123456789;

            for (i = 1; i <= 5; i++)
            {
                Equiv.equiv_random(n, ref seed, ref npart, ref a);

                Equiv.equiv_print2(n, a, "  The partition:");
            }
        }

        public static void equiv_random_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EQUIV_RANDOM_TEST tests EQUIV_RANDOM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 May 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] a = new int[4];
            int i;
            int n = 4;
            int npart = 0;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("EQUIV_RANDOM_TEST");
            Console.WriteLine("  EQUIV_RANDOM selects a random set partition.");

            seed = 123456789;

            for (i = 1; i <= 5; i++)
            {
                Equiv.equiv_random(n, ref seed, ref npart, ref a);

                Equiv.equiv_print(n, a, "  The partition:");
            }

        }

    }
}