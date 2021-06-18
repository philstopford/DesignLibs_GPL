using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest
{
    partial class Program
    {
        static void part_enum_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PART_ENUM_TEST tests PART_ENUM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 December 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;

            Console.WriteLine("");
            Console.WriteLine("PART_ENUM_TEST");
            Console.WriteLine("  PART_ENUM enumerates the partitions of N.");
            Console.WriteLine("");

            for (n = 0; n <= 10; n++)
            {
                Console.WriteLine("  " + n.ToString().PadLeft(2)
                    + "  " + Ranking.part_enum(n).ToString().PadLeft(6) + "");
            }
        }

        static void part_rsf_check_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PART_RSF_CHECK_TEST tests PART_RSF_CHECK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] a = new int[1];
            int[] a1 =  {
                1, 4, 4, 6
            }
            ;
            int[] a2 =  {
                1, 4, 4, 6
            }
            ;
            int[] a3 =  {
                -9, 4, 4, 16
            }
            ;
            int[] a4 =  {
                6, 4, 4, 1
            }
            ;
            int[] a5 =  {
                1, 4, 5, 6
            }
            ;
            int[] a6 =  {
                1, 4, 4, 6
            }
            ;
            bool check;
            int n = 0;
            int npart = 0;
            int test;

            Console.WriteLine("");
            Console.WriteLine("PART_RSF_CHECK TEST");
            Console.WriteLine("  PART_RSF_CHECK checks a reverse standard form partition.");

            for (test = 1; test <= 6; test++)
            {
                if (test == 1)
                {
                    n = 0;
                    npart = 4;
                    a = typeMethods.i4vec_copy_new(npart, a1);
                }
                else if (test == 2)
                {
                    n = 15;
                    npart = 0;
                    a = typeMethods.i4vec_copy_new(4, a2);
                }
                else if (test == 3)
                {
                    n = 15;
                    npart = 4;
                    a = typeMethods.i4vec_copy_new(npart, a3);
                }
                else if (test == 4)
                {
                    n = 15;
                    npart = 4;
                    a = typeMethods.i4vec_copy_new(npart, a4);
                }
                else if (test == 5)
                {
                    n = 15;
                    npart = 4;
                    a = typeMethods.i4vec_copy_new(npart, a5);
                }
                else if (test == 6)
                {
                    n = 15;
                    npart = 4;
                    a = typeMethods.i4vec_copy_new(npart, a6);
                }

                Console.WriteLine("");
                Console.WriteLine("  Partition in RSF form.");
                Console.WriteLine("  Partition of N = " + n + "");
                Console.WriteLine("  Number of parts NPART = " + npart + "");
                typeMethods.i4vec_transpose_print(npart, a, "");
                check = Ranking.part_rsf_check(n, npart, a);
                Console.WriteLine("  Check = " + check + "");
            }
        }

        static void part_sf_check_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PART_SF_CHECK_TEST tests PART_SF_CHECK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] a = new int[1];
            int[] a1 =  {
                6, 4, 4, 1
            }
            ;
            int[] a2 =  {
                6, 4, 4, 1
            }
            ;
            int[] a3 =  {
                16, 4, 4, -9
            }
            ;
            int[] a4 =  {
                1, 4, 4, 6
            }
            ;
            int[] a5 =  {
                6, 4, 5, 1
            }
            ;
            int[] a6 =  {
                6, 4, 4, 1
            }
            ;
            bool check;
            int n = 0;
            int npart = 0;
            int test;

            Console.WriteLine("");
            Console.WriteLine("PART_SF_CHECK TEST");
            Console.WriteLine("  PART_SF_CHECK checks a standard form partition.");

            for (test = 1; test <= 6; test++)
            {
                if (test == 1)
                {
                    n = 0;
                    npart = 4;
                    a = typeMethods.i4vec_copy_new(npart, a1);
                }
                else if (test == 2)
                {
                    n = 15;
                    npart = 0;
                    a = typeMethods.i4vec_copy_new(4, a2);
                }
                else if (test == 3)
                {
                    n = 15;
                    npart = 4;
                    a = typeMethods.i4vec_copy_new(npart, a3);
                }
                else if (test == 4)
                {
                    n = 15;
                    npart = 4;
                    a = typeMethods.i4vec_copy_new(npart, a4);
                }
                else if (test == 5)
                {
                    n = 15;
                    npart = 4;
                    a = typeMethods.i4vec_copy_new(npart, a5);
                }
                else if (test == 6)
                {
                    n = 15;
                    npart = 4;
                    a = typeMethods.i4vec_copy_new(npart, a6);
                }

                Console.WriteLine("");
                Console.WriteLine("  Partition in SF form.");
                Console.WriteLine("  Partition of N = " + n + "");
                Console.WriteLine("  Number of parts NPART = " + npart + "");
                typeMethods.i4vec_transpose_print(npart, a, "");
                check = Ranking.part_sf_check(n, npart, a);
                Console.WriteLine("  Check = " + check + "");
            }
        }

        static void part_successor_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PART_SUCCESSOR_TEST tests PART_SUCCESSOR.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 December 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int n = 8;
            int npart = 0;
            int rank;
            int rank_old;
            int[] t;

            Console.WriteLine("");
            Console.WriteLine("PART_SUCCESSOR_TEST");
            Console.WriteLine("  PART_SUCCESSOR produces partitions of N,");

            t = new int[n];

            rank = -1;

            for (;;)
            {
                rank_old = rank;

                Ranking.part_successor(n, ref npart, ref t, ref rank);

                if (rank <= rank_old)
                {
                    break;
                }

                string cout = "  " + rank.ToString().PadLeft(4);
                for (i = 0; i < npart; i++)
                {
                    cout += "  " + t[i].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

        static void part_sf_conjugate_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PART_SF_CONJUGATE_TEST tests PART_SF_CONJUGATE.
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
            int[] b;
            int i;
            int n = 8;
            int npart = 0;
            int npartb = 0;
            int rank = 0;
            int rank_old;
            int[] t;

            Console.WriteLine("");
            Console.WriteLine("PART_SF_CONJUGATE_TEST");
            Console.WriteLine("  PART_SF_CONJUGATE produces the conjugate of a partition.");
            Console.WriteLine("");
            Console.WriteLine("  Partitions of N = " + n + "");
            //
            //  List.
            //
            t = new int[n];

            rank = -1;

            for (;;)
            {
                rank_old = rank;

                Ranking.part_successor(n, ref npart, ref t, ref rank);

                if (rank <= rank_old)
                {
                    break;
                }

                string cout = "  " + rank.ToString().PadLeft(4);
                for (i = 0; i < npart; i++)
                {
                    cout += "  " + t[i].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);

                b = Ranking.part_sf_conjugate(n, npart, t, ref npartb);
                typeMethods.i4vec_transpose_print(npartb, b, "  Con:");
            }
        }

        static void part_sf_majorize_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PART_SF_MAJORIZE_TEST tests PART_SF_MAJORIZE.
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
            int N = 8;

            int[] a =  {
                2, 2, 2, 1, 1, 0, 0, 0
            }
            ;
            int[] b =  {
                3, 1, 1, 1, 1, 1, 0, 0
            }
            ;
            int[] c =  {
                2, 2, 1, 1, 1, 1, 0, 0
            }
            ;
            int n = N;
            int nparta = 5;
            int npartb = 6;
            int npartc = 6;
            int result;

            Console.WriteLine("");
            Console.WriteLine("PART_SF_MAJORIZE_TEST");
            Console.WriteLine("  PART_SF_MAJORIZE determines if one partition");
            Console.WriteLine("  majorizes another.");
            Console.WriteLine("");
            Console.WriteLine("  Partitions of N = " + n + "");
            Console.WriteLine("");
            typeMethods.i4vec_transpose_print(nparta, a, "  A: ");
            typeMethods.i4vec_transpose_print(npartb, b, "  B: ");
            typeMethods.i4vec_transpose_print(npartc, c, "  C: ");

            result = Ranking.part_sf_majorize(n, nparta, a, npartb, b);
            Console.WriteLine("");
            Console.WriteLine("  A compare B: " + result + "");
            result = Ranking.part_sf_majorize(n, npartb, b, npartc, c);
            Console.WriteLine("  B compare C: " + result + "");
            result = Ranking.part_sf_majorize(n, npartc, c, nparta, a);
            Console.WriteLine("  C compare A: " + result + "");
            result = Ranking.part_sf_majorize(n, npartc, c, npartc, c);
            Console.WriteLine("  C compare C: " + result + "");
        }

        static void part_table_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PART_TABLE_TEST tests PART_TABLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 December 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int maxn = 10;
            int[] p;

            Console.WriteLine("");
            Console.WriteLine("PART_TABLE_TEST");
            Console.WriteLine("  PART_TABLE tabulates partitions of N.");

            p = Ranking.part_table(maxn);

            Console.WriteLine("");

            for (i = 0; i <= maxn; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                    + "  " + p[i].ToString().PadLeft(4) + "");
            }
        }
    }
}