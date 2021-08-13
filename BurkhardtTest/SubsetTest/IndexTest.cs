using System;
using Burkardt.Types;

namespace SubsetTestNS
{
    using Index = Burkardt.IndexNS.Index;

    public static class IndexTest
    {
        public static void index_box_next_2d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_BOX_NEXT_2D_TEST tests INDEX_BOX_NEXT_2D.
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
            int i = 0;
            int j = 0;
            bool more = false;
            int n1 = 5;
            int n2 = 3;
            int n;

            Console.WriteLine("");
            Console.WriteLine("INDEX_BOX_NEXT_2D_TEST");
            Console.WriteLine("  INDEX_BOX_NEXT_2D produces IJ indices that");
            Console.WriteLine("  lie on the surface of a box in 2D.");
            Console.WriteLine("");
            Console.WriteLine("  The box has logical dimensions:");
            Console.WriteLine(n1.ToString().PadLeft(3) + "  "
                                                       + n2.ToString().PadLeft(3) + "");
            Console.WriteLine("");
            Console.WriteLine("   #    I   J");
            Console.WriteLine("");

            more = false;
            n = 0;

            for (;;)
            {
                Index.index_box_next_2d(n1, n2, ref i, ref j, ref more);

                if (!more)
                {
                    break;
                }

                n = n + 1;
                Console.WriteLine(n.ToString().PadLeft(3) + "  "
                                                          + i.ToString().PadLeft(3) + "  "
                                                          + j.ToString().PadLeft(3) + "");
            }
        }

        public static void index_box_next_3d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_BOX_NEXT_3D_TEST tests INDEX_BOX_NEXT_3D.
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
            int i = 0;
            int j = 0;
            int k = 0;
            bool more = false;
            int n1 = 5;
            int n2 = 3;
            int n3 = 4;
            int n = 0;

            Console.WriteLine("");
            Console.WriteLine("INDEX_BOX_NEXT_3D_TEST");
            Console.WriteLine("  INDEX_BOX_NEXT_3D produces IJK indices that");
            Console.WriteLine("  lie on the surface of a box.");
            Console.WriteLine("");
            Console.WriteLine("  The box has logical dimensions:");
            Console.WriteLine(n1.ToString().PadLeft(3) + "  "
                                                       + n2.ToString().PadLeft(3) + "  "
                                                       + n3.ToString().PadLeft(3) + "");
            Console.WriteLine("");
            Console.WriteLine("   #    I   J   K");
            Console.WriteLine("");

            more = false;
            n = 0;

            for (;;)
            {
                Index.index_box_next_3d(n1, n2, n3, ref i, ref j, ref k, ref more);

                if (!more)
                {
                    break;
                }

                n = n + 1;
                Console.WriteLine(n.ToString().PadLeft(3) + "  "
                                                          + i.ToString().PadLeft(3) + "  "
                                                          + j.ToString().PadLeft(3) + "  "
                                                          + k.ToString().PadLeft(3) + "");

            }
        }

        public static void index_box2_next_2d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_BOX2_NEXT_2D_TEST tests INDEX_BOX2_NEXT_2D.
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
            int i = 0;
            int ic = 10;
            int j = 0;
            int jc = 20;
            bool more = false;
            int n1 = 4;
            int n2 = 3;
            int n = 0;

            Console.WriteLine("");
            Console.WriteLine("INDEX_BOX2_NEXT_2D_TEST");
            Console.WriteLine("  INDEX_BOX2_NEXT_2D produces IJ indices that");
            Console.WriteLine("  lie on the surface of a box2 in 2D.");
            Console.WriteLine("");
            Console.WriteLine("  The box has half-widths:");
            Console.WriteLine(n1.ToString().PadLeft(3) + "  "
                                                       + n2.ToString().PadLeft(3) + "");
            Console.WriteLine("");
            Console.WriteLine("  and has center cell:");
            Console.WriteLine(ic.ToString().PadLeft(3) + "  "
                                                       + jc.ToString().PadLeft(3) + "");
            Console.WriteLine("");
            Console.WriteLine("   #    I   J");
            Console.WriteLine("");

            more = false;
            n = 0;

            for (;;)
            {
                Index.index_box2_next_2d(n1, n2, ic, jc, ref i, ref j, ref more);

                if (!more)
                {
                    break;
                }

                n = n + 1;
                Console.WriteLine(n.ToString().PadLeft(3) + "  "
                                                          + i.ToString().PadLeft(3) + "  "
                                                          + j.ToString().PadLeft(3) + "");
            }
        }

        public static void index_box2_next_3d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_BOX2_NEXT_3D_TEST tests INDEX_BOX2_NEXT_3D.
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
            int i = 0;
            int ic = 10;
            int j = 0;
            int jc = 20;
            int k = 0;
            int kc = 30;
            bool more;
            int n1 = 5;
            int n2 = 3;
            int n3 = 4;
            int n = 0;

            Console.WriteLine("");
            Console.WriteLine("INDEX_BOX2_NEXT_3D_TEST");
            Console.WriteLine("  INDEX_BOX2_NEXT_3D produces IJK indices that");
            Console.WriteLine("  lie on the surface of a box.");
            Console.WriteLine("");
            Console.WriteLine("  The box has half widths:");
            Console.WriteLine(n1.ToString().PadLeft(3) + "  "
                                                       + n2.ToString().PadLeft(3) + "  "
                                                       + n3.ToString().PadLeft(3) + "");
            Console.WriteLine("");
            Console.WriteLine("  and central cell:");
            Console.WriteLine(ic.ToString().PadLeft(3) + "  "
                                                       + jc.ToString().PadLeft(3) + "  "
                                                       + kc.ToString().PadLeft(3) + "");
            Console.WriteLine("");
            Console.WriteLine("  We will only print a PORTION of the data!");
            Console.WriteLine("");
            Console.WriteLine("   #    I   J   K");
            Console.WriteLine("");

            more = false;
            n = 0;

            for (;;)
            {
                Index.index_box2_next_3d(n1, n2, n3, ic, jc, kc, ref i, ref j, ref k, ref more);

                if (!more)
                {
                    break;
                }

                n = n + 1;

                if (n <= 10 || 370 <= n)
                {
                    Console.WriteLine(n.ToString().PadLeft(3) + "  "
                                                              + i.ToString().PadLeft(3) + "  "
                                                              + j.ToString().PadLeft(3) + "  "
                                                              + k.ToString().PadLeft(3) + "");
                }

            }
        }

        public static void index_next0_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_NEXT0_TEST tests INDEX_NEXT0.
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
            int N = 3;

            int[] a = new int[N];
            int hi = 3;
            int i;
            bool more = false;

            Console.WriteLine("");
            Console.WriteLine("INDEX_NEXT0_TEST");
            Console.WriteLine("  INDEX_NEXT0 generates all indices of an");
            Console.WriteLine("  array of given shape, with");
            Console.WriteLine("  lower limit 1 and given upper limit.");
            Console.WriteLine("");
            Console.WriteLine("  Number of index entries = " + N + "");
            Console.WriteLine("  Coordinate maximum HI =   " + hi + "");

            Console.WriteLine("");
            Console.WriteLine("  Index arrays:");
            Console.WriteLine("");

            more = false;

            for (;;)
            {
                Index.index_next0(N, hi, ref a, ref more);

                string cout = "";
                for (i = 0; i < N; i++)
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

        public static void index_next1_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_NEXT1_TEST tests INDEX_NEXT1.
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
            int N = 3;

            int[] a = new int[N];
            int[] hi = { 4, 2, 3 };
            int i;
            bool more = false;

            Console.WriteLine("");
            Console.WriteLine("INDEX_NEXT1_TEST");
            Console.WriteLine("  INDEX_NEXT1 generates all indices of an");
            Console.WriteLine("  array of given shape, with");
            Console.WriteLine("  lower limit 1 and given upper limits.");
            Console.WriteLine("");
            Console.WriteLine("  Number of index entries = " + N + "");

            typeMethods.i4vec1_print(N, hi, "  Coordinate maximum indices:");

            Console.WriteLine("");
            Console.WriteLine("  Index arrays:");
            Console.WriteLine("");

            more = false;

            for (;;)
            {
                Index.index_next1(N, hi, ref a, ref more);

                string cout = "";
                for (i = 0; i < N; i++)
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

        public static void index_next2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_NEXT2_TEST tests INDEX_NEXT2.
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
            int N = 3;

            int[] a = new int[N];
            int[] hi = { 11, -3, 1 };
            int i;
            int[] lo = { 10, -5, 0 };
            bool more = false;

            Console.WriteLine("");
            Console.WriteLine("INDEX_NEXT2_TEST");
            Console.WriteLine("  INDEX_NEXT2 generates all indices of an");
            Console.WriteLine("  array of given shape with given");
            Console.WriteLine("  lower and upper limits.");
            Console.WriteLine("");
            Console.WriteLine("  Number of index entries = " + N + "");
            Console.WriteLine("");
            Console.WriteLine("  Coordinate, Maximum Index");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine((i + 1).ToString().PadLeft(8) + "  "
                                                                + lo[i].ToString().PadLeft(8) + "  "
                                                                + hi[i].ToString().PadLeft(8) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("Index arrays:");
            Console.WriteLine("");

            more = false;

            for (;;)
            {
                Index.index_next2(N, lo, hi, ref a, ref more);

                string cout = "";
                for (i = 0; i < N; i++)
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

        public static void index_rank0_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_RANK0_TEST tests INDEX_RANK0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            int[] a = { 3, 1, 2 };
            int hi = 3;
            int rank;

            Console.WriteLine("");
            Console.WriteLine("INDEX_RANK0_TEST");
            Console.WriteLine("  INDEX_RANK0 ranks an index with");
            Console.WriteLine("  lower limit 1 and given upper limit.");
            Console.WriteLine("");
            Console.WriteLine("  Number of index entries = " + N + "");
            Console.WriteLine("");
            Console.WriteLine("  Coordinate maximum Index = " + hi + "");
            Console.WriteLine("");

            typeMethods.i4vec1_print(N, a, "  The index array:");

            rank = Index.index_rank0(N, hi, a);

            Console.WriteLine("");
            Console.WriteLine("  The rank of this object is " + rank + "");

        }

        public static void index_rank1_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_RANK1_TEST tests INDEX_RANK1.
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
            int N = 3;

            int[] a = { 4, 1, 2 };
            int[] hi = { 4, 2, 3 };
            int i;
            int rank;

            Console.WriteLine("");
            Console.WriteLine("INDEX_RANK1_TEST");
            Console.WriteLine("  INDEX_RANK1 ranks an index with");
            Console.WriteLine("  lower limit 1 and given upper limits.");
            Console.WriteLine("");
            Console.WriteLine("  Number of index entries = " + N + "");
            Console.WriteLine("");
            Console.WriteLine("  Coordinate, Maximum Index");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine((i + 1).ToString().PadLeft(10) + "  "
                                                                 + hi[i].ToString().PadLeft(10) + "");
            }

            typeMethods.i4vec1_print(N, a, "  The index array:");

            rank = Index.index_rank1(N, hi, a);

            Console.WriteLine("");
            Console.WriteLine("  The rank of this object is " + rank + "");

        }

        public static void index_rank2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_RANK2_TEST tests INDEX_RANK2.
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
            int N = 3;

            int[] a = { 1, 11, 5 };
            int[] hi = { 2, 11, 6 };
            int i;
            int[] lo = { 1, 10, 4 };
            int rank = 0;

            Console.WriteLine("");
            Console.WriteLine("INDEX_RANK2_TEST");
            Console.WriteLine("  INDEX_RANK2 ranks an index with given");
            Console.WriteLine("  lower and upper limits.");
            Console.WriteLine("");
            Console.WriteLine("  Number of index entries = " + N + "");
            Console.WriteLine("");
            Console.WriteLine("  Coordinate, Minimum index, Maximum Index");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine((i + 1).ToString().PadLeft(10) + "  "
                                                                 + lo[i].ToString().PadLeft(10) + "  "
                                                                 + hi[i].ToString().PadLeft(10) + "");
            }

            typeMethods.i4vec1_print(N, a, "  The index array:");

            rank = Index.index_rank2(N, lo, hi, a);

            Console.WriteLine("");
            Console.WriteLine("  The rank of this object is " + rank + "");

        }

        public static void index_unrank0_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_UNRANK0_TEST tests INDEX_UNRANK0.
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
            int N = 3;

            int[] a = new int[N];
            int hi = 3;
            int i;
            int maxrank;
            int rank;

            Console.WriteLine("");
            Console.WriteLine("INDEX_UNRANK0_TEST");
            Console.WriteLine("  INDEX_UNRANK0 unranks a multi-index.");
            Console.WriteLine("");
            Console.WriteLine("  The multi-index has dimension " + N + "");
            Console.WriteLine("");
            Console.WriteLine("  The upper limit is HI = " + hi + "");
            Console.WriteLine("");
            Console.WriteLine("  Rank, Multi-Index:");
            Console.WriteLine("");

            maxrank = (int)Math.Pow(hi, N);

            for (rank = 1; rank <= maxrank; rank++)
            {
                Index.index_unrank0(N, hi, rank, ref a);
                string cout = rank.ToString().PadLeft(3) + "  ";
                for (i = 0; i < N; i++)
                {
                    cout += a[i].ToString().PadLeft(6) + "  ";
                }

                Console.WriteLine(cout);
            }
        }

        public static void index_unrank1_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_UNRANK1_TEST tests INDEX_UNRANK1.
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
            int N = 3;

            int[] a = new int[N];
            int[] hi = { 4, 2, 3 };
            int i;
            int maxrank;
            int rank;

            Console.WriteLine("");
            Console.WriteLine("INDEX_UNRANK1_TEST");
            Console.WriteLine("  INDEX_UNRANK1 unranks a multi-index.");
            Console.WriteLine("");
            Console.WriteLine("  The multi-index has dimension " + N + "");

            typeMethods.i4vec1_print(N, hi, "  The upper limits:");

            Console.WriteLine("");
            Console.WriteLine("  Rank, Multi-Index:");
            Console.WriteLine("");

            maxrank = typeMethods.i4vec_product(N, hi);

            for (rank = 1; rank <= maxrank; rank++)
            {
                Index.index_unrank1(N, hi, rank, ref a);
                string cout = rank.ToString().PadLeft(3) + "  ";
                for (i = 0; i < N; i++)
                {
                    cout += a[i].ToString().PadLeft(6) + "  ";
                }

                Console.WriteLine(cout);
            }
        }

        public static void index_unrank2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INDEX_UNRANK2_TEST tests INDEX_UNRANK2.
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
            int N = 3;

            int[] a = new int[N];
            int[] hi = { 2, 11, 6 };
            int i;
            int[] lo = { 1, 10, 4 };
            int rank;

            Console.WriteLine("");
            Console.WriteLine("INDEX_UNRANK2_TEST");
            Console.WriteLine("  INDEX_UNRANK2 unranks a multi-index.");
            Console.WriteLine("");
            Console.WriteLine("  The multi-index has dimension " + N + "");
            Console.WriteLine("");
            Console.WriteLine("  The lower and upper limits are:");
            Console.WriteLine("");

            for (i = 0; i < N; i++)
            {
                Console.WriteLine(i.ToString().PadLeft(10) + "  "
                                                           + lo[i].ToString().PadLeft(10) + "  "
                                                           + hi[i].ToString().PadLeft(10) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Rank, Multi-Index:");
            Console.WriteLine("");

            rank = 7;

            Index.index_unrank2(N, lo, hi, rank, ref a);
            string cout = rank.ToString().PadLeft(3) + "  ";
            for (i = 0; i < N; i++)
            {
                cout += a[i].ToString().PadLeft(6) + "  ";
            }

            Console.WriteLine(cout);

        }

    }
}