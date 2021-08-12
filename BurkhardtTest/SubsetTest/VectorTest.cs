using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace SubsetTest
{
    public static class VectorTest
    {
        public static void vec_colex_next_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VEC_COLEX_NEXT_TEST tests VEC_COLEX_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;

            int[] a = new int[DIM_NUM];
            int base_ = 3;
            int i;
            bool more;

            Console.WriteLine("");
            Console.WriteLine("VEC_COLEX_NEXT_TEST");
            Console.WriteLine("  VEC_COLEX_NEXT generates all DIM_NUM-vectors");
            Console.WriteLine("  in colex order in a given base_ BASE.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + DIM_NUM + "");
            Console.WriteLine("  The base_ BASE =                 " + base_ + "");

            Console.WriteLine("");

            more = false;

            for (;;)
            {
                typeMethods.vec_colex_next(DIM_NUM, base_, ref a, ref more);

                if (!more)
                {
                    break;
                }

                string cout = "";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += a[i].ToString().PadLeft(4) + "  ";
                }

                Console.WriteLine(cout);
            }

        }

        public static void vec_colex_next2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VEC_COLEX_NEXT2_TEST tests VEC_COLEX_NEXT2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;

            int[] a = new int[DIM_NUM];
            int[] base_ = { 2, 1, 3 };
            int i;
            bool more;

            Console.WriteLine("");
            Console.WriteLine("VEC_COLEX_NEXT2_TEST");
            Console.WriteLine("  VEC_COLEX_NEXT2 generates all DIM_NUM-vectors");
            Console.WriteLine("  in colex order in a given base_ BASE.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + DIM_NUM + "");
            Console.WriteLine("");
            Console.WriteLine("  The base_ vector:");
            Console.WriteLine("");
            string cout = "";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += base_[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("");

            more = false;

            for (;;)
            {
                typeMethods.vec_colex_next2(DIM_NUM, base_, ref a, ref more);

                if (!more)
                {
                    break;
                }

                cout = "";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += a[i].ToString().PadLeft(4) + "  ";
                }

                Console.WriteLine(cout);
            }

        }

        public static void vec_colex_next3_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VEC_COLEX_NEXT3_TEST tests VEC_COLEX_NEXT3.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 August 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;

            int[] a = new int[DIM_NUM];
            int[] base_ = { 2, 1, 3 };
            int i;
            bool more;

            Console.WriteLine("");
            Console.WriteLine("VEC_COLEX_NEXT3_TEST");
            Console.WriteLine("  VEC_COLEX_NEXT3 generates all DIM_NUM-vectors");
            Console.WriteLine("  in colex order in a given base_ BASE.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + DIM_NUM + "");
            Console.WriteLine("");
            Console.WriteLine("  The base_ vector:");
            Console.WriteLine("");
            string cout = "";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += base_[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("");

            more = false;

            for (;;)
            {
                typeMethods.vec_colex_next3(DIM_NUM, base_, ref a, ref more);

                if (!more)
                {
                    break;
                }

                cout = "";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += a[i].ToString().PadLeft(4) + "  ";
                }

                Console.WriteLine(cout);
            }

        }

        public static void vec_gray_next_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VEC_GRAY_NEXT_TEST tests VEC_GRAY_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 4;

            int[] a = new int[N];
            int[] active = new int [N];
            int[] base_ = { 2, 2, 1, 4 };
            int change = 0;
            int[] dir = new int[N];
            bool done;
            int i;
            int prod;
            int rank;

            prod = 1;
            for (i = 0; i < N; i++)
            {
                prod = prod * base_[i];
            }

            Console.WriteLine("");
            Console.WriteLine("VEC_GRAY_NEXT_TEST");
            Console.WriteLine("  VEC_GRAY_NEXT generates product space elements.");
            Console.WriteLine("");
            Console.WriteLine("  The number of components is " + N + "");
            Console.WriteLine("  The number of elements is " + prod + "");
            Console.WriteLine("  Each component has its own number of degrees of");
            Console.WriteLine("  freedom.");
            Console.WriteLine("");
            string cout = "  Rank Change     ";
            for (i = 0; i < N; i++)
            {
                cout += base_[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("");

            rank = 0;
            done = true;

            for (;;)
            {
                rank = rank + 1;

                typeMethods.vec_gray_next(N, base_, a, ref done, active, dir, ref change);

                if (done)
                {
                    break;
                }

                cout = rank.ToString().PadLeft(4)
                       + change.ToString().PadLeft(4);
                for (i = 0; i < N; i++)
                {
                    cout += a[i].ToString().PadLeft(4) + "  ";
                }

                Console.WriteLine(cout);

            }

        }

        public static void vec_gray_rank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VEC_GRAY_RANK_TEST tests VEC_GRAY_RANK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 4;

            int[] a = new int[N];
            int[] base_ = { 2, 2, 1, 4 };
            int i;
            int prod;
            int rank;

            prod = 1;
            for (i = 0; i < N; i++)
            {
                prod = prod * base_[i];
            }

            Console.WriteLine("");
            Console.WriteLine("VEC_GRAY_RANK_TEST");
            Console.WriteLine("  VEC_GRAY_RANK ranks product space elements.");
            Console.WriteLine("");
            Console.WriteLine("  The number of components is " + N + "");
            Console.WriteLine("  The number of elements is " + prod + "");
            Console.WriteLine("  Each component has its own number of degrees of");
            Console.WriteLine("  freedom, which, for this example, are:");

            for (i = 0; i < N; i++)
            {
                a[i] = base_[i] / 2;
            }

            rank = Ranking.vec_gray_rank(N, base_, ref a);

            Console.WriteLine("");
            Console.WriteLine("  VEC_GRAY_RANK reports the element");
            Console.WriteLine("");
            string cout = "";
            for (i = 0; i < N; i++)
            {
                cout += a[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("");
            Console.WriteLine("  has rank " + rank + "");


            return;

        }

        public static void vec_gray_unrank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VEC_GRAY_UNRANK_TEST tests VEC_GRAY_UNRANK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 4;

            int[] a = new int[N];
            int[] base_ = { 2, 2, 1, 4 };
            int i;
            int prod;
            int rank;

            prod = 1;
            for (i = 0; i < N; i++)
            {
                prod = prod * base_[i];
            }

            Console.WriteLine("");
            Console.WriteLine("VEC_GRAY_UNRANK_TEST");
            Console.WriteLine("  VEC_GRAY_UNRANK unranks product space elements.");
            Console.WriteLine("");
            Console.WriteLine("  The number of components is " + N + "");
            Console.WriteLine("  The number of elements is " + prod + "");

            rank = 7;
            Ranking.vec_gray_unrank(N, base_, rank, ref a);

            Console.WriteLine("");
            Console.WriteLine("  VEC_GRAY_UNRANK reports the element of rank " + rank + "  is:");
            Console.WriteLine("");
            string cout = "";
            for (i = 0; i < N; i++)
            {
                cout += a[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);

        }

        public static void vec_lex_next_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VEC_LEX_NEXT_TEST tests VEC_LEX_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;

            int[] a = new int[DIM_NUM];
            int base_ = 3;
            int i;
            bool more;

            Console.WriteLine("");
            Console.WriteLine("VEC_LEX_NEXT_TEST");
            Console.WriteLine("  VEC_LEX_NEXT generates all DIM_NUM-vectors");
            Console.WriteLine("  in a given base.  Here we use base_ " + base_ + "");
            Console.WriteLine("");

            more = false;

            for (;;)
            {
                typeMethods.vec_lex_next(DIM_NUM, base_, ref a, ref more);

                if (!more)
                {
                    break;
                }

                string cout = "";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += a[i].ToString().PadLeft(4) + "  ";
                }

                Console.WriteLine(cout);

            }

        }

        public static void vec_random_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VEC_RANDOM_TEST tests VEC_RANDOM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 October 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            int[] a = new int[N];
            int base_;
            int i;
            int j;
            int seed;

            base_ = 3;
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("VEC_RANDOM_TEST");
            Console.WriteLine("  VEC_RANDOM generates a random N-vector");
            Console.WriteLine("  in a given base.");
            Console.WriteLine("  Here, we use base_ " + base_ + "");
            Console.WriteLine("");

            for (i = 1; i <= 5; i++)
            {
                typeMethods.vec_random(N, base_, ref seed, ref a);
                string cout = i.ToString().PadLeft(4) + "    ";
                for (j = 0; j < N; j++)
                {
                    cout += a[j].ToString().PadLeft(4) + "  ";
                }

                Console.WriteLine(cout);
            }

        }

        public static void vector_constrained_next_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_CONSTRAINED_NEXT_TEST tests VECTOR_CONSTRAINED_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 March 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            int constraint = 0;
            int i;
            int j;
            bool more;
            int[] x = new int[N];
            int[] x_max = { 4, 5, 3 };
            int[] x_min = { 2, 2, 1 };
            int x_prod;

            Console.WriteLine("");
            Console.WriteLine("VECTOR_CONSTRAINED_NEXT_TEST");
            Console.WriteLine("  VECTOR_CONSTRAINED_NEXT:");
            Console.WriteLine("  Consider vectors:");
            Console.WriteLine("    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),");
            Console.WriteLine("  Set");
            Console.WriteLine("    P = Product X_MAX(1:N)");
            Console.WriteLine("  Accept only vectors for which:");
            Console.WriteLine("    sum ( (X(1:N)-1) * P / X_MAX(1:N) ) <= P");

            more = false;

            Console.WriteLine("");
            Console.WriteLine("  X_MIN:");
            string cout = "";
            for (j = 0; j < N; j++)
            {
                cout += "  " + x_min[j].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
            Console.WriteLine("");
            Console.WriteLine("  X_MAX:");
            cout = "";
            for (j = 0; j < N; j++)
            {
                cout += x_max[j].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);

            i = 0;

            x_prod = 1;
            for (j = 0; j < N; j++)
            {
                x_prod = x_prod * x_max[j];
            }

            Console.WriteLine("");
            Console.WriteLine("  Maximum allowed CONSTRAINT = P = " + x_prod + "");
            Console.WriteLine("");

            for (;;)
            {
                typeMethods.vector_constrained_next(N, x_min, x_max, ref x, ref constraint, ref more);

                if (!more)
                {
                    break;
                }

                i = i + 1;
                cout = "  " + i.ToString().PadLeft(8);
                cout += "  " + constraint.ToString().PadLeft(12);
                for (j = 0; j < N; j++)
                {
                    cout += "  " + x[j].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);
            }

        }

        public static void vector_constrained_next2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_CONSTRAINED_NEXT2_TEST tests VECTOR_CONSTRAINED_NEXT2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 March 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N_MAX = 3;

            int constraint = 0;
            int i;
            int j;
            bool more;
            int n;
            int[] x = new int[N_MAX];
            int[] x_max = { 5, 6, 4 };
            int[] x_min = { 1, 1, 1 };
            int x_prod;
            string cout = "";
            Console.WriteLine("");
            Console.WriteLine("VECTOR_CONSTRAINED_NEXT2_TEST");
            Console.WriteLine("  VECTOR_CONSTRAINED_NEXT2:");
            Console.WriteLine("  Consider vectors:");
            Console.WriteLine("    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),");
            Console.WriteLine("  Set");
            Console.WriteLine("    P = Product X_MAX(1:N)");
            Console.WriteLine("  Accept only vectors for which:");
            Console.WriteLine("    sum ( X(1:N) * P / X_MAX(1:N) ) <= P");

            for (n = 2; n <= N_MAX; n++)
            {
                more = false;

                Console.WriteLine("");
                Console.WriteLine("  X_MIN:");
                cout = "";
                for (j = 0; j < n; j++)
                {
                    cout += "  " + x_min[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
                Console.WriteLine("");
                Console.WriteLine("  X_MAX:");
                cout = "";
                for (j = 0; j < n; j++)
                {
                    cout += "  " + x_max[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);

                i = 0;

                x_prod = 1;
                for (j = 0; j < n; j++)
                {
                    x_prod = x_prod * x_max[j];
                }

                Console.WriteLine("");
                Console.WriteLine("  Maximum allowed CONSTRAINT = P = " + x_prod + "");
                Console.WriteLine("");

                for (;;)
                {
                    typeMethods.vector_constrained_next2(n, x_min, x_max, ref x, ref constraint, ref more);

                    if (!more)
                    {
                        break;
                    }

                    i = i + 1;
                    cout = "  " + i.ToString().PadLeft(8);
                    cout += "  " + constraint.ToString().PadLeft(12);
                    for (j = 0; j < n; j++)
                    {
                        cout += "  " + x[j].ToString().PadLeft(8);
                    }

                    Console.WriteLine(cout);
                }
            }

        }

        public static void vector_constrained_next3_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_CONSTRAINED_NEXT3_TEST tests VECTOR_CONSTRAINED_NEXT3.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 April 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N_MAX = 3;

            double constraint = 0;
            int i;
            int j;
            bool more;
            int n;
            int[] x = new int[N_MAX];
            int[] x_max = { 5, 6, 4 };
            int[] x_min = { 1, 1, 1 };
            string cout;

            Console.WriteLine("");
            Console.WriteLine("VECTOR_CONSTRAINED_NEXT3_TEST");
            Console.WriteLine("  VECTOR_CONSTRAINED_NEXT3:");
            Console.WriteLine("  Consider vectors:");
            Console.WriteLine("    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),");
            Console.WriteLine("  Set");
            Console.WriteLine("    CONSTRAINT = sum ( X(1:N) / X_MAX(1:N) )");
            Console.WriteLine("  Accept only vectors for which:");
            Console.WriteLine("    CONSTRAINT <= 1");

            for (n = 2; n <= N_MAX; n++)
            {
                more = false;

                Console.WriteLine("");
                Console.WriteLine("  X_MIN:");
                cout = "";
                for (j = 0; j < n; j++)
                {
                    cout += "  " + x_min[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
                Console.WriteLine("");
                Console.WriteLine("  X_MAX:");
                cout = "";
                for (j = 0; j < n; j++)
                {
                    cout += "  " + x_max[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
                Console.WriteLine("");

                i = 0;

                for (;;)
                {
                    typeMethods.vector_constrained_next3(n, x_min, x_max, ref x, ref constraint, ref more);

                    if (!more)
                    {
                        break;
                    }

                    i = i + 1;
                    cout = "  " + i.ToString().PadLeft(8);
                    cout += "  " + constraint.ToString().PadLeft(12);
                    for (j = 0; j < n; j++)
                    {
                        cout += "  " + x[j].ToString().PadLeft(8);
                    }

                    Console.WriteLine(cout);
                }
            }

        }

        public static void vector_constrained_next4_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_CONSTRAINED_NEXT4_TEST tests VECTOR_CONSTRAINED_NEXT4.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 May 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N_MAX = 3;

            double[] alpha = { 4.0, 3.0, 5.0 };
            int i;
            int j;
            bool more;
            int n;
            double q = 20.0;
            double total;
            int[] x = new int[N_MAX];
            int[] x_max = { 2, 6, 4 };
            int[] x_min = { 1, 0, 1 };
            string cout;

            Console.WriteLine("");
            Console.WriteLine("VECTOR_CONSTRAINED_NEXT4_TEST");
            Console.WriteLine("  VECTOR_CONSTRAINED_NEXT4:");
            Console.WriteLine("  Consider vectors:");
            Console.WriteLine("    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),");
            Console.WriteLine("  Set");
            Console.WriteLine("    TOTAL = sum ( ALPHA(1:N) * X(1:N) )");
            Console.WriteLine("  Accept only vectors for which:");
            Console.WriteLine("    TOTAL <= Q");

            for (n = 2; n <= N_MAX; n++)
            {
                more = false;

                Console.WriteLine("");
                cout = "  ALPHA:";
                for (j = 0; j < n; j++)
                {
                    cout += "  " + alpha[j].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);
                cout = "  Q:    ";
                cout += "  " + q.ToString().PadLeft(8);
                Console.WriteLine(cout);
                cout = "  X_MIN:";
                for (j = 0; j < n; j++)
                {
                    cout += "  " + x_min[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
                cout = "  X_MAX:";
                for (j = 0; j < n; j++)
                {
                    cout += "  " + x_max[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
                Console.WriteLine("");

                i = 0;

                for (;;)
                {
                    typeMethods.vector_constrained_next4(n, alpha, x_min, x_max, ref x, q, ref more);

                    if (!more)
                    {
                        break;
                    }

                    total = 0.0;
                    for (j = 0; j < n; j++)
                    {
                        total = total + alpha[j] * (double)x[j];
                    }

                    i = i + 1;
                    cout = "  " + i.ToString().PadLeft(8);
                    cout += "  " + total.ToString().PadLeft(14);
                    for (j = 0; j < n; j++)
                    {
                        cout += "  " + x[j].ToString().PadLeft(8);
                    }

                    Console.WriteLine(cout);
                }
            }

        }

        public static void vector_constrained_next5_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_CONSTRAINED_NEXT5_TEST tests VECTOR_CONSTRAINED_NEXT5
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 May 2015subset.sh
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int base_;
            int i;
            int j;
            bool more;
            int n;
            int sum_max;
            int sum_min;
            int[] x = new int[3];

            Console.WriteLine("");
            Console.WriteLine("VECTOR_CONSTRAINED_NEXT5_TEST");
            Console.WriteLine("  VECTOR_CONSTRAINED_NEXT5:");
            Console.WriteLine("  Generate integer vectors X such that:");
            Console.WriteLine("    SUM_MIN <= sum ( X(1:N) ) <= SUM_MAX,");
            Console.WriteLine("  We require every X(I) to be at least 1.");

            n = 3;
            sum_min = 5;
            sum_max = 7;
            base_ = 0;
            more = false;

            Console.WriteLine("");
            Console.WriteLine("  N =       " + n + "");
            Console.WriteLine("  SUM_MIN = " + sum_min + "");
            Console.WriteLine("  SUM_MAX = " + sum_max + "");
            Console.WriteLine("");
            Console.WriteLine("         #        X(1)      X(2)      X(3)");
            Console.WriteLine("");

            i = 0;

            for (;;)
            {
                typeMethods.vector_constrained_next5(n, ref x, sum_min, sum_max, ref base_, ref more);

                if (!more)
                {
                    break;
                }

                i = i + 1;
                string cout = "  " + i.ToString().PadLeft(8);
                for (j = 0; j < n; j++)
                {
                    cout += "  " + x[j].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);
            }
        }

        public static void vector_constrained_next6_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_CONSTRAINED_NEXT6_TEST tests VECTOR_CONSTRAINED_NEXT6.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N_MAX = 3;

            double[] alpha = { 4.0, 3.0, 5.0 };
            int i;
            int j;
            bool more;
            int n;
            double q_max = 20.0;
            double q_min = 16.0;
            double total;
            int[] x = new int[N_MAX];
            int[] x_max = { 2, 6, 4 };
            int[] x_min = { 1, 0, 1 };
            string cout;

            Console.WriteLine("");
            Console.WriteLine("VECTOR_CONSTRAINED_NEXT6_TEST");
            Console.WriteLine("  VECTOR_CONSTRAINED_NEXT6:");
            Console.WriteLine("  Consider vectors:");
            Console.WriteLine("    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),");
            Console.WriteLine("  Set");
            Console.WriteLine("    TOTAL = sum ( ALPHA(1:N) * X(1:N) )");
            Console.WriteLine("  Accept only vectors for which:");
            Console.WriteLine("    Q_MIN <= TOTAL <= Q_MAX");

            for (n = 2; n <= N_MAX; n++)
            {
                more = false;

                Console.WriteLine("");
                cout = "  ALPHA:";
                for (j = 0; j < n; j++)
                {
                    cout += "  " + alpha[j].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);
                cout = "  Q_MIN:";
                cout += "  " + q_min.ToString().PadLeft(8);
                Console.WriteLine(cout);
                cout = "  Q_MAX:";
                cout += "  " + q_max.ToString().PadLeft(8);
                Console.WriteLine(cout);
                cout = "  X_MIN:";
                for (j = 0; j < n; j++)
                {
                    cout += "  " + x_min[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
                cout = "  X_MAX:";
                for (j = 0; j < n; j++)
                {
                    cout += "  " + x_max[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
                Console.WriteLine("");

                i = 0;

                for (;;)
                {
                    typeMethods.vector_constrained_next6(n, alpha, x_min, x_max, ref x, q_min,
                        q_max, ref more);

                    if (!more)
                    {
                        break;
                    }

                    total = 0.0;
                    for (j = 0; j < n; j++)
                    {
                        total = total + alpha[j] * (double)x[j];
                    }

                    i = i + 1;
                    cout = "  " + i.ToString().PadLeft(8);
                    cout += "  " + total.ToString().PadLeft(14);
                    for (j = 0; j < n; j++)
                    {
                        cout += "  " + x[j].ToString().PadLeft(8);
                    }

                    Console.WriteLine(cout);
                }
            }

        }

        public static void vector_constrained_next7_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_CONSTRAINED_NEXT7_TEST tests VECTOR_CONSTRAINED_NEXT7.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 July 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N_MAX = 3;

            double[] alpha = { 4.0, 3.0, 5.0 };
            int i;
            int j;
            bool more;
            int n;
            double q_max = 20.0;
            double q_min = 16.0;
            double total;
            int[] x = new int[N_MAX];
            int[] x_max = { 2, 6, 4 };
            string cout;

            Console.WriteLine("");
            Console.WriteLine("VECTOR_CONSTRAINED_NEXT7_TEST");
            Console.WriteLine("  VECTOR_CONSTRAINED_NEXT7:");
            Console.WriteLine("  Consider vectors:");
            Console.WriteLine("    0 <= X(1:N) <= X_MAX(1:N),");
            Console.WriteLine("  Set");
            Console.WriteLine("    TOTAL = sum ( ALPHA(1:N) * X(1:N) )");
            Console.WriteLine("  Accept only vectors for which:");
            Console.WriteLine("    Q_MIN <= TOTAL <= Q_MAX");

            for (n = 2; n <= N_MAX; n++)
            {
                more = false;

                Console.WriteLine("");
                cout = "  ALPHA:";
                for (j = 0; j < n; j++)
                {
                    cout += "  " + alpha[j].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);
                cout = "  Q_MIN:";
                cout += "  " + q_min.ToString().PadLeft(8);
                Console.WriteLine(cout);
                cout = "  Q_MAX:";
                cout += "  " + q_max.ToString().PadLeft(8);
                Console.WriteLine(cout);
                cout = "  X_MAX:";
                for (j = 0; j < n; j++)
                {
                    cout += "  " + x_max[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
                Console.WriteLine("");

                i = 0;

                for (;;)
                {
                    typeMethods.vector_constrained_next7(n, alpha, x_max, ref x, q_min,
                        q_max, ref more);

                    if (!more)
                    {
                        break;
                    }

                    total = 0.0;
                    for (j = 0; j < n; j++)
                    {
                        total = total + alpha[j] * (double)x[j];
                    }

                    i = i + 1;
                    cout = "  " + i.ToString().PadLeft(8);
                    cout += "  " + total.ToString().PadLeft(14);
                    for (j = 0; j < n; j++)
                    {
                        cout += "  " + x[j].ToString().PadLeft(8);
                    }

                    Console.WriteLine(cout);
                }
            }

        }

        public static void vector_next_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_NEXT_TEST tests VECTOR_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 July 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N_MAX = 3;

            int i;
            int j;
            bool more;
            int n;
            int[] x = new int[N_MAX];
            int[] x_max = { 2, 6, 4 };
            int[] x_min = { 1, 4, 3 };
            string cout;

            Console.WriteLine("");
            Console.WriteLine("VECTOR_NEXT_TEST");
            Console.WriteLine("  VECTOR_NEXT:");
            Console.WriteLine("  Generate all vectors X such that:");
            Console.WriteLine("    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),");

            for (n = 2; n <= N_MAX; n++)
            {
                more = false;

                Console.WriteLine("");
                cout = "    X_MIN:";
                for (j = 0; j < n; j++)
                {
                    cout += "  " + x_min[j].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);

                i = 0;

                for (;;)
                {
                    typeMethods.vector_next(n, x_min, x_max, ref x, ref more);

                    if (!more)
                    {
                        break;
                    }

                    i = i + 1;
                    cout = "  " + i.ToString().PadLeft(8);
                    for (j = 0; j < n; j++)
                    {
                        cout += "  " + x[j].ToString().PadLeft(8);
                    }

                    Console.WriteLine(cout);
                }

                cout = "    X_MAX:";
                for (j = 0; j < n; j++)
                {
                    cout += "  " + x_max[j].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);
            }

        }

    }
}