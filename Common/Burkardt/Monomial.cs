using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt
{
    public static class Monomial
    {
        public static void mono_print ( int m, int[] f, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONO_PRINT prints a monomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int F[M], the exponents.
        //
        //    Input, string TITLE, a title.
        //
        {
            int i;

            string cout = title;
            cout += "x^(";
            for ( i = 0; i < m; i++ )
            {
                cout += f[i];
                if ( i < m - 1 )
                {
                    cout += ",";
                }
                else
                {
                    cout += ").";
                    Console.WriteLine(cout);
                }
            }
        }
        
        public static double[] monomial_value(int m, int n, int[] e, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_VALUE evaluates a monomial.
            //
            //  Discussion:
            //
            //    This routine evaluates a monomial of the form
            //
            //      product ( 1 <= i <= m ) x(i)^e(i)
            //
            //    where the exponents are nonnegative integers.  Note that
            //    if the combination 0^0 is encountered, it should be treated
            //    as 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of points at which the
            //    monomial is to be evaluated.
            //
            //    Input, int E[M], the exponents.
            //
            //    Input, double X[M*N], the point coordinates.
            //
            //    Output, double MONOMIAL_VALUE[N], the value of the monomial.
            //
        {
            int i;
            int j;
            double[] v;

            v = new double[n];

            for (j = 0; j < n; j++)
            {
                v[j] = 1.0;
            }

            for (i = 0; i < m; i++)
            {
                if (0 != e[i])
                {
                    for (j = 0; j < n; j++)
                    {
                        v[j] = v[j] * Math.Pow(x[i + j * m], e[i]);
                    }
                }
            }

            return v;
        }

        public static int mono_between_enum(int d, int n1, int n2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONO_BETWEEN_ENUM enumerates monomials in D dimensions of degrees in a range.
            //
            //  Discussion:
            //
            //    For D = 3, we have the following table:
            //
            //     N2 0  1  2  3  4  5  6   7   8
            //   N1 +----------------------------
            //    0 | 1  4 10 20 35 56 84 120 165
            //    1 | 0  3  9 19 34 55 83 119 164
            //    2 | 0  0  6 16 31 52 80 116 161
            //    3 | 0  0  0 10 25 46 74 110 155
            //    4 | 0  0  0  0 15 36 64 100 145
            //    5 | 0  0  0  0  0 21 49  85 130
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 November 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int D, the spatial dimension.
            //
            //    Input, int N1, N2, the minimum and maximum degrees.
            //    0 <= N1 <= N2.
            //
            //    Output, int MONO_BETWEEN_ENUM, the number of monomials 
            //    in D variables, of total degree between N1 and N2 inclusive.
            //
        {
            int n0;
            int n1_copy;
            int value;

            n1_copy = Math.Max(n1, 0);

            if (n2 < n1_copy)
            {
                value = 0;
                return value;
            }

            if (n1_copy == 0)
            {
                value = typeMethods.i4_choose(n2 + d, n2);
            }
            else if (n1_copy == n2)
            {
                value = typeMethods.i4_choose(n2 + d - 1, n2);
            }
            else
            {
                n0 = n1_copy - 1;
                value = typeMethods.i4_choose(n2 + d, n2) - typeMethods.i4_choose(n0 + d, n0);
            }

            return value;
        }

        public static void mono_between_next_grlex(int d, int n1, int n2, ref int[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONO_BETWEEN_NEXT_GRLEX: grlex next monomial, degree between N1 and N2.
            //
            //  Discussion:
            //
            //    We consider all monomials in a D dimensional space, with total
            //    degree N between N1 and N2, inclusive.
            //
            //    For example:
            //
            //    D = 3
            //    N1 = 2
            //    N2 = 3
            //
            //    #  X(1)  X(2)  X(3)  Degree
            //      +------------------------
            //    1 |  0     0     2        2
            //    2 |  0     1     1        2
            //    3 |  0     2     0        2
            //    4 |  1     0     1        2
            //    5 |  1     1     0        2
            //    6 |  2     0     0        2
            //      |
            //    7 |  0     0     3        3
            //    8 |  0     1     2        3
            //    9 |  0     2     1        3
            //   10 |  0     3     0        3
            //   11 |  1     0     2        3
            //   12 |  1     1     1        3
            //   13 |  1     2     0        3
            //   14 |  2     0     1        3
            //   15 |  2     1     0        3
            //   16 |  3     0     0        3
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 December 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int D, the spatial dimension.
            //
            //    Input, int N1, N2, the minimum and maximum degrees.
            //    0 <= N1 <= N2.
            //
            //    Input/output, int X[D], the current monomial.
            //    To start the sequence, set X = [ 0, 0, ..., 0, N1 ].
            //    The last value in the sequence is X = [ N2, 0, ..., 0, 0 ].
            //
        {
            if (n1 < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("MONO_BETWEEN_NEXT_GRLEX - Fatal error!");
                Console.WriteLine("  N1 < 0.");
                return;
            }

            if (n2 < n1)
            {
                Console.WriteLine("");
                Console.WriteLine("MONO_BETWEEN_NEXT_GRLEX - Fatal error!");
                Console.WriteLine("  N2 < N1.");
                return;
            }

            if (typeMethods.i4vec_sum(d, x) < n1)
            {
                Console.WriteLine("");
                Console.WriteLine("MONO_BETWEEN_NEXT_GRLEX - Fatal error!");
                Console.WriteLine("  Input X sums to less than N1.");
                return;
            }

            if (n2 < typeMethods.i4vec_sum(d, x))
            {
                Console.WriteLine("");
                Console.WriteLine("MONO_BETWEEN_NEXT_GRLEX - Fatal error!");
                Console.WriteLine("  Input X sums to more than N2.");
                return;
            }

            if (n2 == 0)
            {
                return;
            }

            if (x[0] == n2)
            {
                x[0] = 0;
                x[d - 1] = n1;
            }
            else
            {
                mono_next_grlex(d, ref x);
            }
        }

        public static void mono_next_grlex(int d, ref int[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONO_NEXT_GRLEX returns the next monomial in grlex order.
            //
            //  Discussion:
            //
            //    Example:
            //
            //    D = 3
            //
            //    #  X(1)  X(2)  X(3)  Degree
            //      +------------------------
            //    1 |  0     0     0        0
            //      |
            //    2 |  0     0     1        1
            //    3 |  0     1     0        1
            //    4 |  1     0     0        1
            //      |
            //    5 |  0     0     2        2
            //    6 |  0     1     1        2
            //    7 |  0     2     0        2
            //    8 |  1     0     1        2
            //    9 |  1     1     0        2
            //   10 |  2     0     0        2
            //      |
            //   11 |  0     0     3        3
            //   12 |  0     1     2        3
            //   13 |  0     2     1        3
            //   14 |  0     3     0        3
            //   15 |  1     0     2        3
            //   16 |  1     1     1        3
            //   17 |  1     2     0        3
            //   18 |  2     0     1        3
            //   19 |  2     1     0        3
            //   20 |  3     0     0        3
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 December 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int D, the spatial dimension.
            //
            //    Input/output, int X[D], the current monomial.
            //    The first element is X = [ 0, 0, ..., 0, 0 ].
            //
        {
            int i;
            int im1 = 0;
            int j;
            int t = 0;
            //
            //  Ensure that 1 <= D.
            //
            if (d < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("MONO_NEXT_GRLEX - Fatal error!");
                Console.WriteLine("  D < 1");
                return;
            }

            //
            //  Ensure that 0 <= X(I).
            //
            for (i = 0; i < d; i++)
            {
                if (x[i] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("MONO_NEXT_GRLEX - Fatal error!");
                    Console.WriteLine("  X[I] < 0");
                    return;
                }
            }

            //
            //  Find I, the index of the rightmost nonzero entry of X.
            //
            i = 0;
            for (j = d; 1 <= j; j--)
            {
                if (0 < x[j - 1])
                {
                    i = j;
                    break;
                }
            }

            //
            //  set T = X(I)
            //  set X(I) to zero,
            //  increase X(I-1) by 1,
            //  increment X(D) by T-1.
            //
            if (i == 0)
            {
                x[d - 1] = 1;
                return;
            }
            else if (i == 1)
            {
                t = x[0] + 1;
                im1 = d;
            }
            else if (1 < i)
            {
                t = x[i - 1];
                im1 = i - 1;
            }

            x[i - 1] = 0;
            x[im1 - 1] = x[im1 - 1] + 1;
            x[d - 1] = x[d - 1] + t - 1;

            return;
        }

        public static int mono_rank_grlex(int m, int[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONO_RANK_GRLEX computes the graded lexicographic rank of a monomial.
            //
            //  Discussion:
            //
            //    The graded lexicographic ordering is used, over all monomials of
            //    dimension M, for degree NM = 0, 1, 2, ...
            //
            //    For example, if M = 3, the ranking begins:
            //
            //    Rank  Sum    1  2  3
            //    ----  ---   -- -- --
            //       1    0    0  0  0
            //
            //       2    1    0  0  1
            //       3    1    0  1  0
            //       4    1    1  0  1
            //
            //       5    2    0  0  2
            //       6    2    0  1  1
            //       7    2    0  2  0
            //       8    2    1  0  1
            //       9    2    1  1  0
            //      10    2    2  0  0
            //
            //      11    3    0  0  3
            //      12    3    0  1  2
            //      13    3    0  2  1
            //      14    3    0  3  0
            //      15    3    1  0  2
            //      16    3    1  1  1
            //      17    3    1  2  0
            //      18    3    2  0  1
            //      19    3    2  1  0
            //      20    3    3  0  0
            //
            //      21    4    0  0  4
            //      ..   ..   .. .. ..
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 December 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //    1 <= D.
            //
            //    Input, int XC[M], the monomial.
            //    For each 1 <= I <= M, we have 0 <= XC(I).
            //
            //    Output, int MONO_RANK_GRLEX, the rank.
            //
        {
            int i;
            int j;
            int ks;
            int n;
            int nm;
            int ns;
            int rank;
            int tim1;
            int[] xs;
            //
            //  Ensure that 1 <= M.
            //
            if (m < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("MONO_RANK_GRLEX - Fatal error!");
                Console.WriteLine("  M < 1");
                return (1);
            }

            //
            //  Ensure that 0 <= X(I).
            //
            for (i = 0; i < m; i++)
            {
                if (x[i] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("MONO_RANK_GRLEX - Fatal error!");
                    Console.WriteLine("  X[I] < 0");
                    return (1);
                }
            }

            //
            //  NM = sum ( X )
            //
            nm = typeMethods.i4vec_sum(m, x);
            //
            //  Convert to KSUBSET format.
            //
            ns = nm + m - 1;
            ks = m - 1;
            xs = new int[ks];
            xs[0] = x[0] + 1;
            for (i = 2; i < m; i++)
            {
                xs[i - 1] = xs[i - 2] + x[i - 1] + 1;
            }

            //
            //  Compute the rank.
            //
            rank = 1;

            for (i = 1; i <= ks; i++)
            {
                if (i == 1)
                {
                    tim1 = 0;
                }
                else
                {
                    tim1 = xs[i - 2];
                }

                if (tim1 + 1 <= xs[i - 1] - 1)
                {
                    for (j = tim1 + 1; j <= xs[i - 1] - 1; j++)
                    {
                        rank = rank + typeMethods.i4_choose(ns - j, ks - i);
                    }
                }
            }

            for (n = 0; n < nm; n++)
            {
                rank = rank + typeMethods.i4_choose(n + m - 1, n);
            }
            
            return rank;
        }

        public static int mono_total_enum(int d, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONO_TOTAL_ENUM enumerates monomials in D dimensions of degree equal to N.
            //
            //  Discussion:
            //
            //    For D = 3, we have the following values:
            //
            //    N  VALUE
            //
            //    0    1
            //    1    3
            //    2    6
            //    3   10
            //    4   15
            //    5   21
            //
            //    In particular, VALUE(3,3) = 10 because we have the 10 monomials:
            //
            //      x^3, x^2y, x^2z, xy^2, xyz, xz^3, y^3, y^2z, yz^2, z^3.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 November 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int D, the spatial dimension.
            //
            //    Input, int N, the maximum degree.
            //
            //    Output, int MONO_TOTAL_ENUM, the number of monomials in D variables,
            //    of total degree N.
            //
        {
            int value = typeMethods.i4_choose(n + d - 1, n);

            return value;
        }

        public static void mono_total_next_grlex(int d, int n, ref int[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONO_TOTAL_NEXT_GRLEX: grlex next monomial with total degree equal to N.
            //
            //  Discussion:
            //
            //    We consider all monomials in a D dimensional space, with total degree N.
            //
            //    For example:
            //
            //    D = 3
            //    N = 3
            //
            //    #  X(1)  X(2)  X(3)  Degree
            //      +------------------------
            //    1 |  0     0     3        3
            //    2 |  0     1     2        3
            //    3 |  0     2     1        3
            //    4 |  0     3     0        3
            //    5 |  1     0     2        3
            //    6 |  1     1     1        3
            //    7 |  1     2     0        3
            //    8 |  2     0     1        3
            //    9 |  2     1     0        3
            //   10 |  3     0     0        3
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 December 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int D, the spatial dimension.
            //
            //    Input, int N, the degree.
            //    0 <= N.
            //
            //    Input/output, int X[D], the current monomial.
            //    To start the sequence, set X = [ 0, 0, ..., 0, N ].
            //    The last value in the sequence is X = [ N, 0, ..., 0, 0 ].
            //
        {
            if (n < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("MONO_TOTAL_NEXT_GRLEX - Fatal error!");
                Console.WriteLine("  N < 0.");
                return;
            }

            if (typeMethods.i4vec_sum(d, x) != n)
            {
                Console.WriteLine("");
                Console.WriteLine("MONO_TOTAL_NEXT_GRLEX - Fatal error!");
                Console.WriteLine("  Input X does not sum to N.");
                return;
            }

            if (n == 0)
            {
                return;
            }

            if (x[0] == n)
            {
                x[0] = 0;
                x[d - 1] = n;
            }
            else
            {
                mono_next_grlex(d, ref x);
            }
        }

        public static int[] mono_unrank_grlex(int d, int rank)

            /******************************************************************************/
            //
            //  Purpose:
            //
            //    MONO_UNRANK_GRLEX computes the composition of given grlex rank.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int D, the spatial dimension.
            //    1 <= D.
            //
            //    Input, int RANK, the rank.
            //    1 <= RANK.
            //
            //    Output, int MONO_UNRANK_GRLEX[D], the monomial X of the given rank.
            //    For each I, 0 <= XC[I] <= NM, and 
            //    sum ( 1 <= I <= D ) XC[I] = NM.
            //
        {
            int i;
            int j;
            int ks;
            int nm;
            int ns;
            int r;
            int rank1;
            int rank2;
            int[] x;
            int[] xs;
            //
            //  Ensure that 1 <= D.
            //
            if (d < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("MONO_UNRANK_GRLEX - Fatal error!");
                Console.WriteLine("  D < 1");
                Console.WriteLine("  D = " + d + "");
                return null;
            }

            //
            //  Ensure that 1 <= RANK.
            //
            if (rank < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("MONO_UNRANK_GRLEX - Fatal error!");
                Console.WriteLine("  RANK < 1");
                Console.WriteLine("  RANK = " + rank + "");
                return null;
            }

            //
            //  Special case D == 1.
            //
            if (d == 1)
            {
                x = new int[d];
                x[0] = rank - 1;
                return x;
            }

            //
            //  Determine the appropriate value of NM.
            //  Do this by adding up the number of compositions of sum 0, 1, 2, 
            //  ..., without exceeding RANK.  Moreover, RANK - this sum essentially
            //  gives you the rank of the composition within the set of compositions
            //  of sum NM.  And that's the number you need in order to do the
            //  unranking.
            //
            rank1 = 1;
            nm = -1;
            for (;;)
            {
                nm = nm + 1;
                r = typeMethods.i4_choose(nm + d - 1, nm);
                if (rank < rank1 + r)
                {
                    break;
                }

                rank1 = rank1 + r;
            }

            rank2 = rank - rank1;
            //
            //  Convert to KSUBSET format.
            //  Apology: an unranking algorithm was available for KSUBSETS,
            //  but not immediately for compositions.  One day we will come back
            //  and simplify all this.
            //
            ks = d - 1;
            ns = nm + d - 1;
            xs = new int[ks];

            j = 1;

            for (i = 1; i <= ks; i++)
            {
                r = typeMethods.i4_choose(ns - j, ks - i);

                while (r <= rank2 && 0 < r)
                {
                    rank2 = rank2 - r;
                    j = j + 1;
                    r = typeMethods.i4_choose(ns - j, ks - i);
                }

                xs[i - 1] = j;
                j = j + 1;
            }

            //
            //  Convert from KSUBSET format to COMP format.
            //
            x = new int[d];
            x[0] = xs[0] - 1;
            for (i = 2; i < d; i++)
            {
                x[i - 1] = xs[i - 1] - xs[i - 2] - 1;
            }

            x[d - 1] = ns - xs[ks - 1];
            
            return x;
        }

        public static int mono_upto_enum(int d, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONO_UPTO_ENUM enumerates monomials in D dimensions of degree up to N.
            //
            //  Discussion:
            //
            //    For D = 2, we have the following values:
            //
            //    N  VALUE
            //
            //    0    1
            //    1    3
            //    2    6
            //    3   10
            //    4   15
            //    5   21
            //
            //    In particular, VALUE(2,3) = 10 because we have the 10 monomials:
            //
            //      1, x, y, x^2, xy, y^2, x^3, x^2y, xy^2, y^3.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 November 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int D, the spatial dimension.
            //
            //    Input, int N, the maximum degree.
            //
            //    Output, int MONO_UPTO_ENUM, the number of monomials in
            //    D variables, of total degree N or less.
            //
        {
            int value = typeMethods.i4_choose(n + d, n);

            return value;
        }

        public static void mono_upto_next_grlex(int m, int n, ref int[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONO_UPTO_NEXT_GRLEX: grlex next monomial with total degree up to N.
            //
            //  Discussion:
            //
            //    We consider all monomials in a M dimensional space, with total
            //    degree up to N.
            //
            //    For example:
            //
            //    M = 3
            //    N = 3
            //
            //    #  X(1)  X(2)  X(3)  Degree
            //      +------------------------
            //    1 |  0     0     0        0
            //      |
            //    2 |  0     0     1        1
            //    3 |  0     1     0        1
            //    4 |  1     0     0        1
            //      |
            //    5 |  0     0     2        2
            //    6 |  0     1     1        2
            //    7 |  0     2     0        2
            //    8 |  1     0     1        2
            //    9 |  1     1     0        2
            //   10 |  2     0     0        2
            //      |
            //   11 |  0     0     3        3
            //   12 |  0     1     2        3
            //   13 |  0     2     1        3
            //   14 |  0     3     0        3
            //   15 |  1     0     2        3
            //   16 |  1     1     1        3
            //   17 |  1     2     0        3
            //   18 |  2     0     1        3
            //   19 |  2     1     0        3
            //   20 |  3     0     0        3
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 December 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the maximum degree.
            //    0 <= N.
            //
            //    Input/output, int X[M], the current monomial.
            //    To start the sequence, set X = [ 0, 0, ..., 0, 0 ].
            //    The last value in the sequence is X = [ N, 0, ..., 0, 0 ].
            //
        {
            if (n < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("MONO_UPTO_NEXT_GRLEX - Fatal error!");
                Console.WriteLine("  N < 0.");
                return;
            }

            if (typeMethods.i4vec_sum(m, x) < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("MONO_UPTO_NEXT_GRLEX - Fatal error!");
                Console.WriteLine("  Input X sums to less than 0.");
                return;
            }

            if (n < typeMethods.i4vec_sum(m, x))
            {
                Console.WriteLine("");
                Console.WriteLine("MONO_UPTO_NEXT_GRLEX - Fatal error!");
                Console.WriteLine("  Input X sums to more than N.");
                return;
            }

            if (n == 0)
            {
                return;
            }

            if (x[0] == n)
            {
                x[0] = 0;
                x[m - 1] = n;
            }
            else
            {
                mono_next_grlex(m, ref x);
            }
        }

        public static int[] mono_upto_random(int m, int n, ref int seed, ref int rank)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONO_UPTO_RANDOM: random monomial with total degree less than or equal to N.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 November 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the degree.
            //    0 <= N.
            //
            //    Input/output, int &SEED, the random number seed.
            //
            //    Output, int &RANK, the rank of the monomial.
            //
            //    Output, int MONO_UPTO_RANDOM[M], the random monomial.
            //
        {
            int rank_max;
            int rank_min;
            int[] x;

            rank_min = 1;
            rank_max = mono_upto_enum(m, n);
            rank = UniformRNG.i4_uniform_ab(rank_min, rank_max, ref seed);
            x = mono_unrank_grlex(m, rank);

            return x;
        }

        public static double[] mono_value(int d, int nx, int[] f, double[] x, int xIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONO_VALUE evaluates a monomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int D, the spatial dimension.
        //
        //    Input, int NX, the number of evaluation points.
        //
        //    Input, int F[D], the exponents of the monomial.
        //
        //    Input, double X[D*NX], the coordinates of the evaluation points.
        //
        //    Output, double MONO_VALUE[NX], the value of the monomial at X.
        //
        {
            int i;
            int j;
            double[] v;

            v = new double[nx];

            for (j = 0; j < nx; j++)
            {
                v[j] = 1.0;
                for (i = 0; i < d; i++)
                {
                    v[j] = v[j] * Math.Pow(x[xIndex + (i + j * d)], f[i]);
                }
            }

            return v;
        }

        public static double[] monomial_value_1d ( int n, int e, double[] x )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_VALUE_1D evaluates a monomial in 1D.
            //
            //  Discussion:
            //
            //    This routine evaluates a monomial of the form
            //
            //      product ( 1 <= i <= m ) x(i)^e(i)
            //
            //    where the exponents are nonnegative integers.  Note that
            //    if the combination 0^0 is encountered, it should be treated
            //    as 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points at which the
            //    monomial is to be evaluated.
            //
            //    Input, int E, the exponents.
            //
            //    Input, double X[N], the point coordinates.
            //
            //    Output, double MONOMIAL_VALUE_1D[N], the value of the monomial.
            //
        {
            int j;
            double[] v;

            v = new double[n];

            for ( j = 0; j < n; j++ )
            {
                v[j] = Math.Pow ( x[j], e );
            }

            return v;
        }

    }
}