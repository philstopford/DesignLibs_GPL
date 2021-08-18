using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8col_variance(int m, int n, double[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8COL_VARIANCE returns the variances of an R8COL.
            //
            //  Discussion:
            //
            //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
            //    each of length M.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, double A[M*N], the array whose variances are desired.
            //
            //    Output, double R8COL_VARIANCE[N], the variances of the rows.
            //
        {
            double[] variance = new double[n];

            for (int j = 0; j < n; j++)
            {
                double mean = 0.0;
                for (int i = 0; i < m; i++)
                {
                    mean = mean + a[i + j * m];
                }

                mean = mean / (double) (m);

                variance[j] = 0.0;
                for (int i = 0; i < m; i++)
                {
                    variance[j] = variance[j] + Math.Pow(a[i + j * m] - mean, 2);
                }

                if (1 < m)
                {
                    variance[j] = variance[j] / (double) (m - 1);
                }
                else
                {
                    variance[j] = 0.0;
                }
            }

            return variance;
        }

        public static double[] r8col_mean(int m, int n, double[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8COL_MEAN returns the column means of an R8COL.
            //
            //  Discussion:
            //
            //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
            //    each of length M.
            //
            //  Example:
            //
            //    A =
            //      1  2  3
            //      2  6  7
            //
            //    R8COL_MEAN =
            //      1.5  4.0  5.0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns.
            //
            //    Input, double A[M*N], the array to be examined.
            //
            //    Output, double R8COL_MEAN[N], the means, or averages, of the columns.
            //
        {
            double[] mean = new double[n];

            for (int j = 0; j < n; j++)
            {
                mean[j] = 0.0;
                for (int i = 0; i < m; i++)
                {
                    mean[j] = mean[j] + a[i + j * m];
                }

                mean[j] = mean[j] / (double) (m);
            }

            return mean;
        }

        public static void r8col_permute ( int m, int n, int[] p, int base_, double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_PERMUTE permutes an R8COL in place.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    This routine permutes an array of real "objects", but the same
        //    logic can be used to permute an array of objects of any arithmetic
        //    type, or an array of objects of any complexity.  The only temporary
        //    storage required is enough to store a single object.  The number
        //    of data movements made is N + the number of cycles of order 2 or more,
        //    which is never more than N + N/2.
        //
        //  Example:
        //
        //    Input:
        //
        //      M = 2
        //      N = 5
        //      P = (   2,    4,    5,    1,    3 )
        //      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
        //          (11.0, 22.0, 33.0, 44.0, 55.0 )
        //      BASE = 1
        //
        //    Output:
        //
        //      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
        //             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the length of objects.
        //
        //    Input, int N, the number of objects.
        //
        //    Input, int P[N], the permutation.  P(I) = J means
        //    that the I-th element of the output array should be the J-th
        //    element of the input array.
        //
        //    Input, int BASE, is 0 for a 0-based permutation and 1 for a
        //    1-based permutation.
        //
        //    Input/output, double A[M*N], the array to be permuted.
        //
        {
            double[] a_temp;
            int i;
            int iget;
            int iput;
            int istart;
            int j;

            if (!perm_check(n, p, base_))
            {
                Console.WriteLine("");
                Console.WriteLine("R8COL_PERMUTE - Fatal error!");
                Console.WriteLine("  PERM_CHECK rejects this permutation.");
                return;
            }

            //
            //  In order for the sign negation trick to work, we need to assume that the
            //  entries of P are strictly positive.  Presumably, the lowest number is BASE.
            //  So temporarily add 1-BASE to each entry to force positivity.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = p[i] + 1 - base_;
            }

            a_temp = new double[m];
            //
            //  Search for the next element of the permutation that has not been used.
            //
            for (istart = 1; istart <= n; istart++)
            {
                if (p[istart - 1] < 0)
                {
                    continue;
                }
                else if (p[istart - 1] == istart)
                {
                    p[istart - 1] = -p[istart - 1];
                    continue;
                }
                else
                {
                    for (i = 0; i < m; i++)
                    {
                        a_temp[i] = a[i + (istart - 1) * m];
                    }

                    iget = istart;
                    //
                    //  Copy the new value into the vacated entry.
                    //
                    for (;;)
                    {
                        iput = iget;
                        iget = p[iget - 1];

                        p[iput - 1] = -p[iput - 1];

                        if (iget < 1 || n < iget)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("R8COL_PERMUTE - Fatal error!");
                            Console.WriteLine("  Entry IPUT = " + iput + " of the permutation has");
                            Console.WriteLine("  an illegal value IGET = " + iget + ".");
                            return;
                        }

                        if (iget == istart)
                        {
                            for (i = 0; i < m; i++)
                            {
                                a[i + (iput - 1) * m] = a_temp[i];
                            }

                            break;
                        }

                        for (i = 0; i < m; i++)
                        {
                            a[i + (iput - 1) * m] = a[i + (iget - 1) * m];
                        }
                    }
                }
            }

            //
            //  Restore the signs of the entries.
            //
            for (j = 0; j < n; j++)
            {
                p[j] = -p[j];
            }

            //
            //  Restore the base of the entries.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = p[i] - 1 + base_;
            }
        }

        public static void r8col_permute(int m, int n, int[] p, ref double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8COL_PERMUTE permutes an R8COL in place.
            //
            //  Discussion:
            //
            //    An R8COL is an M by N array of double precision values, regarded
            //    as an array of N columns of length M.
            //
            //    This routine permutes an array of real "objects", but the same
            //    logic can be used to permute an array of objects of any arithmetic
            //    type, or an array of objects of any complexity.  The only temporary
            //    storage required is enough to store a single object.  The number
            //    of data movements made is N + the number of cycles of order 2 or more,
            //    which is never more than N + N/2.
            //
            //  Example:
            //
            //    Input:
            //
            //      M = 2
            //      N = 5
            //      P = (   2,    4,    5,    1,    3 )
            //      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
            //          (11.0, 22.0, 33.0, 44.0, 55.0 )
            //
            //    Output:
            //
            //      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
            //             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the length of objects.
            //
            //    Input, int N, the number of objects.
            //
            //    Input/output, double A[M*N], the array to be permuted.
            //
            //    Input, int P[N], the permutation.  P(I) = J means
            //    that the I-th element of the output array should be the J-th
            //    element of the input array.  P must be a legal permutation
            //    of the integers from 1 to N, otherwise the algorithm will
            //    fail catastrophically.
            //
        {
            double[] a_temp;
            int i;
            int iget;
            int iput;
            int istart;
            int j;

            a_temp = new double[m];
            //
            //  Need to increment the entries by 1 in order to use the sign trick.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = p[i] + 1;
            }

            //
            //  Search for the next element of the permutation that has not been used.
            //
            for (istart = 1; istart <= n; istart++)
            {
                if (p[istart - 1] < 0)
                {
                    continue;
                }
                else if (p[istart - 1] == istart)
                {
                    p[istart - 1] = -p[istart - 1];
                    continue;
                }
                else
                {
                    for (i = 0; i < m; i++)
                    {
                        a_temp[i] = a[i + (istart - 1) * m];
                    }

                    iget = istart;
                    //
                    //  Copy the new value into the vacated entry.
                    //
                    for (;;)
                    {
                        iput = iget;
                        iget = p[iget - 1];

                        p[iput - 1] = -p[iput - 1];

                        if (iget < 1 || n < iget)
                        {
                            Console.WriteLine();
                            Console.WriteLine("R8COL_PERMUTE - Fatal error!");
                            Console.WriteLine("  Entry IPUT = " + iput + " of the permutation has");
                            Console.WriteLine("  an illegal value IGET = " + iget + "");
                            return;
                        }

                        if (iget == istart)
                        {
                            for (i = 0; i < m; i++)
                            {
                                a[i + (iput - 1) * m] = a_temp[i];
                            }

                            break;
                        }

                        for (i = 0; i < m; i++)
                        {
                            a[i + (iput - 1) * m] = a[i + (iget - 1) * m];
                        }
                    }
                }
            }

            //
            //  Restore the signs of the entries.
            //
            for (j = 0; j < n; j++)
            {
                p[j] = -p[j];
            }

            //
            //  Need to unincrement the entries by 1.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = p[i] - 1;
            }
        }

        public static void r8col_separation(int m, int n, double[] a, ref double d_min, ref double d_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8COL_SEPARATION returns the "separation" of an R8COL.
            //
            //  Discussion:
            //
            //    D_MIN is the minimum distance between two columns,
            //    D_MAX is the maximum distance between two columns.
            //
            //    The distances are measured using the Loo norm.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns 
            //    in the array.  If N < 2, it does not make sense to call this routine.
            //
            //    Input, double A[M*N], the array whose variances are desired.
            //
            //    Output, double &D_MIN, &D_MAX, the minimum and maximum distances.
            //
        {
            double d;
            int i;
            int j1;
            int j2;

            d_min = r8_huge();
            d_max = 0.0;

            for (j1 = 0; j1 < n; j1++)
            {
                for (j2 = j1 + 1; j2 < n; j2++)
                {
                    d = 0.0;
                    for (i = 0; i < m; i++)
                    {
                        d = Math.Max(d, Math.Abs(a[i + j1 * m] - a[i + j2 * m]));
                    }

                    d_min = Math.Min(d_min, d);
                    d_max = Math.Max(d_max, d);
                }
            }
        }
    }
}