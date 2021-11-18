using System;
using System.Globalization;
using Burkardt.SortNS;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static int r8row_compare(int m, int n, double[] a, int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_COMPARE compares two rows in an R8ROW.
        //
        //  Discussion:
        //
        //    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
        //    each of length N.
        //
        //  Example:
        //
        //    Input:
        //
        //      M = 4, N = 3, I = 2, J = 4
        //
        //      A = (
        //        1. 5. 9.
        //        2. 6. 10.
        //        3. 7. 11.
        //        4. 8. 12. )
        //
        //    Output:
        //
        //      R8ROW_COMPARE = -1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the M by N array.
        //
        //    Input, int I, J, the rows to be compared.
        //    I and J must be between 0 and M - 1.
        //
        //    Output, int R8ROW_COMPARE, the results of the comparison:
        //    -1, row I < row J,
        //     0, row I = row J,
        //    +1, row J < row I.
        //
    {
        //
        //  Check.
        //
        if (i < 0 || m <= i)
        {
            Console.WriteLine("");
            Console.WriteLine("R8ROW_COMPARE - Fatal error!");
            Console.WriteLine("  Row index I is out of bounds.");
            Console.WriteLine("  I = " + i + "");
            return 1;
        }

        if (j < 0 || m <= j)
        {
            Console.WriteLine("");
            Console.WriteLine("R8ROW_COMPARE - Fatal error!");
            Console.WriteLine("  Row index J is out of bounds.");
            Console.WriteLine("  J = " + j + "");
            return 1;
        }

        int value = 0;

        if (i == j)
        {
            return value;
        }

        int k = 0;

        while (k < n)
        {
            if (a[i + k * m] < a[j + k * m])
            {
                value = -1;
                return value;
            }

            if (a[j + k * m] < a[i + k * m])
            {
                value = +1;
                return value;
            }

            k += 1;
        }

        return value;
    }

    public static double[] r8row_indicator_new(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_INDICATOR_NEW sets up an "indicator" R8ROW.
        //
        //  Discussion:
        //
        //    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
        //    each of length N.
        //
        //    The value of each entry suggests its location, as in:
        //
        //      11  12  13  14
        //      21  22  23  24
        //      31  32  33  34
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Output, double R8ROW_INDICATOR_NEW[M*N], the table.
        //
    {
        int i;

        double[] a = new double[m * n];

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (i = 1; i <= m; i++)
        {
            int j;
            for (j = 1; j <= n; j++)
            {
                a[i - 1 + (j - 1) * m] = fac * i + j;
            }
        }

        return a;
    }

    public static double[] r8row_max(int m, int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_MAX returns the maximums of an R8ROW.
        //
        //  Example:
        //
        //    A =
        //      1  2  3
        //      2  6  7
        //
        //    MAX =
        //      3
        //      7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the array.
        //
        //    Input, double A[M*N], the array to be examined.
        //
        //    Output, double R8ROW_MAX[M], the maximums of the rows.
        //
    {
        int i;
        double[] amax = new double[m];

        for (i = 0; i < m; i++)
        {
            amax[i] = a[i + 0 * m];

            int j;
            for (j = 1; j < n; j++)
            {
                if (amax[i] < a[i + j * m])
                {
                    amax[i] = a[i + j * m];
                }
            }
        }

        return amax;
    }

    public static double[] r8row_mean(int m, int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_MEAN returns the means of an R8ROW.
        //
        //  Example:
        //
        //    A =
        //      1  2  3
        //      2  6  7
        //
        //    MEAN =
        //      2
        //      5
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the array.
        //
        //    Input, double A[M*N], the array to be examined.
        //
        //    Output, double R8ROW_MEAN[M], the means, or averages, of the rows.
        //
    {
        int i;
        double[] mean = new double[m];

        for (i = 0; i < m; i++)
        {
            mean[i] = 0.0;
            int j;
            for (j = 0; j < n; j++)
            {
                mean[i] += a[i + j * m];
            }

            mean[i] /= n;
        }

        return mean;
    }

    public static double[] r8row_min(int m, int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_MIN returns the minimums of an R8ROW.
        //
        //  Example:
        //
        //    A =
        //      1  2  3
        //      2  6  7
        //
        //    MIN =
        //      1
        //      2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the array.
        //
        //    Input, double A[M*N], the array to be examined.
        //
        //    Output, double R8ROW_MIN[M], the minimums of the rows.
        //
    {
        int i;
        double[] amin = new double[m];

        for (i = 0; i < m; i++)
        {
            amin[i] = a[i + 0 * m];
            int j;
            for (j = 1; j < n; j++)
            {
                if (a[i + j * m] < amin[i])
                {
                    amin[i] = a[i + j * m];
                }
            }
        }

        return amin;
    }

    public static void r8row_print(int m, int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_PRINT prints an R8ROW.
        //
        //  Discussion:
        //
        //    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
        //    each of length N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows in A.
        //
        //    Input, int N, the number of columns in A.
        //
        //    Input, double A[M*N], the M by N matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8row_print_some(m, n, a, 1, 1, m, n, title);

    }

    public static void r8row_print_some(int m, int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_PRINT_SOME prints some of an R8ROW.
        //
        //  Discussion:
        //
        //    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
        //    each of length N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, double A[M*N], the matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

        int j2lo;

        Console.WriteLine("");
        Console.WriteLine(title + "");

        if (m <= 0 || n <= 0)
        {
            Console.WriteLine("");
            Console.WriteLine("  (None)");
            return;
        }

        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            int j2hi = j2lo + INCX - 1;
            if (n < j2hi)
            {
                j2hi = n;
            }

            if (jhi < j2hi)
            {
                j2hi = jhi;
            }

            Console.WriteLine("");
            //
            //  For each column J in the current range...
            //
            //  Write the header.
            //
            string cout = "  Col:    ";
            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                cout += (j - 1).ToString(CultureInfo.InvariantCulture).PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("");
            int i2lo = ilo switch
            {
                //
                //  Determine the range of the rows in this strip.
                //
                > 1 => ilo,
                _ => 1
            };

            int i2hi;
            i2hi = ihi < m ? ihi : m;

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = (i - 1).ToString(CultureInfo.InvariantCulture).PadLeft(5) + ": ";
                for (j = j2lo; j <= j2hi; j++)
                {
                    cout += a[i - 1 + (j - 1) * m].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void r8row_reverse(int m, int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_REVERSE reverses the order of the rows of an R8MAT.
        //
        //  Discussion:
        //
        //    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
        //    each of length N.
        //
        //    To reverse the rows is to start with something like
        //
        //      11 12 13 14 15
        //      21 22 23 24 25
        //      31 32 33 34 35
        //      41 42 43 44 45
        //      51 52 53 54 55
        //
        //    and return
        //
        //      51 52 53 54 55
        //      41 42 43 44 45
        //      31 32 33 34 35
        //      21 22 23 24 25
        //      11 12 13 14 15
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, double A[M*N], the matrix whose rows are to be flipped.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m / 2; i++)
            {
                (a[i + j * m], a[m - 1 - i + j * m]) = (a[m - 1 - i + j * m], a[i + j * m]);
            }
        }
    }

    public static double[] r8row_running_average(int m, int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_RUNNING_AVERAGE computes the running averages of an R8ROW.
        //
        //  Discussion:
        //
        //    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
        //    each of length N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of items in each row.
        //
        //    Input, double V[M*N], the data.
        //
        //    Output, double R8ROW_RUNNING_AVERAGE[M*(N+1)], the running average.  
        //    A(I,J) is the average value of V(I,1:J-1).
        //
    {
        int i;
        int j;

        double[] a = new double[m * (n + 1)];
        //
        //  Sum.
        //
        for (i = 0; i < m; i++)
        {
            a[i + 0 * m] = 0.0;
            for (j = 1; j < n + 1; j++)
            {
                a[i + j * m] = a[i + (j - 1) * m] + v[i + (j - 1) * m];
            }
        }

        //
        //  Average.
        //
        for (i = 0; i < m; i++)
        {
            for (j = 1; j < n + 1; j++)
            {
                a[i + j * m] /= j;
            }
        }

        return a;
    }

    public static double[] r8row_running_sum(int m, int n, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_RUNNING_SUM computes the running sum of an R8ROW.
        //
        //  Discussion:
        //
        //    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
        //    each of length N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of items in each row.
        //
        //    Input, double V[M,N], the data.
        //
        //    Output, double R8ROW_RUNNING_SUM[M*(N+1)], the running sums.  
        //    S(I,J) is the sum of V(i,1:J-1).
        //
    {
        int i;

        double[] s = new double[m * (n + 1)];
        //
        //  Sum.
        //
        for (i = 0; i < m; i++)
        {
            s[i + 0 * m] = 0.0;
            int j;
            for (j = 1; j < n + 1; j++)
            {
                s[i + j * m] = s[i + (j - 1) * m] + v[i + (j - 1) * m];
            }
        }

        return s;
    }

    public static void r8row_sort_heap_a(int m, int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_SORT_HEAP_A ascending heapsorts an R8ROW.
        //
        //  Discussion:
        //
        //    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
        //    each of length N.
        //
        //    In lexicographic order, the statement "X < Y", applied to two real
        //    vectors X and Y of length M, means that there is some index I, with
        //    1 <= I <= M, with the property that
        //
        //      X(J) = Y(J) for J < I,
        //    and
        //      X(I) < Y(I).
        //
        //    In other words, the first time they differ, X is smaller.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, double A[M*N].
        //    On input, the array of M rows of N-vectors.
        //    On output, the rows of A have been sorted in lexicographic order.
        //
    {
        SortHeapExternalData data = new();

        switch (m)
        {
            case <= 0:
                return;
        }

        switch (n)
        {
            case <= 1:
                return;
        }

        //
        //  Initialize.
        //
        int i = 0;
        int indx = 0;
        int isgn = 0;
        int j = 0;
        //
        //  Call the external heap sorter.
        //
        for (;;)
        {
            Sort.sort_heap_external(ref data, m, ref indx, ref i, ref j, isgn);
            //
            //  Interchange the I and J objects.
            //
            if (0 < indx)
            {
                r8row_swap(m, n, ref a, i - 1, j - 1);
            }
            //
            //  Compare the I and J objects.
            //
            else if (indx < 0)
            {
                isgn = r8row_compare(m, n, a, i - 1, j - 1);
            }
            else
            {
                break;
            }
        }

    }

    public static double[] r8row_sum(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_SUM returns the sums of the rows of an R8ROW.
        //
        //  Discussion:
        //
        //    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
        //    each of length N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the M by N array.
        //
        //    Output, double ROWSUM[M], the sum of the entries of
        //    each row.
        //
    {
        int i;

        double[] rowsum = new double[m];

        for (i = 0; i < m; i++)
        {
            rowsum[i] = 0.0;
            int j;
            for (j = 0; j < n; j++)
            {
                rowsum[i] += a[i + j * m];
            }
        }

        return rowsum;
    }

    public static void r8row_swap(int m, int n, ref double[] a, int irow1, int irow2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_SWAP swaps two rows of an R8ROW.
        //
        //  Discussion:
        //
        //    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
        //    each of length N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, double A[M*N], an array of data.
        //
        //    Input, int IROW1, IROW2, the two rows to swap.
        //    These indices should be between 0 and M - 1.
        //
    {
        int j;
        //
        //  Check.
        //
        if (irow1 < 0 || m <= irow1)
        {
            Console.WriteLine("");
            Console.WriteLine("R8ROW_SWAP - Fatal error!");
            Console.WriteLine("  IROW1 is out of range.");
            return;
        }

        if (irow2 < 0 || m < irow2)
        {
            Console.WriteLine("");
            Console.WriteLine("R8ROW_SWAP - Fatal error!");
            Console.WriteLine("  IROW2 is out of range.");
            return;
        }

        if (irow1 == irow2)
        {
            return;
        }

        for (j = 0; j < n; j++)
        {
            (a[irow1 + j * m], a[irow2 + j * m]) = (a[irow2 + j * m], a[irow1 + j * m]);
        }

    }

    public static double[] r8row_to_r8vec(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_TO_R8VEC converts an R8ROW into an R8VEC.
        //
        //  Discussion:
        //
        //    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
        //    each of length N.
        //
        //  Example:
        //
        //    M = 3, N = 4
        //
        //    A =
        //      11 12 13 14
        //      21 22 23 24
        //      31 32 33 34
        //
        //    R8ROW_TO_R8VEC = ( 11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the M by N array.
        //
        //    Output, double R8ROW_TO_R8VEC[M*N], a vector containing the M rows of A.
        //
    {
        int j;

        double[] x = new double[m * n];

        int k = 0;
        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                x[k] = a[i + j * m];
                k += 1;
            }
        }

        return x;
    }

    public static void r8row_transpose_print(int m, int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_TRANSPOSE_PRINT prints an R8ROW, transposed.
        //
        //  Discussion:
        //
        //    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
        //    each of length N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], an M by N matrix to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8row_transpose_print_some(m, n, a, 1, 1, m, n, title);

    }

    public static void r8row_transpose_print_some(int m, int n, double[] a, int ilo, int jlo,
            int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_TRANSPOSE_PRINT_SOME prints some of an R8ROW, transposed.
        //
        //  Discussion:
        //
        //    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
        //    each of length N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], an M by N matrix to be printed.
        //
        //    Input, int ILO, JLO, the first row and column to print.
        //
        //    Input, int IHI, JHI, the last row and column to print.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

        int i2lo;
        int i2lo_hi;

        Console.WriteLine("");
        Console.WriteLine(title + "");

        if (m <= 0 || n <= 0)
        {
            Console.WriteLine("");
            Console.WriteLine("  (None)");
            return;
        }

        int i2lo_lo = ilo switch
        {
            < 1 => 1,
            _ => ilo
        };

        i2lo_hi = ihi < m ? m : ihi;

        for (i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo += INCX)
        {
            int i2hi = i2lo + INCX - 1;

            if (m < i2hi)
            {
                i2hi = m;
            }

            if (ihi < i2hi)
            {
                i2hi = ihi;
            }

            int inc = i2hi + 1 - i2lo;

            Console.WriteLine("");
            string cout = "  Row: ";
            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                cout += (i - 1).ToString(CultureInfo.InvariantCulture).PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Col");
            Console.WriteLine("");

            int j2lo = jlo switch
            {
                < 1 => 1,
                _ => jlo
            };

            int j2hi;
            j2hi = n < jhi ? n : jhi;

            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                cout = (j - 1).ToString(CultureInfo.InvariantCulture).PadLeft(5) + ":";
                int i2;
                for (i2 = 1; i2 <= inc; i2++)
                {
                    i = i2lo - 1 + i2;
                    cout += a[i - 1 + (j - 1) * m].ToString(CultureInfo.InvariantCulture).PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8row_variance(int m, int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_VARIANCE returns the variances of an R8ROW.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 October 2004
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
        //    Output, double R8ROW_VARIANCE[M], the variances of the rows.
        //
    {
        int i;
        double[] variance = new double[m];

        for (i = 0; i < m; i++)
        {
            double mean = 0.0;
            int j;
            for (j = 0; j < n; j++)
            {
                mean += a[i + j * m];
            }

            mean /= n;

            variance[i] = 0.0;
            for (j = 0; j < n; j++)
            {
                variance[i] += Math.Pow(a[i + j * m] - mean, 2);
            }

            switch (n)
            {
                case > 1:
                    variance[i] /= n - 1;
                    break;
                default:
                    variance[i] = 0.0;
                    break;
            }

        }

        return variance;
    }

    public static double[] r8rows_to_r8mat(int m, int n, double[] r8rows)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROWS_TO_R8MAT converts a row-major vector to an R8MAT.
        //
        //  Discussion:
        //
        //    This function allows me to declare a vector of the right type and length,
        //    fill it with data that I can display row-wise, and then have the
        //    data copied into a column-wise doubly dimensioned array array.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the array.
        //
        //    Input, double R8ROWS[M*N], the rowwise data.
        //
        //    Output, double R8ROWS_TO_R8MAT[M*N], the doubly-dimensioned columnwise data.
        //
    {
        int i;

        double[] r8mat = new double[m * n];

        int k = 0;
        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                r8mat[i + j * m] = r8rows[k];
                k += 1;
            }
        }

        return r8mat;
    }
}