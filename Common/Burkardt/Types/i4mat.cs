using System;
using System.IO;
using System.Linq;
using Burkardt.Table;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void i4mat_01_rowcolsum(int m, int n, int[] r, int[] c, ref int[] a, ref bool error)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_01_ROWCOLSUM creates a 0/1 I4MAT with given row and column sums.
//
//  Discussion:
//
//    Given an M vector R and N vector C, there may exist one or more
//    M by N matrices with entries that are 0 or 1, whose row sums are R
//    and column sums are C.
//
//    For convenience, this routine requires that the entries of R and C
//    be given in nonincreasing order.
//
//    There are several requirements on R and C.  The simple requirements
//    are that the entries of R and C must be nonnegative, that the entries
//    of R must each be no greater than N, and those of C no greater than M,
//    and that the sum of the entries of R must equal the sum of the entries 
//    of C.
//
//    The final technical requirement is that if we form R*, the conjugate
//    partition of R, then C is majorized by R*, that is, that every partial
//    sum from 1 to K of the entries of C is no bigger than the sum of the same
//    entries of R*, for every K from 1 to N.
//
//    Given these conditions on R and C, there is at least one 0/1 matrix
//    with the given row and column sums.
//
//    The conjugate partition of R is constructed as follows:
//      R*(1) is the number of entries of R that are 1 or greater.
//      R*(2) is the number of entries of R that are 2 or greater.
//      ...
//      R*(N) is the number of entries of R that are N (can't be greater).
//
//  Example:
//
//    M = N = 5
//    R = ( 3, 2, 2, 1, 1 )
//    C = ( 2, 2, 2, 2, 1 )
//
//    A =
//      1 0 1 0 1
//      1 0 0 1 0
//      0 1 0 1 0
//      0 1 0 0 0
//      0 0 1 0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack van Lint, Richard Wilson,
//    A Course in Combinatorics,
//    Oxford, 1992, pages 148-156.
//
//    James Sandeson,
//    Testing Ecobool Patterns,
//    American Scientist,
//    Volume 88, July-August 2000, pages 332-339.
//
//    Ian Saunders,
//    Algorithm AS 205,
//    Enumeration of R x C Tables with Repeated Row Totals,
//    Applied Statistics,
//    Volume 33, Number 3, pages 340-352, 1984.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, int R[M], C[N], the row and column sums desired for the array.
//    Both vectors must be arranged in descending order.
//    The elements of R must be between 0 and N.
//    The elements of C must be between 0 and M.
//
//    Output, int A[M*N], the M by N matrix with the given row and
//    column sums.
//    Each entry of A is 0 or 1.
//
//    Output, bool &ERROR, is true if an error occurred.
//
    {
        int i;
        int j;
//
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                a[i + j * m] = 0;
            }
        }

//
//  Check conditions.
//
        error = false;

        if (i4vec_sum(m, r) != i4vec_sum(n, c))
        {
            Console.WriteLine("");
            Console.WriteLine("I4MAT_01_ROWCOLSUM - Fatal error!");
            Console.WriteLine("  Row sums R and column sums C don't have the same sum!");
            error = true;
            return;
        }

        if (!i4vec_descends(m, ref r))
        {
            Console.WriteLine("");
            Console.WriteLine("I4MAT_01_ROWCOLSUM - Fatal error!");
            Console.WriteLine("  Row sum vector R is not descending!");
            error = true;
            return;
        }

        if (n < r[0] || r[m - 1] < 0)
        {
            error = true;
            return;
        }

        if (!i4vec_descends(n, ref c))
        {
            Console.WriteLine("");
            Console.WriteLine("I4MAT_01_ROWCOLSUM - Fatal error!");
            Console.WriteLine("  Column sum vector C is not descending!");
            error = true;
            return;
        }

        if (m < c[0] || c[n - 1] < 0)
        {
            error = true;
            return;
        }

//
//  Compute the conjugate of R.
//
        int[] r_conj = new int[n];

        for (i = 0; i < n; i++)
        {
            r_conj[i] = 0;
        }

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < r[i]; j++)
            {
                r_conj[j] += 1;
            }
        }

//
//  C must be majorized by R_CONJ.
//
        int r_sum = 0;
        int c_sum = 0;
        for (i = 0; i < n; i++)
        {
            r_sum += r_conj[i];
            c_sum += c[i];
            if (r_sum >= c_sum)
            {
                continue;
            }

            error = true;
            return;
        }

        switch (error)
        {
            case true:
                return;
        }

        int[] r2 = new int[m];

//
//  We need a temporary copy of R that we can decrement.
//
        for (i = 0; i < m; i++)
        {
            r2[i] = r[i];
        }

        for (j = n - 1; 0 <= j; j--)
        {
            i = i4vec_maxloc_last(m, r2);

            int k;
            for (k = 1; k <= c[j]; k++)
            {
//
//  By adding 1 rather than setting A(I,J) to 1, we were able to spot
//  an error where the index was "sticking".
//
                a[i + j * m] += 1;

                r2[i] -= 1;

                switch (i)
                {
                    case > 0:
                        i -= 1;
                        break;
//
                    default:
                    {
                        i = i4vec_maxloc_last(m, r2);
                        switch (i)
                        {
                            case 0 when k < c[j]:
                                i = 1 + i4vec_maxloc_last(m - 1, r2, xIndex: +1);
                                break;
                        }

                        break;
                    }
                }
            }
        }

    }

    public static void i4mat_01_rowcolsum2(int m, int n, int[] r, int[] c, ref int[] a,
        ref bool error)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_01_ROWCOLSUM2 creates a 0/1 I4MAT with given row and column sums.
//
//  Discussion:
//
//    This routine uses network flow optimization to compute the results.
//
//    Given an M vector R and N vector C, there may exist one or more
//    M by N matrices with entries that are 0 or 1, whose row sums are R
//    and column sums are C.
//
//    For convenience, this routine requires that the entries of R and C
//    be given in nonincreasing order.
//
//    There are several requirements on R and C.  The simple requirements
//    are that the entries of R and C must be nonnegative, that the entries
//    of R must each no greater than N, and those of C no greater than M,
//    and that the sum of the entries of R must equal the sum of the 
//    entries of C.
//
//    The final technical requirement is that if we form R*, the conjugate
//    partition of R, then C is majorized by R*, that is, that every partial
//    sum from 1 to K of the entries of C is no bigger than the sum of the same
//    entries of R*, for every K from 1 to N.
//
//    Given these conditions on R and C, there is at least one 0/1 matrix
//    with the given row and column sums.
//
//    The conjugate partition of R is constructed as follows:
//      R*(1) is the number of entries of R that are 1 or greater.
//      R*(2) is the number of entries of R that are 2 or greater.
//      ...
//      R*(N) is the number of entries of R that are N (can't be greater).
//
//  Example:
//
//    M = N = 5
//    R = ( 3, 2, 2, 1, 1 )
//    C = ( 2, 2, 2, 2, 1 )
//
//    A =
//      1 0 1 0 1
//      1 0 0 1 0
//      0 1 0 1 0
//      0 1 0 0 0
//      0 0 1 0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//    Jack van Lint, Richard Wilson,
//    A Course in Combinatorics,
//    Oxford, 1992, pages 148-156.
//
//    James Sandeson,
//    Testing Ecobool Patterns,
//    American Scientist,
//    Volume 88, July-August 2000, pages 332-339.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//    These values do not have to be equal.
//
//    Input, int R[M], C[N], the row and column sums desired for the array.
//    Both vectors must be arranged in descending order.
//    The elements of R must be between 0 and N.
//    The elements of C must be between 0 and M.
//    One of the conditions for a solution to exist is that the sum of the
//    elements in R equal the sum of the elements in C.
//
//    Output, int A[M*N], the matrix with the given row and column sums.
//    Each entry of A is 0 or 1.
//
//    Output, bool &ERROR, is true if an error occurred.
//
    {
        int i;
        int j;

        error = false;

        int[] capflo = new int[2 * 2 * (m + m * n + n)];
        int[] iendpt = new int[2 * 2 * (m + m * n + n)];
//
//  There are M + N + 2 nodes.  The last two are the special source and sink.
//
        int source = m + n + 1;
        int isink = m + n + 2;
//
//  The source is connected to each of the R nodes.
//
        int k = 0;

        for (i = 0; i < m; i++)
        {
            iendpt[0 + 2 * k] = source;
            iendpt[1 + 2 * k] = i + 1;
            capflo[0 + 2 * k] = r[i];
            capflo[1 + 2 * k] = 0;
            k += 1;

            iendpt[0 + 2 * k] = i + 1;
            iendpt[1 + 2 * k] = source;
            capflo[0 + 2 * k] = r[i];
            capflo[1 + 2 * k] = 0;
            k += 1;
        }

//
//  Every R node is connected to every C node, with capacity 1.
//
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                iendpt[0 + 2 * k] = i + 1;
                iendpt[1 + 2 * k] = j + 1 + m;
                capflo[0 + 2 * k] = 1;
                capflo[1 + 2 * k] = 0;
                k += 1;

                iendpt[0 + 2 * k] = j + 1 + m;
                iendpt[1 + 2 * k] = i + 1;
                capflo[0 + 2 * k] = 1;
                capflo[1 + 2 * k] = 0;
                k += 1;
            }
        }

//
//  Every C node is connected to the sink.
//
        for (j = 0; j < n; j++)
        {
            iendpt[0 + 2 * k] = j + 1 + m;
            iendpt[1 + 2 * k] = isink;
            capflo[0 + 2 * k] = c[j];
            capflo[1 + 2 * k] = 0;
            k += 1;

            iendpt[0 + 2 * k] = isink;
            iendpt[1 + 2 * k] = j + 1 + m;
            capflo[0 + 2 * k] = c[j];
            capflo[1 + 2 * k] = 0;
            k += 1;
        }

//
//  Determine the maximum flow on the network.
//
        int nedge = k;

//network_flow_max(nnode, nedge, iendpt, capflo, source, isink,
//    icut, node_flow);
//
//  We have a perfect solution if, and only if, the edges leading from the
//  source, and the edges leading to the sink, are all saturated.
//
        for (k = 0; k < nedge; k++)
        {
            i = iendpt[0 + 2 * k];
            j = iendpt[1 + 2 * k] - m;

            if (i <= m && 1 <= j && j <= n)
            {
                if (capflo[1 + 2 * k] != 0 && capflo[1 + 2 * k] != 1)
                {
                    error = true;
                }
            }

            if (iendpt[0 + 2 * k] == source)
            {
                if (capflo[0 + 2 * k] != capflo[1 + 2 * k])
                {
                    error = true;
                }
            }

            if (iendpt[1 + 2 * k] != isink)
            {
                continue;
            }

            if (capflo[0 + 2 * k] != capflo[1 + 2 * k])
            {
                error = true;
            }

        }

//
//  If we have a solution, then A(I,J) = the flow on the edge from
//  R node I to C node J.
//
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                a[i + j * m] = 0;
            }
        }

        for (k = 0; k < nedge; k++)
        {
            i = iendpt[0 + 2 * k];
            j = iendpt[1 + 2 * k] - m;

            if (i <= m && 1 <= j && j <= n)
            {
                a[i + j * m] = capflo[1 + 2 * k];
            }
        }
    }

    public static int i4mat_max(int m, int n, int[] a)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_MAX returns the maximum of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
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
//    Input, int A[M*N], the M by N matrix.
//
//    Output, int I4MAT_MAX, the maximum entry of A.
//
    {
        int j;

        int value = -i4_huge();

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                if (value < a[i + j * m])
                {
                    value = a[i + j * m];
                }
            }
        }

        return value;
    }

    public static int i4mat_min(int m, int n, int[] a)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_MIN returns the minimum of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2009
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
//    Input, int A[M*N], the M by N matrix.
//
//    Output, int I4MAT_MIN, the minimum entry of A.
//
    {
        int j;

        int value = i4_huge();

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                if (a[i + j * m] < value)
                {
                    value = a[i + j * m];
                }
            }
        }

        return value;
    }

    public static int[] i4mat_histogram(int m, int n, int[] a, int histo_num)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_HISTOGRAM computes a histogram of the elements of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an array of I4's.
//
//    It is assumed that the entries in the vector A are nonnegative.
//    Only values between 0 and HISTO_NUM will be histogrammed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input, int A[M*N], the array to examine.
//
//    Input, int HISTO_NUM, the maximum value for which a
//    histogram entry will be computed.
//
//    Output, int I4MAT_HISTOGRAM[HISTO_NUM+1], contains the number of
//    entries of A with the values of 0 through HISTO_NUM.
//
    {
        int i;
        int j;

        int[] histo_gram = new int[histo_num + 1];

        for (i = 0; i <= histo_num; i++)
        {
            histo_gram[i] = 0;
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                switch (a[i + j * m])
                {
                    case >= 0 when a[i + j * m] <= histo_num:
                        histo_gram[a[i + j * m]] += 1;
                        break;
                }
            }
        }

        return histo_gram;
    }

    public static int[] i4mat_border_add(int m, int n, int[] table)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_BORDER_ADD adds a "border" to an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an array of I4's.
//
//    We suppose the input data gives values of a quantity on nodes
//    in the interior of a 2D grid, and we wish to create a new table
//    with additional positions for the nodes that would be on the
//    border of the 2D grid.
//
//                  0 0 0 0 0 0
//      * * * *     0 * * * * 0
//      * * * * --> 0 * * * * 0
//      * * * *     0 * * * * 0
//                  0 0 0 0 0 0
//
//    The illustration suggests the situation in which a 3 by 4 array
//    is input, and a 5 by 6 array is to be output.
//
//    The old data is shifted to its correct positions in the new array.
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
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, int TABLE[M*N], the data.
//
//    Output, int TABLE2[(M+2)*(N+2)], the augmented data.
//
    {
        int[] table2 = new int[(m + 2) * (n + 2)];

        for (int j = 0; j < n + 2; j++)
        {
            for (int i = 0; i < m + 2; i++)
            {
                if (i == 0 || i == m + 1 || j == 0 || j == n + 1)
                {
                    table2[i + j * (m + 2)] = 0;
                }
                else
                {
                    table2[i + j * (m + 2)] = table[i - 1 + (j - 1) * m];
                }
            }
        }

        return table2;
    }

    public static void i4mat_u1_inverse(int n, int[] a, ref int[] b)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_U1_INVERSE inverts a unit upper triangular I4MAT.
//
//  Discussion:
//
//    A unit upper triangular matrix is a matrix with only 1's on the main
//    diagonal, and only 0's below the main diagonal.  Above the main
//    diagonal, the entries may be assigned any value.
//
//    It may be surprising to note that the inverse of an integer unit upper
//    triangular matrix is also an integer unit upper triangular matrix.
//
//    Note that this routine can invert a matrix in place, that is, with no
//    extra storage.  If the matrix is stored in A, then the call
//
//      i4mat_u1_inverse ( n, a, a )
//
//    will result in A being overwritten by its inverse, which can
//    save storage if the original value of A is not needed later.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of rows and columns in the matrix.
//
//    Input, int A[N*N], the unit upper triangular matrix
//    to be inverted.
//
//    Output, int B[N*N], the inverse matrix.
//
    {
        int j;

        for (j = n; 1 <= j; j--)
        {
            int i;
            for (i = n; 1 <= i; i--)
            {
                int isum = i == j ? 1 : 0;

                int k;
                for (k = i + 1; k <= j; k++)
                {
                    isum -= a[i - 1 + (k - 1) * n] * b[k - 1 + (j - 1) * n];
                }

                b[i - 1 + (j - 1) * n] = isum;
            }
        }
    }

    public static void i4mat_write(string output_filename, int m, int n, int[] table)
//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_WRITE writes an I4MAT file with no header.
//
//  Discussion:
//
//    An I4MAT is an array of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, int TABLE[M*N], the data.
//
    {
        string[] outData = new string[n];
        for (int j = 0; j < n; j++)
        {
            string line = "";
            for (int i = 0; i < m; i++)
            {
                line += "  " + table[i + j * m].ToString().PadLeft(10);
            }

            outData[j] = line;
        }

        try
        {
            File.WriteAllLines(output_filename, outData);
        }
        catch (Exception)
        {
            Console.WriteLine();
            Console.WriteLine("I4MAT_WRITE - Fatal error!");
            Console.WriteLine("  Could not open the output file: \"" + output_filename + "\"");
            throw;
        }
    }

    public static void i4mat_flip_cols(int m, int n, ref int[] a)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_FLIP_COLS swaps the columns of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//    To "flip" the columns of an I4MAT is to start with something like
//
//      11 12 13 14 15
//      21 22 23 24 25
//      31 32 33 34 35
//      41 42 43 44 45
//      51 52 53 54 55
//
//    and return
//
//      15 14 13 12 11
//      25 24 23 22 21
//      35 34 33 32 31
//      45 44 43 42 41
//      55 54 53 52 51
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int A[M*N], the matrix whose columns are to be flipped.
//
    {
        int i;

        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 0; j < n / 2; j++)
            {
                (a[i + j * m], a[i + (n - 1 - j) * m]) = (a[i + (n - 1 - j) * m], a[i + j * m]);
            }
        }
    }

    public static void i4mat_flip_rows(int m, int n, ref int[] a)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_FLIP_ROWS swaps the rows of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//    To "flip" the rows of an I4MAT is to start with something like
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
//    22 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int A[M*N], the matrix whose rows are to be flipped.
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

    public static TableHeader i4mat_header_read(string input_filename)
//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_HEADER_READ reads the header from an I4MAT file.
//
//  Discussion:
//
//    An I4MAT is an array of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int &M, the number of spatial dimensions.
//
//    Output, int &N, the number of points
//
    {
        TableHeader ret = TableMisc.readHeader(input_filename);

        switch (ret.m)
        {
            case <= 0:
                Console.WriteLine();
                Console.WriteLine("I4MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_COLUMN_COUNT failed.");
                ret.code = 1;
                break;
        }

        switch (ret.n)
        {
            case <= 0:
                Console.WriteLine();
                Console.WriteLine("I4MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_ROW_COUNT failed.");
                ret.code = 1;
                break;
        }

        return ret;
    }


    public static int[] i4mat_data_read(string input_filename, int m, int n)
//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_DATA_READ reads data from an I4MAT file.
//
//  Discussion:
//
//    An I4MAT is an array of I4's.
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with '#' are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly (or at least)
//    M real numbers, representing the coordinates of a point.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, int I4MAT_DATA_READ[M*N], the data.
//
    {
        string[] lines;

        try
        {
            lines = File.ReadLines(input_filename).ToArray();
        }
        catch (Exception)
        {
            Console.WriteLine();
            Console.WriteLine("R8MAT_DATA_READ - Fatal error!");
            Console.WriteLine("  Could not open the input file: \"" + input_filename + "\"");
            throw;
        }

        int[] table = new int[m * n];

        int j = 0;
        int l = 0;

        while (j < n)
        {
            string line = lines[l];
            l++;
            if (line[0] == '#' || s_len_trim(line) == 0)
            {
                continue;
            }

            i4vec res = s_to_i4vec(line, m);

            bool error = res.error;
            int[] x = res.ivec;

            switch (error)
            {
                case false:
                {
                    int i;
                    for (i = 0; i < m; i++)
                    {
                        table[i + j * m] = x[i];
                    }

                    break;
                }
            }

            j += 1;
        }

        return table;
    }

    public static int[] i4mat_border_cut(int m, int n, int[] table)
//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_BORDER_CUT cuts the "border" of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an array of I4's.
//
//    We suppose the input data gives values of a quantity on nodes
//    on a 2D grid, and we wish to create a new table corresponding only
//    to those nodes in the interior of the 2D grid.
//
//      0 0 0 0 0 0
//      0 * * * * 0    * * * *
//      0 * * * * 0 -> * * * *
//      0 * * * * 0    * * * *
//      0 0 0 0 0 0
//
//    The illustration suggests the situation in which a 5 by 6 array
//    is input, and a 3 by 4 array is to be output.
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
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, int TABLE[M*N], the data.
//
//    Output, int TABLE2[(M-2)*(N-2)], the "interior" data.
//
    {
        if (m <= 2 || n <= 2)
        {
            return Array.Empty<int>();
        }

        int[] table2 = new int[(m - 2) * (n - 2)];

        for (int j = 0; j < n - 2; j++)
        {
            for (int i = 0; i < m - 2; i++)
            {
                table2[i + j * (m - 2)] = table[i + 1 + (j + 1) * m];
            }
        }

        return table2;
    }

    public static void i4mat_shortest_path(int n, ref int[] m)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_SHORTEST_PATH computes the shortest distance between all pairs of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Floyd,
//    Algorithm 97, Shortest Path,
//    Communications of the ACM,
//    Volume 5, Number 6, June 1962, page 345.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input/output, int M[N*N].
//    On input, M(I,J) contains the length of the direct link between 
//    nodes I and J, or HUGE if there is no direct link.
//    On output, M(I,J) contains the distance between nodes I and J,
//    that is, the length of the shortest path between them.  If there
//    is no such path, then M(I,J) will remain HUGE.
//
    {
        int i;
        const int i4_inf = 2147483647;

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                switch (m[j + i * n])
                {
                    case < i4_inf:
                    {
                        int k;
                        for (k = 0; k < n; k++)
                        {
                            switch (m[i + k * n])
                            {
                                case < i4_inf:
                                {
                                    int s = m[j + i * n] + m[i + k * n];
                                    if (s < m[j + k * n])
                                    {
                                        m[j + k * n] = s;
                                    }

                                    break;
                                }
                            }
                        }

                        break;
                    }
                }
            }
        }
    }

    public static int i4mat_sum(int m, int n, ref int[] a)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_SUM returns the sum of the entries of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2018
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
//    Input, int A[M*N], the M by N matrix.
//
//    Output, int I4MAT_SUM, the sum of the entries.
//
    {
        int j;

        int value = 0;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                value += a[i + j * m];
            }
        }

        return value;
    }

    public static void i4mat_transpose(int m, int n, ref int[] a)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE transposes an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2018
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
//    Input, int A[M*N], the M by N matrix.
//
//    Output, int A[N*M], the transposed matrix.
//
    {
        int i;
        int j;

        int[] b = new int[n * m];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                b[j + i * n] = a[i + j * m];
            }
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                a[j + i * n] = b[j + i * n];
            }
        }
    }

    public static void i4mat_transpose_print(int m, int n, int[] a, string title)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 January 2005
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
//    Input, int A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
    {
        i4mat_transpose_print_some(m, n, a, 1, 1, m, n, title);
    }

    public static void i4mat_transpose_print_some(int m, int n, int[] a, int ilo, int jlo,
        int ihi, int jhi, string title)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2014
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
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
    {
        const int INCX = 10;

        Console.WriteLine();
        Console.WriteLine(title);

        if (m <= 0 || n <= 0)
        {
            Console.WriteLine();
            Console.WriteLine("  (None)");
            return;
        }

//
//  Print the columns of the matrix, in strips of INCX.
//
        for (int i2lo = ilo; i2lo <= ihi; i2lo += INCX)
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

            Console.WriteLine();
//
//  For each row I in the current range...
//
//  Write the header.
//
            string line = "  Row: ";
            for (int i = i2lo; i <= i2hi; i++)
            {
                string t = i - 1 + "  ";
                line += t.PadLeft(6);
            }

            Console.WriteLine(line);
            Console.WriteLine();
            Console.WriteLine("  Col");
            Console.WriteLine();
//
//  Determine the range of the rows in this strip.
//
            int j2lo = jlo;
            j2lo = j2lo switch
            {
                < 1 => 1,
                _ => j2lo
            };

            int j2hi = jhi;
            if (n < j2hi)
            {
                j2hi = n;
            }

            for (int j = j2lo; j <= j2hi; j++)
            {
//
//  Print out (up to INCX) entries in column J, that lie in the current strip.
//
                string t = j - 1 + ":";
                line = "";
                line += t.PadLeft(5);

                for (int i = i2lo; i <= i2hi; i++)
                {
                    t = a[i - 1 + (j - 1) * m] + "  ";
                    line += t.PadLeft(6);
                }

                Console.WriteLine(line);

                Console.WriteLine();
            }
        }
    }

    public static void i4mat_perm0(int n, ref int[] a, int[] p)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PERM0 permutes the rows and columns of a square I4MAT.
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
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, int A[N*N].
//    On input, the matrix to be permuted.
//    On output, the permuted matrix.
//
//    Input, int P[N], the permutation.  P(I) is the new number of row
//    and column I.
//
    {
        int i;
        int is_ = 0;
        int nc = 0;

        Permutation.perm0_cycle(n, p, ref is_, ref nc, 1);
//
//  Temporarily increment P by 1.
//
        for (i = 0; i < n; i++)
        {
            p[i] += 1;
        }

        for (i = 1; i <= n; i++)
        {
            int i1 = -p[i - 1];

            switch (i1)
            {
                case > 0:
                {
                    int lc = 0;

                    for (;;)
                    {
                        i1 = p[i1 - 1];
                        lc += 1;

                        if (i1 <= 0)
                        {
                            break;
                        }
                    }

                    i1 = i;

                    int j;
                    for (j = 1; j <= n; j++)
                    {
                        switch (p[j - 1])
                        {
                            case <= 0:
                            {
                                int j2 = j;
                                int k = lc;
                                for (;;)
                                {
                                    int j1 = j2;
                                    int it = a[i1 - 1 + (j1 - 1) * n];

                                    for (;;)
                                    {
                                        i1 = Math.Abs(p[i1 - 1]);
                                        j1 = Math.Abs(p[j1 - 1]);

                                        (a[i1 - 1 + (j1 - 1) * n], it) = (it, a[i1 - 1 + (j1 - 1) * n]);

                                        if (j1 != j2)
                                        {
                                            continue;
                                        }

                                        k -= 1;

                                        if (i1 == i)
                                        {
                                            break;
                                        }
                                    }

                                    j2 = Math.Abs(p[j2 - 1]);

                                    if (k == 0)
                                    {
                                        break;
                                    }
                                }

                                break;
                            }
                        }
                    }

                    break;
                }
            }
        }

//
//  Restore the positive signs of the data.
//
        for (i = 0; i < n; i++)
        {
            p[i] = Math.Abs(p[i]);
        }

//
//  Decrement P by 1.
//
        for (i = 0; i < n; i++)
        {
            p[i] -= 1;
        }

    }

    public static void i4mat_2perm0(int m, int n, int[] a, int[] p, int[] q)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_2PERM0 permutes the rows and columns of a rectangular I4MAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int M, number of rows in the matrix.
//
//    Input, int N, number of columns in the matrix.
//
//    Input/output, int A[M*N].
//    On input, the matrix to be permuted.
//    On output, the permuted matrix.
//
//    Input, int P[M], the row permutation.  P(I) is the new number of row I.
//
//    Input, int Q[N].  The column permutation.  Q(I) is the new number
//    of column I.  Note that this routine allows you to pass a single array
//    as both P and Q.
//
    {
        int i;
        int is_ = 0;
        int j;
        int nc = 0;
//
//  Wretched maneuvers to deal with necessity of 1-based values,
//  and to handle case where P and Q are same vector.
//
        int[] p1 = i4vec_copy_new(m, p);
        Permutation.perm0_cycle(m, p1, ref is_, ref nc, 1);
        for (i = 0; i < m; i++)
        {
            p1[i] += 1;
        }

        int[] q1 = i4vec_copy_new(n, q);
        Permutation.perm0_cycle(n, q1, ref is_, ref nc, 1);
        for (j = 0; j < n; j++)
        {
            q1[j] += 1;
        }

        for (i = 1; i <= m; i++)
        {
            int i1 = -p1[i - 1];

            switch (i1)
            {
                case > 0:
                {
                    int lc = 0;

                    for (;;)
                    {
                        i1 = p1[i1 - 1];
                        lc += 1;

                        if (i1 <= 0)
                        {
                            break;
                        }
                    }

                    i1 = i;

                    for (j = 1; j <= n; j++)
                    {
                        switch (q1[j - 1])
                        {
                            case <= 0:
                            {
                                int j2 = j;
                                int k = lc;

                                for (;;)
                                {
                                    int j1 = j2;
                                    int it = a[i1 - 1 + (j1 - 1) * m];

                                    for (;;)
                                    {
                                        i1 = Math.Abs(p1[i1 - 1]);
                                        j1 = Math.Abs(q1[j1 - 1]);

                                        (it, a[i1 - 1 + (j1 - 1) * m]) = (a[i1 - 1 + (j1 - 1) * m], it);

                                        if (j1 != j2)
                                        {
                                            continue;
                                        }

                                        k -= 1;

                                        if (i1 == i)
                                        {
                                            break;
                                        }
                                    }

                                    j2 = Math.Abs(q1[j2 - 1]);

                                    if (k == 0)
                                    {
                                        break;
                                    }
                                }

                                break;
                            }
                        }
                    }

                    break;
                }
            }
        }
    }

    public static void i4mat_print(int m, int n, int[] a, string title)
//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT prints an I4MAT, with an optional title.
//
//  Discussion:
//
//    An I4MAT is an array of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 2003
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
//    Input, int A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
    {
        i4mat_print_some(m, n, a, 1, 1, m, n, title);
    }

    public static void i4mat_print_some(int m, int n, int[] a, int ilo, int jlo, int ihi,
        int jhi, string title)
//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT_SOME prints some of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an array of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 2004
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
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
    {
        const int INCX = 10;


        Console.WriteLine();
        Console.WriteLine(title);
//
//  Print the columns of the matrix, in strips of INCX.
//
        for (int j2lo = jlo; j2lo <= jhi; j2lo += INCX)
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

            Console.WriteLine();
//
//  For each column J in the current range...
//
//  Write the header.
//
            string cout = "  Col: ";
            for (int j = j2lo; j <= j2hi; j++)
            {
                cout += j.ToString().PadLeft(6) + "  ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
//
//  Determine the range of the rows in this strip.
//
            int i2lo = ilo switch
            {
                > 1 => 1,
                _ => ilo
            };

            int i2hi = ihi < m ? ihi : m;

            for (int i = i2lo; i <= i2hi; i++)
            {
//
//  Print out (up to INCX) entries in row I, that lie in the current strip.
//
                cout = i.ToString().PadLeft(5) + "  ";
                for (int j = j2lo; j <= j2hi; j++)
                {
                    cout += a[i - 1 + (j - 1) * m].ToString().PadLeft(6) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }


    public static int[] i4mat_indicator_new(int m, int n)
//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_INDICATOR_NEW sets up an "indicator" I4MAT.
//
//  Discussion:
//
//    An I4MAT is an array of I4's.
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
//    Output, int I4MAT_INDICATOR_NEW[M*N], the indicator matrix.
//
    {
        int[] table = new int[m * n];

        int fac = (int) Math.Pow(10.0, Math.Floor(Math.Log10(n) + 1));

        for (int i = 1; i <= m; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                table[i - 1 + (j - 1) * m] = fac * i + j;
            }
        }

        return table;
    }

    public static void i4mat_copy(int m, int n, int[] a1, ref int[] a2)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_COPY copies one I4MAT to another.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A1[M*N], the matrix to be copied.
//
//    Output, int A2[M*N], the copy of A1.
//
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a2[i + j * m] = a1[i + j * m];
            }
        }
    }

    public static int[] i4mat_copy_new(int m, int n, int[] a1)

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_COPY_NEW copies an I4MAT to a "new" I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A1[M*N], the matrix to be copied.
//
//    Output, int I4MAT_COPY_NEW[M*N], the copy of A1.
//
    {
        int j;

        int[] a2 = new int[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a2[i + j * m] = a1[i + j * m];
            }
        }

        return a2;
    }
}
