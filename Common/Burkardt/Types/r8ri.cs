using System;
using System.Globalization;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8ri_dif2(int n, int nz, ref int[] ija, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RI_DIF2 stores the second difference matrix in R8RI format.
        //
        //  Discussion:
        //
        //    An R8RI matrix is in row indexed sparse storage form, using an index
        //    array IJA and a value array A.  The first N entries of A store the
        //    diagonal elements in order.  The first N entries of IJA store the index
        //    of the first off-diagonal element of the corresponding row; if there is
        //    no off-diagonal element in that row, it is one greater than the index
        //    in A of the most recently stored element in the previous row.
        //    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
        //    greater than the index in A of the last off-diagonal element of the
        //    last row.  Location N+1 of A is not used.  Entries in A with index
        //    N+2 or greater contain the off-diagonal values, ordered by row, and
        //    then by column.  Entries in IJA with index N+2 or greater contain the
        //    column number of the corresponding element in A.
        //
        //  Example:
        //
        //    A:
        //      3 0 1 0 0
        //      0 4 0 0 0
        //      0 7 5 9 0
        //      0 0 0 0 2
        //      0 0 0 6 8
        //
        //    NZ = 11
        //
        //    IJA:
        //      7  8  8 10 11 12  3  2  4  5  4
        //
        //    A:
        //      3  4  5  0  8  *  1  7  9  2  6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
        //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
        //    Third Edition,
        //    Cambridge University Press, 2007,
        //    ISBN13: 978-0-521-88068-8,
        //    LC: QA297.N866.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the size required for the RI
        //    or "row indexed" sparse storage.  NZ = 3*N-1.
        //
        //    Output, int IJA[NZ], the index vector.
        //
        //    Output, double A[NZ], the value vector.
        //
    {
        int i;
        //
        //  Diagonal elements of A.
        //
        for (i = 0; i < n; i++)
        {
            a[i] = 2.0;
        }

        //
        //  First N entries of IJA store first offdiagonal of each row.
        //
        int k = n + 1;

        for (i = 0; i < n; i++)
        {
            ija[i] = k;
            if (i == 0 || i == n - 1)
            {
                k += 1;
            }
            else
            {
                k += 2;
            }
        }

        //
        //  IJA(N+1) stores one beyond last element of A.
        //
        ija[n] = k;
        a[n] = 0.0;
        //
        //  IJA(N+2), A(N+2) and beyond store column and value.
        //
        k = n;

        for (i = 0; i < n; i++)
        {
            switch (i)
            {
                case 0:
                    k += 1;
                    ija[k] = i + 1;
                    a[k] = -1.0;
                    break;
                default:
                {
                    if (i < n - 1)
                    {
                        k += 1;
                        ija[k] = i - 1;
                        a[k] = -1.0;
                        k += 1;
                        ija[k] = i + 1;
                        a[k] = -1.0;
                    }
                    else if (i == n - 1)
                    {
                        k += 1;
                        ija[k] = i - 1;
                        a[k] = -1.0;
                    }

                    break;
                }
            }
        }

        return a;
    }

    public static double[] r8ri_indicator(int n, int nz, int[] ija)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RI_INDICATOR returns the R8RI indicator matrix for given sparsity.
        //
        //  Discussion:
        //
        //    An R8RI matrix is in row indexed sparse storage form, using an index
        //    array IJA and a value array A.  The first N entries of A store the
        //    diagonal elements in order.  The first N entries of IJA store the index
        //    of the first off-diagonal element of the corresponding row; if there is
        //    no off-diagonal element in that row, it is one greater than the index
        //    in A of the most recently stored element in the previous row.
        //    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
        //    greater than the index in A of the last off-diagonal element of the
        //    last row.  Location N+1 of A is not used.  Entries in A with index
        //    N+2 or greater contain the off-diagonal values, ordered by row, and
        //    then by column.  Entries in IJA with index N+2 or greater contain the
        //    column number of the corresponding element in A.
        //
        //  Example:
        //
        //    A:
        //      3 0 1 0 0
        //      0 4 0 0 0
        //      0 7 5 9 0
        //      0 0 0 0 2
        //      0 0 0 6 8
        //
        //    NZ = 11
        //
        //    IJA:
        //      7  8  8 10 11 12  3  2  4  5  4
        //
        //    A:
        //      3  4  5  0  8  *  1  7  9  2  6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
        //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
        //    Third Edition,
        //    Cambridge University Press, 2007,
        //    ISBN13: 978-0-521-88068-8,
        //    LC: QA297.N866.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the size required for the RI
        //    or "row indexed" sparse storage.  NZ = 3*N-1.
        //
        //    Input, int IJA[NZ], the index vector.
        //
        //    Output, double R8RI_INDICATOR[NZ], the value vector.
        //
    {
        int i;

        double[] a = r8vec_zeros_new(nz);

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);
        //
        //  Diagonal elements of A.
        //
        for (i = 0; i < n; i++)
        {
            a[i] = fac * (i + 1) + i + 1;
        }

        for (i = 0; i < n; i++)
        {
            int k;
            for (k = ija[i]; k < ija[i + 1]; k++)
            {
                int j = ija[k];
                a[k] = fac * (i + 1) + j + 1;
            }
        }

        return a;
    }

    public static double[] r8ri_mtv(int n, int nz, int[] ija, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RI_MTV multiplies the transpose of an R8RI matrix times a vector.
        //
        //  Discussion:
        //
        //    An R8RI matrix is in row indexed sparse storage form, using an index
        //    array IJA and a value array A.  The first N entries of A store the
        //    diagonal elements in order.  The first N entries of IJA store the index
        //    of the first off-diagonal element of the corresponding row; if there is
        //    no off-diagonal element in that row, it is one greater than the index
        //    in A of the most recently stored element in the previous row.
        //    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
        //    greater than the index in A of the last off-diagonal element of the
        //    last row.  Location N+1 of A is not used.  Entries in A with index
        //    N+2 or greater contain the off-diagonal values, ordered by row, and
        //    then by column.  Entries in IJA with index N+2 or greater contain the
        //    column number of the corresponding element in A.
        //
        //  Example:
        //
        //    A:
        //      3 0 1 0 0
        //      0 4 0 0 0
        //      0 7 5 9 0
        //      0 0 0 0 2
        //      0 0 0 6 8
        //
        //    NZ = 11
        //
        //    IJA:
        //      7  8  8 10 11 12  3  2  4  5  4
        //
        //    A:
        //      3  4  5  0  8  *  1  7  9  2  6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
        //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
        //    Third Edition,
        //    Cambridge University Press, 2007,
        //    ISBN13: 978-0-521-88068-8,
        //    LC: QA297.N866.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the size required for the RI
        //    or "row indexed" sparse storage.
        //
        //    Input, int IJA[NZ], the index vector.
        //
        //    Input, double A[NZ], the value vector.
        //
        //    Input, double X[N], the vector to be multiplied.
        //
        //    Output, double R8RI_MTV[N], the product A'*X.
        //
    {
        int i;

        if (ija[0] != n + 1)
        {
            Console.WriteLine("");
            Console.WriteLine("R8RI_MTV - Fatal error!");
            Console.WriteLine("  The values IJA[0] and N are inconsistent.");
            return null;
        }

        double[] b = r8vec_zeros_new(n);

        for (i = 0; i < n; i++)
        {
            b[i] = a[i] * x[i];
        }

        for (i = 0; i < n; i++)
        {
            int k;
            for (k = ija[i]; k < ija[i + 1]; k++)
            {
                int j = ija[k];
                b[j] += a[k] * x[i];
            }
        }

        return b;
    }

    public static double[] r8ri_mv(int n, int nz, int[] ija, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RI_MV multiplies an R8RI matrix times a vector.
        //
        //  Discussion:
        //
        //    An R8RI matrix is in row indexed sparse storage form, using an index
        //    array IJA and a value array A.  The first N entries of A store the
        //    diagonal elements in order.  The first N entries of IJA store the index
        //    of the first off-diagonal element of the corresponding row; if there is
        //    no off-diagonal element in that row, it is one greater than the index
        //    in A of the most recently stored element in the previous row.
        //    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
        //    greater than the index in A of the last off-diagonal element of the
        //    last row.  Location N+1 of A is not used.  Entries in A with index
        //    N+2 or greater contain the off-diagonal values, ordered by row, and
        //    then by column.  Entries in IJA with index N+2 or greater contain the
        //    column number of the corresponding element in A.
        //
        //  Example:
        //
        //    A:
        //      3 0 1 0 0
        //      0 4 0 0 0
        //      0 7 5 9 0
        //      0 0 0 0 2
        //      0 0 0 6 8
        //
        //    NZ = 11
        //
        //    IJA:
        //      7  8  8 10 11 12  3  2  4  5  4
        //
        //    A:
        //      3  4  5  0  8  *  1  7  9  2  6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
        //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
        //    Third Edition,
        //    Cambridge University Press, 2007,
        //    ISBN13: 978-0-521-88068-8,
        //    LC: QA297.N866.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the size required for the RI
        //    or "row indexed" sparse storage.
        //
        //    Input, int IJA[NZ], the index vector.
        //
        //    Input, double A[NZ], the value vector.
        //
        //    Input, double X[N], the vector to be multiplied.
        //
        //    Output, double R8RI_MTV[N], the product A*X.
        //
    {
        int i;

        if (ija[0] != n + 1)
        {
            Console.WriteLine("");
            Console.WriteLine("R8RI_MV - Fatal error!");
            Console.WriteLine("  The values IJA[0] and N are inconsistent.");
            return null;
        }

        double[] b = r8vec_zeros_new(n);

        for (i = 0; i < n; i++)
        {
            b[i] = a[i] * x[i];
            int k;
            for (k = ija[i]; k < ija[i + 1]; k++)
            {
                b[i] += a[k] * x[ija[k]];
            }
        }

        return b;
    }

    public static void r8ri_print(int n, int nz, int[] ija, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RI_PRINT prints an R8RI matrix.
        //
        //  Discussion:
        //
        //    An R8RI matrix is in row indexed sparse storage form, using an index
        //    array IJA and a value array A.  The first N entries of A store the
        //    diagonal elements in order.  The first N entries of IJA store the index
        //    of the first off-diagonal element of the corresponding row; if there is
        //    no off-diagonal element in that row, it is one greater than the index
        //    in A of the most recently stored element in the previous row.
        //    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
        //    greater than the index in A of the last off-diagonal element of the
        //    last row.  Location N+1 of A is not used.  Entries in A with index
        //    N+2 or greater contain the off-diagonal values, ordered by row, and
        //    then by column.  Entries in IJA with index N+2 or greater contain the
        //    column number of the corresponding element in A.
        //
        //  Example:
        //
        //    A:
        //      3 0 1 0 0
        //      0 4 0 0 0
        //      0 7 5 9 0
        //      0 0 0 0 2
        //      0 0 0 6 8
        //
        //    NZ = 11
        //
        //    IJA:
        //      7  8  8 10 11 12  3  2  4  5  4
        //
        //    A:
        //      3  4  5  0  8  *  1  7  9  2  6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the size required for the RI
        //    or "row indexed" sparse storage.
        //
        //    Input, int IJA[NZ], the index vector.
        //
        //    Input, double A[NZ], the value vector.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8ri_print_some(n, nz, ija, a, 0, 0, n - 1, n - 1, title);
    }

    public static void r8ri_print_some(int n, int nz, int[] ija, double[] a, int ilo, int jlo,
            int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RI_PRINT_SOME prints some of an R8RI matrix.
        //
        //  Discussion:
        //
        //    An R8RI matrix is in row indexed sparse storage form, using an index
        //    array IJA and a value array A.  The first N entries of A store the
        //    diagonal elements in order.  The first N entries of IJA store the index
        //    of the first off-diagonal element of the corresponding row; if there is
        //    no off-diagonal element in that row, it is one greater than the index
        //    in A of the most recently stored element in the previous row.
        //    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
        //    greater than the index in A of the last off-diagonal element of the
        //    last row.  Location N+1 of A is not used.  Entries in A with index
        //    N+2 or greater contain the off-diagonal values, ordered by row, and
        //    then by column.  Entries in IJA with index N+2 or greater contain the
        //    column number of the corresponding element in A.
        //
        //  Example:
        //
        //    A:
        //      3 0 1 0 0
        //      0 4 0 0 0
        //      0 7 5 9 0
        //      0 0 0 0 2
        //      0 0 0 6 8
        //
        //    NZ = 11
        //
        //    IJA:
        //      7  8  8 10 11 12  3  2  4  5  4
        //
        //    A:
        //      3  4  5  0  8  *  1  7  9  2  6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the size required for the RI
        //    or "row indexed" sparse storage.
        //
        //    Input, int IJA[NZ], the index vector.
        //
        //    Input, double A[NZ], the value vector.
        //
        //    Input, int ILO, JLO, IHI, JHI, the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int incx = 5;
        int j2lo;

        double[] arow = r8vec_zeros_new(n);

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += incx)
        {
            int j2hi = j2lo + incx - 1;
            j2hi = Math.Min(j2hi, n);
            j2hi = Math.Min(j2hi, jhi);

            Console.WriteLine("");
            string cout = "  Col: ";
            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                cout += j.ToString(CultureInfo.InvariantCulture).PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            int i2lo = Math.Max(ilo, 0);
            int i2hi = Math.Min(ihi, n - 1);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  1) Assume everything is zero.
                //
                for (j = j2lo; j <= j2hi; j++)
                {
                    arow[j] = 0.0;
                }

                //
                //  2) Diagonal entry?
                //
                if (j2lo <= i && i <= j2hi)
                {
                    arow[i] = a[i];
                }

                //
                //  3) Now examine all the offdiagonal entries.
                //
                int k;
                for (k = ija[i]; k < ija[i + 1]; k++)
                {
                    j = ija[k];
                    if (j2lo <= j && j <= j2hi)
                    {
                        arow[j] = a[k];
                    }
                }

                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";
                for (j = j2lo; j <= j2hi; j++)
                {
                    cout += arow[j].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }

    }

    public static double[] r8ri_random(int n, int nz, int[] ija, ref int seed)

        //*****************************************************************************/
        //
        //  Purpose:
        //
        //    R8RI_RANDOM randomizes an R8RI matrix for given sparsity.
        //
        //  Discussion:
        //
        //    An R8RI matrix is in row indexed sparse storage form, using an index
        //    array IJA and a value array A.  The first N entries of A store the
        //    diagonal elements in order.  The first N entries of IJA store the index
        //    of the first off-diagonal element of the corresponding row; if there is
        //    no off-diagonal element in that row, it is one greater than the index
        //    in A of the most recently stored element in the previous row.
        //    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
        //    greater than the index in A of the last off-diagonal element of the
        //    last row.  Location N+1 of A is not used.  Entries in A with index
        //    N+2 or greater contain the off-diagonal values, ordered by row, and
        //    then by column.  Entries in IJA with index N+2 or greater contain the
        //    column number of the corresponding element in A.
        //
        //  Example:
        //
        //    A:
        //      3 0 1 0 0
        //      0 4 0 0 0
        //      0 7 5 9 0
        //      0 0 0 0 2
        //      0 0 0 6 8
        //
        //    NZ = 11
        //
        //    IJA:
        //      7  8  8 10 11 12  3  2  4  5  4
        //
        //    A:
        //      3  4  5  0  8  *  1  7  9  2  6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
        //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
        //    Third Edition,
        //    Cambridge University Press, 2007,
        //    ISBN13: 978-0-521-88068-8,
        //    LC: QA297.N866.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the size required for the RI
        //    or "row indexed" sparse storage.  NZ = 3*N-1.
        //
        //    Input, int IJA[NZ], the index vector.
        //
        //    Input/output, int *SEED, a seed for the random number
        //    generator.
        //
        //    Output, double A[NZ], the value vector.
        //
    {
        int i;

        double[] a = r8vec_zeros_new(nz);
        //
        //  Diagonal elements of A.
        //
        for (i = 0; i < n; i++)
        {
            a[i] = UniformRNG.r8_uniform_01(ref seed);
        }

        for (i = 0; i < n; i++)
        {
            int k;
            for (k = ija[i]; k < ija[i + 1]; k++)
            {
                a[k] = UniformRNG.r8_uniform_01(ref seed);
            }
        }

        return a;
    }

    public static double[] r8ri_to_r8ge(int n, int nz, int[] ija, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RI_TO_R8GE converts an R8RI matrix to R8GE form.
        //
        //  Discussion:
        //
        //    An R8RI matrix is in row indexed sparse storage form.
        //
        //    A R8GE matrix is in general storage.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
        //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
        //    Third Edition,
        //    Cambridge University Press, 2007,
        //    ISBN13: 978-0-521-88068-8,
        //    LC: QA297.N866.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the size required for the RI
        //    or "row indexed" sparse storage.
        //
        //    Input, int IJA[NZ], the index vector.
        //
        //    Input, double A[NZ], the value vector.
        //
        //    Output, double R8RI_TO_R8GE[N*N], the matrix stored in GE 
        //    or "general" format.
        //
    {
        int i;
        int j;
        int k;

        double[] a_r8ge = r8vec_zeros_new(n * n);

        for (k = 0; k < n; k++)
        {
            i = k;
            j = k;
            a_r8ge[i + j * n] = a[k];
        }

        for (i = 0; i < n; i++)
        {
            for (k = ija[i]; k < ija[i + 1]; k++)
            {
                j = ija[k];
                a_r8ge[i + j * n] = a[k];
            }
        }

        return a_r8ge;
    }

    public static double[] r8ri_zeros(int n, int nz, int[] ija)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RI_ZEROS zeros an R8RI matrix.
        //
        //  Discussion:
        //
        //    An R8RI matrix is in row indexed sparse storage form.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the size required for the RI
        //    or "row indexed" sparse storage.
        //
        //    Input, int IJA[NZ], the index vector.
        //
        //    Output, double R8RI_ZEROS[NZ], the value vector.
        //
    {
        double[] a = r8vec_zeros_new(nz);

        return a;
    }
}