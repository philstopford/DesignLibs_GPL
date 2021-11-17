using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8s3_diagonal(int m, int n, int nz_num, int sym, int[] row, int[] col,
            ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_DIAGONAL reorders an R8S3 matrix so diagonal entries are first.
        //
        //  Discussion:
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //    There is a symmetry option for square matrices.  If the symmetric storage
        //    option is used, the format specifies that only nonzeroes on the diagonal
        //    and lower triangle are stored.  However, this routine makes no attempt
        //    to enforce this.  The only thing it does is to "reflect" any nonzero
        //    offdiagonal value.  Moreover, no check is made for the erroneous case
        //    in which both A(I,J) and A(J,I) are specified, but with different values.
        //
        //    This routine reorders the entries of A so that the first N entries
        //    are exactly the diagonal entries of the matrix, in order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in 
        //    the matrix.
        //
        //    Input, int SYM, is 0 if the matrix is not symmetric, 
        //    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
        //    only the nonzeroes on the diagonal and in the lower triangle are stored.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and 
        //    column indices of the nonzero elements.
        //
        //    Input/output, double A[NZ_NUM], the nonzero elements 
        //    of the matrix.
        //
    {
        int k;

        int found = 0;

        for (k = 0; k < nz_num; k++)
        {
            while (row[k] == col[k])
            {
                if (row[k] == k)
                {
                    found += 1;
                    break;
                }

                int i = row[k];

                int j = row[i];
                row[i] = row[k];
                row[k] = j;

                j = col[i];
                col[i] = col[k];
                col[k] = j;

                (a[i], a[k]) = (a[k], a[i]);

                found += 1;

                if (Math.Min(m, n) <= found)
                {
                    break;
                }
            }

            if (Math.Min(m, n) <= found)
            {
                break;
            }
        }

        if (found >= Math.Min(m, n))
        {
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("R8S3_DIAGONAL - Warning!");
        Console.WriteLine("  Number of diagonal entries expected: " + Math.Min(m, n) + "");
        Console.WriteLine("  Number found was " + found + "");

    }

    public static void r8s3_dif2(int m, int n, int nz_num, int sym, int[] row, int[] col,
            double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_DIF2 sets up an R8S3 second difference matrix.
        //
        //  Discussion:
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int SYM, is 0 if the matrix is not symmetric, 
        //    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
        //    only the nonzeroes on the diagonal and in the lower triangle are stored.
        //
        //    Output, int ROW[NZ_NUM], COL[NZ_NUM], the row and column
        //    indices of the nonzero elements.
        //
        //    Output, double A[NZ_NUM], the indicator matrix.
        //
    {
        int i;
        int j;

        int k = 0;
        //
        //  Diagonal entries.
        //
        for (j = 0; j < n; j++)
        {
            i = j;
            row[k] = i;
            col[k] = j;
            a[k] = 2.0;
            k += 1;
        }

        //
        //  Offdiagonal nonzeros, by column.
        //
        for (j = 0; j < n; j++)
        {
            if (sym != 1)
            {
                switch (j)
                {
                    case > 0:
                        i = j - 1;
                        row[k] = i;
                        col[k] = j;
                        a[k] = -1.0;
                        k += 1;
                        break;
                }
            }

            if (j + 1 > m - 1)
            {
                continue;
            }

            i = j + 1;
            row[k] = i;
            col[k] = j;
            a[k] = -1.0;
            k += 1;
        }
    }

    public static double[] r8s3_indicator(int m, int n, int nz_num, int sym, int[] row,
            int[] col)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_INDICATOR sets up an R8S3 indicator matrix.
        //
        //  Discussion:
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int SYM, is 0 if the matrix is not symmetric, and 1
        //    if the matrix is symmetric.  If the matrix is symmetric, then
        //    only the nonzeroes on the diagonal and in the lower triangle are stored.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
        //    of the nonzero elements.
        //
        //    Output, double R8S3_INDICATOR[NZ_NUM], the indicator matrix.
        //
    {
        int k;

        double[] a = new double[nz_num];

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (k = 0; k < nz_num; k++)
        {
            int i = row[k];
            int j = col[k];
            a[k] = fac * (i + 1) + j + 1;
        }

        return a;
    }

    public static void r8s3_jac_sl(int n, int nz_num, int sym, int[] row, int[] col,
            double[] a, double[] b, ref double[] x, int it_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_JAC_SL solves an R8S3 system using Jacobi iteration.
        //
        //  Discussion:
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //    There is a symmetry option for square matrices.  If the symmetric storage
        //    option is used, the format specifies that only nonzeroes on the diagonal
        //    and lower triangle are stored.  However, this routine makes no attempt
        //    to enforce this.  The only thing it does is to "reflect" any nonzero
        //    offdiagonal value.  Moreover, no check is made for the erroneous case
        //    in which both A(I,J) and A(J,I) are specified, but with different values.
        //
        //    This routine REQUIRES that the matrix be square, that the matrix
        //    have nonzero diagonal entries, and that the first N entries of
        //    the array A be exactly the diagonal entries of the matrix, in order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in 
        //    the matrix.
        //
        //    Input, int SYM, is 0 if the matrix is not symmetric, 
        //    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
        //    only the nonzeroes on the diagonal and in the lower triangle are stored.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column 
        //    indices of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements of the matrix.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Input/output, double X[N], an approximate solution 
        //    to the system.
        //
        //    Input, int IT_MAX, the maximum number of iterations.
        //
    {
        int it_num;

        double[] x_new = new double[n];

        for (it_num = 1; it_num <= it_max; it_num++)
        {
            //
            //  Initialize to right hand side.
            //
            int j;
            for (j = 0; j < n; j++)
            {
                x_new[j] = b[j];
            }

            //
            //  Subtract off-diagonal terms.
            //
            int k;
            for (k = n; k < nz_num; k++)
            {
                int i = row[k];
                j = col[k];
                x_new[i] -= a[k] * x[j];
                switch (sym)
                {
                    case 1:
                        x_new[j] -= a[k] * x[i];
                        break;
                }
            }

            //
            //  Divide by diagonal terms and update.
            //
            for (j = 0; j < n; j++)
            {
                x[j] = x_new[j] / a[j];
            }
        }
    }

    public static double[] r8s3_mtv(int m, int n, int nz_num, int sym, int[] row, int[] col,
            double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_MTV multiplies an R8VEC times an R8S3 matrix.
        //
        //  Discussion:
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //    There is a symmetry option for square matrices.  If the symmetric storage
        //    option is used, the format specifies that only nonzeroes on the diagonal
        //    and lower triangle are stored.  However, this routine makes no attempt
        //    to enforce this.  The only thing it does is to "reflect" any nonzero
        //    offdiagonal value.  Moreover, no check is made for the erroneous case
        //    in which both A(I,J) and A(J,I) are specified, but with different values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in
        //    the matrix.
        //
        //    Input, int SYM, is 0 if the matrix is not symmetric, 
        //    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
        //    only the nonzeroes on the diagonal and in the lower triangle are stored.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column
        //    indices of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements of the matrix.
        //
        //    Input, double X[M], the vector to be multiplied by A'.
        //
        //    Output, double B[N], the product A' * x.
        //
    {
        int i;
        int j;
        int k;

        double[] b = r8vec_zeros_new(n);

        for (k = 0; k < nz_num; k++)
        {
            i = col[k];
            j = row[k];
            b[i] += a[k] * x[j];
        }

        switch (sym)
        {
            //
            //  Handle the symmetric option.
            //
            case 1 when m == n:
            {
                for (k = 0; k < nz_num; k++)
                {
                    i = row[k];
                    j = col[k];
                    if (i != j)
                    {
                        b[i] += a[k] * x[j];
                    }
                }

                break;
            }
        }

        return b;
    }

    public static double[] r8s3_mv(int m, int n, int nz_num, int sym, int[] row, int[] col,
            double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_MV multiplies an R8S3 matrix by an R8VEC.
        //
        //  Discussion:
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //    There is a symmetry option for square matrices.  If the symmetric storage
        //    option is used, the format specifies that only nonzeroes on the diagonal
        //    and lower triangle are stored.  However, this routine makes no attempt
        //    to enforce this.  The only thing it does is to "reflect" any nonzero
        //    offdiagonal value.  Moreover, no check is made for the erroneous case
        //    in which both A(I,J) and A(J,I) are specified, but with different values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in 
        //    the matrix.
        //
        //    Input, int SYM, is 0 if the matrix is not symmetric, 
        //    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
        //    only the nonzeroes on the diagonal and in the lower triangle are stored.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column 
        //    indices of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements of the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double B[M], the product A * x.
        //
    {
        int i;
        int j;
        int k;

        double[] b = r8vec_zeros_new(m);

        for (k = 0; k < nz_num; k++)
        {
            i = row[k];
            j = col[k];
            b[i] += a[k] * x[j];
        }

        switch (sym)
        {
            //
            //  Handle the symmetric option.
            //
            case 1 when m == n:
            {
                for (k = 0; k < nz_num; k++)
                {
                    i = col[k];
                    j = row[k];
                    if (i != j)
                    {
                        b[i] += a[k] * x[j];
                    }
                }

                break;
            }
        }

        return b;
    }

    public static void r8s3_print(int m, int n, int nz_num, int sym, int[] row, int[] col,
            double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_PRINT prints an R8S3 matrix.
        //
        //  Discussion:
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //    There is a symmetry option for square matrices.  If the symmetric storage
        //    option is used, the format specifies that only nonzeroes on the diagonal
        //    and lower triangle are stored.  However, this routine makes no attempt
        //    to enforce this.  The only thing it does is to "reflect" any nonzero
        //    offdiagonal value.  Moreover, no check is made for the erroneous case
        //    in which both A(I,J) and A(J,I) are specified, but with different values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in the matrix.
        //
        //    Input, int SYM, is 0 if the matrix is not symmetric, and 1
        //    if the matrix is symmetric.  If the matrix is symmetric, then
        //    only the nonzeroes on the diagonal and in the lower triangle are stored.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
        //    of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements 
        //    of the matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8s3_print_some(m, n, nz_num, sym, row, col, a, 0, 0, m - 1, n - 1, title);
    }

    public static void r8s3_print_some(int m, int n, int nz_num, int sym, int[] row, int[] col,
            double[] a, int ilo, int jlo, int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_PRINT_SOME prints some of an R8S3 matrix.
        //
        //  Discussion:
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //    There is a symmetry option for square matrices.  If the symmetric storage
        //    option is used, the format specifies that only nonzeroes on the diagonal
        //    and lower triangle are stored.  However, this routine makes no attempt
        //    to enforce this.  The only thing it does is to "reflect" any nonzero
        //    offdiagonal value.  Moreover, no check is made for the erroneous case
        //    in which both A(I,J) and A(J,I) are specified, but with different values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in the matrix.
        //
        //    Input, int SYM, is 0 if the matrix is not symmetric, and 1
        //    if the matrix is symmetric.  If the matrix is symmetric, then
        //    only the nonzeroes on the diagonal and in the lower triangle are stored.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
        //    of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements 
        //    of the matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

        int[] index = new int[INCX];
        int j2lo;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            int j2hi = j2lo + INCX - 1;
            j2hi = Math.Min(j2hi, n - 1);
            j2hi = Math.Min(j2hi, jhi);

            int inc = j2hi + 1 - j2lo;

            Console.WriteLine("");

            string cout = "  Col: ";
            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                cout += j.ToString().PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            int i2lo = Math.Max(ilo, 0);
            int i2hi = Math.Min(ihi, m - 1);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                bool nonzero = false;

                int j2;
                for (j2 = 0; j2 < inc; j2++)
                {
                    index[j2] = -1;
                }

                int k;
                for (k = 0; k < nz_num; k++)
                {
                    if (i == row[k] && j2lo <= col[k] && col[k] <= j2hi)
                    {
                        j2 = col[k] - j2lo + 1;

                        if (a[k] != 0.0)
                        {
                            index[j2 - 1] = k;
                            nonzero = true;
                        }
                    }
                    else
                    {
                        switch (sym)
                        {
                            case 1 when m == n && i == col[k] && j2lo <= row[k] && row[k] <= j2hi:
                            {
                                j2 = row[k] - j2lo + 1;

                                if (a[k] != 0.0)
                                {
                                    index[j2 - 1] = k;
                                    nonzero = true;
                                }

                                break;
                            }
                        }
                    }
                }

                switch (nonzero)
                {
                    case true:
                    {
                        cout = i.ToString().PadLeft(5) + " ";
                        for (j2 = 0; j2 < inc; j2++)
                        {
                            double aij = index[j2] switch
                            {
                                >= 0 => a[index[j2]],
                                _ => 0.0
                            };

                            cout += aij.ToString().PadLeft(14);
                        }

                        Console.WriteLine(cout);
                        break;
                    }
                }
            }
        }

        Console.WriteLine("");
    }

    public static double[] r8s3_random(int m, int n, int nz_num, int sym, int[] row, int[] col,
            ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_RANDOM randomizes an R8S3 matrix.
        //
        //  Discussion:
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int SYM, is 0 if the matrix is not symmetric, and 1
        //    if the matrix is symmetric.  If the matrix is symmetric, then
        //    only the nonzeroes on the diagonal and in the lower triangle are stored.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
        //    of the nonzero elements.
        //
        //    Input/output, integer &SEED, a seed for the random number generator.
        //
        //    Output, double R8S3_RANDOM[NZ_NUM], the matrix.
        //
    {
        int k;

        double[] a = new double[nz_num];

        for (k = 0; k < nz_num; k++)
        {
            a[k] = UniformRNG.r8_uniform_01(ref seed);
        }

        return a;
    }

    public static void r8s3_read(string input_file, int m, int n, int nz_num, ref int[] row,
            ref int[] col, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_READ reads a square R8S3 matrix from a file.
        //
        //  Discussion:
        //
        //    This routine needs the value of NZ_NUM, which can be determined
        //    by a call to R8S3_READ_SIZE.
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //    There is a symmetry option for square matrices.  If the symmetric storage
        //    option is used, the format specifies that only nonzeroes on the diagonal
        //    and lower triangle are stored.  However, this routine makes no attempt
        //    to enforce this.  The only thing it does is to "reflect" any nonzero
        //    offdiagonal value.  Moreover, no check is made for the erroneous case
        //    in which both A(I,J) and A(J,I) are specified, but with different values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string INPUT_FILE, the name of the file to be read.
        //
        //    Unused, int M, N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in the matrix.
        //
        //    Output, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
        //    of the nonzero elements.
        //
        //    Output, double A[NZ_NUM], the nonzero elements of the matrix.
        //
    {
        string[] input;
        int k = 0;

        try
        {
            input = File.ReadAllLines(input_file);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("R8S3_READ - Fatal error!");
            Console.WriteLine("  Could not open the input file: \"" + input_file + "\"");
            return;
        }

        foreach (string line in input)
        {
            string[] tokens = Helpers.splitStringByWhitespace(line);
            row[k] = Convert.ToInt32(tokens[0]);
            col[k] = Convert.ToInt32(tokens[1]);
            a[k] = Convert.ToInt32(tokens[2]);
            k++;
        }
    }

    public static void r8s3_read_size(string input_file, ref int m, ref int n, ref int nz_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_READ_SIZE reads the size of a square R8S3 matrix from a file.
        //
        //  Discussion:
        //
        //    The value of NZ_NUM is simply the number of records in the input file.
        //
        //    The value of N is determined as the maximum entry in the row and column
        //    vectors.
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //    There is a symmetry option for square matrices.  If the symmetric storage
        //    option is used, the format specifies that only nonzeroes on the diagonal
        //    and lower triangle are stored.  However, this routine makes no attempt
        //    to enforce this.  The only thing it does is to "reflect" any nonzero
        //    offdiagonal value.  Moreover, no check is made for the erroneous case
        //    in which both A(I,J) and A(J,I) are specified, but with different values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string INPUT_FILE, the name of the file to 
        //    be read.
        //
        //    Output, int &M, &N, the order of the matrix.
        //
        //    Output, int &NZ_NUM, the number of nonzero elements in the matrix.
        //
    {
        string[] input;

        m = 0;
        n = 0;
        nz_num = 0;

        try
        {
            input = File.ReadAllLines(input_file);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("R8S3_READ_SIZE - Fatal error!");
            Console.WriteLine("  Could not open the input file: \"" + input_file + "\"");
            return;
        }

        foreach (string line in input)
        {
            string[] tokens = Helpers.splitStringByWhitespace(line);
            int row_k = Convert.ToInt32(tokens[0]);
            int col_k = Convert.ToInt32(tokens[0]);

            nz_num += 1;
            m = Math.Max(m, row_k + 1);
            n = Math.Max(n, col_k + 1);
        }
    }

    public static double[] r8s3_res(int m, int n, int nz_num, int sym, int[] row, int[] col,
            double[] a, double[] x, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_RES computes the residual r=b-A*x for an R8S3 matrix.
        //
        //  Discussion:
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //    There is a symmetry option for square matrices.  If the symmetric storage
        //    option is used, the format specifies that only nonzeroes on the diagonal
        //    and lower triangle are stored.  However, this routine makes no attempt
        //    to enforce this.  The only thing it does is to "reflect" any nonzero
        //    offdiagonal value.  Moreover, no check is made for the erroneous case
        //    in which both A(I,J) and A(J,I) are specified, but with different values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in 
        //    the matrix.
        //
        //    Input, int SYM, is 0 if the matrix is not symmetric, 
        //    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
        //    only the nonzeroes on the diagonal and in the lower triangle are stored.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column 
        //    indices of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements of the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Input, double B[M], the right hand side.
        //
        //    Output, double R[M], the residual b-A*x.
        //
    {
        int i;

        double[] r = r8s3_mv(m, n, nz_num, sym, row, col, a, x);

        for (i = 0; i < m; i++)
        {
            r[i] = b[i] - r[i];
        }

        return r;
    }

    public static double[] r8s3_to_r8ge(int m, int n, int nz_num, int sym, int[] row,
            int[] col, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_TO_R8GE copies an R8S3 matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //    There is a symmetry option for square matrices.  If the symmetric storage
        //    option is used, the format specifies that only nonzeroes on the diagonal
        //    and lower triangle are stored.  However, this routine makes no attempt
        //    to enforce this.  The only thing it does is to "reflect" any nonzero
        //    offdiagonal value.  Moreover, no check is made for the erroneous case
        //    in which both A(I,J) and A(J,I) are specified, but with different values.
        //
        //    The R8GE storage format is used for a general M by N matrix.  A storage 
        //    space is made for each entry.  The two dimensional logical
        //    array can be thought of as a vector of M*N entries, starting with
        //    the M entries in the column 1, then the M entries in column 2
        //    and so on.  Considered as a vector, the entry A(I,J) is then stored
        //    in vector location I+(J-1)*M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in 
        //    the matrix.
        //
        //    Input, int SYM, is 0 if the matrix is not symmetric, 
        //    and 1 if the matrix is symmetric.  The symmetric case only makes sense
        //    if the matrix is also square, that is, M = N.  In this case, only
        //    the nonzeroes on the diagonal and in the lower triangle are stored.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column
        //    indices of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements of the matrix.
        //
        //    Output, double B[M*N], the R8GE matrix.
        //
    {
        int i;
        int j;
        int k;

        double[] b = new double[m * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                b[i + j * m] = 0.0;
            }
        }

        for (k = 0; k < nz_num; k++)
        {
            i = row[k];
            j = col[k];
            b[i + j * m] += a[k];
            switch (sym)
            {
                case 1 when m == n && i != j:
                    b[j + i * m] += a[k];
                    break;
            }
        }

        return b;
    }

    public static void r8s3_write(int m, int n, int nz_num, int sym, int[] row, int[] col,
            double[] a, string output_file)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_WRITE writes a square R8S3 matrix to a file.
        //
        //  Discussion:
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //    There is a symmetry option for square matrices.  If the symmetric storage
        //    option is used, the format specifies that only nonzeroes on the diagonal
        //    and lower triangle are stored.  However, this routine makes no attempt
        //    to enforce this.  The only thing it does is to "reflect" any nonzero
        //    offdiagonal value.  Moreover, no check is made for the erroneous case
        //    in which both A(I,J) and A(J,I) are specified, but with different values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in the matrix.
        //
        //    Input, int SYM, is 0 if the matrix is not symmetric, and 1
        //    if the matrix is symmetric.  If the matrix is symmetric, then
        //    only the nonzeroes on the diagonal and in the lower triangle are stored.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
        //    of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements 
        //    of the matrix.
        //
        //    Input, string OUTPUT_FILE, the name of the file to which
        //    the information is to be written.
        //
    {
        int k;
        List<string> output = new();

        for (k = 0; k < nz_num; k++)
        {
            output.Add("  " + row[k].ToString().PadLeft(8)
                            + "  " + col[k].ToString().PadLeft(8)
                            + "  " + a[k].ToString().PadLeft(16) + "");
        }

        try
        {
            File.WriteAllLines(output_file, output);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("R8S3_WRITE - Fatal error!");
            Console.WriteLine("  Could not open the output file.");
        }

    }

    public static double[] r8s3_zeros(int m, int n, int nz_num, int sym, int[] row, int[] col)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8S3_ZEROS zeros an R8S3 indicator matrix.
        //
        //  Discussion:
        //
        //    The R8S3 storage format corresponds to the SLAP Triad format.
        //
        //    The R8S3 storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
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
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int SYM, is 0 if the matrix is not symmetric, and 1
        //    if the matrix is symmetric.  If the matrix is symmetric, then
        //    only the nonzeroes on the diagonal and in the lower triangle are stored.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
        //    of the nonzero elements.
        //
        //    Output, double R8S3_ZEROS[NZ_NUM], the matrix.
        //
    {
        double[] a = r8vec_zeros_new(nz_num);

        return a;
    }
}