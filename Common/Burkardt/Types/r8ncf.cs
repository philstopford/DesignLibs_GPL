using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8ncf_dif2(int m, int n, int nz_num, int[] rowcol)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8NCF_DIF2 sets up an R8NCF second difference matrix.
        //
        //  Discussion:
        //
        //    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
        //    a real array containing the nonzero values, a 2 by NZ_NUM integer 
        //    array storing the row and column of each nonzero entry.
        //
        //    The R8NCF format is used by NSPCG.  NSPCG requires that the information
        //    for the diagonal entries of the matrix must come first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in
        //    the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int ROWCOL[2,NZ_NUM], the coordinates of 
        //    the nonzero entries.
        //
        //    Output, double A[NZ_NUM], the matrix.
        //
    {
        int k;

        double[] a = r8vec_zeros_new(nz_num);

        for (k = 0; k < nz_num; k++)
        {
            int i = rowcol[0 + k * 2];
            int j = rowcol[1 + k * 2];

            if (j == i - 1)
            {
                a[k] = -1.0;
            }
            else if (j == i)
            {
                a[k] = 2.0;
            }
            else if (j == i + 1)
            {
                a[k] = -1.0;
            }
        }

        return a;
    }

    public static int r8ncf_dif2_nz_num(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8NCF_DIF2_NZ_NUM counts nonzeros in an R8NCF second difference matrix.
        //
        //  Discussion:
        //
        //    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
        //    a real array containing the nonzero values, a 2 by NZ_NUM integer 
        //    array storing the row and column of each nonzero entry.
        //
        //    The R8NCF format is used by NSPCG.  NSPCG requires that the information
        //    for the diagonal entries of the matrix must come first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in
        //    the matrix.
        //
        //    Output, int NZ_NUM, the number of nonzero entries.
        //
    {
        int nz_num;

        if (m < n)
        {
            nz_num = 3 * m - 1;
        }
        else if (m == n)
        {
            nz_num = 3 * n - 2;
        }
        else
        {
            nz_num = 3 * n - 1;
        }

        return nz_num;
    }

    public static int[] r8ncf_dif2_rowcol(int m, int n, int nz_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8NCF_DIF2_ROWCOL sets indexing for an R8NCF second difference matrix.
        //
        //  Discussion:
        //
        //    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
        //    a real array containing the nonzero values, a 2 by NZ_NUM integer 
        //    array storing the row and column of each nonzero entry.
        //
        //    The R8NCF format is used by NSPCG.  NSPCG requires that the information
        //    for the diagonal entries of the matrix must come first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in
        //    the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Output, int ROWCOL[2*NZ_NUM], the coordinates of 
        //    the nonzero entries.
        //
    {
        int i;

        int[] rowcol = i4vec_zeros_new(2 * nz_num);

        int k = 0;

        for (i = 0; i < m; i++)
        {
            int j = i - 1;
            switch (j)
            {
                case >= 0 when j < n:
                    rowcol[0 + k * 2] = i;
                    rowcol[1 + k * 2] = j;
                    k += 1;
                    break;
            }

            j = i;
            if (j < n)
            {
                rowcol[0 + k * 2] = i;
                rowcol[1 + k * 2] = j;
                k += 1;
            }

            j = i + 1;
            if (j < n)
            {
                rowcol[0 + k * 2] = i;
                rowcol[1 + k * 2] = j;
                k += 1;
            }

        }

        return rowcol;
    }

    public static double[] r8ncf_indicator(int m, int n, int nz_num, int[] rowcol)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8NCF_INDICATOR sets up an R8NCF indicator matrix.
        //
        //  Discussion:
        //
        //    The R8NCF storage format stores NZ_NUM, the number of nonzeros,
        //    a real array containing the nonzero values, a 2 by NZ_NUM integer
        //    array storing the row and column of each nonzero entry.
        //
        //    The R8NCF format is used by NSPCG.  NSPCG requires that the information
        //    for the diagonal entries of the matrix must come first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int ROWCOL[2*NZ_NUM], the coordinates of the nonzero entries.
        //
        //    Output, double A[NZ_NUM], the indicator matrix.
        //
    {
        int k;

        double[] a = r8vec_zeros_new(nz_num);

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (k = 0; k < nz_num; k++)
        {
            int i = rowcol[0 + k * 2];
            int j = rowcol[1 + k * 2];
            a[k] = fac * (i + 1) + j + 1;
        }

        return a;
    }

    public static double[] r8ncf_mtv(int m, int n, int nz_num, int[] rowcol, double[] a,
            double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8NCF_MTV multiplies an R8VEC times an R8NCF matrix.
        //
        //  Discussion:
        //
        //    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
        //    a real array containing the nonzero values, a 2 by NZ_NUM integer 
        //    array storing the row and column of each nonzero entry.
        //
        //    The R8NCF format is used by NSPCG.  NSPCG requires that the information
        //    for the diagonal entries of the matrix must come first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2016
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
        //    Input, int ROWCOL[2*NZ_NUM], the row and column
        //    indices of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements of the matrix.
        //
        //    Input, double X[M], the vector to be multiplied by A'.
        //
        //    Output, double R8NCF_MTV[N], the product A' * x.
        //
    {
        int k;

        double[] b = r8vec_zeros_new(n);

        for (k = 0; k < nz_num; k++)
        {
            int i = rowcol[0 + k * 2];
            int j = rowcol[1 + k * 2];
            b[j] += a[k] * x[i];
        }

        return b;
    }

    public static double[] r8ncf_mv(int m, int n, int nz_num, int[] rowcol, double[] a,
            double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8NCF_MV multiplies an R8NCF matrix by an R8VEC.
        //
        //  Discussion:
        //
        //    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
        //    a real array containing the nonzero values, a 2 by NZ_NUM integer 
        //    array storing the row and column of each nonzero entry.
        //
        //    The R8NCF format is used by NSPCG.  NSPCG requires that the information
        //    for the diagonal entries of the matrix must come first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2016
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
        //    Input, int ROWCOL[2*NZ_NUM], the row and column 
        //    indices of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements of the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8NCF_MV[M], the product A * x.
        //
    {
        int k;

        double[] b = r8vec_zeros_new(m);

        for (k = 0; k < nz_num; k++)
        {
            int i = rowcol[0 + k * 2];
            int j = rowcol[1 + k * 2];
            b[i] += a[k] * x[j];
        }

        return b;
    }

    public static void r8ncf_print(int m, int n, int nz_num, int[] rowcol,
            double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8NCF_PRINT prints an R8NCF matrix.
        //
        //  Discussion:
        //
        //    The R8NCF storage format stores NZ_NUM, the number of nonzeros,
        //    a real array containing the nonzero values, a 2 by NZ_NUM integer
        //    array storing the row and column of each nonzero entry.
        //
        //    The R8NCF format is used by NSPCG.  NSPCG requires that the information
        //    for the diagonal entries of the matrix must come first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in the matrix.
        //
        //    Input, int ROWCOL[2*NZ_NUM], the row and column indices
        //    of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements of the matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8ncf_print_some(m, n, nz_num, rowcol, a, 0, 0, m - 1, n - 1, title);
    }

    public static void r8ncf_print_some(int m, int n, int nz_num, int[] rowcol,
            double[] a, int ilo, int jlo, int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8NCF_PRINT_SOME prints some of an R8NCF matrix.
        //
        //  Discussion:
        //
        //    The R8NCF storage format stores NZ_NUM, the number of nonzeros,
        //    a real array containing the nonzero values, a 2 by NZ_NUM integer
        //    array storing the row and column of each nonzero entry.
        //
        //    The R8NCF format is used by NSPCG.  NSPCG requires that the information
        //    for the diagonal entries of the matrix must come first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in the matrix.
        //
        //    Input, int ROWCOL[2*NZ_NUM], the row and column indices
        //    of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements of the matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

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
                cout = i.ToString().PadLeft(5);
                for (j = j2lo; j <= j2hi; j++)
                {
                    double aij = 0.0;
                    int k;
                    for (k = 0; k < nz_num; k++)
                    {
                        if (rowcol[0 + k * 2] == i && rowcol[1 + k * 2] == j)
                        {
                            aij = a[k];
                            break;
                        }
                    }

                    cout += "  " + aij.ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8ncf_random(int m, int n, int nz_num, int[] rowcol, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8NCF_RANDOM randomizes an R8NCF matrix.
        //
        //  Discussion:
        //
        //    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
        //    a real array containing the nonzero values, a 2 by NZ_NUM integer 
        //    array storing the row and column of each nonzero entry.
        //
        //    The R8NCF format is used by NSPCG.  NSPCG requires that the information
        //    for the diagonal entries of the matrix must come first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in
        //    the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int ROWCOL[2*NZ_NUM], the coordinates of 
        //    the nonzero entries.
        //
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double R8NCF_RANDOM[NZ_NUM], the indicator matrix.
        //
    {
        int k;

        double[] a = r8vec_zeros_new(nz_num);

        for (k = 0; k < nz_num; k++)
        {
            a[k] = UniformRNG.r8_uniform_01(ref seed);
        }

        return a;
    }

    public static double[] r8ncf_to_r8ge(int m, int n, int nz_num, int[] rowcol, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8NCF_TO_R8GE converts an R8NCF matrix to R8GE format.
        //
        //  Discussion:
        //
        //    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
        //    a real array containing the nonzero values, a 2 by NZ_NUM integer 
        //    array storing the row and column of each nonzero entry.
        //
        //    The R8NCF format is used by NSPCG.  NSPCG requires that the information
        //    for the diagonal entries of the matrix must come first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in
        //    the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int ROWCOL[2*NZ_NUM], the coordinates of 
        //    the nonzero entries.
        //
        //    Input, double A[NZ_NUM], the matrix.
        //
        //    Output, double R8NCF_TO_R8GE[M*N], the R8GE matrix.
        //
    {
        int k;

        double[] a_r8ge = r8vec_zeros_new(m * n);

        for (k = 0; k < nz_num; k++)
        {
            int i = rowcol[0 + k * 2];
            int j = rowcol[1 + k * 2];
            a_r8ge[i + j * m] += a[k];
        }

        return a_r8ge;
    }

    public static double[] r8ncf_zeros(int m, int n, int nz_num, int[] rowcol)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8NCF_ZEROS zeros an R8NCF matrix.
        //
        //  Discussion:
        //
        //    The R8NCF storage format stores NZ_NUM, the number of nonzeros,
        //    a real array containing the nonzero values, a 2 by NZ_NUM integer
        //    array storing the row and column of each nonzero entry.
        //
        //    The R8NCF format is used by NSPCG.  NSPCG requires that the information
        //    for the diagonal entries of the matrix must come first.
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
        //    Input, int M, N, the number of rows and columns in the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int ROWCOL[2*NZ_NUM], the coordinates of the nonzero entries.
        //
        //    Output, double R8NCF_ZEROS[NZ_NUM], the matrix.
        //
    {
        double[] a = r8vec_zeros_new(nz_num);

        return a;
    }

}