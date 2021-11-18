using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8cc_dif2(int m, int n, int nz_num, ref int[] col, ref int[] row, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_DIF2 returns the second difference as an R8CC matrix.
        //
        //  Discussion:
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in A.
        //
        //    Output, int COL[N+1], points to the first element of each column.
        //
        //    Output, int ROW[NZ_NUM], contains the row indices of the elements.
        //
        //    Output, double A[NZ_NUM], the R8CC matrix.
        //
    {
        int i;
        int j;
        //
        //  Column pointers
        //
        col[0] = 0;
        col[1] = 2;
        for (j = 2; j < n; j++)
        {
            col[j] = col[j - 1] + 3;
        }

        col[n] = col[n - 1] + 2;
        //
        //  Row indices
        //
        int k = 0;
        row[k] = 0;
        k += 1;
        row[k] = 1;
        k += 1;
        for (j = 1; j < n - 1; j++)
        {
            for (i = j - 1; i <= j + 1; i++)
            {
                row[k] = i;
                k += 1;
            }
        }

        row[k] = m - 2;
        k += 1;
        row[k] = m - 1;
        k += 1;
        //
        //  Values
        //
        k = 0;

        j = 0;

        i = 0;
        a[k] = 2.0;
        k += 1;
        i = 1;
        a[k] = -1.0;
        k += 1;

        for (j = 1; j < n - 1; j++)
        {
            i = j - 1;
            a[k] = -1.0;
            k += 1;
            i = j;
            a[k] = 2.0;
            k += 1;
            i = j + 1;
            a[k] = -1.0;
            k += 1;
        }

        j = n - 1;
        i = m - 2;
        a[k] = -1.0;
        k += 1;
        i = m - 1;
        a[k] = 2.0;
        k += 1;
    }

    public static double r8cc_get(int m, int n, int nz_num, int[] col, int[] row,
            double[] a, int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_GET gets a value of an R8CC matrix.
        //
        //  Discussion:
        //
        //    It is legal to request entries of the matrix for which no storage
        //    was set aside.  In that case, a zero value will be returned.
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int COL[N+1], indicate where each column's data begins.
        //
        //    Input, int ROW[NZ_NUM], the row indices.
        //
        //    Input, double A[NZ_NUM], the nonzero entries.
        //
        //    Input, int I, J, the indices of the value to retrieve.
        //
        //    Output, double R8CC_GET, the value of A(I,J).
        //
    {
        //
        //  Seek sparse index K corresponding to full index (I,J).
        //
        int k = r8cc_ijk(m, n, nz_num, col, row, i, j);
        double aij = k switch
        {
            //
            //  If no K was found, then be merciful, and simply return 0.
            //
            -1 => 0.0,
            _ => a[k]
        };

        return aij;
    }

    public static int r8cc_ijk(int m, int n, int nz_num, int[] col, int[] row, int i,
            int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_IJK seeks K, the sparse index of (I,J), the full index of an R8CC matrix.
        //
        //  Discussion:
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int COL[N+1], indicate where each column's data begins.
        //
        //    Input, int ROW[NZ_NUM], the row indices.
        //
        //    Input, int I, J, the indices of the value to retrieve.
        //
        //    Output, int R8CC_IJK, the index of the sparse matrix in which entry
        //    (I,J) is stored, or -1 if no such entry exists.
        //
    {
        //
        //  Determine the part of ROW containing row indices of entries
        //  in column J.
        //
        int k1 = col[j];
        int k2 = col[j + 1] - 1;
        //
        //  Seek the location K for which ROW(K) = I.
        //  
        int k = i4vec_search_binary_a(k2 + 1 - k1, row, i, aIndex: +k1);

        if (k != -1)
        {
            k += k1;
        }

        return k;
    }

    public static void r8cc_inc(int m, int n, int nz_num, int[] col, int[] row, ref double[] a,
            int i, int j, double aij)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_INC increments a value of an R8CC matrix.
        //
        //  Discussion:
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int COL[N+1], indicate where each column's data begins.
        //
        //    Input, int ROW[NZ_NUM], the row indices.
        //
        //    Input/output, double A[NZ_NUM], the nonzero entries.
        //    On output, entry (I,J) has been incremented.
        //
        //    Input, int I, J, the indices of the value to retrieve.
        //
        //    Input, double AIJ, the value to be added to A(I,J).
        //
    {
        //
        //  Seek sparse index K corresponding to full index (I,J).
        //
        int k = r8cc_ijk(m, n, nz_num, col, row, i, j);
        switch (k)
        {
            //
            //  If no K was found, we fail.
            //
            case -1:
                Console.WriteLine("");
                Console.WriteLine("R8CC_INC - Fatal error!");
                Console.WriteLine("  R8CC_IJK could not find the entry.");
                Console.WriteLine("  Row I = " + i + "");
                Console.WriteLine("  Col J = " + j + "");
                return;
            default:
                a[k] += aij;
                break;
        }
    }

    public static double[] r8cc_indicator(int m, int n, int nz_num, int[] col, int[] row)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_INDICATOR sets up an R8CC indicator matrix.
        //
        //  Discussion:
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in A.
        //
        //    Input, int COL[N+1], points to the first element of each column.
        //
        //    Input, int ROW[NZ_NUM], contains the row indices of the elements.
        //
        //    Output, double R8CC_INDICATOR[NZ_NUM], the R8CC matrix.
        //
    {
        int j;

        double[] a = new double[nz_num];

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (j = 0; j < n; j++)
        {
            int k;
            for (k = col[j]; k <= col[j + 1] - 1; k++)
            {
                int i = row[k];
                a[k] = fac * (i + 1) + j + 1;
            }
        }

        return a;
    }

    public static void r8cc_kij(int m, int n, int nz_num, int[] col, int[] row, int k,
            ref int i, ref int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_KIJ seeks (I,J), the full index of K, the sparse index of an R8CC matrix.
        //
        //  Discussion:
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int COL[N+1], indicate where each column's data begins.
        //
        //    Input, int ROW[NZ_NUM], the row indices.
        //
        //    Input, int K, the sparse index of an entry of the matrix.
        //    1 <= K <= NZ_NUM.
        //
        //    Output, int &I, &J, the full indices corresponding to the sparse
        //    index K.
        //
    {
        int jj;

        i = -1;
        j = -1;

        if (k < 0 || nz_num <= k)
        {
            return;
        }

        //
        //  The row index is easy.
        //
        i = row[k];
        //
        //  Determine the column by bracketing in COl.
        //
        for (jj = 0; jj < n; jj++)
        {
            int k1 = col[jj];
            int k2 = col[jj + 1] - 1;
            if (k1 > k || k > k2)
            {
                continue;
            }

            j = jj;
            break;
        }

        switch (j)
        {
            case -1:
                return;
        }
    }

    public static double[] r8cc_mtv(int m, int n, int nz_num, int[] col, int[] row,
            double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_MTV multiplies a vector times an R8CC matrix.
        //
        //  Discussion:
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in A.
        //
        //    Input, int COL[N+1], points to the first element of each column.
        //
        //    Input, int ROW[NZ_NUM], contains the row indices of the elements.
        //
        //    Input, double A[NZ_NUM], the R8CC matrix.
        //
        //    Input, double X[M], the vector to be multiplied.
        //
        //    Output, double R8CC_MTV[N], the product A' * X.
        //
    {
        int j;

        double[] b = r8vec_zeros_new(n);

        for (j = 0; j < n; j++)
        {
            int k;
            for (k = col[j]; k <= col[j + 1] - 1; k++)
            {
                int i = row[k];
                b[j] += a[k] * x[i];
            }
        }

        return b;
    }

    public static double[] r8cc_mv(int m, int n, int nz_num, int[] col, int[] row,
            double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_MV multiplies an R8CC matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in A.
        //
        //    Input, int COL[N+1], points to the first element of each column.
        //
        //    Input, int ROW[NZ_NUM], contains the row indices of the elements.
        //
        //    Input, double A[NZ_NUM], the R8CC matrix.
        //
        //    Input, double X[N], the vector to be multiplied.
        //
        //    Output, double R8CC_MV[M], the product A * X.
        //
    {
        int j;

        double[] b = r8vec_zeros_new(m);

        for (j = 0; j < n; j++)
        {
            int k;
            for (k = col[j]; k <= col[j + 1] - 1; k++)
            {
                int i = row[k];
                b[i] += a[k] * x[j];
            }
        }

        return b;
    }

    public static void r8cc_print(int m, int n, int nz_num, int[] col, int[] row,
            double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_PRINT prints an R8CC matrix.
        //
        //  Discussion:
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in A.
        //
        //    Input, int COL[N+1], points to the first element of each column.
        //
        //    Input, int ROW[NZ_NUM], contains the row indices of the elements.
        //
        //    Input, double A[NZ_NUM], the R8CC matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8cc_print_some(m, n, nz_num, col, row, a, 0, 0, m - 1, n - 1, title);
    }

    public static void r8cc_print_some(int m, int n, int nz_num, int[] col, int[] row,
            double[] a, int ilo, int jlo, int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_PRINT_SOME prints some of an R8CC matrix.
        //
        //  Discussion:
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in A.
        //
        //    Input, int COL[N+1], points to the first element of each column.
        //
        //    Input, int ROW[NZ_NUM], contains the row indices of the elements.
        //
        //    Input, double A[NZ_NUM], the R8CC matrix.
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
            string cout = "  Col:  ";

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
            int i2hi = Math.Min(ihi, m - 1);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                //
                //  Now consider each column J in J2LO to J2HI,
                //  and look at every nonzero, and check if it occurs in row I.
                //
                for (j = j2lo; j <= j2hi; j++)
                {
                    double value = 0.0;
                    int k;
                    for (k = col[j]; k <= col[j + 1] - 1; k++)
                    {
                        if (row[k] == i)
                        {
                            value = a[k];
                        }
                    }

                    cout += value.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8cc_random(int m, int n, int nz_num, int[] col, int[] row,
            ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_RANDOM randomizes an R8CC matrix.
        //
        //  Discussion:
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in A.
        //
        //    Input, int COL[N+1], points to the first element of each column.
        //
        //    Input, int ROW[NZ_NUM], contains the row indices of the elements.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8CC_RANDOM[NZ_NUM], the R8CC matrix.
        //
    {
        double[] a = UniformRNG.r8vec_uniform_01_new(nz_num, ref seed);

        return a;
    }

    public static void r8cc_read(string col_file, string row_file, string a_file, int m,
            int n, int nz_num, ref int[] col, ref int[] row, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_READ reads an R8CC matrix from three files.
        //
        //  Discussion:
        //
        //    This routine needs the values of M, N, and NZ_NUM, which can be 
        //    determined by a call to R8CC_READ_SIZE.
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, string COL_FILE, ROW_FILE, A_FILE, the names of the 
        //    files containing the column pointers, row indices, and matrix entries.
        //
        //    Input, int M, N, the number of rows and columns in the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in the matrix.
        //
        //    Output, int COL[N+1], the column pointers.
        //
        //    Output, int ROW[NZ_NUM], the row indices.
        //
        //    Output, double A[NZ_NUM], the nonzero elements of the matrix.
        //
    {
        string[] input;
        int k;

        try
        {
            input = File.ReadAllLines(col_file);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("R8CC_READ - Fatal error!");
            Console.WriteLine("  Could not open the file \"" + col_file + "\".");
            return;
        }

        for (k = 0; k < n + 1; k++)
        {
            col[k] = Convert.ToInt32(input[k]);
        }

        //
        //  Read the row information.
        //
        try
        {
            input = File.ReadAllLines(row_file);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("R8CC_READ - Fatal error!");
            Console.WriteLine("  Could not open the file \"" + row_file + "\".");
            return;
        }

        for (k = 0; k < nz_num; k++)
        {
            row[k] = Convert.ToInt32(input[k]);
        }

        //
        //  Read the value information.
        //
        try
        {
            input = File.ReadAllLines(a_file);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("R8CC_READ - Fatal error!");
            Console.WriteLine("  Could not open the file \"" + a_file + "\".");
            return;
        }

        for (k = 0; k < nz_num; k++)
        {
            a[k] = Convert.ToInt32(input[k]);
        }
    }

    public static void r8cc_read_size(string col_file, string row_file, ref int m, ref int n,
            ref int nz_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_READ_SIZE reads the sizes of an R8CC sparse matrix from a file.
        //
        //  Discussion:
        //
        //    The value of M is "guessed" to be the largest value that occurs in
        //    the ROW file.  However, if a row index of 0 is encountered, then
        //    the value of M is incremented by 1.
        //
        //    The value of N is the number of records in the COL file minus 1.
        //
        //    The value of NZ_NUM is simply the number of records in the ROW file.
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, string COL_FILE, *ROW_FILE, the names of the 
        //    column and row files that describe the structure of the matrix.
        //
        //    Output, int &M, &N, the inferred number of rows and columns 
        //    in the sparse matrix.
        //
        //    Output, int &NZ_NUM, the number of nonzero entries in the
        //    sparse matrix.
        //
    {
        string[] input;
        string[] input2;
        //
        //  Check the COL file first.
        //
        try
        {
            input = File.ReadAllLines(col_file);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("R8CC_READ_SIZE - Fatal error!");
            Console.WriteLine("  Could not open the file \"" + col_file + "\".");
            return;
        }

        n = input.Length;


        //
        //  Check the ROW file.
        //
        //  For unfathomable reasons, if I use "INPUT" for this file,
        //  I can get a file open failure.  Rather than make right the
        //  world, I gave up and accessed "INPUT2".
        //
        m = 0;
        nz_num = 0;

        try
        {
            input2 = File.ReadAllLines(row_file);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("R8CC_READ_SIZE - Fatal error!");
            Console.WriteLine("  Could not open the file \"" + row_file + "\".");
            return;
        }

        foreach (string line in input2)
        {
            int row = Convert.ToInt32(line);

            nz_num += 1;
            m = Math.Max(m, row);
        }

        m += 1;
    }

    public static void r8cc_set(int m, int n, int nz_num, int[] col, int[] row, ref double[] a,
            int i, int j, double aij)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_SET sets a value of an R8CC matrix.
        //
        //  Discussion:
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int COL[N+1], indicate where each column's data begins.
        //
        //    Input, int ROW[NZ_NUM], the row indices.
        //
        //    Input/output, double A[NZ_NUM], the nonzero entries.
        //    On output, the entry of A corresponding to (I,J) has been reset.
        //
        //    Input, int I, J, the indices of the value to retrieve.
        //
        //    Input, double AIJ, the new value of A(I,J).
        //
    {
        //
        //  Seek sparse index K corresponding to full index (I,J).
        //
        int k = r8cc_ijk(m, n, nz_num, col, row, i, j);
        switch (k)
        {
            //
            //  If no K was found, we fail.
            //
            case -1:
                Console.WriteLine("");
                Console.WriteLine("R8CC_SET - Fatal error!");
                Console.WriteLine("  R8CC_IJK could not find the entry.");
                Console.WriteLine("  Row I = " + i + "");
                Console.WriteLine("  Col J = " + j + "");
                return;
            default:
                a[k] = aij;
                break;
        }
    }

    public static double[] r8cc_to_r8ge(int m, int n, int nz_num, int[] col, int[] row,
            double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_TO_R8GE converts an R8CC matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in A.
        //
        //    Input, int COL[N+1], points to the first element of each column.
        //
        //    Input, int ROW[NZ_NUM], contains the row indices of the elements.
        //
        //    Input, double A[NZ_NUM], the R8CC matrix.
        //
        //    Input, double R8CC_TO_R8GE[M*N], the R8GE matrix.
        //
    {
        int j;

        double[] b = r8vec_zeros_new(m * n);

        if (col[0] < 0 || nz_num < col[0])
        {
            Console.WriteLine("");
            Console.WriteLine("R8CC_TO_R8GE - Fatal error!");
            Console.WriteLine("  COL[" + 0 + "] = " + col[0] + "");
            return null;
        }

        for (j = 0; j < n; j++)
        {
            if (col[j + 1] < 0 || nz_num < col[j + 1] - 1)
            {
                Console.WriteLine("");
                Console.WriteLine("R8CC_TO_R8GE - Fatal error!");
                Console.WriteLine("  COL[" + j + "] = " + col[j + 1] + "");
                return null;
            }

            int k;
            for (k = col[j]; k <= col[j + 1] - 1; k++)
            {
                int i = row[k];
                if (i < 0 || m <= i)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8CC_TO_R8GE - Fatal error!");
                    Console.WriteLine("  ROW[" + k + "] = " + i + "");
                    return null;
                }

                b[i + j * m] = a[k];
            }
        }

        return b;
    }

    public static void r8cc_write(string col_file, string row_file, string a_file, int m, int n,
            int nz_num, int[] col, int[] row, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_WRITE writes an R8CC matrix to three files.
        //
        //  Discussion:
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, string COL_FILE, ROW_FILE, A_FILE, the names of the 
        //    files containing the column pointers, row entries, and matrix entries.
        //
        //    Input, int M, N, the number of rows and columns in the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in the matrix.
        //
        //    Input, int COL[N+1], the column pointers.
        //
        //    Input, int ROW[NZ_NUM], the row indices.
        //
        //    Input, double A[NZ_NUM], the nonzero elements of the matrix.
        //
    {
        List<string> output = new();
        int k;


        for (k = 0; k < n + 1; k++)
        {
            output.Add(col[k] + "");
        }

        try
        {
            File.WriteAllLines(col_file, output);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("R8CC_WRITE - Fatal error!");
            Console.WriteLine("  Could not open the file \"" + col_file + "\".");
            return;
        }

        output.Clear();

        for (k = 0; k < nz_num; k++)
        {
            output.Add(row[k] + "");
        }

        try
        {
            File.WriteAllLines(row_file, output);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("R8CC_WRITE - Fatal error!");
            Console.WriteLine("  Could not open the file \"" + row_file + "\".");
            return;
        }

        output.Clear();


        for (k = 0;
             k < nz_num;
             k++)
        {
            output.Add(a[k] + "");
        }

        try
        {
            File.WriteAllLines(row_file, output);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("R8CC_WRITE - Fatal error!");
            Console.WriteLine("  Could not open the file \"" + a_file + "\".");
        }
    }

    public static double[] r8cc_zeros(int m, int n, int nz_num, int[] col, int[] row)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CC_ZEROS zeros an R8CC matrix.
        //
        //  Discussion:
        //
        //    The R8CC format is the double precision sparse compressed column
        //    format.  Associated with this format, we have an M by N matrix
        //    with NZ_NUM nonzero entries.  We construct the column pointer
        //    vector COL of length N+1, such that entries of column J will be
        //    stored in positions COL(J) through COL(J+1)-1.  This indexing
        //    refers to both the ROW and A vectors, which store the row indices
        //    and the values of the nonzero entries.  The entries of the
        //    ROW vector corresponding to each column are assumed to be
        //    ascending sorted.
        //
        //    The R8CC format is equivalent to the MATLAB "sparse" format,
        //    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in A.
        //
        //    Input, int COL[N+1], points to the first element of each column.
        //
        //    Input, int ROW[NZ_NUM], contains the row indices of the elements.
        //
        //    Output, double R8CC_ZERO[NZ_NUM], the R8CC matrix.
        //
    {
        double[] a = r8vec_zeros_new(nz_num);

        return a;
    }
}