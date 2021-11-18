using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8cmat_delete ( int m, int n, ref double[][] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CMAT_DELETE frees memory associated with an R8CMAT.
        //
        //  Discussion:
        //
        //    This function releases the memory associated with an R8CMAT.
        //
        //    An R8CMAT is a column-major array, storing element (I,J)
        //    as A[J][I], and can be created by a command like:
        //      double **a;
        //      a = r8cmat_new ( m, n );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the array.
        //
        //    Input, double **A, the pointer to the array.
        //
    {
        a = null;
    }

    public static double[][] r8cmat_new(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CMAT_NEW allocates a new R8CMAT.
        //
        //  Discussion:
        //
        //    An R8CMAT is a column-major array, storing element (I,J)
        //    as A[J][I], and can be created by a command like:
        //      double **a;
        //      a = r8cmat_new ( m, n );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix.
        //
        //    Output, double **R8CMAT_NEW, a new matrix.
        //
    {
        int j;

        double[][] a = new double[n][];

        for (j = 0; j < n; j++)
        {
            a[j] = new double[m];
        }

        return a;
    }

    public static void r8cmat_print(int m, int n, double[][] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CMAT_PRINT prints an R8CMAT.
        //
        //  Discussion:
        //
        //    An R8CMAT is a column-major array, storing element (I,J)
        //    as A[J][I], and can be created by a command like:
        //      double **a;
        //      a = r8cmat_new ( m, n );
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
        //    Input, double **A = A[M][N], the M by N matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8cmat_print_some(m, n, a, 1, 1, m, n, title);
    }

    public static void r8cmat_print_some(int m, int n, double[][] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CMAT_PRINT_SOME prints some of an R8CMAT.
        //
        //  Discussion:
        //
        //    An R8CMAT is a column-major array, storing element (I,J)
        //    as A[J][I], and can be created by a command like:
        //      double **a;
        //      a = r8cmat_new ( m, n );
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
        //    Input, double **A = A[M][N], the matrix.
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
                cout += (j - 1).ToString().PadLeft(7) + "       ";
            }

            Console.WriteLine("");
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
            if (ihi < m)
            {
                i2hi = ihi;
            }
            else
            {
                i2hi = m;
            }

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = (i - 1).ToString().PadLeft(5) + ": ";
                for (j = j2lo; j <= j2hi; j++)
                {
                    cout += a[j - 1][i - 1].ToString().PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8cmat_to_r8mat_new(int m, int n, double[][] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CMAT_TO_R8MAT_NEW copies data from an R8CMAT to an R8MAT.
        //
        //  Discussion:
        //
        //    An R8CMAT is a column-major array, storing element (I,J)
        //    as A[J][I], and can be created by a command like:
        //      double **a;
        //      a = r8cmat_new ( m, n );
        //
        //    An R8MAT is a column-major array stored as a vector, so
        //    that element (I,J) of the M by N array is stored in location
        //    I+J*M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double **A = double A[M][N], the data, stored as an R8CMAT.
        //
        //    Output, double R8CMAT_TO_R8MAT_NEW[M*N], the data, stored as an R8MAT.
        //
    {
        int j;

        double[] b = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                b[i + j * m] = a[j][i];
            }
        }

        return b;
    }

    public static double[][] r8cmat_zeros_new(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CMAT_ZEROS_NEW allocates and zeros a new R8CMAT.
        //
        //  Discussion:
        //
        //    An R8CMAT is a column-major array, storing element (I,J)
        //    as A[J][I], and can be created by a command like:
        //      double **a;
        //      a = r8cmat_new ( m, n );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix.
        //
        //    Output, double **R8CMAT_ZEROS_NEW, a new matrix.
        //
    {
        int j;

        double[][] a = new double[n][];

        for (j = 0; j < n; j++)
        {
            a[j] = new double[m];
        }

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a[j][i] = 0.0;
            }
        }

        return a;
    }
}