using System;
using System.Globalization;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[][] r8rmat_copy_new(int m, int n, double[][] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RMAT_COPY_NEW makes a new copy of an R8RMAT .
        //
        //  Discussion:
        //
        //    An R8RMAT is a matrix stored in row major form, using M pointers
        //    to the beginnings of rows.
        //
        //    A declaration of the form
        //      double **a;
        //    is necesary.  Then an assignment of the form:
        //      a = r8rmat_new ( m, n );
        //    allows the user to assign entries to the matrix using typical
        //    2D array notation:
        //      a[2][3] = 17.0;
        //      y = a[1][0];
        //    and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double **A, the array to copy.
        //
        //    Output, double **R8RMAT_COPY_NEW, the copied array.
        //
    {
        int i;

        double[][] b = r8rmat_new(m, n);

        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                b[i][j] = a[i][j];
            }
        }

        return b;
    }
        
    public static void r8rmat_delete ( int m, int n, ref double[][] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RMAT_DELETE frees memory associated with an R8RMAT.
        //
        //  Discussion:
        //
        //    This function releases the memory associated with an R8RMAT.
        // 
        //    An R8RMAT is a row-major array that was created by a 
        //    command like:
        //
        //      double **a;
        //      a = r8rmat_new ( m, n );
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

    public static double[] r8rmat_fs_new(int n, double[][] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RMAT_FS_NEW factors and solves an R8RMAT system with one right hand side.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, double **A, the coefficient matrix of the linear system.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Output, double R8RMAT_FS_NEW[N], the solution of the linear system.
        //
    {
        int i;
        int j;
        int k;

        double[][] a2 = r8rmat_copy_new(n, n, a);
        double[] x = r8vec_copy_new(n, b);

        for (k = 0; k < n; k++)
        {
            //
            //  Find the maximum element in column I.
            //
            int p = k;

            for (i = k + 1; i < n; i++)
            {
                if (Math.Abs(a2[p][k]) < Math.Abs(a2[i][k]))
                {
                    p = i;
                }
            }

            switch (a2[p][k])
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R8RMAT_FS_NEW - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + k + "");
                    return null;
            }

            //
            //  Switch rows K and P.
            //
            double t;
            if (k != p)
            {
                for (j = 0; j < n; j++)
                {
                    t = a2[k][j];
                    a2[k][j] = a2[p][j];
                    a2[p][j] = t;
                }

                t = x[k];
                x[k] = x[p];
                x[p] = t;
            }

            //
            //  Scale the pivot row.
            //
            t = a2[k][k];
            a2[k][k] = 1.0;
            for (j = k + 1; j < n; j++)
            {
                a2[k][j] /= t;
            }

            x[k] /= t;
            //
            //  Use the pivot row to eliminate lower entries in that column.
            //
            for (i = k + 1; i < n; i++)
            {
                if (a2[i][k] == 0.0)
                {
                    continue;
                }

                t = -a2[i][k];
                a2[i][k] = 0.0;
                for (j = k + 1; j < n; j++)
                {
                    a2[i][j] += t * a2[k][j];
                }

                x[i] += t * x[k];
            }
        }

        //
        //  Back solve.
        //
        for (j = n - 1; 1 <= j; j--)
        {
            for (i = 0; i < j; i++)
            {
                x[i] -= a2[i][j] * x[j];
            }
        }
            
        return x;
    }

    public static double[][] r8rmat_new(int m, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RMAT_NEW allocates a new R8RMAT.
        //
        //  Discussion:
        //
        //    An R8RMAT is a row-major array that was created by a 
        //    command like:
        //
        //      double **a;
        //      a = r8rmat_new ( m, n );
        //
        //    The user assigns entries to the matrix using typical
        //    2D array notation:
        //      a[2][3] = 17.0;
        //      y = a[1][0];
        //    and so on.
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
        //    Output, double **R8RMAT_NEW, a new matrix.
        //
    {
        double[][] a = new double[m][];

        for (int i = 0; i < m; i++)
        {
            a[i] = new double[n];
        }

        return a;
    }


    public static void r8rmat_print(int m, int n, double[][] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RMAT_PRINT prints an R8RMAT.
        //
        //  Discussion:
        //
        //    An R8RMAT is a row-major array that was created by a 
        //    command like:
        //
        //      double **a;
        //      a = r8rmat_new ( m, n );
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
        r8rmat_print_some(m, n, a, 1, 1, m, n, title);
    }

    public static void r8rmat_print_some(int m, int n, double[][] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RMAT_PRINT_SOME prints some of an R8RMAT.
        //
        //  Discussion:
        //
        //    An R8RMAT is a row-major array that was created by a 
        //    command like:
        //
        //      double **a;
        //      a = r8rmat_new ( m, n );
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

            int i2hi = ihi < m ? ihi : m;

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = (i - 1).ToString(CultureInfo.InvariantCulture).PadLeft(5) + ": ";
                for (j = j2lo; j <= j2hi; j++)
                {
                    cout += a[i - 1][j - 1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8rmat_to_r8mat(int m, int n, double[][] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RMAT_TO_R8MAT copies data from an R8RMAT to an R8MAT.
        //
        //  Discussion:
        //
        //    An R8RMAT is a row-major array that was created by a 
        //    command like:
        //
        //    double **a;
        //    a = r8rmat_new ( m, n );
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
        //    Input, double **A = double A[M][N], the data, stored as an R8RMAT.
        //
        //    Output, double R8RMAT_TO_R8MAT[M*N], the data, stored as an R8MAT.
        //
    {
        int j;

        double[] b = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                b[i + j * m] = a[i][j];
            }
        }

        return b;
    }

    public static double[][] r8rmat_zeros(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8RMAT_ZEROS allocates and zeroes a new R8RMAT.
        //
        //  Discussion:
        //
        //    An R8RMAT is a row-major array that was created by a 
        //    command like:
        //
        //      double **a;
        //      a = r8rmat_new ( m, n );
        //
        //    The user assigns entries to the matrix using typical
        //    2D array notation:
        //      a[2][3] = 17.0;
        //      y = a[1][0];
        //    and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix.
        //
        //    Output, double **R8RMAT_ZEROS, a new matrix.
        //
    {
        int i;

        double[][] a = new double[m][];

        for (i = 0; i < m; i++)
        {
            a[i] = new double[n];
        }

        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                a[i][j] = 0.0;
            }
        }

        return a;
    }
        
}