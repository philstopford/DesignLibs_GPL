using System;
using System.Globalization;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8pp_delete ( int m, int n, ref double[][] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PP_DELETE frees the memory set aside by R8PP_NEW.
        //
        //  Discussion:
        //
        //    An R8PP is a pointer to pointers to R8's, and is a sort of
        //    variably-dimensioned matrix.
        //
        //    This function releases the memory associated with an array that was 
        //    created by a command like:
        //
        //      double **a;
        //      a = r8pp_new ( m, n );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the array.
        //
        //    Input, double **A, the pointer to the pointers.
        //
    {
        a = null;
    }

    public static void r8pp_print(int n, double[] a, string title)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PP_PRINT prints a R8PP matrix.
        //
        //  Discussion:
        //
        //    The R8PP storage format is appropriate for a symmetric positive
        //    definite matrix.  Only the upper triangle of the matrix is stored,
        //    by successive partial columns, in an array of length (N*(N+1))/2,
        //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 April 2006
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
        //    Input, double A[(N*(N+1))/2], the R8PP matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8pp_print_some(n, a, 1, 1, n, n, title);
    }

    public static void r8pp_print_some(int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PP_PRINT_SOME prints some of a R8PP matrix.
        //
        //  Discussion:
        //
        //    The R8PP storage format is appropriate for a symmetric positive
        //    definite matrix.  Only the upper triangle of the matrix is stored,
        //    by successive partial columns, in an array of length (N*(N+1))/2,
        //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 April 2006
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
        //    Input, double A[(N*(N+1))/2], the R8PP matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (int j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            int j2hi = j2lo + INCX - 1;
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
            int i2lo = Math.Max(ilo, 1);
            int i2hi = Math.Min(ihi, n);

            for (int i = i2lo; i <= i2hi; i++)
            {
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  ";
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                for (j = j2lo; j <= j2hi; j++)
                {
                    double aij = i <= j ? a[i - 1 + j * (j - 1) / 2] : a[j - 1 + i * (i - 1) / 2];

                    cout += aij.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double r8pp_det(int n, double[] a_lu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PP_DET computes the determinant of a matrix factored by R8PP_FA.
        //
        //  Discussion:
        //
        //    The R8PP storage format is appropriate for a symmetric positive
        //    definite matrix.  Only the upper triangle of the matrix is stored,
        //    by successive partial columns, in an array of length (N*(N+1))/2,
        //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A_LU[(N*(N+1))/2], the LU factors from R8PP_FA.
        //
        //    Output, double R8PP_DET, the determinant of A.
        //
    {
        int i;

        double det = 1.0;

        int k = 0;
        for (i = 0; i < n; i++)
        {
            det *= a_lu[k];
            k = k + i + 2;
        }

        det *= det;

        return det;
    }

    public static double[] r8pp_dif2(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PP_DIF2 sets up an R8PP second difference matrix.
        //
        //  Discussion:
        //
        //    The R8PP storage format is appropriate for a symmetric positive
        //    definite matrix.  Only the upper triangle of the matrix is stored,
        //    by successive partial columns, in an array of length (N*(N+1))/2,
        //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double R8PP_DIF2((N*(N+1))/2), the R8PP matrix.
        //
    {
        int j;

        double[] a = new double[n * (n + 1) / 2];

        int k = 0;
        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < j - 1; i++)
            {
                a[k] = 0.0;
                k += 1;
            }

            switch (j)
            {
                case > 0:
                    a[k] = -1.0;
                    k += 1;
                    break;
            }

            a[k] = 2.0;
            k += 1;
        }

        return a;
    }

    public static double[] r8pp_fa(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PP_FA factors an R8PP matrix.
        //
        //  Discussion:
        //
        //    The R8PP storage format is appropriate for a symmetric positive
        //    definite matrix.  Only the upper triangle of the matrix is stored,
        //    by successive partial columns, in an array of length (N*(N+1))/2,
        //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 February 2004
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //    ISBN 0-89871-172-X
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[(N*(N+1))/2], the R8PP matrix.
        //
        //    Output, double R8PP_FA[(N*(N+1))/2], an upper triangular matrix R, stored 
        //    in packed form, so that A = R'*R.
        //
    {
        int i;
        int j;

        double[] b = new double[n * (n + 1) / 2];

        for (i = 0; i < n * (n + 1) / 2; i++)
        {
            b[i] = a[i];
        }

        int jj = 0;

        for (j = 1; j <= n; j++)
        {
            double s = 0.0;
            int kj = jj;
            int kk = 0;

            int k;
            for (k = 1; k <= j - 1; k++)
            {
                kj += 1;
                double t = b[kj - 1];
                for (i = 1; i <= k - 1; i++)
                {
                    t -= b[kk + i - 1] * b[jj + i - 1];
                }

                kk += k;
                t /= b[kk - 1];
                b[kj - 1] = t;
                s += t * t;
            }

            jj += j;
            s = b[jj - 1] - s;

            switch (s)
            {
                case <= 0.0:
                    return null;
                default:
                    b[jj - 1] = Math.Sqrt(s);
                    break;
            }
        }

        return b;
    }

    public static double[] r8pp_indicator(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PP_INDICATOR sets up an R8PP indicator matrix.
        //
        //  Discussion:
        //
        //    The R8PP storage format is appropriate for a symmetric positive
        //    definite matrix.  Only the upper triangle of the matrix is stored,
        //    by successive partial columns, in an array of length (N*(N+1))/2,
        //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 February 2013
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
        //    Output, double R8PP_INDICATOR((N*(N+1))/2), the R8PP matrix.
        //
    {
        int j;

        double[] a = new double[n * (n + 1) / 2];

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        int k = 0;
        for (j = 1; j <= n; j++)
        {
            int i;
            for (i = 1; i <= j; i++)
            {
                a[k] = fac * i + j;
                k += 1;
            }
        }

        return a;
    }

    public static double[] r8pp_mv(int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PP_MV multiplies an R8PP matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8PP storage format is appropriate for a symmetric positive
        //    definite matrix.  Only the upper triangle of the matrix is stored,
        //    by successive partial columns, in an array of length (N*(N+1))/2,
        //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[(N*(N+1))/2], the R8PP matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8PP_MV[N], the product A * x.
        //
    {
        int i;

        double[] b = r8vec_zeros_new(n);

        for (i = 0; i < n; i++)
        {
            int j;
            int k;
            for (j = 0; j < i; j++)
            {
                k = j + i * (i + 1) / 2;
                b[i] += a[k] * x[j];
            }

            for (j = i; j < n; j++)
            {
                k = i + j * (j + 1) / 2;
                b[i] += a[k] * x[j];
            }
        }

        return b;
    }

    public static double[][] r8pp_new(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PP_NEW allocates a new R8PP.
        //
        //  Discussion:
        //
        //    An R8PP is a pointer to pointers to R8's, and is a sort of
        //    variably-dimensioned matrix.
        //
        //    A declaration of the form
        //      double **a;
        //    is necesary.  Then an assignment of the form:
        //      a = r8pp_new ( m, n );
        //    allows the user to assign entries to the matrix using typical
        //    2D array notation:
        //      a[2][3] = 17;
        //      y = a[1][0];
        //    and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix.
        //
        //    Output, double **R8PP_NEW, a pointer to the pointers to the M by N array.
        //
    {
        int i;

        double[][] a = new double [m][];

        for (i = 0; i < m; i++)
        {
            a[i] = new double[n];
        }

        return a;
    }

    public static double[] r8pp_random(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PP_RANDOM randomizes an R8PP matrix.
        //
        //  Discussion:
        //
        //    The R8PP storage format is appropriate for a symmetric positive
        //    definite matrix.  Only the upper triangle of the matrix is stored,
        //    by successive partial columns, in an array of length (N*(N+1))/2,
        //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
        //
        //    The matrix is computed by setting a "random" upper triangular
        //    Cholesky factor R, and then computing A = R'*R.
        //    The randomness is limited by the fact that all the entries of
        //    R will be between 0 and 1.  A truly random R is only required
        //    to have positive entries on the diagonal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2003
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
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8PP_RANDOM[(N*(N+1))/2], the R8PP matrix.
        //
    {
        int i;

        double[] a = new double[n * (n + 1) / 2];

        for (i = 0; i < n * (n + 1) / 2; i++)
        {
            a[i] = 0.0;
        }

        for (i = n; 1 <= i; i--)
        {
            //
            //  Set row I of R.
            //
            int ij;
            int j;
            for (j = i; j <= n; j++)
            {
                ij = i + j * (j - 1) / 2;
                a[ij - 1] = UniformRNG.r8_uniform_01(ref seed);
            }

            //
            //  Consider element J of row I, last to first.
            //
            for (j = n; i <= j; j--)
            {
                //
                //  Add multiples of row I to lower elements of column J.
                //
                ij = i + j * (j - 1) / 2;

                int k;
                for (k = i + 1; k <= j; k++)
                {
                    int kj = k + j * (j - 1) / 2;
                    int ik = i + k * (k - 1) / 2;
                    a[kj - 1] += a[ik - 1] * a[ij - 1];
                }

                //
                //  Reset element J.
                //
                int ii = i + i * (i - 1) / 2;
                a[ij - 1] = a[ii - 1] * a[ij - 1];
            }
        }

        return a;
    }

    public static double[] r8pp_sl(int n, double[] a_lu, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PP_SL solves an R8PP system factored by R8PP_FA.
        //
        //  Discussion:
        //
        //    The R8PP storage format is appropriate for a symmetric positive
        //    definite matrix.  Only the upper triangle of the matrix is stored,
        //    by successive partial columns, in an array of length (N*(N+1))/2,
        //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 October 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //    ISBN 0-89871-172-X
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A_LU[(N*(N+1))/2], the LU factors from R8PP_FA.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R8PP_SL[N], the solution.
        //
    {
        int i;
        int k;
        double t;

        double[] x = new double[n];

        int kk = 0;

        for (k = 1; k <= n; k++)
        {
            t = 0.0;
            for (i = 0; i < k - 1; i++)
            {
                t += a_lu[kk + i] * x[i];
            }

            kk += k;
            x[k - 1] = (b[k - 1] - t) / a_lu[kk - 1];
        }

        for (k = n; 1 <= k; k--)
        {
            x[k - 1] /= a_lu[kk - 1];
            kk -= k;
            t = -x[k - 1];
            for (i = 0; i < k - 1; i++)
            {
                x[i] += t * a_lu[kk + i];
            }
        }

        return x;
    }

    public static double[] r8pp_to_r8ge(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PP_TO_R8GE copies an R8PP matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8PP storage format is appropriate for a symmetric positive
        //    definite matrix.  Only the upper triangle of the matrix is stored,
        //    by successive partial columns, in an array of length (N*(N+1))/2,
        //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[(N*(N+1))/2], the R8PP matrix.
        //
        //    Output, double R8PP_TO_R8GE[N*N], the R8GE matrix.
        //
    {
        int i;
        int j;

        double[] b = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                b[i + j * n] = 0.0;
            }
        }

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (i <= j)
                {
                    b[i + j * n] = a[i + j * (j + 1) / 2];
                }
                else
                {
                    b[i + j * n] = a[j + i * (i + 1) / 2];
                }
            }
        }

        return b;
    }

    public static double[] r8pp_zeros(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PP_ZEROS zeros an R8PP matrix.
        //
        //  Discussion:
        //
        //    The R8PP storage format is appropriate for a symmetric positive
        //    definite matrix.  Only the upper triangle of the matrix is stored,
        //    by successive partial columns, in an array of length (N*(N+1))/2,
        //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //    N must be positive.
        //
        //    Output, double R8PP_ZERO[(N*(N+1))/2], the R8PP matrix.
        //
    {
        int k;

        double[] a = new double[n * (n + 1) / 2];

        for (k = 0; k < n * (n + 1) / 2; k++)
        {
            a[k] = 0.0;
        }

        return a;
    }

}