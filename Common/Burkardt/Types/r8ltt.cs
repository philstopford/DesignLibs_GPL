using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8ltt_det(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LTT_DET computes the determinant of a R8LTT matrix.
        //
        //  Discussion:
        //
        //    The R8LTT storage format is used for an N by N lower triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the matrix.
        //
        //    Output, double R8LTT_DET, the determinant of the matrix.
        //
    {
        double det;

        det = Math.Pow(a[0], n);

        return det;
    }

    public static double[] r8ltt_indicator(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LTT_INDICATOR sets up a R8LTT indicator matrix.
        //
        //  Discussion:
        //
        //    The R8LTT storage format is used for an N by N lower triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double A[N], the matrix.
        //
    {
        double[] a;
        int j;

        a = new double[n];

        for (j = 0; j < n; j++)
        {
            a[j] = j + 1;
        }

        return a;
    }

    public static double[] r8ltt_inverse(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LTT_INVERSE computes the inverse of a R8LTT matrix.
        //
        //  Discussion:
        //
        //    The R8LTT storage format is used for an N by N lower triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the matrix to be inverted.
        //
        //    Output, double R8LTT_INVERSE[N], the inverse matrix.
        //
    {
        double[] b;
        double d;
        int i;
        int j;
        double[] p;
        double[] pn;
        double[] pnn;
        //
        //  Initialize B.
        //
        d = 1.0 / a[0];
        b = new double[n];
        b[0] = d;
        for (i = 1; i < n; i++)
        {
            b[i] = 0.0;
        }

        //
        //  Set the strict lower triangle.
        //
        p = new double[n];
        p[0] = 0.0;
        for (i = 1; i < n; i++)
        {
            p[i] = a[i];
        }

        //
        //  PN will hold powers of P.
        //
        pn = new double[n];
        pn[0] = 1.0;
        for (i = 1; i < n; i++)
        {
            pn[i] = 0.0;
        }

        //
        //  Add N-1 powers of strict lower triangle.
        //
        for (j = 1; j < n; j++)
        {
            d = -d / a[0];
            pnn = r8ltt_mm(n, p, pn);
            for (i = 0; i < n; i++)
            {
                b[i] += d * pnn[i];
                pn[i] = pnn[i];
            }
        }

        return b;
    }

    public static double[] r8ltt_mm(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LTT_MM computes C = A * B, where A and B are R8LTT matrices.
        //
        //  Discussion:
        //
        //    The R8LTT storage format is used for an N by N lower triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrices.
        //
        //    Input, double A[N], the first factor.
        //
        //    Input, double B[N], the second factor.
        //
        //    Output, double R8LTT_MM[N], the product.
        //
    {
        double[] c;
        double[] d;
        double[] e;
        int k;

        d = new double[n];

        for (k = 0; k < n; k++)
        {
            d[k] = b[n - 1 - k];
        }

        e = r8ltt_mtv(n, a, d);

        c = new double[n];

        for (k = 0; k < n; k++)
        {
            c[k] = e[n - 1 - k];
        }

        return c;
    }

    public static double[] r8ltt_mtm(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LTT_MTM computes C = A' * B, where A and B are R8LTT matrices.
        //
        //  Discussion:
        //
        //    The R8LTT storage format is used for an N by N lower triangular Toeplitz
        //    matrix.
        //
        //    Note that the result C is a dense matrix, of type R8GE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrices.
        //
        //    Input, double A[N], B[N], the factors.
        //
        //    Output, double R8LTT_MTM[N*N], the product.
        //
    {
        double[] c;
        int i;
        int j;
        int k;

        c = new double[n * n];

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                c[i + j * n] = 0.0;
                for (k = Math.Max(i, j); k < n; k++)
                {
                    c[i + j * n] += a[k - i] * b[k - j];
                }
            }
        }

        return c;
    }

    public static double[] r8ltt_mtv(int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LTT_MTV computes b = A'*x, where A is an R8LTT matrix.
        //
        //  Discussion:
        //
        //    The R8LTT storage format is used for an N by N lower triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8LTT_MTV[N], the product A' * x.
        //
    {
        double[] b;
        int d;
        int i;
        int j;

        b = new double[n];

        for (j = 0; j < n; j++)
        {
            b[j] = 0.0;
        }

        for (d = 0; d < n; d++)
        {
            for (i = d; i < n; i++)
            {
                j = i - d;
                b[j] += a[i - j] * x[i];
            }
        }

        return b;
    }

    public static double[] r8ltt_mv(int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LTT_MV computes b=A*x, where A is an R8LTT matrix.
        //
        //  Discussion:
        //
        //    The R8LTT storage format is used for an N by N lower triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8LTT_MV[N], the product A * x.
        //
    {
        double[] b;
        int d;
        int i;
        int j;

        b = new double[n];

        for (i = 0; i < n; i++)
        {
            b[i] = 0.0;
        }

        for (d = 0; d < n; d++)
        {
            for (i = d; i < n; i++)
            {
                j = i - d;
                b[i] += a[i - j] * x[j];
            }
        }

        return b;
    }

    public static void r8ltt_print(int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LTT_PRINT prints an R8LTT matrix.
        //
        //  Discussion:
        //
        //    The R8LTT storage format is used for an N by N lower triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8LTT matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8ltt_print_some(n, a, 0, 0, n - 1, n - 1, title);
    }

    public static void r8ltt_print_some(int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LTT_PRINT_SOME prints some of an R8LTT matrix.
        //
        //  Discussion:
        //
        //    The R8LTT storage format is used for an N by N lower triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8LTT matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //    0 <= ILO <= IHI < M.
        //    0 <= JLO <= JHI < N.
        //
        //    Input, string TITLE, a title.
        //
    {
        int INCX = 5;

        int i;
        int i2hi;
        int i2lo;
        int j;
        int j2hi;
        int j2lo;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            j2hi = j2lo + INCX - 1;
            j2hi = Math.Min(j2hi, n - 1);
            j2hi = Math.Min(j2hi, jhi);

            Console.WriteLine("");
            cout = "  Col: ";

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
            i2lo = Math.Max(ilo, 0);
            i2hi = Math.Min(ihi, n - 1);

            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = i.ToString().PadLeft(4) + "  ";

                for (j = j2lo; j <= j2hi; j++)
                {
                    if (i < j)
                    {
                        cout += "              ";
                    }
                    else
                    {
                        cout += a[i - j].ToString().PadLeft(12) + "  ";
                    }
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8ltt_random(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LTT_RANDOM randomizes an R8LTT matrix.
        //
        //  Discussion:
        //
        //    The R8LTT storage format is used for an N by N lower triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8LTT_RANDOM[N], the R8LTT matrix.
        //
    {
        double[] a;

        a = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        return a;
    }

    public static double[] r8ltt_sl(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LTT_SL solves a linear system A*x=b with an R8LTT matrix.
        //
        //  Discussion:
        //
        //    The R8LTT storage format is used for an N by N lower triangular Toeplitz
        //    matrix.
        //
        //    No factorization of the lower triangular matrix is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the R8LTT matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R8LTT_SL[N], the solution vector.
        //
    {
        int i;
        int j;
        double[] x;

        x = new double[n];

        for (j = 0; j < n; j++)
        {
            x[j] = b[j];
        }

        for (j = 0; j < n; j++)
        {
            x[j] /= a[0];
            for (i = j + 1; i < n; i++)
            {
                x[i] -= a[i - j] * x[j];
            }
        }

        return x;
    }

    public static double[] r8ltt_slt(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LTT_SLT solves a linear system A'*x=b with an R8LTT matrix.
        //
        //  Discussion:
        //
        //    The R8LTT storage format is used for an N by N lower triangular Toeplitz
        //    matrix.
        //
        //    No factorization of the lower triangular matrix is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the R8LTT matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R8LTT_SLT[N], the solution vector.
        //
    {
        int i;
        int j;
        double[] x;

        x = new double[n];
        for (j = 0; j < n; j++)
        {
            x[j] = b[j];
        }

        for (i = n - 1; 0 <= i; i--)
        {
            x[i] /= a[0];
            for (j = 0; j < i; j++)
            {
                x[j] -= a[i - j] * x[i];
            }
        }

        return x;
    }

    public static double[] r8ltt_to_r8ge(int n, double[] a_utt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LTT_TO_R8GE copies an R8LTT matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a general M by N matrix.  A storage 
        //    space is made for each entry.  The two dimensional logical
        //    array can be thought of as a vector of M*N entries, starting with
        //    the M entries in the column 1, then the M entries in column 2
        //    and so on.  Considered as a vector, the entry A(I,J) is then stored
        //    in vector location I+(J-1)*M.
        //
        //    The R8LTT storage format is used for an N by N lower triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A_UTT[N], the R8LTT matrix.
        //
        //    Output, double R8LTT_TO_R8GE[N*N], the R8GE matrix.
        //
    {
        double[] a_ge;
        int d;
        int i;
        int j;

        a_ge = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                a_ge[i + j * n] = 0.0;
            }
        }

        for (d = 0; d < n; d++)
        {
            for (i = d; i < n; i++)
            {
                j = i - d;
                a_ge[i + j * n] = a_utt[i - j];
            }
        }

        return a_ge;
    }

    public static double[] r8ltt_zeros(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LTT_ZEROS zeros an R8LTT matrix.
        //
        //  Discussion:
        //
        //    The R8LTT storage format is used for an N by N lower triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Output, double R8LTT_ZEROS[M*N], the R8LTT matrix.
        //
    {
        double[] a;
        int i;

        a = new double[n];

        for (i = 0; i < n; i++)
        {
            a[i] = 0.0;
        }

        return a;
    }
}