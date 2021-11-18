using System;
using System.Globalization;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8utt_det(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTT_DET computes the determinant of a R8UTT matrix.
        //
        //  Discussion:
        //
        //    The R8UTT storage format is used for an N by N upper triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2015
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
        //    Output, double R8UTT_DET, the determinant of the matrix.
        //
    {
        double det = Math.Pow(a[0], n);

        return det;
    }

    public static double[] r8utt_indicator(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTT_INDICATOR sets up a R8UTT indicator matrix.
        //
        //  Discussion:
        //
        //    The R8UTT storage format is used for an N by N upper triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2015
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
        int j;

        double[] a = new double[n];

        for (j = 0; j < n; j++)
        {
            a[j] = j + 1;
        }

        return a;
    }

    public static double[] r8utt_inverse(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTT_INVERSE computes the inverse of a R8UTT matrix.
        //
        //  Discussion:
        //
        //    The R8UTT storage format is used for an N by N upper triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2015
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
        //    Output, double R8UTT_INVERSE[N], the inverse matrix.
        //
    {
        int i;
        int j;
        //
        //  Initialize B.
        //
        double d = 1.0 / a[0];
        double[] b = new double[n];
        b[0] = d;
        for (i = 1; i < n; i++)
        {
            b[i] = 0.0;
        }

        //
        //  Set the strict upper triangle.
        //
        double[] p = new double[n];
        p[0] = 0.0;
        for (i = 1; i < n; i++)
        {
            p[i] = a[i];
        }

        //
        //  PN will hold powers of P.
        //
        double[] pn = new double[n];
        pn[0] = 1.0;
        for (i = 1; i < n; i++)
        {
            pn[i] = 0.0;
        }

        //
        //  Add N-1 powers of strict upper triangle.
        //
        for (j = 1; j < n; j++)
        {
            d = -d / a[0];
            double[] pnn = r8utt_mm(n, p, pn);
            for (i = 0; i < n; i++)
            {
                b[i] += d * pnn[i];
                pn[i] = pnn[i];
            }
        }

        return b;
    }

    public static double[] r8utt_mm(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTT_MM computes C = A * B, where A and B are R8UTT matrices.
        //
        //  Discussion:
        //
        //    The R8UTT storage format is used for an N by N upper triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2015
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
        //    Output, double R8UTT_MM[N], the product.
        //
    {
        int k;

        double[] d = new double[n];

        for (k = 0; k < n; k++)
        {
            d[k] = b[n - 1 - k];
        }

        double[] e = r8utt_mv(n, a, d);

        double[] c = new double[n];

        for (k = 0; k < n; k++)
        {
            c[k] = e[n - 1 - k];
        }

        return c;
    }

    public static double[] r8utt_mtm(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTT_MTM computes C = A' * B, where A and B are R8UTT matrices.
        //
        //  Discussion:
        //
        //    The R8UTT storage format is used for an N by N upper triangular Toeplitz
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
        //    16 November 2015
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
        //    Output, double R8UTT_MTM[N*N], the product.
        //
    {
        int i;

        double[] c = new double[n * n];

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                c[i + j * n] = 0.0;
                int k;
                for (k = 0; k <= Math.Min(i, j); k++)
                {
                    c[i + j * n] += a[i - k] * b[j - k];
                }
            }
        }

        return c;
    }

    public static double[] r8utt_mtv(int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTT_MTV computes b = A'*x, where A is an R8UTT matrix.
        //
        //  Discussion:
        //
        //    The R8UTT storage format is used for an N by N upper triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2015
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
        //    Output, double R8UTT_MTV[N], the product A' * x.
        //
    {
        int d;
        int j;

        double[] b = new double[n];

        for (j = 0; j < n; j++)
        {
            b[j] = 0.0;
        }

        for (d = 0; d < n; d++)
        {
            for (j = d; j < n; j++)
            {
                int i = j - d;
                b[j] += a[j - i] * x[i];
            }
        }

        return b;
    }

    public static double[] r8utt_mv(int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTT_MV computes b=A*x, where A is an R8UTT matrix.
        //
        //  Discussion:
        //
        //    The R8UTT storage format is used for an N by N upper triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2015
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
        //    Output, double R8UTT_MV[N], the product A * x.
        //
    {
        int d;
        int i;

        double[] b = new double[n];

        for (i = 0; i < n; i++)
        {
            b[i] = 0.0;
        }

        for (d = 0; d < n; d++)
        {
            int j;
            for (j = d; j < n; j++)
            {
                i = j - d;
                b[i] += a[j - i] * x[j];
            }
        }

        return b;
    }

    public static void r8utt_print(int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTT_PRINT prints an R8UTT matrix.
        //
        //  Discussion:
        //
        //    The R8UTT storage format is used for an N by N upper triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8UTT matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8utt_print_some(n, a, 0, 0, n - 1, n - 1, title);

    }

    public static void r8utt_print_some(int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTT_PRINT_SOME prints some of an R8UTT matrix.
        //
        //  Discussion:
        //
        //    The R8UTT storage format is used for an N by N upper triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8UTT matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //    0 <= ILO <= IHI < M.
        //    0 <= JLO <= JHI < N.
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
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";

                for (j = j2lo; j <= j2hi; j++)
                {
                    if (j < i)
                    {
                        cout += "              ";
                    }
                    else
                    {
                        cout += a[j - i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                    }
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8utt_random(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTT_RANDOM randomizes an R8UTT matrix.
        //
        //  Discussion:
        //
        //    The R8UTT storage format is used for an N by N upper triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2015
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
        //    Output, double R8UTT_RANDOM[N], the R8UTT matrix.
        //
    {
        double[] a = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        return a;
    }

    public static double[] r8utt_sl(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTT_SL solves a linear system A*x=b with an R8UTT matrix.
        //
        //  Discussion:
        //
        //    The R8UTT storage format is used for an N by N upper triangular Toeplitz
        //    matrix.
        //
        //    No factorization of the upper triangular matrix is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the R8UTT matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R8UTT_SL[N], the solution vector.
        //
    {
        int j;

        double[] x = new double[n];

        for (j = 0; j < n; j++)
        {
            x[j] = b[j];
        }

        for (j = n - 1; 0 <= j; j--)
        {
            x[j] /= a[0];
            int i;
            for (i = 0; i < j; i++)
            {
                x[i] -= a[j - i] * x[j];
            }
        }

        return x;
    }

    public static double[] r8utt_slt(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTT_SLT solves a linear system A'*x=b with an R8UTT matrix.
        //
        //  Discussion:
        //
        //    The R8UTT storage format is used for an N by N upper triangular Toeplitz
        //    matrix.
        //
        //    No factorization of the upper triangular matrix is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the R8UTT matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R8UTT_SLT[N], the solution vector.
        //
    {
        int j;

        double[] x = new double[n];
        for (j = 0; j < n; j++)
        {
            x[j] = b[j];
        }

        for (j = 0; j < n; j++)
        {
            x[j] /= a[0];
            int i;
            for (i = j + 1; i < n; i++)
            {
                x[i] -= x[j] * a[i - j];
            }
        }

        return x;
    }

    public static double[] r8utt_to_r8ge(int n, double[] a_utt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTT_TO_R8GE copies an R8UTT matrix to an R8GE matrix.
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
        //    The R8UTT storage format is used for an N by N upper triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A_UTT[N], the R8UTT matrix.
        //
        //    Output, double R8UTT_TO_R8GE[N*N], the R8GE matrix.
        //
    {
        int d;
        int i;
        int j;

        double[] a_ge = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                a_ge[i + j * n] = 0.0;
            }
        }

        for (d = 0; d < n; d++)
        {
            for (j = d; j < n; j++)
            {
                i = j - d;
                a_ge[i + j * n] = a_utt[j - i];
            }
        }

        return a_ge;
    }

    public static double[] r8utt_zeros(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTT_ZEROS zeros an R8UTT matrix.
        //
        //  Discussion:
        //
        //    The R8UTT storage format is used for an N by N upper triangular Toeplitz
        //    matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Output, double R8UTT_ZEROS[M*N], the R8UTT matrix.
        //
    {
        int i;

        double[] a = new double[n];

        for (i = 0; i < n; i++)
        {
            a[i] = 0.0;
        }

        return a;
    }
}