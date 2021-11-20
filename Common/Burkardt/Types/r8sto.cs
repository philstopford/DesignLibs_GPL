using System;
using System.Globalization;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8sto_dif2(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8STO_DIF2 sets the second difference as an R8STO matrix.
        //
        //  Discussion:
        //
        //    The R8STO storage format is used for a symmetric Toeplitz matrix.
        //    It stores the N elements of the first row.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2015
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
        //    Output, double R8STO_DIF2[N], the R8STO matrix.
        //
    {
        double[] a = r8vec_zeros_new(n);
        a[0] = 2.0;
        a[1] = -1.0;

        return a;
    }

    public static double[] r8sto_indicator(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8STO_INDICATOR sets up an R8STO indicator matrix.
        //
        //  Discussion:
        //
        //    The R8STO storage format is used for a symmetric Toeplitz matrix.
        //    It stores the N elements of the first row.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 January 2004
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
        //    Output, double R8STO_INDICATOR[N], the R8STO matrix.
        //
    {
        int j;

        double[] a = new double[n];

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        const int i = 1;
        int k = 0;
        for (j = 1; j <= n; j++)
        {
            a[k] = fac * i + j;
            k += 1;
        }

        return a;
    }

    public static double[] r8sto_inverse(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8STO_INVERSE computes the inverse of an R8STO matrix.
        //
        //  Discussion:
        //
        //    The R8STO storage format is used for a symmetric Toeplitz matrix.
        //    It stores the N elements of the first row.
        //
        //    For this routine, the matrix is also required to be positive definite.
        //
        //    The original implementation of the algorithm assumed that the
        //    diagonal element was 1.  The algorithm has been modified so that
        //    this is no longer necessary.
        //
        //    The inverse matrix is NOT guaranteed to be a Toeplitz matrix.  
        //    It is guaranteed to be symmetric and persymmetric.
        //    The inverse matrix is returned in general storage, that is,
        //    as an "SGE" matrix.
        //
        //  Example:
        //
        //    To compute the inverse of
        //
        //     1.0 0.5 0.2
        //     0.5 1.0 0.5
        //     0.2 0.5 1.0
        //
        //    we input:
        //
        //      N = 3
        //      A = { 1.0, 0.5, 0.2 }
        //
        //    with output:
        //
        //      B = ( 1/56) * [ 75, -40,   5,
        //                     -40,  96, -40,
        //                       5, -40,  75 ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Gene Golub, Charles Van Loan,
        //    Section 4.7.3, "Computing the Inverse",
        //    Matrix Computations,
        //    Third Edition,
        //    Johns Hopkins, 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the system.
        //
        //    Input, double A[N], the R8STO matrix.
        //
        //    Output, double R8STO_INVERSE[N*N], the inverse of the matrix.
        //
    {
        int i;
        int j;

        double[] a2 = new double[n - 1];
        double[] b = new double[n * n];

        for (i = 0; i < n - 1; i++)
        {
            a2[i] = a[i + 1] / a[0];
        }

        double[] v = r8sto_yw_sl(n - 1, a2);
        //
        //  Compute the N-th entry of V.
        //
        double t = 0.0;
        for (i = 0; i < n - 1; i++)
        {
            t += a2[i] * v[i];
        }

        double vn = 1.0 / (1.0 + t);
        //
        //  Reverse the first N-1 entries of V.
        //
        for (i = 0; i < (n - 1) / 2; i++)
        {
            j = n - 2 - i;
            t = v[i];
            v[i] = v[j];
            v[j] = t;
        }

        //
        //  Scale the entries.
        //
        for (i = 0; i < n - 1; i++)
        {
            v[i] = vn * v[i];
        }

        //
        //  Set the boundaries of B.
        //
        b[0 + 0 * n] = vn;
        for (j = 1; j < n; j++)
        {
            b[0 + j * n] = v[n - j - 1];
        }

        for (j = 0; j < n - 1; j++)
        {
            b[n - 1 + j * n] = v[j];
        }

        b[n - 1 + (n - 1) * n] = vn;

        for (i = 1; i < n - 1; i++)
        {
            b[i + 0 * n] = v[n - 1 - i];
            b[i + (n - 1) * n] = v[i];
        }

        //
        //  Fill the interior.
        //
        for (i = 2; i <= 1 + (n - 1) / 2; i++)
        {
            for (j = i; j <= n - i + 1; j++)
            {
                t = b[i - 2 + (j - 2) * n] + (v[n - j] * v[n - i] - v[i - 2] * v[j - 2]) / vn;
                b[i - 1 + (j - 1) * n] = t;
                b[j - 1 + (i - 1) * n] = t;
                b[n - i + (n - j) * n] = t;
                b[n - j + (n - i) * n] = t;
            }
        }

        //
        //  Scale B.
        //
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                b[i + j * n] /= a[0];
            }
        }

        return b;
    }

    public static double[] r8sto_mv(int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8STO_MV multiplies an R8STO matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8STO storage format is used for a symmetric Toeplitz matrix.
        //    It stores the N elements of the first row.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the R8STO matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8STO_MV[N], the product A * x.
        //
    {
        int i;

        double[] b = r8vec_zeros_new(n);

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j <= i - 1; j++)
            {
                b[i] += a[i - j] * x[j];
            }

            for (j = i; j < n; j++)
            {
                b[i] += a[j - i] * x[j];
            }

        }

        return b;
    }

    public static void r8sto_print(int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8STO_PRINT prints an R8STO matrix.
        //
        //  Discussion:
        //
        //    The R8STO storage format is used for a symmetric Toeplitz matrix.
        //    It stores the N elements of the first row.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2015
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
        //    Input, double A[N], the R8STO matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8sto_print_some(n, a, 0, 0, n - 1, n - 1, title);
    }

    public static void r8sto_print_some(int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8STO_PRINT_SOME prints some of am R8STO matrix.
        //
        //  Discussion:
        //
        //    The R8STO storage format is used for a symmetric Toeplitz matrix.
        //    It stores the N elements of the first row.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2015
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
        //    Input, double A[N], the R8STO matrix.
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
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                for (j = j2lo; j <= j2hi; j++)
                {
                    double aij = i <= j ? a[j - i] : a[i - j];

                    cout += aij.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }

        Console.WriteLine("");

    }

    public static double[] r8sto_random(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8STO_RANDOM randomizes an R8STO matrix.
        //
        //  Discussion:
        //
        //    The R8STO storage format is used for a symmetric Toeplitz matrix.
        //    It stores the N elements of the first row.
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
        //    Output, double R8STO_RANDOM[N], the R8STO matrix.
        //
    {
        double[] r = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        return r;
    }

    public static double[] r8sto_sl(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8STO_SL solves an R8STO system.
        //
        //  Discussion:
        //
        //    The R8STO storage format is used for a symmetric Toeplitz matrix.
        //    It stores the N elements of the first row.
        //
        //    The matrix is also required to be positive definite.
        //
        //    This implementation of the algorithm assumes that the diagonal element
        //    (the first element of A) is 1.
        //
        //    Note that there is a typographical error in the presentation
        //    of this algorithm in the reference, and another in the presentation
        //    of a sample problem.  Both involve sign errors.  A minor error
        //    makes the algorithm incorrect for the case N = 1.
        //
        //  Example:
        //
        //    To solve
        //
        //     1.0 0.5 0.2    x1    4.0
        //     0.5 1.0 0.5 *  x2 = -1.0
        //     0.2 0.5 1.0    x3    3.0
        //
        //    we input:
        //
        //      N = 3
        //      A = (/ 1.0, 0.5, 0.2 /)
        //      B = (/ 4.0, -1.0, 3.0 /)
        //
        //    with output:
        //
        //      X = (/ 355, -376, 285 /) / 56
        //        = (/ 6.339, -6.714, 5.089 /)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Gene Golub, Charles Van Loan,
        //    Section 4.7.3, "The General Right Hand Side Problem",
        //    Matrix Computations,
        //    Third Edition,
        //    Johns Hopkins, 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the system.
        //
        //    Input, double A[N], the R8STO matrix, with the EXTRA CONDITION
        //    that the first entry is 1.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Output, double R8STO_SL[N], the solution of the linear system.
        //
    {
        double[] x = new double[n];
        double[] y = new double[n];

        int k = 0;
        double beta = 1.0;
        x[k] = b[k] / beta;

        if (k < n - 1)
        {
            y[k] = -a[1] / beta;
        }

        for (k = 1; k <= n - 1; k++)
        {
            beta = (1.0 - y[k - 1] * y[k - 1]) * beta;

            x[k] = b[k];
            int i;
            for (i = 1; i <= k; i++)
            {
                x[k] -= a[i] * x[k - i];
            }

            x[k] /= beta;

            for (i = 1; i <= k; i++)
            {
                x[i - 1] += x[k] * y[k - i];
            }

            if (k >= n - 1)
            {
                continue;
            }

            y[k] = -a[k + 1];
            for (i = 1; i <= k; i++)
            {
                y[k] -= a[i] * y[k - i];
            }

            y[k] /= beta;

            for (i = 1; i <= k; i++)
            {
                y[i - 1] += y[k] * y[k - i];
            }
        }

        return x;
    }

    public static double[] r8sto_to_r8ge(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8STO_TO_R8GE copies an R8STO matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8STO storage format is used for a symmetric Toeplitz matrix.
        //    It stores the N elements of the first row.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double R8STO_TO_R8GE[N], the R8STO matrix.
        //
        //    Output, double R8STO_TO_R8GE[N*N], the R8GE matrix.
        //
    {
        int i;

        double[] b = new double[n * n];

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j < i; j++)
            {
                b[i + j * n] = a[i - j];
            }

            for (j = i; j < n; j++)
            {
                b[i + j * n] = a[j - i];
            }
        }

        return b;
    }

    public static double[] r8sto_yw_sl(int n, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8STO_YW_SL solves the Yule-Walker equations for an R8STO matrix.
        //
        //  Discussion:
        //
        //    The R8STO storage format is used for a symmetric Toeplitz matrix.
        //    It stores the N elements of the first row.
        //
        //    The matrix is also required to be positive definite.
        //
        //    This implementation of the algorithm assumes that the diagonal element
        //    is 1.
        //
        //    The real symmetric Toeplitz matrix can be described by N numbers, which,
        //    for convenience, we will label B(0:N-1).  We assume there is one more
        //    number, B(N).  If we let A be the symmetric Toeplitz matrix whose first
        //    row is B(0:N-1), then the Yule-Walker equations are:
        //
        //      A * X = -B(1:N)
        //
        //  Example:
        //
        //    To solve
        //
        //     1.0 0.5 0.2    x1   0.5
        //     0.5 1.0 0.5 *  x2 = 0.2
        //     0.2 0.5 1.0    x3   0.1
        //
        //    we input:
        //
        //      N = 3
        //      B = (/ 0.5, 0.2, 0.1 /)
        //
        //    with output:
        //
        //      X = (/ -75, 12, -5 /) / 140
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Gene Golub, Charles Van Loan,
        //    Section 4.7.2, "Solving the Yule-Walker Equations",
        //    Matrix Computations,
        //    Third Edition,
        //    Johns Hopkins, 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the system.
        //
        //    Input, double B[N], defines the linear system.  The first entry of the
        //    symmetric Toeplitz matrix is assumed to be a 1, which is NOT stored.  The N-1
        //    remaining elements of the first row of are stored in B, followed by
        //    the remaining scalar that defines the linear system.
        //
        //    Output, double R8STO_YW_SL[N], the solution of the linear system.
        //
    {
        int i;

        double[] x = new double[n];
        double[] x2 = new double[n];

        x[0] = -b[0];
        double beta = 1.0;
        double alpha = -b[0];

        for (i = 1; i <= n - 1; i++)
        {
            beta = (1.0 - alpha * alpha) * beta;

            alpha = b[i];
            int j;
            for (j = 1; j <= i; j++)
            {
                alpha += b[i - j] * x[j - 1];
            }

            alpha = -alpha / beta;

            for (j = 1; j <= i; j++)
            {
                x2[j - 1] = x[j - 1];
            }

            for (j = 1; j <= i; j++)
            {
                x[j - 1] += alpha * x2[i - j];
            }

            x[i] = alpha;
        }

        return x;
    }

    public static double[] r8sto_zeros(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8STO_ZEROS zeros an R8STO matrix.
        //
        //  Discussion:
        //
        //    The R8STO storage format is used for a symmetric Toeplitz matrix.
        //    It stores the N elements of the first row.
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
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Output, double R8STO_ZERO[N], the R8STO matrix.
        //
    {
        double[] a = r8vec_zeros_new(n);

        return a;
    }
}