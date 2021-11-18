using System;
using System.Numerics;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8ci_det(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CI_DET returns the determinant of an R8CI matrix.
        //
        //  Discussion:
        //
        //    The R8CI storage format is used for an N by N circulant matrix.
        //    An N by N circulant matrix A has the property that the entries on
        //    row I appear again on row I+1, shifted one position to the right,
        //    with the final entry of row I appearing as the first of row I+1.
        //
        //    A circulant matrix data structure simply records the first row.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis,
        //    Circulant Matrices,
        //    Wiley, 1979.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the R8CI matrix.
        //
        //    Output, double R8CI_DET, the complex eigenvalues.
        //
    {
        int i;

        Complex[] lambda = r8ci_eval(n, a);

        Complex detc = 1.0;
        for (i = 0; i < n; i++)
        {
            detc *= lambda[i];
        }

        double det = detc.Real;

        return det;
    }

    public static double[] r8ci_dif2(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CI_DIF2 sets up an R8CI second difference matrix.
        //
        //  Discussion:
        //
        //    This is actually a periodic second difference matrix.
        //
        //    The R8CI storage format is used for an N by N circulant matrix.
        //    An N by N circulant matrix A has the property that the entries on
        //    row I appear again on row I+1, shifted one position to the right,
        //    with the final entry of row I appearing as the first of row I+1.
        //    The R8CI format simply records the first row of the matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 3.
        //
        //    Output, double R8CI_DIF2[N], the R8CI matrix.
        //
    {
        double[] a = r8vec_zeros_new(n);

        a[0] = 2.0;
        a[1] = -1.0;
        a[n - 1] = -1.0;

        return a;
    }

    public static Complex[] r8ci_eval(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CI_EVAL returns the eigenvalues of an R8CI matrix.
        //
        //  Discussion:
        //
        //    The R8CI storage format is used for an N by N circulant matrix.
        //    An N by N circulant matrix A has the property that the entries on
        //    row I appear again on row I+1, shifted one position to the right,
        //    with the final entry of row I appearing as the first of row I+1.
        //
        //    A circulant matrix data structure simply records the first row.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis,
        //    Circulant Matrices,
        //    Wiley, 1979.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the R8CI matrix.
        //
        //    Output, Complex R8CI_EVAL[N], the complex eigenvalues.
        //
    {
        int i;

        Complex[] lambda = new Complex[n];

        Complex[] w = c8vec_unity_new(n);

        for (i = 0; i < n; i++)
        {
            lambda[i] = a[n - 1];
        }

        for (i = n - 2; 0 <= i; i--)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                lambda[j] = lambda[j] * w[j] + a[i];
            }
        }

        c8vec_sort_a_l2(n, ref lambda);

        return lambda;
    }

    public static double[] r8ci_indicator(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CI_INDICATOR sets up an R8CI indicator matrix.
        //
        //  Discussion:
        //
        //    The R8CI storage format is used for an N by N circulant matrix.
        //    An N by N circulant matrix A has the property that the entries on
        //    row I appear again on row I+1, shifted one position to the right,
        //    with the final entry of row I appearing as the first of row I+1.
        //    The R8CI format simply records the first row of the matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 January 2004
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
        //    Output, double R8CI_INDICATOR[N], the R8CI matrix.
        //
    {
        int j;

        double[] a = new double[n];

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        int i = 1;

        for (j = 1; j <= n; j++)
        {
            a[j - 1] = fac * i + j;
        }

        return a;
    }

    public static double[] r8ci_mtv(int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CI_MTV multiplies a vector times an R8CI matrix.
        //
        //  Discussion:
        //
        //    The R8CI storage format is used for an N by N circulant matrix.
        //    An N by N circulant matrix A has the property that the entries on
        //    row I appear again on row I+1, shifted one position to the right,
        //    with the final entry of row I appearing as the first of row I+1.
        //
        //    A circulant matrix data structure simply records the first row.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 December 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the R8CI matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8CI_MTV[N], the product A' * X.
        //
    {
        int i;

        double[] b = r8vec_zeros_new(n);

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j <= i; j++)
            {
                b[i] += a[i - j] * x[j];
            }

            for (j = i + 1; j < n; j++)
            {
                b[i] += a[n + i - j] * x[j];
            }
        }

        return b;
    }

    public static double[] r8ci_mv(int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CI_MV multiplies an R8CI matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8CI storage format is used for an N by N circulant matrix.
        //    An N by N circulant matrix A has the property that the entries on
        //    row I appear again on row I+1, shifted one position to the right,
        //    with the final entry of row I appearing as the first of row I+1.
        //
        //    A circulant matrix data structure simply records the first row.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the R8CI matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8CI_MV[N], the product A * x.
        //
    {
        int i;

        double[] b = r8vec_zeros_new(n);

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j <= i - 1; j++)
            {
                b[i] += a[j - i + n] * x[j];
            }

            for (j = i; j < n; j++)
            {
                b[i] += a[j - i] * x[j];
            }
        }

        return b;
    }

    public static void r8ci_print(int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CI_PRINT prints an R8CI matrix.
        //
        //  Discussion:
        //
        //    The R8CI storage format is used for an N by N circulant matrix.
        //    An N by N circulant matrix A has the property that the entries on
        //    row I appear again on row I+1, shifted one position to the right,
        //    with the final entry of row I appearing as the first of row I+1.
        //
        //    A circulant matrix data structure simply records the first row.
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
        //    Input, double A[N], the R8CI matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8ci_print_some(n, a, 1, 1, n, n, title);
    }

    public static void r8ci_print_some(int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CI_PRINT_SOME prints some of an R8CI matrix.
        //
        //  Discussion:
        //
        //    The R8CI storage format is used for an N by N circulant matrix.
        //    An N by N circulant matrix A has the property that the entries on
        //    row I appear again on row I+1, shifted one position to the right,
        //    with the final entry of row I appearing as the first of row I+1.
        //
        //    A circulant matrix data structure simply records the first row.
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
        //    Input, double A[N], the R8CI matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

        int j;
        int j2lo;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            int j2hi = j2lo + INCX - 1;
            j2hi = Math.Min(j2hi, n);
            j2hi = Math.Min(j2hi, jhi);

            Console.WriteLine("");
            string cout = "  Col: ";
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
            int i2lo = Math.Max(ilo, 1);
            int i2hi = Math.Min(ihi, n);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = i.ToString().PadLeft(6) + "  ";

                for (j = j2lo; j <= j2hi; j++)
                {
                    if (i <= j)
                    {
                        cout += a[j - i].ToString().PadLeft(12) + "  ";
                    }
                    else
                    {
                        cout += a[n + j - i].ToString().PadLeft(12) + "  ";
                    }
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8ci_random(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CI_RANDOM randomizes an R8CI matrix.
        //
        //  Discussion:
        //
        //    The R8CI storage format is used for an N by N circulant matrix.
        //    An N by N circulant matrix A has the property that the entries on
        //    row I appear again on row I+1, shifted one position to the right,
        //    with the final entry of row I appearing as the first of row I+1.
        //
        //    A circulant matrix data structure simply records the first row.
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
        //    Output, double R8CI_RANDOM[N], the R8CI matrix.
        //
    {
        int i;

        double[] a = new double[n];

        for (i = 0; i < n; i++)
        {
            a[i] = UniformRNG.r8_uniform_01(ref seed);
        }

        return a;
    }

    public static double[] r8ci_sl(int n, double[] a, double[] b, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CI_SL solves an R8CI system.
        //
        //  Discussion:
        //
        //    The R8CI storage format is used for an N by N circulant matrix.
        //    An N by N circulant matrix A has the property that the entries on
        //    row I appear again on row I+1, shifted one position to the right,
        //    with the final entry of row I appearing as the first of row I+1.
        //
        //    A circulant matrix data structure simply records the first row.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 September 2003
        //
        //  Author:
        //
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the R8CI matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Input, int JOB, specifies the system to solve.
        //    0, solve A * x = b.
        //    nonzero, solve A' * x = b.
        //
        //    Output, double R8CI_SL[N], the solution of the linear system.
        //
    {
        int i;
        int nsub;
        double r1;
        double r2;
        double r3;
        double r5;
        double r6;

        double[] work = new double[2 * n - 2];
        double[] x = new double[n];

        switch (job)
        {
            case 0:
            {
                //
                //  Solve the system with the principal minor of order 1.
                //
                r1 = a[0];
                x[0] = b[0] / r1;

                r2 = 0.0;
                //
                //  Recurrent process for solving the system.
                //
                for (nsub = 2; nsub <= n; nsub++)
                {
                    //
                    //  Compute multiples of the first and last columns of
                    //  the inverse of the principal minor of order N.
                    //
                    r5 = a[n + 2 - nsub - 1];
                    r6 = a[nsub - 1];

                    switch (nsub)
                    {
                        case > 2:
                        {
                            work[nsub - 2] = r2;

                            for (i = 1; i <= nsub - 2; i++)
                            {
                                r5 += a[n - i] * work[nsub - i - 1];
                                r6 += a[i] * work[n - 2 + i];
                            }

                            break;
                        }
                    }

                    r2 = -r5 / r1;
                    r3 = -r6 / r1;
                    r1 += r5 * r3;

                    switch (nsub)
                    {
                        case > 2:
                        {
                            r6 = work[n - 1];
                            work[n + nsub - 3] = 0.0;
                            for (i = 2; i <= nsub - 1; i++)
                            {
                                r5 = work[n - 2 + i];
                                work[n - 2 + i] = work[i - 1] * r3 + r6;
                                work[i - 1] += r6 * r2;
                                r6 = r5;
                            }

                            break;
                        }
                    }

                    work[n - 1] = r3;
                    //
                    //  Compute the solution of the system with the principal minor of order NSUB.
                    //
                    r5 = 0.0;
                    for (i = 1; i <= nsub - 1; i++)
                    {
                        r5 += a[n - i] * x[nsub - i - 1];
                    }

                    r6 = (b[nsub - 1] - r5) / r1;
                    for (i = 1; i <= nsub - 1; i++)
                    {
                        x[i - 1] += work[n + i - 2] * r6;
                    }

                    x[nsub - 1] = r6;
                }

                break;
            }
            default:
            {
                //
                //  Solve the system with the principal minor of order 1.
                //
                r1 = a[0];
                x[0] = b[0] / r1;

                r2 = 0.0;
                //
                //  Recurrent process for solving the system.
                //
                for (nsub = 2; nsub <= n; nsub++)
                {
                    //
                    //  Compute multiples of the first and last columns of
                    //  the inverse of the principal minor of order N.
                    //
                    r5 = a[nsub - 1];
                    r6 = a[n + 1 - nsub];

                    switch (nsub)
                    {
                        case > 2:
                        {
                            work[nsub - 2] = r2;
                            for (i = 1; i <= nsub - 2; i++)
                            {
                                r5 += a[i] * work[nsub - i - 1];
                                r6 += a[n - i] * work[n - 2 + i];
                            }

                            break;
                        }
                    }

                    r2 = -r5 / r1;
                    r3 = -r6 / r1;
                    r1 += r5 * r3;

                    switch (nsub)
                    {
                        case > 2:
                        {
                            r6 = work[n - 1];
                            work[n + nsub - 3] = 0.0;
                            for (i = 2; i <= nsub - 1; i++)
                            {
                                r5 = work[n - 2 + i];
                                work[n - 2 + i] = work[i - 1] * r3 + r6;
                                work[i - 1] += r6 * r2;
                                r6 = r5;
                            }

                            break;
                        }
                    }

                    work[n - 1] = r3;
                    //
                    //  Compute the solution of the system with the principal minor of order NSUB.
                    //
                    r5 = 0.0;
                    for (i = 1; i <= nsub - 1; i++)
                    {
                        r5 += a[i] * x[nsub - i - 1];
                    }

                    r6 = (b[nsub - 1] - r5) / r1;
                    for (i = 1; i <= nsub - 1; i++)
                    {
                        x[i - 1] += work[n - 2 + i] * r6;
                    }

                    x[nsub - 1] = r6;
                }

                break;
            }
        }

        return x;
    }

    public static double[] r8ci_to_r8ge(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CI_TO_R8GE copies an R8CI matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8CI storage format is used for an N by N circulant matrix.
        //    An N by N circulant matrix A has the property that the entries on
        //    row I appear again on row I+1, shifted one position to the right,
        //    with the final entry of row I appearing as the first of row I+1.
        //
        //    A circulant matrix data structure simply records the first row.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N], the R8CI matrix.
        //
        //    Output, double R8CI_TO_R8GE[N*N], the R8GE matrix.
        //
    {
        int j;

        double[] b = new double[n * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                int k = (j - i) % n;
                b[i + j * n] = a[k];
            }
        }

        return b;
    }

    public static double[] r8ci_zeros(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CI_ZEROS zeros an R8CI matrix.
        //
        //  Discussion:
        //
        //    The R8CI storage format is used for an N by N circulant matrix.
        //    An N by N circulant matrix A has the property that the entries on
        //    row I appear again on row I+1, shifted one position to the right,
        //    with the final entry of row I appearing as the first of row I+1.
        //
        //    A circulant matrix data structure simply records the first row.
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
        //    Output, double R8CI_ZERO[N], the R8CI matrix.
        //
    {
        double[] a = r8vec_zeros_new(n);

        return a;
    }
}