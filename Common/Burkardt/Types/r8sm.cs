using System;
using Burkardt.Uniform;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8sm_indicator(int m, int n, ref double[] a, ref double[] u, ref double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8SM_INDICATOR returns the indicator matrix as an R8SM matrix.
            //
            //  Discussion:
            //
            //    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
            //    which is defined by an M by N matrix A, an M vector U, and
            //    an N vector V, by B = A - U * V'
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns 
            //    of the matrix.
            //
            //    Output, double A[M*N], the R8SM matrix.
            //
            //    Output, double U[M], V[N], the R8SM vectors.
            //
        {
            int fac;
            int i;
            int j;

            fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

            for (i = 0; i < m; i++)
            {
                u[i] = -1.0;
            }

            for (j = 0; j < n; j++)
            {
                v[j] = (double) (j + 1);
            }

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    a[i + j * m] = (double) (fac * (i + 1));
                }
            }

        }

        public static double[] r8sm_ml(int n, double[] a_lu, double[] u, double[] v, int[] pivot,
                double[] x, int job)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8SM_ML multiplies a factored square R8SM matrix times a vector.
            //
            //  Discussion:
            //
            //    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
            //    which is defined by an M by N matrix A, an M vector U, and
            //    an N vector V, by B = A - U * V'
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 October 2003
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
            //    Input, double A_LU[N*N], the LU factors from R8GE_FA.
            //
            //    Input, double U[N], V[N], the Sherman Morrison vectors.
            //
            //    Input, int PIVOT[N], the pivot vector computed by R8GE_FA.
            //
            //    Input, double X[N], the vector to be multiplied.
            //
            //    Input, int JOB, specifies the operation to be done:
            //    JOB = 0, compute (A-u*v') * x.
            //    JOB nonzero, compute (A-u*v')' * x.
            //
            //    Output, double R8SM_ML[N], the result of the multiplication.
            //
        {
            double[] b;
            int i;
            double ux;
            double vx;

            b = r8ge_ml(n, a_lu, pivot, x, job);

            if (job == 0)
            {
                vx = 0.0;
                for (i = 0; i < n; i++)
                {
                    vx = vx + v[i] * x[i];
                }

                for (i = 0; i < n; i++)
                {
                    b[i] = b[i] - u[i] * vx;
                }
            }
            else
            {
                ux = 0.0;
                for (i = 0; i < n; i++)
                {
                    ux = ux + u[i] * x[i];
                }

                for (i = 0; i < n; i++)
                {
                    b[i] = b[i] - v[i] * ux;
                }
            }

            return b;
        }

        public static double[] r8sm_mtv(int m, int n, double[] a, double[] u, double[] v, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8SM_MTV multiplies a vector by an R8SM matrix.
            //
            //  Discussion:
            //
            //    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
            //    which is defined by an M by N matrix A, an M vector U, and
            //    an N vector V, by B = A - U * V'
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 January 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns of the matrix.
            //
            //    Input, double A[M*N], the R8SM matrix.
            //
            //    Input, double U[M], V[N], the R8SM vectors.
            //
            //    Input, double X[M], the vector to be multiplied.
            //
            //    Output, double R8SM_MTV[N], the product (A-u*v')' * X.
            //
        {
            double[] b;
            double dot;
            int i;
            int j;

            b = r8vec_zeros_new(n);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < m; j++)
                {
                    b[i] = b[i] + x[j] * a[j + i * m];
                }

                dot = 0.0;
                for (j = 0; j < m; j++)
                {
                    dot = dot + u[j] * x[j];
                }

                b[i] = b[i] - v[i] * dot;
            }

            return b;
        }

        public static double[] r8sm_mv(int m, int n, double[] a, double[] u, double[] v, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8SM_MV multiplies an R8SM matrix times a vector.
            //
            //  Discussion:
            //
            //    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
            //    which is defined by an M by N matrix A, an M vector U, and
            //    an N vector V, by B = A - U * V'
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 January 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns of the matrix.
            //
            //    Input, double A[M*N], the R8SM matrix.
            //
            //    Input, double U[M], V[N], the R8SM vectors U and V.
            //
            //    Input, double X[N], the vector to be multiplied by (A-u*v').
            //
            //    Output, double R8SM_MV[M], the product (A-u*v') * x.
            //
        {
            double[] b;
            int i;
            int j;
            double vx;

            b = r8vec_zeros_new(m);

            vx = 0.0;
            for (j = 0; j < n; j++)
            {
                vx = vx + v[j] * x[j];
            }

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    b[i] = b[i] + a[i + j * m] * x[j];
                }

                b[i] = b[i] - u[i] * vx;
            }

            return b;
        }

        public static void r8sm_print(int m, int n, double[] a, double[] u, double[] v, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8SM_PRINT prints an R8SM matrix.
            //
            //  Discussion:
            //
            //    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
            //    which is defined by an M by N matrix A, an M vector U, and
            //    an N vector V, by B = A - U * V'
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
            //    Input, int M, N, the number of rows and columns of the matrix.
            //
            //    Input, double A[M*N], the R8SM matrix.
            //
            //    Input, double U[M], V[N], the R8SM vectors.
            //
            //    Input, string TITLE, a title.
            //
        {
            r8sm_print_some(m, n, a, u, v, 1, 1, m, n, title);

        }

        public static void r8sm_print_some(int m, int n, double[] a, double[] u, double[] v, int ilo,
                int jlo, int ihi, int jhi, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8SM_PRINT_SOME prints some of an R8SM matrix.
            //
            //  Discussion:
            //
            //    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
            //    which is defined by an M by N matrix A, an M vector U, and
            //    an N vector V, by B = A - U * V'
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns of the matrix.
            //
            //    Input, double A[M*N], the R8SM matrix.
            //
            //    Input, double U[M], V[N], the R8SM vectors.
            //
            //    Input, int ILO, JLO, IHI, JHI, designate the first row and
            //    column, and the last row and column to be printed.
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
            for (j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX)
            {
                j2hi = j2lo + INCX - 1;
                j2hi = Math.Min(j2hi, n);
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
                i2lo = Math.Max(ilo, 1);
                i2hi = Math.Min(ihi, m);

                for (i = i2lo; i <= i2hi; i++)
                {
                    cout = i.ToString().PadLeft(4) + "  ";
                    //
                    //  Print out (up to) 5 entries in row I, that lie in the current strip.
                    //
                    for (j = j2lo; j <= j2hi; j++)
                    {
                        cout += (a[i - 1 + (j - 1) * m] - u[i - 1] * v[j - 1]).ToString().PadLeft(12) + "  ";
                    }

                    Console.WriteLine(cout);
                }
            }
        }

        public static void r8sm_random(int m, int n, ref int seed, ref double[] a, ref double[] u, ref double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8SM_RANDOM randomizes an R8SM matrix.
            //
            //  Discussion:
            //
            //    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
            //    which is defined by an M by N matrix A, an M vector U, and
            //    an N vector V, by B = A - U * V'
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 May 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns of the matrix.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double A[M*N], the R8SM matrix.
            //
            //    Output, double U[M], V[N], the R8SM vectors.
            //
        {
            int i;
            int j;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    a[i + j * m] = UniformRNG.r8_uniform_01(ref seed);
                }
            }

            for (i = 0; i < m; i++)
            {
                u[i] = UniformRNG.r8_uniform_01(ref seed);
            }

            for (j = 0; j < n; j++)
            {
                v[j] = UniformRNG.r8_uniform_01(ref seed);
            }

        }

        public static double[] r8sm_sl(int n, double[] a_lu, double[] u, double[] v, double[] b,
                int[] pivot, int job)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8SM_SL solves a square R8SM system that has been factored.
            //
            //  Discussion:
            //
            //    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
            //    which is defined by an M by N matrix A, an M vector U, and
            //    an N vector V, by B = A - U * V'
            //
            //    It is assumed that A has been decomposed into its LU factors
            //    by R8GE_FA.  The Sherman Morrison formula allows
            //    us to solve linear systems involving (A-u*v') by solving linear
            //    systems involving A and adjusting the results.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 October 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    David Kahaner, Cleve Moler, Stephen Nash
            //    Numerical Methods and Software,
            //    Prentice Hall, 1989
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //    N must be positive.
            //
            //    Input, double A_LU[N*N], the LU factors from R8GE_FA.
            //
            //    Input, double U[N], V[N], the R8SM vectors U and V.
            //
            //    Input, double B[N], the right hand side vector.
            //
            //    Input, int PIVOT[N], the pivot vector produced by R8GE_FA.
            //
            //    Input, int JOB, specifies the system to solve.
            //    0, solve (A-u*v') * X = B.
            //    nonzero, solve (A-u*v') * X = B.
            //
            //    Output, double R8SM_SL[N], the solution vector, or NULL if
            //    an error occurred.
            //
        {
            double alpha;
            double beta;
            int i;
            int job_local;
            double[] w;
            double[] x;

            x = new double[n];

            if (job == 0)
            {
                //
                //  Solve A' * w = v.
                //
                job_local = 1;
                w = r8ge_sl_new(n, a_lu, pivot, v, job_local);
                //
                //  Set beta = w' * b.
                //
                beta = 0.0;
                for (i = 0; i < n; i++)
                {
                    beta = beta + w[i] * b[i];
                }

                //
                //  Solve A * x = b.
                //
                job_local = 0;
                x = r8ge_sl_new(n, a_lu, pivot, b, job_local);
                //
                //  Solve A * w = u.
                //
                job_local = 0;
                w = r8ge_sl_new(n, a_lu, pivot, u, job_local);
                //
                //  Set alpha = 1 / ( 1 - v' * w ).
                //
                alpha = 1.0;
                for (i = 0; i < n; i++)
                {
                    alpha = alpha - v[i] * w[i];
                }
            }
            else
            {
                //
                //  Solve A * w = u.
                //
                job_local = 0;
                w = r8ge_sl_new(n, a_lu, pivot, u, job_local);
                //
                //  Set beta = w' * b.
                //
                beta = 0.0;
                for (i = 0; i < n; i++)
                {
                    beta = beta + w[i] * b[i];
                }

                //
                //  Solve A' * x = b.
                //
                job_local = 1;
                x = r8ge_sl_new(n, a_lu, pivot, b, job_local);
                //
                //  Solve A' * w = v.
                //
                job_local = 1;
                w = r8ge_sl_new(n, a_lu, pivot, v, job_local);
                //
                //  Set alpha = 1 / ( 1 - u' * w ).
                //
                alpha = 1.0;
                for (i = 0; i < n; i++)
                {
                    alpha = alpha - u[i] * w[i];
                }
            }

            if (alpha == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8SM_SL - Fatal error!");
                Console.WriteLine("  The divisor ALPHA is zero.");
                return null;
            }

            alpha = 1.0 / alpha;
            //
            //  Set b = b + alpha * beta * w.
            //
            for (i = 0; i < n; i++)
            {
                x[i] = x[i] + alpha * beta * w[i];
            }

            return x;
        }

        public static double[] r8sm_to_r8ge(int m, int n, double[] a, double[] u, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8SM_TO_R8GE copies an R8SM matrix to an R8GE matrix.
            //
            //  Discussion:
            //
            //    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
            //    which is defined by an M by N matrix A, an M vector U, and
            //    an N vector V, by B = A - U * V'
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 May 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns of the matrix.
            //
            //    Input, double A[M*N], the R8SM matrix.
            //
            //    Input, double U[M], V[N], the R8SM vectors.
            //
            //    Output, double R8SM_TO_R8GE[M*N], the R8GE matrix.
            //
        {
            double[] b;
            int i;
            int j;

            b = new double[m * n];

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    b[i + j * m] = a[i + j * m] - u[i] * v[j];
                }
            }

            return b;
        }

        public static void r8sm_zeros(int m, int n, ref double[] a, ref double[] u, ref double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8SM_ZEROS zeros an R8SM matrix.
            //
            //  Discussion:
            //
            //    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
            //    which is defined by an M by N matrix A, an M vector U, and
            //    an N vector V, by B = A - U * V'
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 January 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns of the matrix.
            //
            //    Output, double A[M*N], the R8SM matrix.
            //
            //    Output, double U[M], V[N], the R8SM vectors.
            //
        {
            int i;
            int j;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    a[i + j * m] = 0.0;
                }
            }

            for (i = 0; i < m; i++)
            {
                u[i] = 0.0;
            }

            for (i = 0; i < n; i++)
            {
                v[i] = 0.0;
            }

        }
    }
}