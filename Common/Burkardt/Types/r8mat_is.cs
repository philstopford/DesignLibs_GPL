using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static bool r8mat_is_binary(int m, int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_IS_BINARY is true if the entries in an R8MAT are all 0 or 1.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 April 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the dimensions of the array.
            //
            //    Input, double X[M*N], the array to be checked.
            //
            //    Output, bool R8MAT_IS_BINARY is true if are entries are 0 or 1.
            //
        {
            int i;
            int j;
            bool value;

            value = true;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    if (x[i + j * m] != 0.0 && x[i + j * m] != 1.0)
                    {
                        value = false;
                        break;
                    }
                }
            }

            return value;
        }


        public static bool r8mat_is_in_01(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_IS_IN_01 is TRUE if the entries of an R8MAT are in the range [0,1].
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 October 2004
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
            //    Input, double A[M*N], the matrix.
            //
            //    Output, bool R8MAT_IS_IN_01, is TRUE if every entry of A is
            //    between 0 and 1.
            //
        {
            int i;
            int j;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    if (a[i + j * m] < 0.0 || 1.0 < a[i + j * m])
                    {
                        return false;
                    }
                }
            }

            return true;
        }

        public static bool r8mat_is_insignificant(int m, int n, double[] r, double[] s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_IS_INSIGNIFICANT determines if an R8MAT is relatively insignificant.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the dimension of the matrices.
            //
            //    Input, double R[M*N], the vector to be compared against.
            //
            //    Input, double S[M*N], the vector to be compared.
            //
            //    Output, bool R8MAT_IS_INSIGNIFICANT, is TRUE if S is insignificant
            //    compared to R.
            //
        {
            int i;
            int j;
            double t;
            double tol;
            bool value;

            value = true;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    t = r[i + j * m] + s[i + j * m];
                    tol = double.Epsilon * Math.Abs(r[i + j * m]);

                    if (tol < Math.Abs(r[i + j * m] - t))
                    {
                        value = false;
                        break;
                    }
                }
            }

            return value;
        }

        public static bool r8mat_is_integer(int m, int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_IS_INTEGER is true if an R8MAT only contains integer entries.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 August 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the dimensions of the array.
            //
            //    Input, double X[M*N], the vector to be checked.
            //
            //    Output, bool R8MAT_IS_INTEGER is true if all elements of X
            //    are integers.
            //
        {
            int i;
            int j;
            bool value;

            value = true;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    if (!r8_is_integer(x[i + j * m]))
                    {
                        value = false;
                        break;
                    }
                }
            }

            return value;
        }

        public static bool r8mat_is_significant(int m, int n, double[] r, double[] s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_IS_SIGNIFICANT determines if an R8MAT is relatively significant.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the dimension of the matrices.
            //
            //    Input, double R[M*N], the vector to be compared against.
            //
            //    Input, double S[M*N], the vector to be compared.
            //
            //    Output, bool R8MAT_IS_SIGNIFICANT, is TRUE if S is significant
            //    compared to R.
            //
        {
            int i;
            int j;
            double t;
            double tol;
            bool value;

            value = false;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    t = r[i + j * m] + s[i + j * m];
                    tol = double.Epsilon * Math.Abs(r[i + j * m]);

                    if (tol < Math.Abs(r[i + j * m] - t))
                    {
                        value = true;
                        break;
                    }
                }
            }

            return value;
        }

        public static double r8mat_is_symmetric(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_IS_SYMMETRIC checks an R8MAT for symmetry.
            //
            //  Discussion:
            //
            //    An R8MAT is a matrix of double precision real values.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 July 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the order of the matrix.
            //
            //    Input, double A[M*N], the matrix.
            //
            //    Output, double RMAT_IS_SYMMETRIC, measures the 
            //    Frobenius norm of ( A - A' ), which would be zero if the matrix
            //    were exactly symmetric.
            //
        {
            int i;
            int j;
            const double r8_huge = 1.79769313486231571E+308;
            double value;

            if (m != n)
            {
                value = r8_huge;
                return value;
            }

            value = 0.0;
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    value = value + Math.Pow(a[i + j * m] - a[j + i * m], 2);
                }
            }

            value = Math.Sqrt(value);

            return value;
        }
    }
}