using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8mat_identity(int n, ref double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_IDENTITY sets the square matrix A to the identity.
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
            //    01 December 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of A.
            //
            //    Output, double A[N*N], the N by N identity matrix.
            //
        {
            int i;
            int j;
            int k;

            k = 0;
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    if (i == j)
                    {
                        a[k] = 1.0;
                    }
                    else
                    {
                        a[k] = 0.0;
                    }

                    k = k + 1;
                }
            }

            return;
        }

        public static double[] r8mat_identity_new(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_IDENTITY_NEW returns an identity matrix.
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
            //    06 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of A.
            //
            //    Output, double R8MAT_IDENTITY_NEW[N*N], the N by N identity matrix.
            //
        {
            double[] a;
            int i;
            int j;
            int k;

            a = new double[n * n];

            k = 0;
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    if (i == j)
                    {
                        a[k] = 1.0;
                    }
                    else
                    {
                        a[k] = 0.0;
                    }

                    k = k + 1;
                }
            }

            return a;
        }

        public static double r8mat_is_identity(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_IS_IDENTITY determines if an R8MAT is the identity.
            //
            //  Discussion:
            //
            //    An R8MAT is a matrix of real ( kind = 8 ) values.
            //
            //    The routine returns the Frobenius norm of A - I.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 July 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, double A[N*N], the matrix.
            //
            //    Output, double R8MAT_IS_IDENTITY, the Frobenius norm
            //    of the difference matrix A - I, which would be exactly zero
            //    if A were the identity matrix.
            //
        {
            double error_frobenius;
            int i;
            int j;
            double t;

            error_frobenius = 0.0;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if (i == j)
                    {
                        t = a[i + j * n] - 1.0;
                    }
                    else
                    {
                        t = a[i + j * n];
                    }

                    error_frobenius = error_frobenius + t * t;
                }
            }

            error_frobenius = Math.Sqrt(error_frobenius);

            return error_frobenius;
        }

    }
}