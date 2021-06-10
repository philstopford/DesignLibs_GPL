using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8mat_podet(int n, double[] r)

            /******************************************************************************/
            /*
            Purpose:
            
            R8MAT_PODET computes the determinant of a factored positive definite matrix.
            
            Discussion:
            
            This routine expects to receive R, the upper triangular factor of A,
            computed by R8MAT_POFAC, with the property that A = R' * R.
            
            Licensing:
            
            This code is distributed under the GNU LGPL license.
            
            Modified:
            
            09 June 2013
            
            Author:
            
            C version by John Burkardt.
            
            Reference:
            
            Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
            LINPACK User's Guide,
            SIAM, 1979,
            ISBN13: 978-0-898711-72-1,
            LC: QA214.L56.
            
            Parameters:
            
            Input, int N, the order of the matrix.
            
            Input, double R[N*N], the Cholesky factor of A.
            
            Output, double R8MAT_PODET, the determinant of A.
            */
        {
            double det = 1.0;
            for (int i = 0; i < n; i++)
            {
                det = det * r[i + i * n] * r[i + i * n];
            }

            return det;
        }

        public static double[] r8mat_pofac(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_POFAC factors a real symmetric positive definite matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 June 2013
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, 
            //    Pete Stewart,
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
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
            //    Input, double A[N*N], the matrix to be factored.
            //
            //    Output, double R8MAT_POFAC[N*N], an upper triangular matrix such that
            //    A = R'*R.
            //
        {
            double[] r = new double[n * n];

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i <= j; i++)
                {
                    r[i + j * n] = a[i + j * n];
                }

                for (int i = j + 1; i < n; i++)
                {
                    r[i + j * n] = 0.0;
                }
            }

            for (int j = 0; j < n; j++)
            {
                double s = 0.0;

                for (int k = 0; k < j; k++)
                {
                    double dot = 0.0;
                    for (int i = 0; i < k; i++)
                    {
                        dot = dot + r[i + k * n] * r[i + j * n];
                    }

                    double t = r[k + j * n] - dot;
                    t = t / r[k + k * n];
                    r[k + j * n] = t;
                    s = s + t * t;
                }

                s = r[j + j * n] - s;

                if (s < 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8MAT_POFAC - Fatal error!");
                    Console.WriteLine("  The matrix is not positive definite.");
                    return new double[1];
                }

                if (s == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8MAT_POFAC - Warning!");
                    Console.WriteLine("  The matrix is not strictly positive definite.");
                }

                r[j + j * n] = Math.Sqrt(s);
            }

            return r;
        }

        public static double[] r8mat_poinv(int n, double[] r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_POINV computes the inverse of a factored positive definite matrix.
            //
            //  Discussion:
            //
            //    This routine expects to receive R, the upper triangular factor of A,
            //    computed by R8MAT_POFAC, with the property that A = R' * R.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 June 2013
            //
            //  Author:
            //
            //    C version by John Burkardt.
            //
            //  Reference:
            //
            //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
            //    LINPACK User's Guide,
            //    SIAM, 1979,
            //    ISBN13: 978-0-898711-72-1,
            //    LC: QA214.L56.
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix A.
            //
            //    Input, double R[N*N], the Cholesky factor of A.
            //
            //    Input, double R8MAT_POINV[N*N], the inverse of A.
            //
        {
            double t;

            double[] b = new double[n * n];

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    b[i + j * n] = r[i + j * n];
                }
            }

            for (int k = 0; k < n; k++)
            {
                b[k + k * n] = 1.0 / b[k + k * n];
                for (int i = 0; i < k; i++)
                {
                    b[i + k * n] = -b[i + k * n] * b[k + k * n];
                }

                for (int j = k + 1; j < n; j++)
                {
                    t = b[k + j * n];
                    b[k + j * n] = 0.0;
                    for (int i = 0; i <= k; i++)
                    {
                        b[i + j * n] = b[i + j * n] + t * b[i + k * n];
                    }
                }
            }

            /*
            Form inverse(R) * (inverse(R))'.
            */
            for (int j = 0; j < n; j++)
            {
                for (int k = 0; k < j; k++)
                {
                    t = b[k + j * n];
                    for (int i = 0; i <= k; i++)
                    {
                        b[i + k * n] = b[i + k * n] + t * b[i + j * n];
                    }
                }

                t = b[j + j * n];
                for (int i = 0; i <= j; i++)
                {
                    b[i + j * n] = b[i + j * n] * t;
                }
            }

            return b;
        }
        
    }
}