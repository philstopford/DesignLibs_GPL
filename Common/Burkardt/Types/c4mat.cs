using System;
using System.Numerics;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static Complex[] c4mat_test(int n)

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_TEST returns a test matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, complex <float> C4MAT_TEST[N*N], the matrix.
//
        {
            Complex I = new Complex(0.0, 1.0);

            Complex[] a = new Complex[n * n];

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    float angle = (float)(2.0 * Math.PI * (float) (i * j) / (float) (n));
                    
                    a[i + j * n] = Complex.Exp(I * angle) / Math.Sqrt((float) (n));
                }
            }

            return a;
        }

        public static Complex[] c4mat_test_inverse(int n)

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_TEST_INVERSE returns the inverse of a test matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, complex <float> C4MAT_TEST_INVERSE[N*N], the matrix.
//
        {
            Complex[] a;
            int i;
            int j;
            Complex t;

            a = c4mat_test(n);

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < j; i++)
                {
                    t = Complex.Conjugate(a[i + j * n]);
                    a[i + j * n] = Complex.Conjugate(a[j + i * n]);
                    a[j + i * n] = t;
                }
            }

            return a;
        }
    }
}