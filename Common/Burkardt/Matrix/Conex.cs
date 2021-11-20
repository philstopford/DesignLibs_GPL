using System;
using Burkardt.Types;

namespace Burkardt.MatrixNS;

public static partial class Matrix
{
    public static double[] conex1(double alpha)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONEX1 returns the CONEX1 matrix.
        //
        //  Discussion:
        //
        //    The CONEX1 matrix is a counterexample to the LINPACK condition
        //    number estimator RCOND available in the LINPACK routine DGECO.
        //
        //  Formula:
        //
        //    1  -1 -2*ALPHA   0
        //    0   1    ALPHA    -ALPHA
        //    0   1  1+ALPHA  -1-ALPHA
        //    0   0  0           ALPHA
        //
        //  Example:
        //
        //    ALPHA = 100
        //
        //    1  -1  -200     0
        //    0   1   100  -100
        //    0   1   101  -101
        //    0   0     0   100
        //
        //  Properties:
        //
        //    A is generally not symmetric: A' /= A.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Alan Cline, RK Rew,
        //    A set of counterexamples to three condition number estimators,
        //    SIAM Journal on Scientific and Statistical Computing,
        //    Volume 4, 1983, pages 602-611.
        //
        //  Parameters:
        //
        //    Input, double ALPHA, the scalar defining A.  
        //    A common value is 100.0.
        //
        //    Output, double CONEX1[4*4], the matrix.
        //
    {
        const int n = 4;

        double[] a = new double[n * n];

        a[0] = 1.0;
        a[1] = 0.0;
        a[2] = 0.0;
        a[3] = 0.0;

        a[0 + 1 * n] = -1.0;
        a[1 + 1 * n] = 1.0;
        a[2 + 1 * n] = 1.0;
        a[3 + 1 * n] = 0.0;

        a[0 + 2 * n] = -2.0 * alpha;
        a[1 + 2 * n] = alpha;
        a[2 + 2 * n] = 1.0 + alpha;
        a[3 + 2 * n] = 0.0;

        a[0 + 3 * n] = 0.0;
        a[1 + 3 * n] = -alpha;
        a[2 + 3 * n] = -1.0 - alpha;
        a[3 + 3 * n] = alpha;

        return a;
    }

    public static double[] conex1_inverse(double alpha)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONEX1_INVERSE returns the inverse of the CONEX1 matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double ALPHA, the scalar defining A.  
        //
        //    Output, double CONEX1_INVERSE[4*4], the matrix.
        //
    {
        const int n = 4;

        double[] a = new double[n * n];

        a[0] = 1.0;
        a[1] = 0.0;
        a[2] = 0.0;
        a[3] = 0.0;

        a[0 + 1 * n] = 1.0 - alpha;
        a[1 + 1 * n] = 1.0 + alpha;
        a[2 + 1 * n] = -1.0;
        a[3 + 1 * n] = 0.0;

        a[0 + 2 * n] = alpha;
        a[1 + 2 * n] = -alpha;
        a[2 + 2 * n] = 1.0;
        a[3 + 2 * n] = 0.0;

        a[0 + 3 * n] = 2.0;
        a[1 + 3 * n] = 0.0;
        a[2 + 3 * n] = 1.0 / alpha;
        a[3 + 3 * n] = 1.0 / alpha;

        return a;
    }

    public static double[] conex2(double alpha)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONEX2 returns the CONEX2 matrix.
        //
        //  Formula:
        //
        //    1   1-1/ALPHA^2  -2
        //    0   1/ALPHA      -1/ALPHA
        //    0   0             1
        //
        //  Example:
        //
        //    ALPHA = 100
        //
        //    1  0.9999  -2
        //    0  0.01    -0.01
        //    0  0        1
        //
        //  Properties:
        //
        //    A is generally not symmetric: A' /= A.
        //
        //    A is upper triangular.
        //
        //    det ( A ) = 1 / ALPHA.
        //
        //    LAMBDA = ( 1, 1/ALPHA, 1 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Alan Cline, RK Rew,
        //    A set of counterexamples to three condition number estimators,
        //    SIAM Journal on Scientific and Statistical Computing,
        //    Volume 4, 1983, pages 602-611.
        //
        //  Parameters:
        //
        //    Input, double ALPHA, the scalar defining A.  
        //    A common value is 100.0.  ALPHA must not be zero.
        //
        //    Output, double CONEX2[3*3], the matrix.
        //
    {
        const int n = 3;

        double[] a = new double[n * n];

        switch (alpha)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("CONEX2 - Fatal error!");
                Console.WriteLine("  The input value of ALPHA was zero.");
                return null;
        }

        a[0] = 1.0;
        a[1] = 0.0;
        a[2] = 0.0;

        a[0 + 1 * n] = (alpha - 1.0) * (alpha + 1.0) / alpha / alpha;
        a[1 + 1 * n] = 1.0 / alpha;
        a[2 + 1 * n] = 0.0;

        a[0 + 2 * n] = -2.0;
        a[1 + 2 * n] = -1.0 / alpha;
        a[2 + 2 * n] = 1.0;

        return a;
    }

    public static double[] conex2_inverse(double alpha)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONEX2_INVERSE returns the inverse of the CONEX2 matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double ALPHA, the scalar defining A.  
        //    A common value is 100.0.  ALPHA must not be zero.
        //
        //    Output, double CONEX2_INVERSE[3*3], the matrix.
        //
    {
        const int n = 3;

        double[] a = new double[n * n];

        switch (alpha)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("CONEX2_INVERSE - Fatal error!");
                Console.WriteLine("  The input value of ALPHA was zero.");
                return null;
        }

        a[0] = 1.0;
        a[1] = 0.0;
        a[2] = 0.0;

        a[0 + 1 * n] = (1.0 - alpha) * (1.0 + alpha) / alpha;
        a[1 + 1 * n] = alpha;
        a[2 + 1 * n] = 0.0;

        a[0 + 2 * n] = (1.0 + alpha * alpha) / alpha / alpha;
        a[1 + 2 * n] = 1.0;
        a[2 + 2 * n] = 1.0;

        return a;
    }

    public static double[] conex3(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONEX3 returns the CONEX3 matrix.
        //
        //  Formula:
        //
        //    if ( I = J and I < N )
        //      A(I,J) =  1.0 for 1<=I<N
        //    else if ( I = J = N )
        //      A(I,J) = -1.0
        //    else if ( J < I )
        //      A(I,J) = -1.0
        //    else
        //      A(I,J) =  0.0
        //
        //  Example:
        //
        //    N = 5
        //
        //     1  0  0  0  0
        //    -1  1  0  0  0
        //    -1 -1  1  0  0
        //    -1 -1 -1  1  0
        //    -1 -1 -1 -1 -1
        //
        //  Properties:
        //
        //    A is generally not symmetric: A' /= A.
        //
        //    A is integral, therefore det ( A ) is integral, and 
        //    det ( A ) * inverse ( A ) is integral.
        //
        //    A is lower triangular.
        //
        //    det ( A ) = -1.
        //
        //    A is unimodular.
        //
        //    LAMBDA = ( 1, 1, 1, 1, -1 )
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
        //  Reference:
        //
        //    Alan Cline, RK Rew,
        //    A set of counterexamples to three condition number estimators,
        //    SIAM Journal on Scientific and Statistical Computing,
        //    Volume 4, 1983, pages 602-611.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double CONEX3[N*N], the matrix.
        //
    {
        int j;

        double[] a = new double[n * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                if (j < i)
                {
                    a[i + j * n] = -1.0;
                }
                else if (j == i && i != n - 1)
                {
                    a[i + j * n] = 1.0;
                }
                else if (j == i && i == n - 1)
                {
                    a[i + j * n] = -1.0;
                }
                else
                {
                    a[i + j * n] = 0.0;
                }
            }
        }

        return a;
    }

    public static double[] conex3_inverse(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONEX3_INVERSE returns the inverse of the CONEX3 matrix.
        //
        //  Example:
        //
        //    N = 5
        //
        //     1  0  0  0  0
        //     1  1  0  0  0
        //     2  1  1  0  0
        //     4  2  1  1  0
        //    -8 -4 -2 -1 -1
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
        //  Reference:
        //
        //    Alan Cline, RK Rew,
        //    A set of counterexamples to three condition number estimators,
        //    SIAM Journal on Scientific and Statistical Computing,
        //    Volume 4, 1983, pages 602-611.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double CONEX3_INVERSE[N*N], the matrix.
        //
    {
        int j;

        double[] a = new double[n * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                if (i < n - 1)
                {
                    if (j < i)
                    {
                        a[i + j * n] = Math.Pow(2.0, i - j - 1);
                    }
                    else if (i == j)
                    {
                        a[i + j * n] = 1.0;
                    }
                    else
                    {
                        a[i + j * n] = 0.0;
                    }
                }
                else if (i == n - 1)
                {
                    if (j < i)
                    {
                        a[i + j * n] = -Math.Pow(2.0, i - j - 1);
                    }
                    else
                    {
                        a[i + j * n] = -1.0;
                    }
                }
            }
        }

        return a;
    }

    public static double[] conex4()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONEX4 returns the CONEX4 matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double CONEX4[4*4], the matrix.
        //
    {
        double[] a_save = {
                7.0, 6.0, 5.0, 5.0,
                10.0, 8.0, 7.0, 7.0,
                8.0, 10.0, 9.0, 6.0,
                7.0, 9.0, 10.0, 5.0
            }
            ;

        double[] a = typeMethods.r8mat_copy_new(4, 4, a_save);

        return a;
    }

    public static double[] conex4_inverse()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONEX4_INVERSE returns the inverse of the CONEX4 matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double CONEX4_INVERSE[4*4], the matrix.
        //
    {
        double[] a_save = {
                -41.0, 25.0, 10.0, -6.0,
                -17.0, 10.0, 5.0, -3.0,
                10.0, -6.0, -3.0, 2.0,
                68.0, -41.0, -17.0, 10.0
            }
            ;

        double[] a = typeMethods.r8mat_copy_new(4, 4, a_save);

        return a;
    }
}