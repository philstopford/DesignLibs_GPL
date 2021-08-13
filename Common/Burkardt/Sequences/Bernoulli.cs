using System;
using Burkardt.Types;

namespace Burkardt.Sequence
{
    public static class Bernoulli
    {
        public static void bernoulli_number(int n, ref double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BERNOULLI_NUMBER computes the value of the Bernoulli numbers B(0) through B(N).
            //
            //  Discussion:
            //
            //    The Bernoulli numbers are rational.
            //
            //    If we define the sum of the M-th powers of the first N integers as:
            //
            //      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
            //
            //    and let C(I,J) be the combinatorial coefficient:
            //
            //      C(I,J) = I! / ( ( I - J )! * J! )
            //
            //    then the Bernoulli numbers B(J) satisfy:
            //
            //      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)^(M+1-J)
            //
            //  First values:
            //
            //   B0  1                   =         1.00000000000
            //   B1 -1/2                 =        -0.50000000000
            //   B2  1/6                 =         1.66666666666
            //   B3  0                   =         0
            //   B4 -1/30                =        -0.03333333333
            //   B5  0                   =         0
            //   B6  1/42                =         0.02380952380
            //   B7  0                   =         0
            //   B8 -1/30                =        -0.03333333333
            //   B9  0                   =         0
            //  B10  5/66                =         0.07575757575
            //  B11  0                   =         0
            //  B12 -691/2730            =        -0.25311355311
            //  B13  0                   =         0
            //  B14  7/6                 =         1.16666666666
            //  B15  0                   =         0
            //  B16 -3617/510            =        -7.09215686274
            //  B17  0                   =         0
            //  B18  43867/798           =        54.97117794486
            //  B19  0                   =         0
            //  B20 -174611/330          =      -529.12424242424
            //  B21  0                   =         0
            //  B22  854,513/138         =      6192.123
            //  B23  0                   =         0
            //  B24 -236364091/2730      =    -86580.257
            //  B25  0                   =         0
            //  B26  8553103/6           =   1425517.16666
            //  B27  0                   =         0
            //  B28 -23749461029/870     = -27298231.0678
            //  B29  0                   =         0
            //  B30  8615841276005/14322 = 601580873.901
            //
            //  Recursion:
            //
            //    With C(N+1,K) denoting the standard binomial coefficient,
            //
            //    B(0) = 1.0
            //    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
            //
            //  Warning:
            //
            //    This recursion, which is used in this routine, rapidly results
            //    in significant errors.
            //
            //  Special Values:
            //
            //    Except for B(1), all Bernoulli numbers of odd index are 0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the highest Bernoulli number to compute.
            //
            //    Output, double B[N+1], B(I) contains the I-th Bernoulli number.
            //
        {
            double b_sum;
            int[] c;
            int i;
            int j;

            if (n < 0)
            {
                return;
            }

            b[0] = 1.0;

            if (n < 1)
            {
                return;
            }

            b[1] = -0.5;

            c = new int[n + 2];
            c[0] = 1;
            c[1] = 2;
            c[2] = 1;

            for (i = 2; i <= n; i++)
            {
                Comb.comb_row_next(i + 1, ref c);
                if ((i % 2) == 1)
                {
                    b[i] = 0.0;
                }
                else
                {
                    b_sum = 0.0;
                    for (j = 0; j <= i - 1; j++)
                    {
                        b_sum = b_sum + b[j] * (double)(c[j]);
                    }

                    b[i] = -b_sum / (double)(c[i]);
                }
            }
        }

        public static void bernoulli_number2(int n, ref double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BERNOULLI_NUMBER2 evaluates the Bernoulli numbers.
            //
            //  Discussion:
            //
            //    The Bernoulli numbers are rational.
            //
            //    If we define the sum of the M-th powers of the first N integers as:
            //
            //      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
            //
            //    and let C(I,J) be the combinatorial coefficient:
            //
            //      C(I,J) = I! / ( ( I - J )! * J! )
            //
            //    then the Bernoulli numbers B(J) satisfy:
            //
            //      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)^(M+1-J)
            //
            //    Note that the Bernoulli numbers grow rapidly.  Bernoulli number
            //    62 is probably the last that can be computed on the VAX without
            //    overflow.
            //
            //    A different method than that used in BERN is employed.
            //
            //  First values:
            //
            //   B0  1                   =         1.00000000000
            //   B1 -1/2                 =        -0.50000000000
            //   B2  1/6                 =         1.66666666666
            //   B3  0                   =         0
            //   B4 -1/30                =        -0.03333333333
            //   B5  0                   =         0
            //   B6  1/42                =         0.02380952380
            //   B7  0                   =         0
            //   B8 -1/30                =        -0.03333333333
            //   B9  0                   =         0
            //  B10  5/66                =         0.07575757575
            //  B11  0                   =         0
            //  B12 -691/2730            =        -0.25311355311
            //  B13  0                   =         0
            //  B14  7/6                 =         1.16666666666
            //  B15  0                   =         0
            //  B16 -3617/510            =        -7.09215686274
            //  B17  0                   =         0
            //  B18  43867/798           =        54.97117794486
            //  B19  0                   =         0
            //  B20 -174611/330          =      -529.12424242424
            //  B21  0                   =         0
            //  B22  854,513/138         =      6192.123
            //  B23  0                   =         0
            //  B24 -236364091/2730      =    -86580.257
            //  B25  0                   =         0
            //  B26  8553103/6           =   1425517.16666
            //  B27  0                   =         0
            //  B28 -23749461029/870     = -27298231.0678
            //  B29  0                   =         0
            //  B30  8615841276005/14322 = 601580873.901
            //
            //  Recursion:
            //
            //    With C(N+1,K) denoting the standard binomial coefficient,
            //
            //    B(0) = 1.0
            //    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
            //
            //  Special Values:
            //
            //    Except for B(1), all Bernoulli numbers of odd index are 0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the highest order Bernoulli number to compute.
            //
            //    Output, double B[N+1], the requested Bernoulli numbers.
            //
        {
            double altpi;
            int i;
            int k;
            int kmax = 400;
            const double r8_pi = 3.141592653589793;
            double sgn;
            double sum2;
            double t;
            double term;
            double tol = 1.0E-06;

            if (n < 0)
            {
                return;
            }

            b[0] = 1.0;

            if (n < 1)
            {
                return;
            }

            b[1] = -0.5;

            if (n < 2)
            {
                return;
            }

            altpi = Math.Log(2.0 * r8_pi);
            //
            //  Initial estimates for B(I), I = 2 to N
            //
            b[2] = Math.Log(2.0);

            for (i = 3; i <= n; i++)
            {
                if ((i % 2) == 1)
                {
                    b[i] = 0.0;
                }
                else
                {
                    b[i] = Math.Log((double)(i * (i - 1))) + b[i - 2];
                }
            }

            b[2] = 1.0 / 6.0;

            if (n <= 3)
            {
                return;
            }

            b[4] = -1.0 / 30.0;

            sgn = -1.0;

            for (i = 6; i <= n; i = i + 2)
            {
                sgn = -sgn;
                t = 2.0 * sgn * Math.Exp(b[i] - (double)(i) * altpi);

                sum2 = 1.0;

                for (k = 2; k <= kmax; k++)
                {
                    term = Math.Pow((double)k, (double)-i);

                    sum2 = sum2 + term;

                    if (term <= tol * sum2)
                    {
                        break;
                    }

                }

                b[i] = t * sum2;

            }

        }

        public static double bernoulli_number3(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BERNOULLI_NUMBER3 computes the value of the Bernoulli number B(N).
            //
            //  Discussion:
            //
            //    The Bernoulli numbers are rational.
            //
            //    If we define the sum of the M-th powers of the first N integers as:
            //
            //      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
            //
            //    and let C(I,J) be the combinatorial coefficient:
            //
            //      C(I,J) = I! / ( ( I - J )! * J! )
            //
            //    then the Bernoulli numbers B(J) satisfy:
            //
            //      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) 
            //        C(M+1,J) B(J) * (N+1)^(M+1-J)
            //
            //  First values:
            //
            //     B0  1                   =         1.00000000000
            //     B1 -1/2                 =        -0.50000000000
            //     B2  1/6                 =         1.66666666666
            //     B3  0                   =         0
            //     B4 -1/30                =        -0.03333333333
            //     B5  0                   =         0
            //     B6  1/42                =         0.02380952380
            //     B7  0                   =         0
            //     B8 -1/30                =        -0.03333333333
            //     B9  0                   =         0
            //    B10  5/66                =         0.07575757575
            //    B11  0                   =         0
            //    B12 -691/2730            =        -0.25311355311
            //    B13  0                   =         0
            //    B14  7/6                 =         1.16666666666
            //    B15  0                   =         0
            //    B16 -3617/510            =        -7.09215686274
            //    B17  0                   =         0
            //    B18  43867/798           =        54.97117794486
            //    B19  0                   =         0
            //    B20 -174611/330          =      -529.12424242424
            //    B21  0                   =         0
            //    B22  854513/138          =      6192.123
            //    B23  0                   =         0
            //    B24 -236364091/2730      =    -86580.257
            //    B25  0                   =         0
            //    B26  8553103/6           =   1425517.16666
            //    B27  0                   =         0
            //    B28 -23749461029/870     = -27298231.0678
            //    B29  0                   =         0
            //    B30  8615841276005/14322 = 601580873.901
            //
            //  Recursion:
            //
            //    With C(N+1,K) denoting the standard binomial coefficient,
            //
            //    B(0) = 1.0
            //    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
            //
            //  Special Values:
            //
            //    Except for B(1), all Bernoulli numbers of odd index are 0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 February 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the Bernoulli number to compute.
            //
            //    Output, double BERNOULLI_NUMBER3, the desired Bernoulli number.
            //
        {
            int i;
            int itmax = 1000;
            const double r8_pi = 3.141592653589793;
            double sum2;
            double term;
            double tol = 5.0E-07;
            double value;

            if (n < 0)
            {
                value = 0.0;
            }
            else if (n == 0)
            {
                value = 1.0;
            }
            else if (n == 1)
            {
                value = -0.5;
            }
            else if (n == 2)
            {
                value = 1.0 / 6.0;
            }
            else if ((n % 2) == 1)
            {
                value = 0.0;
            }
            else
            {
                sum2 = 0.0;

                for (i = 1; i <= itmax; i++)
                {
                    term = 1.0 / Math.Pow((double)i, n);
                    sum2 = sum2 + term;

                    if (Math.Abs(term) < tol || Math.Abs(term) < tol * Math.Abs(sum2))
                    {
                        break;
                    }

                }

                value = 2.0 * sum2 * typeMethods.r8_factorial(n)
                        / Math.Pow((2.0 * r8_pi), n);

                if ((n % 4) == 0)
                {
                    value = -value;
                }

            }

            return value;
        }
    }
}