using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8vec_dif(int n, double h)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_DIF computes coefficients for estimating the N-th derivative.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    The routine computes the N+1 coefficients for a centered finite difference
            //    estimate of the N-th derivative of a function.
            //
            //    The estimate has the form
            //
            //      FDIF(N,X) = Sum (I = 0 to N) COF(I) * F ( X(I) )
            //
            //    To understand the computation of the coefficients, it is enough
            //    to realize that the first difference approximation is
            //
            //      FDIF(1,X) = F(X+DX) - F(X-DX) ) / (2*DX)
            //
            //    and that the second difference approximation can be regarded as
            //    the first difference approximation repeated:
            //
            //      FDIF(2,X) = FDIF(1,X+DX) - FDIF(1,X-DX) / (2*DX)
            //         = F(X+2*DX) - 2 F(X) + F(X-2*DX) / (4*DX)
            //
            //    and so on for higher order differences.
            //
            //    Thus, the next thing to consider is the integer coefficients of
            //    the sampled values of F, which are clearly the Pascal coefficients,
            //    but with an alternating negative sign.  In particular, if we
            //    consider row I of Pascal's triangle to have entries j = 0 through I,
            //    then P(I,J) = P(I-1,J-1) - P(I-1,J), where P(*,-1) is taken to be 0,
            //    and P(0,0) = 1.
            //
            //       1
            //      -1  1
            //       1 -2   1
            //      -1  3  -3   1
            //       1 -4   6  -4   1
            //      -1  5 -10  10  -5  1
            //       1 -6  15 -20  15 -6 1
            //
            //    Next, note that the denominator of the approximation for the
            //    N-th derivative will be (2*DX)^N.
            //
            //    And finally, consider the location of the N+1 sampling
            //    points for F:
            //
            //      X-N*DX, X-(N-2)*DX, X-(N-4)*DX, ..., X+(N-4)*DX, X+(N-2*DX), X+N*DX.
            //
            //    Thus, a formula for evaluating FDIF(N,X) is
            //
            //      fdif = 0.0
            //      do i = 0, n
            //        xi = x + (2*i-n) * h
            //        fdif = fdif + cof(i) * f(xi)
            //      end do
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the derivative to be approximated.
            //    N must be 0 or greater.
            //
            //    Input, double H, the half spacing between points.
            //    H must be positive.
            //
            //    Output, double R8VEC_DIF[N+1], the coefficients needed to approximate
            //    the N-th derivative of a function F.
            //
        {
            double[] cof;
            int i;
            int j;

            if (n < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_DIF - Fatal error!");
                Console.WriteLine("  Derivative order N = " + n + "");
                Console.WriteLine("  but N must be at least 0.");
                return null;
            }

            if (h <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_DIF - Fatal error!");
                Console.WriteLine("  The half sampling spacing is H = " + h + "");
                Console.WriteLine("  but H must be positive.");
                return null;
            }

            cof = new double[n + 1];

            for (i = 0; i <= n; i++)
            {
                cof[i] = 1.0;

                for (j = i - 1; 1 <= j; j--)
                {
                    cof[j] = -cof[j] + cof[j - 1];
                }

                if (0 < i)
                {
                    cof[0] = -cof[0];
                }
            }

            for (i = 0; i <= n; i++)
            {
                cof[i] = cof[i] / Math.Pow(2.0 * h, n);
            }

            return cof;
        }

        public static double r8vec_diff_norm(int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The vector L2 norm is defined as:
        //
        //      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 June 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], B[N], the vectors.
        //
        //    Output, double R8VEC_DIFF_NORM, the L2 norm of A - B.
        //
        {
            int i;
            double value;

            value = 0.0;

            for (i = 0; i < n; i++)
            {
                value = value + (a[i] - b[i]) * (a[i] - b[i]);
            }

            value = Math.Sqrt(value);

            return value;
        }

        public static double r8vec_diff_norm_l1(int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_DIFF_NORM_L1 returns the L1 norm of the difference of R8VEC's.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The vector L1 norm is defined as:
        //
        //      R8VEC_NORM_L1 = sum ( 1 <= I <= N ) abs ( A(I) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], B[N], the vectors.
        //
        //    Output, double R8VEC_DIFF_NORM_L1, the L1 norm of A - B.
        //
        {
            int i;
            double value;

            value = 0.0;

            for (i = 0; i < n; i++)
            {
                value = value + Math.Abs(a[i] - b[i]);
            }

            return value;
        }

        public static double r8vec_diff_norm_l2(int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_DIFF_NORM_L2 returns the L2 norm of the difference of R8VEC's.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The vector L2 norm is defined as:
        //
        //      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 June 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], B[N], the vectors.
        //
        //    Output, double R8VEC_DIFF_NORM_L2, the L2 norm of A - B.
        //
        {
            int i;
            double value;

            value = 0.0;

            for (i = 0; i < n; i++)
            {
                value = value + (a[i] - b[i]) * (a[i] - b[i]);
            }

            value = Math.Sqrt(value);

            return value;
        }

        public static double r8vec_diff_norm_li(int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_DIFF_NORM_LI returns the L-oo norm of the difference of R8VEC's.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The vector L-oo norm is defined as:
        //
        //      R8VEC_NORM_LI = max ( 1 <= I <= N ) abs ( A(I) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], B[N], the vectors.
        //
        //    Output, double R8VEC_DIFF_NORM_LI, the L-oo norm of A - B.
        //
        {
            int i;
            double value;

            value = 0.0;

            for (i = 0; i < n; i++)
            {
                value = Math.Max(value, Math.Abs(a[i] - b[i]));
            }

            return value;
        }

        public static double r8vec_diff_norm_squared(int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_DIFF_NORM_SQUARED: square of the L2 norm of the difference of R8VEC's.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The square of the L2 norm of the difference of A and B is:
        //
        //      R8VEC_DIFF_NORM_SQUARED = sum ( 1 <= I <= N ) ( A[I] - B[I] )^2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 June 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], B[N], the vectors.
        //
        //    Output, double R8VEC_DIFF_NORM_SQUARED, the square of the L2 norm of A - B.
        //
        {
            int i;
            double value;

            value = 0.0;

            for (i = 0; i < n; i++)
            {
                value = value + (a[i] - b[i]) * (a[i] - b[i]);
            }

            return value;
        }
    }
}