using System;

namespace Burkardt.Probability;

public static class Misc
{
    public static double sin_power_int(double a, double b, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIN_POWER_INT evaluates the sine power integral.
        //
        //  Discussion:
        //
        //    The function is defined by
        //
        //      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin ( t ))^n dt
        //
        //    The algorithm uses the following fact:
        //
        //      Integral sin^n ( t ) = (1/n) * (
        //        sin^(n-1)(t) * cos(t) + ( n-1 ) * Integral sin^(n-2) ( t ) dt )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters
        //
        //    Input, double A, B, the limits of integration.
        //
        //    Input, int N, the power of the sine function.
        //
        //    Output, double SIN_POWER_INT, the value of the integral.
        //
    {
        int mlo;
        double value = 0;
        switch (n)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("SIN_POWER_INT - Fatal error!");
                Console.WriteLine("  Power N < 0.");
                return 1.0;
        }

        double sa = Math.Sin(a);
        double sb = Math.Sin(b);
        double ca = Math.Cos(a);
        double cb = Math.Cos(b);

        switch (n % 2)
        {
            case 0:
                value = b - a;
                mlo = 2;
                break;
            default:
                value = ca - cb;
                mlo = 3;
                break;
        }

        for (int m = mlo; m <= n; m += 2)
        {
            value = ((m - 1) * value
                        + Math.Pow(sa, m - 1) * ca - Math.Pow(sb, m - 1) * cb)
                    / m;
        }

        return value;
    }

    public static double euler_constant()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EULER_CONSTANT returns the value of the Euler-Mascheroni constant.
        //
        //  Discussion:
        //
        //    The Euler-Mascheroni constant is often denoted by a lower-case
        //    Gamma.  Gamma is defined as
        //
        //      Gamma = limit ( M -> Infinity )
        //        ( Sum ( 1 <= N <= M ) 1 / N ) - Log ( M )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double EULER_CONSTANT, the value of the
        //    Euler-Mascheroni constant.
        //
    {
        const double value = 0.577215664901532860606512090082402431042;

        return value;
    }

    public static double sphere_unit_area_nd(int dim_num)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_UNIT_AREA_ND computes the surface area of a unit sphere in ND.
        //
        //  Discussion:
        //
        //    The unit sphere in ND satisfies the equation:
        //
        //      Sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
        //
        //    DIM_NUM   Area
        //
        //     2    2        * PI
        //     3    4        * PI
        //     4  ( 2 /   1) * PI^2
        //     5  ( 8 /   3) * PI^2
        //     6  ( 1 /   1) * PI^3
        //     7  (16 /  15) * PI^3
        //     8  ( 1 /   3) * PI^4
        //     9  (32 / 105) * PI^4
        //    10  ( 1 /  12) * PI^5
        //
        //    For the unit sphere, Area(DIM_NUM) = DIM_NUM * Volume(DIM_NUM)
        //
        //    Sphere_Unit_Area ( DIM_NUM ) = 2 * PI^(DIM_NUM/2) / Gamma ( DIM_NUM / 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Output, double SPHERE_UNIT_AREA_ND, the area of the sphere.
        //
    {
        double area;
        int i;
        int m;
            

        switch (dim_num % 2)
        {
            case 0:
            {
                m = dim_num / 2;
                area = 2.0 * Math.Pow(Math.PI, m);
                for (i = 1; i <= m - 1; i++)
                {
                    area /= i;
                }

                break;
            }
            default:
            {
                m = (dim_num - 1) / 2;
                area = Math.Pow(2.0, dim_num) * Math.Pow(Math.PI, m);
                for (i = m + 1; i <= 2 * m; i++)
                {
                    area /= i;
                }

                break;
            }
        }

        return area;
    }

    public static int stirling2_value(int n, int m)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STIRLING2_VALUE computes a Stirling number of the second kind.
        //
        //  Discussion:
        //
        //    S2(N,M) represents the number of distinct partitions of N elements
        //    into M nonempty sets.  For a fixed N, the sum of the Stirling
        //    numbers S2(N,M) is represented by B(N), called "Bell's number",
        //    and represents the number of distinct partitions of N elements.
        //
        //    For example, with 4 objects, there are:
        //
        //    1 partition into 1 set:
        //
        //      (A,B,C,D)
        //
        //    7 partitions into 2 sets:
        //
        //      (A,B,C) (D)
        //      (A,B,D) (C)
        //      (A,C,D) (B)
        //      (A) (B,C,D)
        //      (A,B) (C,D)
        //      (A,C) (B,D)
        //      (A,D) (B,C)
        //
        //    6 partitions into 3 sets:
        //
        //      (A,B) (C) (D)
        //      (A) (B,C) (D)
        //      (A) (B) (C,D)
        //      (A,C) (B) (D)
        //      (A,D) (B) (C)
        //      (A) (B,D) (C)
        //
        //    1 partition into 4 sets:
        //
        //      (A) (B) (C) (D)
        //
        //    So S2(4,1) = 1, S2(4,2) = 7, S2(4,3) = 6, S2(4,4) = 1, and B(4) = 15.
        //
        //
        //  First terms:
        //
        //    N/M: 1    2    3    4    5    6    7    8
        //
        //    1    1    0    0    0    0    0    0    0
        //    2    1    1    0    0    0    0    0    0
        //    3    1    3    1    0    0    0    0    0
        //    4    1    7    6    1    0    0    0    0
        //    5    1   15   25   10    1    0    0    0
        //    6    1   31   90   65   15    1    0    0
        //    7    1   63  301  350  140   21    1    0
        //    8    1  127  966 1701 1050  266   28    1
        //
        //  Recursion:
        //
        //    S2(N,1) = 1 for all N.
        //    S2(I,I) = 1 for all I.
        //    S2(I,J) = 0 if I < J.
        //
        //    S2(N,M) = M * S2(N-1,M) + S2(N-1,M-1)
        //
        //  Properties:
        //
        //    sum ( 1 <= K <= M ) S2(I,K) * S1(K,J) = Delta(I,J)
        //
        //    X^N = sum ( 0 <= K <= N ) S2(N,K) X_K
        //    where X_K is the falling factorial function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows of the table.
        //
        //    Input, int M, the number of columns of the table.
        //
        //    Output, int STIRLING2_VALUE, the value of S2(N,M).
        //
    {
        int value;

        int[] s2 = new int[n * m];

        switch (n)
        {
            case <= 0:
                value = 0;
                return value;
        }

        switch (m)
        {
            case <= 0:
                value = 0;
                return value;
        }

        s2[0 + 0 * n] = 1;
        for (int j = 2; j <= m; j++)
        {
            s2[0 + (j - 1) * n] = 0;
        }

        for (int i = 2; i <= n; i++)
        {
            s2[i - 1 + 0 * n] = 1;
            for (int j = 2; j <= m; j++)
            {
                s2[i - 1 + (j - 1) * n] = j * s2[i - 2 + (j - 1) * n] + s2[i - 2 + (j - 2) * n];
            }
        }

        value = s2[n - 1 + (m - 1) * n];
            
        return value;
    }
        
    public static double trigamma ( double x )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIGAMMA calculates the TriGamma function.
        //
        //  Discussion:
        //
        //    TriGamma(x) = d^2 log ( Gamma ( x ) ) / dx^2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 January 2000
        //
        //  Author:
        //
        //    Original FORTRAN77 version by B Schneider
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    B Schneider,
        //    Trigamma Function,
        //    Algorithm AS 121,
        //    Applied Statistics,
        //    Volume 27, Number 1, page 97-99, 1978.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the trigamma function.
        //    0 < X.
        //
        //    Output, double TRIGAMMA, the value of the
        //    trigamma function at X.
        //
    {
        const double a = 0.0001;
        const double b = 5.0;
        const double b2 =   1.0 / 6.0;
        const double b4 = - 1.0 / 30.0;
        const double b6 =   1.0 / 42.0;
        const double b8 = - 1.0 / 30.0;
        double value;
        switch (x)
        {
            //
            //  1): If X is not positive, fail.
            //
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("TRIGAMMA - Fatal error!");
                Console.WriteLine("  X <= 0.");
                return 1;
        }
        //
        //  2): If X is smaller than A, use a small value approximation.
        //

        if ( x <= a )
        {
            value = 1.0 / x / x;
        }
        //
        //  3): Otherwise, increase the argument to B <= ( X + I ).
        //
        else
        {
            double z = x;
            value = 0.0;

            while ( z < b )
            {
                value += 1.0 / z / z;
                z += 1.0;
            }
            //
            //  ...and then apply an asymptotic formula.
            //
            double y = 1.0 / z / z;

            value = value + 0.5 *
                y + ( 1.0
                      + y * ( b2
                              + y * ( b4
                                      + y * ( b6
                                              + y *   b8 )))) / z;
        }

        return value;
    }

}