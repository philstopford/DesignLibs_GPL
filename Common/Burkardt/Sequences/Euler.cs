using System;
using Burkardt.Types;

namespace Burkardt.Sequence;

public static class Euler
{
    public static void eulerian(int n, ref int[] e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EULERIAN computes the Eulerian number E(N,K).
        //
        //  Definition:
        //
        //    A run in a permutation is a sequence of consecutive ascending values.
        //
        //    E(N,K) is the number of permutations of N objects which contain
        //    exactly K runs.
        //
        //  Examples:
        //
        //     N = 7
        //
        //     1     0     0     0     0     0     0
        //     1     1     0     0     0     0     0
        //     1     4     1     0     0     0     0
        //     1    11    11     1     0     0     0
        //     1    26    66    26     1     0     0
        //     1    57   302   302    57     1     0
        //     1   120  1191  2416  1191   120     1
        //
        //  Recursion:
        //
        //    E(N,K) = K * E(N-1,K) + (N-K+1) * E(N-1,K-1).
        //
        //  Properties:
        //
        //    E(N,1) = E(N,N) = 1.
        //    E(N,K) = 0 if K <= 0 or N < K.
        //    sum ( 1 <= K <= N ) E(N,K) = N!.
        //    X^N = sum ( 0 <= K <= N ) COMB(X+K-1, N ) E(N,K)
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
        //  Reference:
        //
        //    Dennis Stanton and Dennis White,
        //    Constructive Combinatorics,
        //    Springer Verlag, 1986
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows desired.
        //
        //    Output, int E[N*N], the first N rows of Eulerian numbers.
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 1:
                return;
        }

        //
        //  Construct rows 1, 2, ..., N of the Eulerian triangle.
        //
        e[1 - 1 + (1 - 1) * n] = 1;
        for (j = 2; j <= n; j++)
        {
            e[1 - 1 + (j - 1) * n] = 0;
        }

        for (i = 2; i <= n; i++)
        {
            e[i - 1 + (1 - 1) * n] = 1;
            for (j = 2; j <= n; j++)
            {
                e[i - 1 + (j - 1) * n] = j * e[i - 2 + (j - 1) * n] + (i - j + 1) * e[i - 2 + (j - 2) * n];
            }
        }
    }

    public static void euler_number(int n, ref int[] e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EULER_NUMBER computes the Euler numbers.
        //
        //  Discussion:
        //
        //    The Euler numbers can be evaluated in Mathematica with the call
        //
        //      EulerE[n]
        //
        //    These numbers rapidly get too big to store in an ordinary integer!
        //
        //    The terms of odd index are 0.
        //
        //    E(N) = -C(N,N-2) * E(N-2) - C(N,N-4) * E(N-4) - ... - C(N,0) * E(0).
        //
        //  First terms:
        //
        //    E0  = 1
        //    E1  = 0
        //    E2  = -1
        //    E3  = 0
        //    E4  = 5
        //    E5  = 0
        //    E6  = -61
        //    E7  = 0
        //    E8  = 1385
        //    E9  = 0
        //    E10 = -50521
        //    E11 = 0
        //    E12 = 2702765
        //    E13 = 0
        //    E14 = -199360981
        //    E15 = 0
        //    E16 = 19391512145
        //    E17 = 0
        //    E18 = -2404879675441
        //    E19 = 0
        //    E20 = 370371188237525
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input, int N, the index of the last Euler number to compute.
        //
        //    Output, int E[N+1], the Euler numbers from index 0 to N.
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return;
        }

        e[0] = 1;

        switch (n)
        {
            case 0:
                return;
        }

        e[1] = 0;

        switch (n)
        {
            case 1:
                return;
        }

        e[2] = -1;

        for (i = 3; i <= n; i++)
        {
            e[i] = 0;

            switch (i % 2)
            {
                case 0:
                {
                    for (j = 2; j <= i; j += 2)
                    {
                        e[i] -= typeMethods.i4_choose(i, j) * e[i - j];
                    }

                    break;
                }
            }
        }
    }

    public static double euler_number2(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EULER_NUMBER2 computes the Euler numbers.
        //
        //  Discussion:
        //
        //    The Euler numbers can be evaluated in Mathematica with the call
        //
        //      EulerE[n]
        //
        //  First terms:
        //
        //    E0  = 1
        //    E1  = 0
        //    E2  = -1
        //    E3  = 0
        //    E4  = 5
        //    E5  = 0
        //    E6  = -61
        //    E7  = 0
        //    E8  = 1385
        //    E9  = 0
        //    E10 = -50521
        //    E11 = 0
        //    E12 = 2702765
        //    E13 = 0
        //    E14 = -199360981
        //    E15 = 0
        //    E16 = 19391512145
        //    E17 = 0
        //    E18 = -2404879675441
        //    E19 = 0
        //    E20 = 370371188237525
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
        //  Reference:
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input, int N, the index of the Euler number to compute.
        //
        //    Output, double EULER_NUMBER2, the value of E(N).
        //
    {
        double[] e =  {
                1.0, -1.0, 5.0, -61.0, 1385.0,
                -50521.0, 2702765.0
            }
            ;
        int i;
        int itmax = 1000;
            
        double sum1;
        double term;
        double value = 0;

        switch (n)
        {
            case < 0:
                return 0.0;
            case 0:
                return e[0];
        }

        switch (n % 2)
        {
            case 1:
                return 0.0;
        }

        switch (n)
        {
            case <= 12:
                return e[n / 2];
        }

        sum1 = 0.0;

        for (i = 1; i <= itmax; i++)
        {
            term = 1.0 / Math.Pow(2 * i - 1, n + 1);

            switch (i % 2)
            {
                case 1:
                    sum1 += term;
                    break;
                default:
                    sum1 -= term;
                    break;
            }

            if (Math.Abs(term) < 1.0E-10)
            {
                break;
            }

            if (Math.Abs(term) < 1.0E-08 * Math.Abs(sum1))
            {
                break;
            }

        }

        value = Math.Pow(2.0, n + 2) * sum1 * typeMethods.r8_factorial(n)
                / Math.Pow(Math.PI, n + 1);

        if (n % 4 != 0)
        {
            value = -value;
        }

        return value;
    }
}