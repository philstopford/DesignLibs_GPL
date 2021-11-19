using System;
using Burkardt.Types;

namespace Burkardt.SolveNS;

public static class Congruence
{
    public static int congruence(int a, int b, int c, ref bool error)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONGRUENCE solves a congruence of the form ( A * X = C ) mod B.
        //
        //  Discussion:
        //
        //    A, B and C are given integers.  The equation is solvable if and only
        //    if the greatest common divisor of A and B also divides C.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 May 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Eric Weisstein, editor,
        //    CRC Concise Encylopedia of Mathematics,
        //    CRC Press, 2002,
        //    Second edition,
        //    ISBN: 1584883472,
        //    LC: QA5.W45.
        //
        //  Parameters:
        //
        //    Input, int A, B, C, the coefficients of the Diophantine equation.
        //
        //    Output, bool &ERROR, error flag, is TRUE if an error occurred..
        //
        //    Output, int CONGRUENCE, the solution of the Diophantine equation.
        //    X will be between 0 and B-1.
        //
    {
        const int N_MAX = 100;

        int k;
        int[] q = new int[N_MAX];
        bool swap;
        //
        //  Defaults for output parameters.
        //
        error = false;
        int x = 0;
        int y;
        switch (a)
        {
            //
            //  Special cases.
            //
            case 0 when b == 0 && c == 0:
                x = 0;
                return x;
            case 0 when b == 0 && c != 0:
                error = true;
                x = 0;
                return x;
            case 0 when b != 0 && c == 0:
                x = 0;
                return x;
            case 0 when b != 0 && c != 0:
            {
                x = 0;
                if (c % b != 0)
                {
                    error = true;
                }

                return x;
            }
            default:
            {
                if (a != 0 && b == 0 && c == 0)
                {
                    x = 0;
                    return x;
                }

                if (a != 0 && b == 0 && c != 0)
                {
                    x = c / a;
                    if (c % a != 0)
                    {
                        error = true;
                    }

                    return x;
                }
                if (a != 0 && b != 0 && c == 0)
                {
                    //  g = i4_gcd ( a, b );
                    //  x = b / g;
                    x = 0;
                    return x;
                }

                break;
            }
        }

        //
        //  Now handle the "general" case: A, B and C are nonzero.
        //
        //  Step 1: Compute the GCD of A and B, which must also divide C.
        //
        int g = typeMethods.i4_gcd(a, b);

        if (c % g != 0)
        {
            error = true;
            return x;
        }

        int a_copy = a / g;
        int b_copy = b / g;
        int c_copy = c / g;
        //
        //  Step 2: Split A and B into sign and magnitude.
        //
        int a_mag = Math.Abs(a_copy);
        int a_sign = typeMethods.i4_sign(a_copy);
        int b_mag = Math.Abs(b_copy);
        switch (a_mag)
        {
            //
            //  Another special case, A_MAG = 1 or B_MAG = 1.
            //
            case 1:
                x = a_sign * c_copy;
                return x;
        }

        switch (b_mag)
        {
            case 1:
                x = 0;
                return x;
        }

        //
        //  Step 3: Produce the Euclidean remainder sequence.
        //
        if (b_mag <= a_mag)
        {
            swap = false;
            q[0] = a_mag;
            q[1] = b_mag;
        }
        else
        {
            swap = true;
            q[0] = b_mag;
            q[1] = a_mag;
        }

        int n = 3;

        for (;;)
        {
            q[n - 1] = q[n - 3] % q[n - 2];

            if (q[n - 1] == 1)
            {
                break;
            }

            n += 1;

            if (N_MAX >= n)
            {
                continue;
            }

            error = true;
            Console.WriteLine("");
            Console.WriteLine("CONGRUENCE - Fatal error!");
            Console.WriteLine("  Exceeded number of iterations.");
            return 1;
        }

        //
        //  Step 4: Now go backwards to solve X * A_MAG + Y * B_MAG = 1.
        //
        y = 0;
        for (k = n; 2 <= k; k--)
        {
            x = y;
            y = (1 - x * q[k - 2]) / q[k - 1];
        }

        switch (swap)
        {
            //
            //  Step 5: Undo the swapping.
            //
            case true:
                (x, y) = (y, x);
                break;
        }

        //
        //  Step 6: Now apply signs to X and Y so that X * A + Y * B = 1.
        //
        x *= a_sign;
        //
        //  Step 7: Multiply by C, so that X * A + Y * B = C.
        //
        x *= c_copy;
        //
        //  Step 8: Now force 0 <= X < B.
        //
        x %= b;
        switch (x)
        {
            //
            //  Step 9: Force positivity.
            //
            case < 0:
                x += b;
                break;
        }

        return x;
    }
}