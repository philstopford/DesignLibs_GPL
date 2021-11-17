using System;
using Burkardt.Types;

namespace Burkardt.SolveNS;

public static class Diophantine
{
    public static void diophantine(int a, int b, int c, ref bool error, ref int x, ref int y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIOPHANTINE solves a Diophantine equation A * X + B * Y = C.
        //
        //  Discussion:
        //
        //    Given integers A, B and C, produce X and Y so that
        //
        //      A * X + B * Y = C.
        //
        //    In general, the equation is solvable if and only if the
        //    greatest common divisor of A and B also divides C.
        //
        //    A solution (X,Y) of the Diophantine equation also gives the solution
        //    X to the congruence equation:
        //
        //      A * X = C mod ( B ).
        //
        //    Generally, if there is one nontrivial solution, there are an infinite
        //    number of solutions to a Diophantine problem.
        //    If (X0,Y0) is a solution, then so is ( X0+T*B/D, Y0-T*A/D ) where
        //    T is any integer, and D is the greatest common divisor of A and B.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Eric Weisstein, editor,
        //    CRC Concise Encylopedia of Mathematics,
        //    CRC Press, 1998, page 446.
        //
        //  Parameters:
        //
        //    Input, int A, B, C, the coefficients of the Diophantine equation.
        //
        //    Output, bool &ERROR, is TRUE if an error occurred.
        //
        //    Output, int &X, &Y, the solution of the Diophantine equation.
        //    Note that the algorithm will attempt to return a solution with
        //    smallest Euclidean norm.
        //
    {
        const int N_MAX = 100;

        int a_copy;
        int a_mag;
        int a_sign;
        int b_copy;
        int b_mag;
        int b_sign;
        int c_copy;
        int g;
        int k;
        int n;
        int[] q = new int[N_MAX];
        bool swap;
        //
        //  Defaults for output parameters.
        //
        error = false;
        x = 0;
        y = 0;
        switch (a)
        {
            //
            //  Special cases.
            //
            case 0 when b == 0 && c == 0:
                x = 0;
                y = 0;
                return;
            case 0 when b == 0 && c != 0:
                error = true;
                x = 0;
                y = 0;
                return;
            case 0 when b != 0 && c == 0:
                x = 0;
                y = 0;
                return;
            case 0 when b != 0 && c != 0:
            {
                x = 0;
                y = c / b;
                if (c % b != 0)
                {
                    error = true;
                }

                return;
            }
            default:
            {
                if (a != 0 && b == 0 && c == 0)
                {
                    x = 0;
                    y = 0;
                    return;
                }

                if (a != 0 && b == 0 && c != 0)
                {
                    x = c / a;
                    y = 0;
                    if (c % a != 0)
                    {
                        error = true;
                    }

                    return;
                }
                if (a != 0 && b != 0 && c == 0)
                {
                    g = typeMethods.i4_gcd(a, b);
                    x = b / g;
                    y = -a / g;
                    return;
                }

                break;
            }
        }

        //
        //  Now handle the "general" case: A, B and C are nonzero.
        //
        //  Step 1: Compute the GCD of A and B, which must also divide C.
        //
        g = typeMethods.i4_gcd(a, b);

        if (c % g != 0)
        {
            error = true;
            return;
        }

        a_copy = a / g;
        b_copy = b / g;
        c_copy = c / g;
        //
        //  Step 2: Split A and B into sign and magnitude.
        //
        a_mag = Math.Abs(a_copy);
        a_sign = typeMethods.i4_sign(a_copy);
        b_mag = Math.Abs(b_copy);
        b_sign = typeMethods.i4_sign(b_copy);
        switch (a_mag)
        {
            //
            //  Another special case, A_MAG = 1 or B_MAG = 1.
            //
            case 1:
                x = a_sign * c_copy;
                y = 0;
                return;
        }

        switch (b_mag)
        {
            case 1:
                x = 0;
                y = b_sign * c_copy;
                return;
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

        n = 3;

        for (;;)
        {
            q[n - 1] = q[n - 3] % q[n - 2];

            if (q[n - 1] == 1)
            {
                break;
            }

            n += 1;

            if (N_MAX < n)
            {
                error = true;
                Console.WriteLine("");
                Console.WriteLine("DIOPHANTINE - Fatal error!");
                Console.WriteLine("  Exceeded number of iterations.");
                return;
            }
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
                typeMethods.i4_swap(ref x, ref y);
                break;
        }

        //
        //  Step 6: Now apply signs to X and Y so that X * A + Y * B = 1.
        //
        x *= a_sign;
        y *= b_sign;
        //
        //  Step 7: Multiply by C, so that X * A + Y * B = C.
        //
        x *= c_copy;
        y *= c_copy;
        //
        //  Step 8: Given a solution (X,Y), try to find the solution of
        //  minimal magnitude.
        //
        diophantine_solution_minimize(a_copy, b_copy, ref x, ref y);

    }

    public static void diophantine_solution_minimize(int a, int b, ref int x, ref int y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIOPHANTINE_SOLUTION_MINIMIZE seeks a minimal solution of a Diophantine equation.
        //
        //  Discussion:
        //
        //    Given a solution (X,Y) of a Diophantine equation:
        //
        //      A * X + B * Y = C.
        //
        //    then there are an infinite family of solutions of the form
        //
        //      ( X(i), Y(i) ) = ( X + i * B, Y - i * A )
        //
        //    An integral solution of minimal Euclidean norm can be found by
        //    tentatively moving along the vectors (B,-A) and (-B,A) one step
        //    at a time.
        //
        //    When large integer values are input, the real arithmetic used
        //    is essential.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Eric Weisstein, editor,
        //    CRC Concise Encylopedia of Mathematics,
        //    CRC Press, 1998, page 446.
        //
        //  Parameters:
        //
        //    Input, int A, B, the coefficients of the Diophantine equation.
        //    A and B are assumed to be relatively prime.
        //
        //    Input/output, int &X, &Y, on input, a solution of the Diophantine
        //    equation.  On output, a solution of minimal Euclidean norm.
        //
    {
        double fa;
        double fb;
        double fx;
        double fy;
        double norm;
        double norm_new;
        double t;
        int xnew;
        int ynew;
        //
        //  Compute the minimum for T real, and then look nearby.
        //
        fa = a;
        fb = b;
        fx = x;
        fy = y;

        t = (-fb * fx + fa * fy) / (fa * fa + fb * fb);

        x += (int)typeMethods.r8_nint(t) * b;
        y -= (int)typeMethods.r8_nint(t) * a;
        //
        //  Now look nearby.
        //
        norm = fx * fx + fy * fy;

        for (;;)
        {
            xnew = x + b;
            ynew = y - a;

            fx = xnew;
            fy = ynew;

            norm_new = fx * fx + fy * fy;

            if (norm <= norm_new)
            {
                break;
            }

            x = xnew;
            y = ynew;
            norm = norm_new;
        }

        for (;;)
        {
            xnew = x - b;
            ynew = y + a;

            fx = xnew;
            fy = ynew;

            norm_new = fx * fx + fy * fy;

            if (norm <= norm_new)
            {
                break;
            }

            x = xnew;
            y = ynew;
            norm = norm_new;
        }

    }
}