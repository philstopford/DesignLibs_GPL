using System;
using Burkardt.Types;

namespace Burkardt.SolveNS;

public static class Pell
{
    public static void pell_basic(int d, ref int x0, ref int y0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PELL_BASIC returns the fundamental solution for Pell's basic equation.
        //
        //  Discussion:
        //
        //    Pell's equation has the form:
        //
        //      X^2 - D * Y^2 = 1
        //
        //    where D is a given non-square integer, and X and Y may be assumed
        //    to be positive integers.
        //
        //  Example:
        //
        //     D   X0   Y0
        //
        //     2    3    2
        //     3    2    1
        //     5    9    4
        //     6    5    2
        //     7    8    3
        //     8    3    1
        //    10   19    6
        //    11   10    3
        //    12    7    2
        //    13  649  180
        //    14   15    4
        //    15    4    1
        //    17   33    8
        //    18   17    4
        //    19  170   39
        //    20    9    2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2003
        //
        //  Author:
        //
        //   John Burkardt
        //
        //  Reference:
        //
        //    Mark Herkommer,
        //    Number Theory, A Programmer's Guide,
        //    McGraw Hill, 1999, pages 294-307
        //
        //  Parameters:
        //
        //    Input, int D, the coefficient in Pell's equation.  D should be
        //    positive, and not a perfect square.
        //
        //    Output, int &X0, &Y0, the fundamental or 0'th solution.
        //    If X0 = Y0 = 0, then the calculation was canceled because of an error.
        //    Both X0 and Y0 will be nonnegative.
        //
    {
        const int TERM_MAX = 100;

        int[] b = new int[TERM_MAX + 1];
        int i;
        int p = 0;
        int q = 0;
        int r = 0;
        int term_num = 0;
        switch (d)
        {
            //
            //  Check D.
            //
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("PELL_BASIC - Fatal error!");
                Console.WriteLine("  Pell coefficient D <= 0.");
                x0 = 0;
                y0 = 0;
                return;
        }

        typeMethods.i4_sqrt(d, ref q, ref r);

        switch (r)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("PELL_BASIC - Fatal error!");
                Console.WriteLine("  Pell coefficient is a perfect square.");
                x0 = 0;
                y0 = 0;
                return;
        }

        //
        //  Find the continued fraction representation of sqrt ( D ).
        //
        typeMethods.i4_sqrt_cf(d, TERM_MAX, ref term_num, ref b);
        switch (term_num % 2)
        {
            //
            //  If necessary, go for two periods.
            //
            case 1:
            {
                for (i = term_num + 1; i <= 2 * term_num; i++)
                {
                    b[i] = b[i - term_num];
                }

                term_num = 2 * term_num;
                break;
            }
        }

        //
        //  Evaluate the continued fraction using the forward recursion algorithm.
        //
        int pm2 = 0;
        int pm1 = 1;
        int qm2 = 1;
        int qm1 = 0;

        for (i = 0; i < term_num; i++)
        {
            p = b[i] * pm1 + pm2;
            q = b[i] * qm1 + qm2;
            pm2 = pm1;
            pm1 = p;
            qm2 = qm1;
            qm1 = q;
        }

        //
        //  Get the fundamental solution.
        //
        x0 = p;
        y0 = q;
    }

    public static void pell_next(int d, int x0, int y0, int xn, int yn, ref int xnp1, ref int ynp1)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PELL_NEXT returns the next solution of Pell's equation.
        //
        //  Discussion:
        //
        //    Pell's equation has the form:
        //
        //      X^2 - D * Y^2 = 1
        //
        //    where D is a given non-square integer, and X and Y may be assumed
        //    to be positive integers.
        //
        //    To compute X0, Y0, call PELL_BASIC.
        //    To compute X1, Y1, call this routine, with XN and YN set to X0 and Y0.
        //    To compute further solutions, call again with X0, Y0 and the previous
        //    solution.
        //
        //  Example:
        //
        //    ------INPUT--------  --OUTPUT--
        //
        //    D  X0  Y0   XN   YN  XNP1  YNP1
        //
        //    2   3   2    3    2    17    12
        //    2   3   2   17   12    99    70
        //    2   3   2   99   70   577   408
        //    2   3   2  577  408  3363  2378
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2003
        //
        //  Author:
        //
        //   John Burkardt
        //
        //  Reference:
        //
        //    Mark Herkommer,
        //    Number Theory, A Programmer's Guide,
        //    McGraw Hill, 1999, pages 294-307
        //
        //  Parameters:
        //
        //    Input, int D, the coefficient in Pell's equation.
        //
        //    Input, int X0, Y0, the fundamental or 0'th solution.
        //
        //    Input, int XN, YN, the N-th solution.
        //
        //    Output, int &XNP1, &YNP1, the N+1-th solution.
        //
    {
        xnp1 = x0 * xn + d * y0 * yn;
        ynp1 = x0 * yn + y0 * xn;

    }
}