using System;

namespace Burkardt.Bisection;

public static class Integer
{
    public static void bisection_integer(Func<int,int> f, ref int a, ref int b, ref int c, ref int fc )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BISECTION_INTEGER seeks an integer root using bisection.
        //
        //  Discussion:
        //
        //    A function F(X) confined to integer arguments is given, with an
        //    interval [A,B] over which F changes sign.  An integer C is sought
        //    such that A <= C <= B and F(C) = 0.
        //
        //    Because we are restricted to integer arguments, it may the case that
        //   there is no such C.
        //
        //    This routine proceeds by a form of bisection, in which the enclosing
        //    interval is restricted to be defined by integer values.
        //
        //    If the user has given a true change of sign interval [A,B], and if,
        //    in the interval, there is a single integer value C for which F(C) = 0,
        //    with the additional restrictions that F(C-1) and F(C+1) are of opposite
        //    signs, then this procedure should locate and return C.
        //
        //    In particular, if the function F is monotone, and there is an integer
        //    solution C in the interval, then this procedure will find it.
        //
        //    However, in general, even if there is an integer C in the interval,
        //    such that F(C) = 0, this procedure may be unable to find it, particularly
        //    if there are also nonintegral solutions within the same interval.
        //
        //    While any integer function can be used with this program, the bisection
        //    approach is most useful if the integer function is monotone, or
        //    varies slowly, or can be regarded as the restriction to integer arguments
        //    of a continuous (and smoothly varying) function of a real argument.
        //    In such cases, knowing that F is negative at A and positive at B
        //    suggests that F generally increases from A to B, and might attain 
        //    the value 0 at some intermediate argument C.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 August 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int F ( int X ), the name of a user-supplied 
        //    procedure that evaluates the function.
        //
        //   Input, int &A, &B, two arguments that define a change of
        //   sign interval for F.  In other words, F(A) and F(B) must be of opposite
        //   sign.
        //
        //   Output, int &C, &FC, the candidate for the root, as 
        //   determined by the program, and its function value.  If FC is not zero,
        //   then the procedure did not find a root in the interval, and C is only
        //   an "approximate" root.
        //
    {
        int fa;
        int fb;
        int t;
        //
        //  Ensure that F(A) < 0 < F(B).
        //
        fa = f(a);
        fb = f(b);

        switch (fa)
        {
            case 0:
                c = a;
                fc = fa;
                break;
            default:
            {
                switch (fb)
                {
                    case 0:
                        c = b;
                        fc = fb;
                        break;
                    default:
                    {
                        switch (fa)
                        {
                            case < 0 when 0 < fb:
                                break;
                            default:
                            {
                                switch (fb)
                                {
                                    case < 0 when 0 < fa:
                                        t = a;
                                        a = b;
                                        b = t;
                                        t = fa;
                                        fa = fb;
                                        fb = t;
                                        break;
                                    default:
                                        Console.WriteLine("");
                                        Console.WriteLine("BISECTION_INTEGER - Fatal error!");
                                        Console.WriteLine("  No change of sign interval supplied.");
                                        Console.WriteLine("  F(" + a + ") = " + fa + "");
                                        Console.WriteLine("  F(" + b + ") = " + fb + "");
                                        return;
                                }

                                break;
                            }
                        }

                        break;
                    }
                }

                break;
            }
        }

        //
        //  Bisection.
        //
        while (1 < Math.Abs(b - a))
        {
            c = (a + b) / 2;
            fc = f(c);

            switch (fc)
            {
                case 0:
                    return;
                case < 0:
                    a = c;
                    fa = fc;
                    break;
                case > 0:
                    b = c;
                    fb = fc;
                    break;
            }
        }

        //
        //  Interval is empty, with FA < 0 and 0 < FB.
        //  Bisection did not produce an integer solution.
        //  Return the argument with smallest function norm.
        //
        if (-fa < fb)
        {
            c = a;
            fc = fa;
        }
        else
        {
            c = b;
            fc = fb;
        }
    }
}