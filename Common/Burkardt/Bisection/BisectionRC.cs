using System;
using Burkardt.Types;

namespace Burkardt.Bisection;

public class RC_data
{
    public double fa ;
    public double fb ;
    public int state ;
    public double x_local;

}

public static class RC
{
    public static double bisection_rc(ref RC_data data, ref double a, ref double b, double fx, ref int job )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BISECTION_RC seeks a zero of f(x) in a change of sign interval.
        //
        //  Discussion:
        //
        //    The bisection method is used.
        //
        //    This routine uses reverse communication, so that the function is always
        //    evaluated in the calling program.
        //
        //    On the first call, the user sets JOB = 0, and the values of A and B.
        //    Thereafter, the user checks the returned value of JOB and follows 
        //    directions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, double &A, &B, the endpoints of the change of 
        //    sign interval.  These values will change from call to call as the
        //    interval size is reduced.
        //
        //    Input, double FX, the function value at the point X returned on
        //    the previous call.  On the first call, FX is ignored.
        //
        //    Input/output, int &JOB, a communication flag.
        //    The user sets JOB to 0 before the first call.  Thereafter, the program
        //    controls setting the value of JOB, whose output values mean:
        //
        //    Output, double BISECTION_RC, a point X at which the function is to 
        //    be evaluated.
        //
    {
        double x;

        switch (job)
        {
            case 0:
                data.fa = 0.0;
                data.fb = 0.0;
                data.state = 1;
                x = a;
                job = 1;
                break;
            default:
            {
                switch (data.state)
                {
                    case 1:
                        data.fa = fx;
                        x = b;
                        data.state = 2;
                        break;
                    case 2:
                    {
                        data.fb = fx;

                        if (Math.Abs(typeMethods.r8_sign(data.fa) - typeMethods.r8_sign(data.fb)) <= double.Epsilon)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("BISECTION_RC - Fatal error!");
                            Console.WriteLine("  F(A) and F(B) have the same sign.");
                            return 1;
                        }

                        x = (a + b) / 2.0;
                        data.state = 3;
                        break;
                    }
                    default:
                    {
                        if (Math.Abs(typeMethods.r8_sign(fx) - typeMethods.r8_sign(data.fa)) <= double.Epsilon)
                        {
                            a = data.x_local;
                            data.fa = fx;
                        }
                        else
                        {
                            b = data.x_local;
                            data.fb = fx;
                        }

                        x = (a + b) / 2.0;
                        data.state = 3;
                        break;
                    }
                }

                break;
            }
        }

        data.x_local = x;

        return x;
    }
}