using System;
using System.Globalization;
using Burkardt.CDFLib;

namespace Burkardt.Types;

public class r8
{
    public bool error { get; set; }
    public double val { get; set; }
    public int lchar { get; set; }
}

public class r8vec
{
    public bool error { get; set; }
    public double[] rvec { get; set; }
}

public static partial class typeMethods
{
    public static double r8_abs(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ABS returns the absolute value of an R8.
        //
        //  Discussion:
        //
        //    The C++ math library provides the function fabs() which is preferred.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 November 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the quantity whose absolute value is desired.
        //
        //    Output, double R8_ABS, the absolute value of X.
        //
    {
        double value = x switch
        {
            >= 0.0 => +x,
            _ => -x
        };

        return value;
    }

    public static double r8_acos(double c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ACOS computes the arc cosine function, with argument truncation.
        //
        //  Discussion:
        //
        //    If you call your system ACOS routine with an input argument that is
        //    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
        //    This routine truncates arguments outside the range.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 June 2002
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double C, the argument, the cosine of an angle.
        //
        //    Output, double R8_ACOS, an angle whose cosine is C.
        //
    {
        double value = c switch
        {
            <= -1.0 => Math.PI,
            >= 1.0 => 0.0,
            _ => Math.Acos(c)
        };

        return value;
    }

    public static double r8_acosh(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ACOSH returns the inverse hyperbolic cosine of a number.
        //
        //  Discussion:
        //
        //    Applying the inverse function
        //
        //      Y = R8_ACOSH(X)
        //
        //    implies that
        //
        //      X = COSH(Y) = 0.5 * ( EXP(Y) + EXP(-Y) ).
        //
        //    For every X greater than or equal to 1, there are two possible
        //    choices Y such that X = COSH(Y), differing only in sign.  It
        //    is usual to resolve this choice by taking the value of ACOSH(X)
        //    to be nonnegative.
        //
        //  Method:
        //
        //    One formula is:
        //
        //      R8_ACOSH = LOG ( X + SQRT ( X^2 - 1.0 ) )
        //
        //    but this formula suffers from roundoff and overflow problems.
        //    The formula used here was recommended by W Kahan, as discussed
        //    by Moler.
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
        //  Reference:
        //
        //    Cleve Moler,
        //    Trigonometry is a Complex Subject,
        //    MATLAB News and Notes,
        //    Summer 1998.
        //
        //  Parameters:
        //
        //    Input, double X, the number whose inverse hyperbolic cosine is desired.
        //    X should be greater than or equal to 1.
        //
        //    Output, double R8_ACOSH, the inverse hyperbolic cosine of X.  The
        //    principal value (that is, the positive value of the two ) is returned.
        //
    {
        switch (x)
        {
            case < 1.0:
                Console.WriteLine("");
                Console.WriteLine("R8_ACOSH - Fatal error!");
                Console.WriteLine("  Argument X must satisfy 1 <= X.");
                Console.WriteLine("  The input X = " + x + "");
                return 1;
            default:
                double value = 2.0 * Math.Log(
                    Math.Sqrt(0.5 * (x + 1.0)) + Math.Sqrt(0.5 * (x - 1.0)));

                return value;
        }
    }

    public static double r8_add(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ADD adds two R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the numbers to be added.
        //
        //    Output, double R8_ADD, the sum of X and Y.
        //
    {
        double value = x + y;

        return value;
    }

    public static double r8_atanh(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ATANH returns the inverse hyperbolic tangent of a number.
        //
        //  Discussion:
        //
        //    Y = R8_ATANH ( X )
        //
        //    implies that
        //
        //    X = TANH(Y) = ( EXP(Y) - EXP(-Y) ) / ( EXP(Y) + EXP(-Y) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 November 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the number whose inverse hyperbolic 
        //    tangent is desired.  The absolute value of X should be less than 
        //    or equal to 1.
        //
        //    Output, double R8_ATANH, the inverse hyperbolic tangent of X.
        //
    {
        double value = x switch
        {
            <= -1.0 => -r8_huge(),
            >= 1.0 => +r8_huge(),
            _ => 0.5 * Math.Log((1.0 + x) / (1.0 - x))
        };

        return value;
    }


    public static double r8_agm(double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_AGM computes the arithmetic-geometric mean of A and B.
        //
        //  Discussion:
        //
        //    The AGM is defined for nonnegative A and B.
        //
        //    The AGM of numbers A and B is defined by setting
        //
        //      A(0) = A,
        //      B(0) = B
        //
        //      A(N+1) = ( A(N) + B(N) ) / 2
        //      B(N+1) = sqrt ( A(N) * B(N) )
        //
        //    The two sequences both converge to AGM(A,B).
        //
        //    In Mathematica, the AGM can be evaluated by
        //
        //      ArithmeticGeometricMean [ a, b ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 July 2014
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
        //    Input, double A, B, the arguments whose AGM is to be computed.
        //    0 <= A, 0 <= B.
        //
        //    Output, double R8_AGM, the arithmetic-geometric mean of A and B.
        //
    {
        double a2;
        const int it_max = 1000;
        double value;

        switch (a)
        {
            case < 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_AGM - Fatal error!");
                Console.WriteLine("  A < 0.");
                return 1;
        }

        switch (b)
        {
            case < 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_AGM - Fatal error!");
                Console.WriteLine("  B < 0.");
                return 1;
        }

        if (a == 0.0 || b == 0.0)
        {
            value = 0.0;
            return value;
        }

        if (Math.Abs(a - b) <= typeMethods.r8_epsilon())
        {
            value = a;
            return value;
        }

        double a1 = a;
        double b1 = b;

        int it = 0;
        double tol = 100.0 * r8_epsilon();

        for (;;)
        {
            it += 1;

            a2 = (a1 + b1) / 2.0;
            double b2 = Math.Sqrt(a1 * b1);

            if (Math.Abs(a2 - b2) <= tol * (a2 + b2))
            {
                break;
            }

            if (it_max < it)
            {
                break;
            }

            a1 = a2;
            b1 = b2;
        }

        value = a2;

        return value;
    }

    public static double r8_aint(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_AINT truncates an R8 argument to an integer.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    1 September 2011
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_AINT, the truncated version of X.
        //
    {
        double value = x switch
        {
            < 0.0 => -(double) (int) Math.Abs(x),
            _ => Math.Abs(x)
        };

        return value;
    }

    public static double r8_asin(double s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ASIN computes the arc sine function, with argument truncation.
        //
        //  Discussion:
        //
        //    If you call your system ASIN routine with an input argument that is
        //    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
        //    This routine truncates arguments outside the range.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 June 2002
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double S, the argument, the sine of an angle.
        //
        //    Output, double R8_ASIN, an angle whose sine is S.
        //
    {
        double angle = s switch
        {
            <= -1.0 => -Math.PI / 2.0,
            >= 1.0 => Math.PI / 2.0,
            _ => Math.Asin(s)
        };

        return angle;
    }

    public static double r8_asinh(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ASINH returns the inverse hyperbolic sine of a number.
        //
        //  Discussion:
        //
        //    The assertion that:
        //
        //      Y = R8_ASINH ( X )
        //
        //    implies that
        //
        //      X = SINH(Y) = 0.5 * ( EXP(Y) - EXP(-Y) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 November 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the number whose inverse hyperbolic 
        //    sine is desired.
        //
        //    Output, double R8_ASINH, the inverse hyperbolic sine of X.
        //
    {
        double value = Math.Log(x + Math.Sqrt(x * x + 1.0));

        return value;
    }

    public static double r8_radians(double degrees)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_RADIANS converts an angle from degree to radian measure.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double DEGREES, the angle measurement in degrees.
        //
        //    Output, double R8_RADIANS, the angle measurement in radians.
        //
    {
        double value = degrees * Math.PI / 180.0;

        return value;
    }

    public static double r8_relu(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_RELU evaluates the ReLU function of an R8.
        //
        //  Discussion:
        //
        //    An R8 is a double precision real value.
        //
        //    The ReLU function is max(x,0.0).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 January 2019
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double VALUE, the function value.
        //
    {
        double value = x switch
        {
            <= 0.0 => 0.0,
            _ => x
        };

        return value;
    }

    public static double r8_reverse_bytes ( double x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_REVERSE_BYTES reverses the bytes in an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, a value whose bytes are to be reversed.
        //
        //    Output, R8_REVERSE_BYTES, a value with bytes in reverse order;
        //
    {
        byte[] t = BitConverter.GetBytes(x);
        Array.Reverse(t);
        return BitConverter.ToDouble(t);
    }

    public static double r8_secd(double degrees)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SECD returns the secant of an angle given in degrees.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double DEGREES, the angle in degrees.
        //
        //    Output, double R8_SECD, the secant of the angle.
        //
    {
        double radians = Math.PI * (degrees / 180.0);

        double value = 1.0 / Math.Cos(radians);

        return value;
    }

    public static double r8_sech(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SECH evaluates the hyperbolic secant, while avoiding COSH overflow.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the function.
        //
        //    Output, double R8_SECH, the value of the function.
        //
    {
        const double log_huge = 80.0;

        double value = Math.Abs(x) switch
        {
            > log_huge => 0.0,
            _ => 1.0 / Math.Cosh(x)
        };

        return value;
    }

    public static double r8_sigmoid(double l, double b, double m, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    r8_sigmoid evaluates the sigmoid or logistic function.
        //
        //  Discussion:
        //
        //    An R8 is a double value.
        //
        //    The sigmoid function is useful for classification problems in
        //    machine learning.  Its value is always between 0 and 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2019
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double l, the maximum value of the function.  This is often 1.
        //
        //    double b, the cutoff value, where the function equals l/2.
        //    This is often 0.
        //
        //    double m, the slope, which determines the steepness of the curve
        //    and the width of the uncertainty interval.  This is often 1.
        //
        //    double x, the argument.
        //
        //  Output:
        //
        //    double r8_sigmoid, the value.
        //
    {
        double value = l / (1.0 + Math.Exp(-m * (x - b)));

        return value;
    }


    public static void r8_sincos_sum(double a, double b, ref double d, ref double e, ref double f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SINCOS_SUM simplifies a*sin(cx)+b*cos(cx).
        //
        //  Discussion:
        //
        //    The expression
        //      a * sin ( c * x ) + b * cos ( c * x )
        //    can be rewritten as
        //      d * sin ( c * x + e )
        //    or
        //      d * cos ( c * x + f ) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the coefficients in the linear combination.
        //
        //    Output, double &D, &E, &F, the new coefficient, and the shift for
        //    sine or for cosine.
        //
    {
        d = Math.Sqrt(a * a + b * b);
        e = Math.Atan2(b, a);
        f = Math.Atan2(b, a) - Math.PI / 2.0E+00;
        switch (f)
        {
            case < -Math.PI:
                f += 2.0E+00 * Math.PI;
                break;
        }
    }

    public static double r8_sind(double degrees)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SIND returns the sine of an angle given in degrees.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double DEGREES, the angle in degrees.
        //
        //    Output, double R8_SIND, the sine of the angle.
        //
    {
        double radians = Math.PI * (degrees / 180.0);

        double value = Math.Sin(radians);

        return value;
    }

    public static double r8_softplus(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SOFTPLUS evaluates the softplus function of an R8.
        //
        //  Discussion:
        //
        //    An R8 is a double precision real value.
        //
        //    The softplus function is a smoothed (differentiable) version of max(x,0.0).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 September 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double VALUE, the function value.
        //
    {
        double value = x switch
        {
            <= -36.841 => 0.0,
            >= +36.841 => x,
            _ => Math.Log(1.0 + Math.Exp(x))
        };

        return value;
    }

    public static double r8_sqrt_i4(int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SQRT_I4 returns the square root of an I4 as an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the number whose square root is desired.
        //
        //    Output, double R8_SQRT_I4, the value of sqrt(I).
        //
    {
        double value = Math.Sqrt(i);

        return value;
    }


    public static double r8_square(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    r8_square returns the square of an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2019
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double x: the argument.
        //
        //  Output:
        //
        //    double r8_square: the square of x.
        //
    {
        double value = x * x;

        return value;
    }
        
    public static r8 s_to_r8(string s)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_R8 reads an R8 from a string.
        //
        //  Discussion:
        //
        //    This routine will read as many characters as possible until it reaches
        //    the end of the string, or encounters a character which cannot be
        //    part of the real number.
        //
        //    Legal input is:
        //
        //       1 blanks,
        //       2 '+' or '-' sign,
        //       2.5 spaces
        //       3 integer part,
        //       4 decimal point,
        //       5 fraction part,
        //       6 'E' or 'e' or 'D' or 'd', exponent marker,
        //       7 exponent sign,
        //       8 exponent integer part,
        //       9 exponent decimal point,
        //      10 exponent fraction part,
        //      11 blanks,
        //      12 final comma or semicolon.
        //
        //    with most quantities optional.
        //
        //  Example:
        //
        //    S                 R
        //
        //    '1'               1.0
        //    '     1   '       1.0
        //    '1A'              1.0
        //    '12,34,56'        12.0
        //    '  34 7'          34.0
        //    '-1E2ABCD'        -100.0
        //    '-1X2ABCD'        -1.0
        //    ' 2E-1'           0.2
        //    '23.45'           23.45
        //    '-4.2E+2'         -420.0
        //    '17d2'            1700.0
        //    '-14e-2'         -0.14
        //    'e2'              100.0
        //    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string S, the string containing the
        //    data to be read.  Reading will begin at position 1 and
        //    terminate at the end of the string, or when no more
        //    characters can be read to form a legal real.  Blanks,
        //    commas, or other nonnumeric data will, in particular,
        //    cause the conversion to halt.
        //
        //    Output, int *LCHAR, the number of characters read from
        //    the string to form the number, including any terminating
        //    characters such as a trailing comma or blanks.
        //
        //    Output, bool *ERROR, is true if an error occurred.
        //
        //    Output, double S_TO_R8, the real value that was read from the string.
        //
    {
        r8 ret = new() {lchar = -1};
        double rexp;
        const char TAB = (char) 9;

        int nchar = s_len_trim(s);
        int isgn = 1;
        double rtop = 0.0;
        double rbot = 1.0;
        int jsgn = 1;
        int jtop = 0;
        int jbot = 1;
        int ihave = 1;
        int iterm = 0;

        for (;;)
        {
            char c = s[ret.lchar + 1];
            ret.lchar += 1;
            //
            //  Blank or TAB character.
            //
            if (c is ' ' or TAB)
            {
                switch (ihave)
                {
                    case 2:
                        break;
                    case 6:
                    case 7:
                        iterm = 1;
                        break;
                    case > 1:
                        ihave = 11;
                        break;
                }
            }
            else
            {
                switch (c)
                {
                    //
                    //  Comma.
                    //
                    case ',':
                    case ';':
                    {
                        if (ihave != 1)
                        {
                            iterm = 1;
                            ihave = 12;
                            ret.lchar += 1;
                        }

                        break;
                    }
                    //
                    //  Minus sign.
                    //
                    case '-' when ihave == 1:
                        ihave = 2;
                        isgn = -1;
                        break;
                    case '-' when ihave == 6:
                        ihave = 7;
                        jsgn = -1;
                        break;
                    case '-':
                        iterm = 1;
                        break;
                    //
                    //  Plus sign.
                    //
                    case '+' when ihave == 1:
                        ihave = 2;
                        break;
                    case '+' when ihave == 6:
                        ihave = 7;
                        break;
                    case '+':
                        iterm = 1;
                        break;
                    //
                    //  Decimal point.
                    //
                    case '.' when ihave < 4:
                        ihave = 4;
                        break;
                    case '.' when ihave is >= 6 and <= 8:
                        ihave = 9;
                        break;
                    case '.':
                        iterm = 1;
                        break;
                    //
                    default:
                    {
                        switch (char.ToUpper(c))
                        {
                            case 'E':
                            case 'D':
                            {
                                switch (ihave)
                                {
                                    case < 6:
                                        ihave = 6;
                                        break;
                                    default:
                                        iterm = 1;
                                        break;
                                }

                                break;
                            }
                            //
                            default:
                            {
                                switch (ihave)
                                {
                                    case < 11 when c is >= '0' and <= '9':
                                    {
                                        switch (ihave)
                                        {
                                            case <= 2:
                                                ihave = 3;
                                                break;
                                            case 4:
                                                ihave = 5;
                                                break;
                                            case 6:
                                            case 7:
                                                ihave = 8;
                                                break;
                                            case 9:
                                                ihave = 10;
                                                break;
                                        }

                                        int ndig = ch_to_digit(c);

                                        switch (ihave)
                                        {
                                            case 3:
                                                rtop = 10.0 * rtop + ndig;
                                                break;
                                            case 5:
                                                rtop = 10.0 * rtop + ndig;
                                                rbot = 10.0 * rbot;
                                                break;
                                            case 8:
                                                jtop = 10 * jtop + ndig;
                                                break;
                                            case 10:
                                                jtop = 10 * jtop + ndig;
                                                jbot = 10 * jbot;
                                                break;
                                        }

                                        break;
                                    }
                                    //
                                    default:
                                        iterm = 1;
                                        break;
                                }

                                break;
                            }
                        }

                        break;
                    }
                }
            }

            //
            //  If we haven't seen a terminator, and we haven't examined the
            //  entire string, go get the next character.
            //
            if (iterm == 1 || nchar <= ret.lchar + 1)
            {
                break;
            }

        }

        //
        //  If we haven't seen a terminator, and we have examined the
        //  entire string, then we're done, and LCHAR is equal to NCHAR.
        //
        if (iterm != 1 && ret.lchar + 1 == nchar)
        {
            ret.lchar = nchar;
        }

        switch (ihave)
        {
            //
            //  Number seems to have terminated.  Have we got a legal number?
            //  Not if we terminated in states 1, 2, 6 or 7!
            //
            case 1:
            case 2:
            case 6:
            case 7:
                ret.error = true;
                return ret;
        }

        switch (jtop)
        {
            //
            //  Number seems OK.  Form it.
            //
            case 0:
                rexp = 1.0;
                break;
            default:
            {
                switch (jbot)
                {
                    case 1:
                        rexp = Math.Pow(10.0, jsgn * jtop);
                        break;
                    default:
                        rexp = jsgn * jtop;
                        rexp /= jbot;
                        rexp = Math.Pow(10.0, rexp);
                        break;
                }

                break;
            }
        }

        ret.val = isgn * rexp * rtop / rbot;

        return ret;
    }


    public static r8vec s_to_r8vec(string s, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_R8VEC reads an R8VEC from a string.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string S, the string to be read.
        //
        //    Input, int N, the number of values expected.
        //
        //    Output, double RVEC[N], the values read from the string.
        //
        //    Output, bool S_TO_R8VEC, is true if an error occurred.
        //
    {
        r8vec ret = new() {rvec = new double[n]};

        string[] tokens = Helpers.splitStringByWhitespace(s);

        for (int i = 0; i < n; i++)
        {
            ret.rvec[i] = s_to_r8(tokens[i]).val;
        }
            
        return ret;
    }


    public static long r8_nint(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_NINT returns the nearest integer to an R8.
        //
        //  Examples:
        //
        //        X         R8_NINT
        //
        //      1.3         1
        //      1.4         1
        //      1.5         1 or 2
        //      1.6         2
        //      0.0         0
        //     -0.7        -1
        //     -1.1        -1
        //     -1.6        -2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 November 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the value.
        //
        //    Output, int R8_NINT, the nearest integer to X.
        //
    {
        long value = (long) (Math.Abs(x) + 0.5);

        value = x switch
        {
            < 0.0 => -value,
            _ => value
        };

        return value;
    }


    public static void r8utp_print(int n, double[] a, string title)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTP_PRINT prints a R8UTP matrix.
        //
        //  Discussion:
        //
        //    The R8UTP storage format is appropriate for an upper triangular
        //    matrix.  Only the upper triangle of the matrix is stored,
        //    by successive partial columns, in an array of length (N*(N+1))/2,
        //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, double A[(N*(N+1))/2], the matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8utp_print_some(n, a, 1, 1, n, n, title);
    }

    public static void r8utp_print_some(int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UTP_PRINT_SOME prints some of an R8UTP matrix.
        //
        //  Discussion:
        //
        //    The R8UTP storage format is appropriate for an upper triangular
        //    matrix.  Only the upper triangle of the matrix is stored,
        //    by successive partial columns, in an array of length (N*(N+1))/2,
        //    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, double A[(N*(N+1))/2], the matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

        int j2lo;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            int j2hi = j2lo + INCX - 1;
            j2hi = Math.Min(j2hi, n);
            j2hi = Math.Min(j2hi, jhi);

            Console.WriteLine("");
            string cout = "  Col: ";
            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                cout += j.ToString(CultureInfo.InvariantCulture).PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            int i2lo = Math.Max(ilo, 1);
            int i2hi = Math.Min(ihi, n);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  ";
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                for (j = j2lo; j <= j2hi; j++)
                {
                    double aij = i <= j ? a[i - 1 + j * (j - 1) / 2] : 0.0;

                    cout += aij.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double r8_huge()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_HUGE returns a "huge" R8.
        //
        //  Discussion:
        //
        //    The value returned by this function is NOT required to be the
        //    maximum representable R8.  This value varies from machine to machine,
        //    from compiler to compiler, and may cause problems when being printed.
        //    We simply want a "very large" but non-infinite number.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 October 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double R8_HUGE, a "huge" R8 value.
        //
    {
        const double value = 1.0E+30;

        return value;
    }
    public static double r8_zeta(double p)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ZETA estimates the Riemann Zeta function.
        //
        //  Discussion:
        //
        //    For 1 < P, the Riemann Zeta function is defined as:
        //
        //      ZETA ( P ) = Sum ( 1 <= N < oo ) 1 / N^P
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Daniel Zwillinger, editor,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition,
        //    CRC Press, 1996.
        //
        //  Parameters:
        //
        //    Input, double P, the power to which the integers are raised.
        //    P must be greater than 1.  For integral P up to 20, a
        //    precomputed value is returned; otherwise the infinite
        //    sum is approximated.
        //
        //    Output, double R8_ZETA, an approximation to the Riemann
        //    Zeta function.
        //
    {
        double value;

        switch (p)
        {
            case <= 1.0:
                value = r8_huge();
                break;
            case 2.0:
                value = Math.Pow(Math.PI, 2) / 6.0;
                break;
            case 3.0:
                value = 1.2020569032;
                break;
            case 4.0:
                value = Math.Pow(Math.PI, 4) / 90.0;
                break;
            case 5.0:
                value = 1.0369277551;
                break;
            case 6.0:
                value = Math.Pow(Math.PI, 6) / 945.0;
                break;
            case 7.0:
                value = 1.0083492774;
                break;
            case 8.0:
                value = Math.Pow(Math.PI, 8) / 9450.0;
                break;
            case 9.0:
                value = 1.0020083928;
                break;
            case 10.0:
                value = Math.Pow(Math.PI, 10) / 93555.0;
                break;
            case 11.0:
                value = 1.0004941886;
                break;
            case 12.0:
                value = 1.0002460866;
                break;
            case 13.0:
                value = 1.0001227133;
                break;
            case 14.0:
                value = 1.0000612482;
                break;
            case 15.0:
                value = 1.0000305882;
                break;
            case 16.0:
                value = 1.0000152823;
                break;
            case 17.0:
                value = 1.0000076372;
                break;
            case 18.0:
                value = 1.0000038173;
                break;
            case 19.0:
                value = 1.0000019082;
                break;
            case 20.0:
                value = 1.0000009540;
                break;
            default:
            {
                double zsum = 0.0;
                int n = 0;

                for (;;)
                {
                    n += 1;
                    double zsum_old = zsum;
                    zsum += 1.0 / Math.Pow(n, p);
                    if (zsum <= zsum_old)
                    {
                        break;
                    }
                }

                value = zsum;
                break;
            }
        }

        return value;
    }

    public static double r8_beta(double x, double y)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BETA returns the value of the Beta function.
        //
        //  Discussion:
        //
        //    BETA(X,Y) = ( GAMMA(X) * GAMMA(Y) ) / GAMMA(X+Y)
        //
        //    BETA(X,Y) = BETA(Y,X).
        //    BETA(X,Y) = Integral ( 0 <= T <= 1 ) T^(X-1) (1-T)^(Y-1) dT.
        //    BETA(X,Y) = GAMMA(X) * GAMMA(Y) / GAMMA(X+Y)
        //
        //    Both X and Y must be greater than 0.
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
        //    Input, double X, Y, the two parameters that define the Beta function.
        //    X and Y must be greater than 0.
        //
        //    Output, double R8_BETA, the value of the Beta function.
        //
    {
        if (x <= 0.0 || y <= 0.0)
        {
            Console.WriteLine("");
            Console.WriteLine("R8_BETA - Fatal error!");
            Console.WriteLine("  Both X and Y must be greater than 0.");
            return 1.0;
        }

        double value = Math.Exp(
            Helpers.LogGamma(x)
            + Helpers.LogGamma(y)
            - Helpers.LogGamma(x + y));

        return value;
    }



    public static double r8_big()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BIG returns a "big" R8.
        //
        //  Discussion:
        //
        //    The value returned by this function is NOT required to be the
        //    maximum representable R8.
        //    We simply want a "very large" but non-infinite number.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double R8_BIG, a "big" R8 value.
        //
    {
        const double value = 1.0E+30;

        return value;
    }

    public static double r8_chop(int place, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CHOP chops an R8 to a given number of binary places.
        //
        //  Example:
        //
        //    3.875 = 2 + 1 + 1/2 + 1/4 + 1/8.
        //
        //    The following values would be returned for the 'chopped' value of
        //    3.875:
        //
        //    PLACE  Value
        //
        //       1      2
        //       2      3     = 2 + 1
        //       3      3.5   = 2 + 1 + 1/2
        //       4      3.75  = 2 + 1 + 1/2 + 1/4
        //       5+     3.875 = 2 + 1 + 1/2 + 1/4 + 1/8
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PLACE, the number of binary places to preserve.
        //    PLACE = 0 means return the integer part of X.
        //    PLACE = 1 means return the value of X, correct to 1/2.
        //    PLACE = 2 means return the value of X, correct to 1/4.
        //    PLACE = -1 means return the value of X, correct to 2.
        //
        //    Input, double X, the number to be chopped.
        //
        //    Output, double R8_CHOP, the chopped number.
        //
    {
        int temp = (int) r8_log_2(x);
        double fac = Math.Pow(2.0, temp - place + 1);
        double value = (int) (x / fac) * fac;

        return value;
    }

    public static double r8_cosd(double degrees)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_COSD returns the cosine of an angle given in degrees.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double DEGREES, the angle in degrees.
        //
        //    Output, double R8_COSD, the cosine of the angle.
        //
    {
        double radians = Math.PI * (degrees / 180.0);

        double value = Math.Cos(radians);

        return value;
    }

    public static double r8_cot(double angle)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_COT returns the cotangent of an angle.
        //
        //  Discussion:
        //
        //    R8_COT ( THETA ) = COS ( THETA ) / SIN ( THETA )
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
        //    Input, double ANGLE, the angle, in radians.
        //
        //    Output, double R8_COT, the cotangent of the angle.
        //
    {
        double value = Math.Cos(angle) / Math.Sin(angle);

        return value;
    }

    public static double r8_cotd(double degrees)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_COTD returns the cotangent of an angle given in degrees.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double DEGREES, the angle in degrees.
        //
        //    Output, double R8_COTD, the cotangent of the angle.
        //
    {
        double radians = Math.PI * (degrees / 180.0);

        double value = Math.Cos(radians) / Math.Sin(radians);

        return value;
    }


    public static double r8_csc(double theta)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CSC returns the cosecant of X.
        //
        //  Discussion:
        //
        //    R8_CSC ( THETA ) = 1.0 / SIN ( THETA )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double THETA, the angle, in radians, whose cosecant is desired.
        //    It must be the case that SIN ( THETA ) is not zero.
        //
        //    Output, double R8_CSC, the cosecant of THETA.
        //
    {
        double value = Math.Sin(theta);

        switch (value)
        {
            case 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("R8_CSC - Fatal error!");
                Console.WriteLine("  Cosecant undefined for THETA = " + theta + "");
                return 1;
            default:
                value = 1.0 / value;

                return value;
        }
    }
        
    public static double r8_modp(double x, double y)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_MODP returns the nonnegative remainder of R8 division.
        //
        //  Discussion:
        //
        //    If
        //      REM = R8_MODP ( X, Y )
        //      RMULT = ( X - REM ) / Y
        //    then
        //      X = Y * RMULT + REM
        //    where REM is always nonnegative.
        //
        //    The MOD function computes a result with the same sign as the
        //    quantity being divided.  Thus, suppose you had an angle A,
        //    and you wanted to ensure that it was between 0 and 360.
        //    Then mod(A,360.0) would do, if A was positive, but if A
        //    was negative, your result would be between -360 and 0.
        //
        //    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.
        //
        //  Example:
        //
        //        I         J     MOD R8_MODP  R8_MODP Factorization
        //
        //      107        50       7       7    107 =  2 *  50 + 7
        //      107       -50       7       7    107 = -2 * -50 + 7
        //     -107        50      -7      43   -107 = -3 *  50 + 43
        //     -107       -50      -7      43   -107 =  3 * -50 + 43
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the number to be divided.
        //
        //    Input, double Y, the number that divides X.
        //
        //    Output, double R8_MODP, the nonnegative remainder when X is divided by Y.
        //
    {
        switch (y)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_MODP - Fatal error!");
                Console.WriteLine("  R8_MODP ( X, Y ) called with Y = " + y + "");
                return 1;
        }

        double value = x - (int) (x / y) * y;

        switch (value)
        {
            case < 0.0:
                value += Math.Abs(y);
                break;
        }

        return value;
    }

    public static double r8_cscd(double degrees)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CSCD returns the cosecant of an angle given in degrees.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double DEGREES, the angle in degrees.
        //
        //    Output, double R8_CSCD, the cosecant of the angle.
        //
    {
        double radians = Math.PI * (degrees / 180.0);

        double value = 1.0 / Math.Sin(radians);

        return value;
    }

    public static double r8_cube_root(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CUBE_ROOT returns the cube root of an R8.
        //
        //  Discussion:
        //
        //    This routine is designed to avoid the possible problems that can occur
        //    when formulas like 0.0^(1/3) or (-1.0)^(1/3) are to be evaluated.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input double X, the number whose cube root is desired.
        //
        //    Output, double R8_CUBE_ROOT, the cube root of X.
        //
    {
        double value = x switch
        {
            > 0.0 => Math.Pow(x, 1.0 / 3.0),
            0.0 => 0.0,
            _ => -Math.Pow(Math.Abs(x), 1.0 / 3.0)
        };

        return value;
    }

    public static double r8_degrees(double radians)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_DEGREES converts an angle from radian to degree measure.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double RADIANS, the angle measurement in radians.
        //
        //    Output, double R8_DEGREES, the angle measurement in degrees.
        //
    {
        double value = radians * 180.0 / Math.PI;

        return value;
    }

    public static double r8_diff(double x, double y, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_DIFF computes (X-Y) to a specified accuracy.
        //
        //  Discussion:
        //
        //    The user controls how many binary digits of accuracy
        //    are to be used.
        //
        //    N determines the accuracy of the value.  If N = 10,
        //    for example, only 11 binary places will be used in the arithmetic.
        //    In general, only N+1 binary places will be used.
        //
        //    N may be zero.  However, a negative value of N should
        //    not be used, since this will cause both X and Y to look like 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the two values whose difference is desired.
        //
        //    Input, int N, the number of binary digits to use.
        //
        //    Output, double R8_DIFF, the value of X-Y.
        //
    {
        double value;

        if (Math.Abs(x - y) <= typeMethods.r8_epsilon())
        {
            value = 0.0;
            return value;
        }

        double pow2 = Math.Pow(2.0, n);
        //
        //  Compute the magnitude of X and Y, and take the larger of the
        //  two.  At least one of the two values is not zero//
        //
        double size = Math.Max(Math.Abs(x), Math.Abs(y));
        //
        //  Make normalized copies of X and Y.  One of the two values will
        //  actually be equal to 1.
        //
        double cx = x / size;
        double cy = y / size;
        //
        //  Here's where rounding comes in.  We know that the larger of the
        //  the two values equals 1.  We multiply both values by 2^N,
        //  where N+1 is the number of binary digits of accuracy we want
        //  to use, truncate the values, and divide back by 2^N.
        //
        cx = (int) (cx * pow2 + 0.5 * r8_sign(cx)) / pow2;
        cy = (int) (cy * pow2 + 0.5 * r8_sign(cy)) / pow2;
        //
        //  Take the difference now.
        //
        value = cx - cy;
        //
        //  Undo the scaling.
        //
        value *= size;

        return value;
    }

    public static int r8_digit(double x, int idigit)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_DIGIT returns a particular decimal digit of an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the number whose IDIGIT-th decimal digit is desired.
        //    Note that if X is zero, all digits will be returned as 0.
        //
        //    Input, int IDIGIT, the position of the desired decimal digit.
        //    A value of 1 means the leading digit, a value of 2 the second digit
        //    and so on.
        //
        //    Output, int R8_DIGIT, the value of the IDIGIT-th decimal digit of X.
        //
    {
        int digit;
        int i;
        int ival = 0;

        switch (x)
        {
            case 0.0:
                digit = 0;
                return digit;
        }

        switch (idigit)
        {
            case <= 0:
                digit = 0;
                return digit;
        }

        //
        //  Force X to lie between 1 and 10.
        //
        x = Math.Abs(x);

        while (x < 1.0)
        {
            x *= 10.0;
        }

        while (10.0 <= x)
        {
            x /= 10.0;
        }

        for (i = 1; i <= idigit; i++)
        {
            ival = (int) x;
            x = (x - ival) * 10.0;
        }

        digit = ival;

        return digit;
    }

    public static double r8_divide_i4(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_DIVIDE_I4 returns an I4 fraction as an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, J, the numerator and denominator.
        //
        //    Output, double R8_DIVIDE_I4, the value of (I/J).
        //
    {
        double value = i / (double) j;

        return value;
    }

    public static double r8_e()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_E returns the value of the base of the natural logarithm system.
        //
        //  Definition:
        //
        //    E = Limit ( N -> +oo ) ( 1 + 1 / N )^N
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
        //    Output, double R8_E, the base of the natural logarithm system.
        //
    {
        const double r8_e_save = 2.718281828459045235360287;

        return r8_e_save;
    }

    public const double _r8_epsilon = 2.220446049250313E-016;
    public static double r8_epsilon()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_EPSILON returns the R8 roundoff unit.
        //
        //  Discussion:
        //
        //    The roundoff unit is a number R which is a power of 2 with the
        //    property that, to the precision of the computer's arithmetic,
        //      1 < 1 + R
        //    but
        //      1 = ( 1 + R / 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double R8_EPSILON, the R8 round-off unit.
        //
    {
        return _r8_epsilon ;//typeMethods.r8_epsilon();
    }

    public class r8EpsilonData
    {
        public double value;

        public r8EpsilonData()
        {
            value = 0;
        }
    }
        
    public static double r8_epsilon_compute(ref r8EpsilonData data)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_EPSILON_COMPUTE computes the R8 roundoff unit.
        //
        //  Discussion:
        //
        //    The roundoff unit is a number R which is a power of 2 with the
        //    property that, to the precision of the computer's arithmetic,
        //      1 < 1 + R
        //    but
        //      1 = ( 1 + R / 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double R8_EPSILON_COMPUTE, the R8 round-off unit.
        //
    {
        switch (data.value)
        {
            case 0.0:
            {
                const double one = 1;

                data.value = one;
                double temp = data.value / 2.0;
                double test = r8_add(one, temp);

                while (one < test)
                {
                    data.value = temp;
                    temp = data.value / 2.0;
                    test = r8_add(one, temp);
                }

                break;
            }
        }

        return data.value;
    }

    public static double r8_erf(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ERF evaluates the error function ERF(X).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 May 2007
        //
        //  Author:
        //
        //    Original FORTRAN77 version by William Cody.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    William Cody,
        //    "Rational Chebyshev approximations for the error function",
        //    Mathematics of Computation, 
        //    1969, pages 631-638.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the error function.
        //
        //    Output, double R8_ERF, the value of the error function.
        //
    {
        double[] a =  {
                3.16112374387056560,
                1.13864154151050156E+02,
                3.77485237685302021E+02,
                3.20937758913846947E+03,
                1.85777706184603153E-01
            }
            ;
        double[] b =  {
                2.36012909523441209E+01,
                2.44024637934444173E+02,
                1.28261652607737228E+03,
                2.84423683343917062E+03
            }
            ;
        double[] c =  {
                5.64188496988670089E-01,
                8.88314979438837594,
                6.61191906371416295E+01,
                2.98635138197400131E+02,
                8.81952221241769090E+02,
                1.71204761263407058E+03,
                2.05107837782607147E+03,
                1.23033935479799725E+03,
                2.15311535474403846E-08
            }
            ;
        double[] d =  {
                1.57449261107098347E+01,
                1.17693950891312499E+02,
                5.37181101862009858E+02,
                1.62138957456669019E+03,
                3.29079923573345963E+03,
                4.36261909014324716E+03,
                3.43936767414372164E+03,
                1.23033935480374942E+03
            }
            ;
        double erfx;
        int i;
        double[] p =  {
                3.05326634961232344E-01,
                3.60344899949804439E-01,
                1.25781726111229246E-01,
                1.60837851487422766E-02,
                6.58749161529837803E-04,
                1.63153871373020978E-02
            }
            ;
        double[] q =  {
                2.56852019228982242,
                1.87295284992346047,
                5.27905102951428412E-01,
                6.05183413124413191E-02,
                2.33520497626869185E-03
            }
            ;
        const double sqrpi = 0.56418958354775628695;
        const double thresh = 0.46875;
        const double xbig = 26.543;
        double xden;
        double xnum;
        const double xsmall = 1.11E-16;
        double xsq;

        double xabs = Math.Abs(x);
        //
        //  Evaluate ERF(X) for |X| <= 0.46875.
        //
        if (xabs <= thresh)
        {
            if (xsmall < xabs)
            {
                xsq = xabs * xabs;
            }
            else
            {
                xsq = 0.0;
            }

            xnum = a[4] * xsq;
            xden = xsq;
            for (i = 0; i < 3; i++)
            {
                xnum = (xnum + a[i]) * xsq;
                xden = (xden + b[i]) * xsq;
            }

            erfx = x * (xnum + a[3]) / (xden + b[3]);
        }
        else
        {
            double del;
            switch (xabs)
            {
                //
                //  Evaluate ERFC(X) for 0.46875 <= |X| <= 4.0.
                //
                case <= 4.0:
                {
                    xnum = c[8] * xabs;
                    xden = xabs;
                    for (i = 0; i < 7; i++)
                    {
                        xnum = (xnum + c[i]) * xabs;
                        xden = (xden + d[i]) * xabs;
                    }

                    erfx = (xnum + c[7]) / (xden + d[7]);
                    xsq = (int)(xabs * 16.0 / 16.0);
                    del = (xabs - xsq) * (xabs + xsq);
                    erfx = Math.Exp(-xsq * xsq) * Math.Exp(-del) * erfx;

                    erfx = x switch
                    {
                        < 0.0 => -erfx,
                        _ => 0.5 - erfx + 0.5
                    };

                    break;
                }
                //
                default:
                {
                    if (xbig <= xabs)
                    {
                        erfx = x switch
                        {
                            > 0.0 => 1.0,
                            _ => -1.0
                        };
                    }
                    else
                    {
                        xsq = 1.0 / (xabs * xabs);

                        xnum = p[5] * xsq;
                        xden = xsq;
                        for (i = 0; i < 4; i++)
                        {
                            xnum = (xnum + p[i]) * xsq;
                            xden = (xden + q[i]) * xsq;
                        }

                        erfx = xsq * (xnum + p[4]) / (xden + q[4]);
                        erfx = (sqrpi - erfx) / xabs;
                        xsq = (int)(xabs * 16.0 / 16.0);
                        del = (xabs - xsq) * (xabs + xsq);
                        erfx = Math.Exp(-xsq * xsq) * Math.Exp(-del) * erfx;

                        erfx = x switch
                        {
                            < 0.0 => -erfx,
                            _ => 0.5 - erfx + 0.5
                        };
                    }

                    break;
                }
            }
        }

        return erfx;
    }

    public static double r8_erf_inverse(double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ERF_INVERSE inverts the error function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double Y, the value of the error function.
        //
        //    Output, double R8_ERF_INVERSE, the value X such that ERF(X) = Y.
        //
    {
        double z = (y + 1.0) / 2.0;

        double x = CDF.normal_01_cdf_inv(z);

        double value = x / Math.Sqrt(2.0);

        return value;
    }

    public static double r8_euler_constant ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_EULER_CONSTANT returns the value of the Euler-Mascheroni constant.
        //
        //  Discussion:
        //
        //    The Euler-Mascheroni constant is often denoted by a lower-case gamma.
        //
        //      gamma = limit ( N -> +oo )
        //        ( sum ( 1 <= I <= N ) 1 / I ) - log ( N )
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
        //    Output, double R8_EULER_CONSTANT, the value of the Euler-Mascheroni constant.
        //
    {
        const double value = 0.577215664901532860606512090082402431042;

        return value;
    }
        
    public static double r8_exp(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_EXP computes the exponential function, avoiding overflow and underflow.
        //
        //  Discussion:
        //
        //    For arguments of very large magnitude, the evaluation of the
        //    exponential function can cause computational problems.  Some languages
        //    and compilers may return an infinite value or a "Not-a-Number".  
        //    An alternative, when dealing with a wide range of inputs, is simply
        //    to truncate the calculation for arguments whose magnitude is too large.
        //    Whether this is the right or convenient approach depends on the problem
        //    you are dealing with, and whether or not you really need accurate
        //    results for large magnitude inputs, or you just want your code to
        //    stop crashing.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the exponential function.
        //
        //    Output, double R8_EXP, the value of exp ( X ).
        //
    {
        const double r8_big = 1.0E+30;
        const double r8_log_max = +69.0776;
        const double r8_log_min = -69.0776;

        double value = x switch
        {
            <= r8_log_min => 0.0,
            < r8_log_max => Math.Exp(x),
            _ => r8_big
        };

        return value;
    }

    public static double r8_fraction(int i, int j)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_FRACTION uses real arithmetic on an integer ratio.
        //
        //  Discussion:
        //
        //    Given integer variables I and J, both FORTRAN and C will evaluate 
        //    an expression such as "I/J" using what is called "integer division",
        //    with the result being an integer.  It is often convenient to express
        //    the parts of a fraction as integers but expect the result to be computed
        //    using real arithmetic.  This function carries out that operation.
        //
        //  Example:
        //
        //       I     J   I/J  R8_FRACTION
        //
        //       1     2     0  0.5
        //       7     4     1  1.75
        //       8     4     2  2.00
        //       9     4     2  2.25
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, J, the arguments.
        //
        //    Output, double R8_FRACTION, the value of the ratio.
        //
    {
        double value = i / (double) j;

        return value;
    }


    public static double r8_mach(int i)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_MACH returns double precision real machine constants.
        //
        //  Discussion:
        //
        //    Assuming that the internal representation of a double precision real
        //    number is in base B, with T the number of base-B digits in the mantissa,
        //    and EMIN the smallest possible exponent and EMAX the largest possible 
        //    exponent, then
        //
        //      R8_MACH(1) = B^(EMIN-1), the smallest positive magnitude.
        //      R8_MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
        //      R8_MACH(3) = B^(-T), the smallest relative spacing.
        //      R8_MACH(4) = B^(1-T), the largest relative spacing.
        //      R8_MACH(5) = log10(B).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 April 2007
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Phyllis Fox, Andrew Hall, Norman Schryer,
        //    Algorithm 528:
        //    Framework for a Portable Library,
        //    ACM Transactions on Mathematical Software,
        //    Volume 4, Number 2, June 1978, page 176-188.
        //
        //  Parameters:
        //
        //    Input, int I, chooses the parameter to be returned.
        //    1 <= I <= 5.
        //
        //    Output, double R8_MACH, the value of the chosen parameter.
        //
    {
        double value;

        switch (i)
        {
            case 1:
                value = 4.450147717014403E-308;
                break;
            case 2:
                value = 8.988465674311579E+307;
                break;
            case 3:
                value = 1.110223024625157E-016;
                break;
            case 4:
                value = 2.220446049250313E-016;
                break;
            case 5:
                value = 0.301029995663981E+000;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("R8_MACH - Fatal error!");
                Console.WriteLine("  The input argument I is out of bounds.");
                Console.WriteLine("  Legal values satisfy 1 <= I <= 5.");
                Console.WriteLine("  I = " + i + "");
                value = 0.0;
                break;
        }

        return value;
    }


    public static double r8_cas(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CAS returns the "casine" of an R8.
        //
        //  Discussion:
        //
        //    The "casine", used in the discrete Hartley transform, is abbreviated
        //    CAS(X), and defined by:
        //
        //      CAS(X) = cos ( X ) + sin( X )
        //             = sqrt ( 2 ) * sin ( X + pi/4 )
        //             = sqrt ( 2 ) * cos ( X - pi/4 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the number whose casine is desired.
        //
        //    Output, double R8_CAS, the casine of X, which will be between
        //    plus or minus the square root of 2.
        //
    {
        double value = Math.Cos(x) + Math.Sin(x);

        return value;
    }

    public static double r8_ceiling(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CEILING rounds an R8 up to the nearest integral R8.
        //
        //  Example:
        //
        //    X        R8_CEILING(X)
        //
        //   -1.1      -1.0
        //   -1.0      -1.0
        //   -0.9       0.0
        //   -0.1       0.0
        //    0.0       0.0
        //    0.1       1.0
        //    0.9       1.0
        //    1.0       1.0
        //    1.1       2.0
        //    2.9       3.0
        //    3.0       3.0
        //    3.14159   4.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the number whose ceiling is desired.
        //
        //    Output, double R8_CEILING, the ceiling of X.
        //
    {
        double value = (int) x;

        if (value < x)
        {
            value += 1.0;
        }

        return value;
    }


    public static double r8_choose(int n, int k)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
        //
        //  Discussion:
        //
        //    The value is calculated in such a way as to avoid overflow and
        //    roundoff.  The calculation is done in R8 arithmetic.
        //
        //    The formula used is:
        //
        //      C(N,K) = N! / ( K! * (N-K)! )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    ML Wolfson, HV Wright,
        //    Algorithm 160:
        //    Combinatorial of M Things Taken N at a Time,
        //    Communications of the ACM,
        //    Volume 6, Number 4, April 1963, page 161.
        //
        //  Parameters:
        //
        //    Input, int N, K, the values of N and K.
        //
        //    Output, double R8_CHOOSE, the number of combinations of N
        //    things taken K at a time.
        //
    {
        int mn;
        int mx;
        double value;

        if (k < n - k)
        {
            mn = k;
            mx = n - k;
        }
        else
        {
            mn = n - k;
            mx = k;
        }

        switch (mn)
        {
            case < 0:
                value = 0.0;
                break;
            case 0:
                value = 1.0;
                break;
            default:
            {
                value = mx + 1;

                int i;
                for (i = 2; i <= mn; i++)
                {
                    value = value * (mx + i) / i;
                }

                break;
            }
        }

        return value;
    }

    public static double r8_mop(int i)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_MOP returns the I-th power of -1 as an R8 value.
        //
        //  Discussion:
        //
        //    An R8 is an double value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the power of -1.
        //
        //    Output, double R8_MOP, the I-th power of -1.
        //
    {
        double value = (i % 2) switch
        {
            0 => 1.0,
            _ => -1.0
        };

        return value;
    }

    public static void r8_to_cfrac(double r, int n, ref int[] a, ref int[] p, ref int[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TO_CFRAC converts a double value to a continued fraction.
        //
        //  Discussion:
        //
        //    The routine is given a double number R.  It computes a sequence of
        //    continued fraction approximations to R, returning the results as
        //    simple fractions of the form P(I) / Q(I).
        //
        //  Example:
        //
        //    X = 2 * PI
        //    N = 7
        //
        //    A = [ *, 6,  3,  1,  1,   7,   2,    146,      3 ]
        //    P = [ 1, 6, 19, 25, 44, 333, 710, 103993, 312689 ]
        //    Q = [ 0, 1,  3,  4,  7,  53, 113,  16551,  49766 ]
        //
        //    (This ignores roundoff error, which will cause later terms to differ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Norman Richert,
        //    Strang's Strange Figures,
        //    American Mathematical Monthly,
        //    Volume 99, Number 2, February 1992, pages 101-107.
        //
        //  Parameters:
        //
        //    Input, double R, the double value.
        //
        //    Input, int N, the number of convergents to compute.
        //
        //    Output, int A[N+1], the partial quotients.
        //
        //    Output, int P[N+2], Q[N+2], the numerators and denominators
        //    of the continued fraction approximations.
        //
    {
        int i;

        switch (r)
        {
            case 0.0:
            {
                for (i = 0; i <= n; i++)
                {
                    a[i] = 0;
                }

                for (i = 0; i <= n + 1; i++)
                {
                    p[i] = 0;
                }

                for (i = 0; i <= n + 1; i++)
                {
                    q[i] = 0;
                }

                return;
            }
        }

        double[] x = new double[n + 1];

        double r_copy = Math.Abs(r);

        p[0] = 1;
        q[0] = 0;

        p[1] = (int)r_copy;
        q[1] = 1;
        x[0] = r_copy;
        a[0] = (int)x[0];

        for (i = 1; i <= n; i++)
        {
            x[i] = 1.0 / (x[i - 1] - a[i - 1]);
            a[i] = (int)x[i];
            p[i + 1] = a[i] * p[i] + p[i - 1];
            q[i + 1] = a[i] * q[i] + q[i - 1];
        }

        switch (r)
        {
            case < 0.0:
            {
                for (i = 0; i <= n + 1; i++)
                {
                    p[i] = -p[i];
                }

                break;
            }
        }
    }

    public static void r8_to_dec(double dval, int dec_digit, ref int mantissa, ref int exponent)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TO_DEC converts a double quantity to a decimal representation.
        //
        //  Discussion:
        //
        //    Given the double value DVAL, the routine computes integers
        //    MANTISSA and EXPONENT so that it is approximatelytruethat:
        //
        //      DVAL = MANTISSA * 10 ** EXPONENT
        //
        //    In particular, only DEC_DIGIT digits of DVAL are used in constructing the
        //    representation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 July 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double DVAL, the value whose decimal representation
        //    is desired.
        //
        //    Input, int DEC_DIGIT, the number of decimal digits to use.
        //
        //    Output, int &MANTISSA, &EXPONENT, the approximate decimal 
        //    representation of DVAL.
        //
    {
        switch (dval)
        {
            //
            //  Special cases.
            //
            case 0.0:
                mantissa = 0;
                exponent = 0;
                return;
        }

        //
        //  Factor DVAL = MANTISSA_DOUBLE * 10^EXPONENT
        //
        double mantissa_double = dval;
        exponent = 0;
        //
        //  Now normalize so that 
        //  10^(DEC_DIGIT-1) <= ABS(MANTISSA_DOUBLE) < 10^(DEC_DIGIT)
        //
        double ten1 = Math.Pow(10.0, dec_digit - 1);
        double ten2 = 10.0 * ten1;

        while (Math.Abs(mantissa_double) < ten1)
        {
            mantissa_double *= 10.0;
            exponent -= 1;
        }

        while (ten2 <= Math.Abs(mantissa_double))
        {
            mantissa_double /= 10.0;
            exponent += 1;
        }

        //
        //  MANTISSA is the integer part of MANTISSA_DOUBLE, rounded.
        //
        mantissa = (int)r8_nint(mantissa_double);
        //
        //  Now divide out any factors of ten from MANTISSA.
        //
        if (mantissa == 0)
        {
            return;
        }

        while (10 * (mantissa / 10) == mantissa)
        {
            mantissa /= 10;
            exponent += 1;
        }

    }

    public static void r8_to_rat(double a, int ndig, ref int iatop, ref int iabot)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TO_RAT converts a real value to a rational value.
        //
        //  Discussion:
        //
        //    The rational value (IATOP/IABOT) is essentially computed by truncating
        //    the decimal representation of the real value after a given number of
        //    decimal digits.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 June 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the real value to be converted.
        //
        //    Input, int NDIG, the number of decimal digits used.
        //
        //    Output, int &IATOP, &IABOT, the numerator and denominator
        //    of the rational value that approximates A.
        //
    {
        int i;

        double factor = Math.Pow(10.0, ndig);

        switch (ndig)
        {
            case > 0:
            {
                iabot = 1;
                for (i = 1; i <= ndig; i++)
                {
                    iabot *= 10;
                }

                iatop = 1;
                break;
            }
            default:
            {
                iabot = 1;
                iatop = 1;
                for (i = 1; i <= -ndig; i++)
                {
                    iatop *= 10;
                }

                break;
            }
        }

        iatop = (int)r8_nint(a * factor) * iatop;
        //
        //  Factor out the greatest common factor.
        //
        int itemp = i4_gcd(iatop, iabot);

        iatop /= itemp;
        iabot /= itemp;
    }

    public static void r8_to_dhms(double r, ref int d, ref int h, ref int m, ref int s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TO_DHMS converts an R8 day value into days, hours, minutes, seconds.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, a real number representing a time period measured in days.
        //
        //    Output, int &D, &H, &M, &S, the equivalent number of days, hours,
        //    minutes and seconds.
        //
    {
        int sign;

        switch (r)
        {
            case >= 0.0:
                sign = 1;
                break;
            default:
                sign = -1;
                r = -r;
                break;
        }

        d = (int) r;

        r -= d;
        r = 24.0 * r;
        h = (int) r;

        r -= h;
        r = 60.0 * r;
        m = (int) r;

        r -= m;
        r = 60.0 * r;
        s = (int) r;

        switch (sign)
        {
            case -1:
                d = -d;
                h = -h;
                m = -m;
                s = -s;
                break;
        }
    }

    public static int r8_to_i4(double xmin, double xmax, double x, int ixmin, int ixmax)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TO_I4 maps real X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].
        //
        //  Discussion:
        //
        //    IX := IXMIN + ( IXMAX - IXMIN ) * ( X - XMIN ) / ( XMAX - XMIN )
        //    IX := min ( IX, max ( IXMIN, IXMAX ) )
        //    IX := max ( IX, min ( IXMIN, IXMAX ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double XMIN, XMAX, the real range.  XMAX and XMIN must not be
        //    equal.  It is not necessary that XMIN be less than XMAX.
        //
        //    Input, double X, the real number to be converted.
        //
        //    Input, int IXMIN, IXMAX, the allowed range of the output
        //    variable.  IXMAX corresponds to XMAX, and IXMIN to XMIN.
        //    It is not necessary that IXMIN be less than IXMAX.
        //
        //    Output, int R8_TO_I4, the value in the range [IXMIN,IXMAX] that
        //    corresponds to X.
        //
    {
        if (Math.Abs(xmax - xmin) <= typeMethods.r8_epsilon())
        {
            Console.WriteLine("");
            Console.WriteLine("R8_TO_I4 - Fatal error!");
            Console.WriteLine("  XMAX = XMIN, making a zero divisor.");
            Console.WriteLine("  XMAX = " + xmax + "");
            Console.WriteLine("  XMIN = " + xmin + "");
            return 1;
        }

        double temp = ((xmax - x) * ixmin
                       + (x - xmin) * ixmax)
                      / (xmax - xmin);

        switch (temp)
        {
            case >= 0.0:
                temp += 0.5;
                break;
            default:
                temp -= 0.5;
                break;
        }

        int ix = (int) temp;

        return ix;
    }

    public static double r8_sum(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SUM returns the sum of two R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the quantities to add.
        //
        //    Output, double R8_SUM, the sum of X and Y.
        //
    {
        double value = x + y;

        return value;
    }

    public static void r8_swap(ref double x, ref double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SWAP switches two R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, double &X, &Y.  On output, the values of X and
        //    Y have been interchanged.
        //
    {
        (x, y) = (y, x);
    }

    public static void r8_swap3(ref double x, ref double y, ref double z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SWAP3 swaps three R8's.
        //
        //  Example:
        //
        //    Input:
        //
        //      X = 1, Y = 2, Z = 3
        //
        //    Output:
        //
        //      X = 2, Y = 3, Z = 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, double &X, &Y, &Z, three values to be swapped.
        //
    {
        double w = x;
        x = y;
        y = z;
        z = w;
    }

    public static double r8_tand(double degrees)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TAND returns the tangent of an angle given in degrees.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double DEGREES, the angle in degrees.
        //
        //    Output, double R8_TAND, the tangent of the angle.
        //
    {
        double radians = Math.PI * (degrees / 180.0);

        double value = Math.Sin(radians) / Math.Cos(radians);

        return value;
    }

    public static double r8_tiny()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TINY returns a "tiny" R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double R8_TINY, a "tiny" R8 value.
        //
    {
        const double value = 0.4450147717014E-307;

        return value;
    }


    public static int r8_to_bin_even(int nbin, double a, double b, double c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TO_BIN_EVEN determines the appropriate "bin" for C in [A,B].
        //
        //  Discussion:
        //
        //    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
        //    An initial bin takes everything less than A, and a final bin takes
        //    everything greater than B.
        //
        //  Example:
        //
        //    NBIN = 7, A = 5, B = 15
        //
        //    C   BIN
        //
        //    1    1
        //    3    1
        //    4.9  1
        //    5    2
        //    6    2
        //    7    3
        //    8    3
        //    9.5  4
        //   13    6
        //   14    6
        //   15    6
        //   15.1  7
        //   99    7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NBIN, the number of bins.  NBIN is normally
        //    at least 3.  If NBIN is 1 or 2, then everything is assigned to bin 1.
        //
        //    Input, double A, B, the lower and upper limits of the bin
        //    interval.  While A is expected to be less than B, the code should
        //    return useful results if A is actually greater than B.
        //
        //    Input, double C, a value to be placed in a bin.
        //
        //    Output, inte R8_TO_BIN_EVEN, the index of the bin to which C is
        //    assigned.
        //
    {
        double a2;
        double b2;
        int bin;
        bool swap;
        switch (nbin)
        {
            //
            //  Take care of special cases.
            //
            case < 1:
                bin = 0;
                return bin;
            case 1:
            case 2:
                bin = 1;
                return bin;
        }

        if (Math.Abs(b - a) <= typeMethods.r8_epsilon())
        {
            bin = 0;
            return bin;
        }

        //
        //  If the limits are descending, then we switch them now, and
        //  unswitch the results at the end.
        //
        if (a < b)
        {
            swap = false;
            a2 = a;
            b2 = b;
        }
        else
        {
            swap = true;
            a2 = b;
            b2 = a;
        }

        //
        //  Compute the bin.
        //
        if (c < a2)
        {
            bin = 1;
        }
        else if (Math.Abs(c - a2) <= typeMethods.r8_epsilon())
        {
            bin = 2;
        }
        else if (Math.Abs(c - b2) <= typeMethods.r8_epsilon())
        {
            bin = nbin - 1;
        }
        else if (b2 < c)
        {
            bin = nbin;
        }
        else
        {
            bin = 2 + (int) ((nbin - 2) * (c - a2) / (b2 - a2));
            bin = Math.Max(bin, 2);
            bin = Math.Min(bin, nbin - 1);
        }

        bin = swap switch
        {
            //
            //  Reverse the switching.
            //
            true => nbin + 1 - bin,
            _ => bin
        };

        return bin;
    }

    public static double r8_to_r8_discrete(double r, double rmin, double rmax, int nr)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TO_R8_DISCRETE maps R to RD in [RMIN, RMAX] with NR possible values.
        //
        //  Discussion:
        //
        //    if ( R < RMIN ) then
        //      RD = RMIN
        //    else if ( RMAX < R ) then
        //      RD = RMAX
        //    else
        //      T = nint ( ( NR - 1 ) * ( R - RMIN ) / ( RMAX - RMIN ) )
        //      RD = RMIN + T * ( RMAX - RMIN ) / real ( NR - 1 )
        //
        //    In the special case where NR = 1, when
        //
        //      XD = 0.5 * ( RMAX + RMIN )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the number to be converted.
        //
        //    Input, double RMAX, RMIN, the maximum and minimum
        //    values for RD.
        //
        //    Input, int NR, the number of allowed values for XD.
        //    NR should be at least 1.
        //
        //    Output, double RD, the corresponding discrete value.
        //
    {
        double rd;
        switch (nr)
        {
            //
            //  Check for errors.
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("R8_TO_R8_DISCRETE - Fatal error!");
                Console.WriteLine("  NR = " + nr + "");
                Console.WriteLine("  but NR must be at least 1.");
                return 1;
            case 1:
                rd = 0.5 * (rmin + rmax);
                return rd;
        }

        if (Math.Abs(rmax - rmin) <= typeMethods.r8_epsilon())
        {
            rd = rmax;
            return rd;
        }

        int f = (int)(nr * (rmax - r) / (rmax - rmin));
        f = Math.Max(f, 0);
        f = Math.Min(f, nr);

        rd = (f * rmin
              + (nr - f) * rmax)
             / nr;

        return rd;
    }

    public static void r8_unswap3(ref double x, ref double y, ref double z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_UNSWAP3 unswaps three R8's.
        //
        //  Example:
        //
        //    Input:
        //
        //      X = 2, Y = 3, Z = 1
        //
        //    Output:
        //
        //      X = 1, Y = 2, Z = 3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, double &X, &Y, &Z, three values to be swapped.
        //
    {
        double w = z;
        z = y;
        y = x;
        x = w;
    }


    public static void r8slmat_print(int m, int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SLMAT_PRINT prints a strict lower triangular R8MAT.
        //
        //  Example:
        //
        //    M = 5, N = 5
        //    A = (/ 21, 31, 41, 51, 32, 42, 52, 43, 53, 54 /)
        //
        //    21
        //    31 32
        //    41 42 43
        //    51 52 53 54
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows in A.
        //
        //    Input, int N, the number of columns in A.
        //
        //    Input, double A[*], the M by N matrix.  Only the strict
        //    lower triangular elements are stored, in column major order.
        //
        //    Input, string TITLE, a title.
        //
    {
        int jlo;

        Console.WriteLine("");
        Console.WriteLine(title + "");

        int jmax = Math.Min(n, m - 1);

        const int nn = 5;

        for (jlo = 1; jlo <= jmax; jlo += nn)
        {
            int jhi = Math.Min(jlo + nn - 1, Math.Min(m - 1, jmax));
            Console.WriteLine("");
            string cout = "  Col   ";
            int j;
            for (j = jlo; j <= jhi; j++)
            {
                cout += j.ToString(CultureInfo.InvariantCulture).PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            int i;
            for (i = jlo + 1; i <= m; i++)
            {
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(5) + ":";
                jhi = Math.Min(jlo + nn - 1, Math.Min(i - 1, jmax));
                for (j = jlo; j <= jhi; j++)
                {
                    int indx = (j - 1) * m + i - j * (j + 1) / 2;
                    cout += " " + a[indx - 1].ToString(CultureInfo.InvariantCulture).PadLeft(12);
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void r8_mant(double x, ref int s, ref double r, ref int l)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_MANT computes the "mantissa" or "fraction part" of an R8.
        //
        //  Discussion:
        //
        //    X = S * R * 2^L
        //
        //    S is +1 or -1,
        //    R is a real between 1.0 and 2.0,
        //    L is an integer.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the real number to be decomposed.
        //
        //    Output, int &S, the "sign" of the number.
        //    S will be -1 if X is less than 0, and +1 if X is greater
        //    than or equal to zero.
        //
        //    Output, double &R, the mantissa of X.  R will be greater
        //    than or equal to 1, and strictly less than 2.  The one
        //    exception occurs if X is zero, in which case R will also
        //    be zero.
        //
        //    Output, int &L, the integer part of the logarithm (base 2) of X.
        //
    {
        s = x switch
        {
            //
            //  Determine the sign.
            //
            < 0.0 => -1,
            _ => 1
        };

        r = x switch
        {
            //
            //  Set R to the absolute value of X, and L to zero.
            //  Then force R to lie between 1 and 2.
            //
            < 0.0 => -x,
            _ => x
        };

        l = 0;
        switch (x)
        {
            //
            //  Time to bail out if X is zero.
            //
            case 0.0:
                return;
        }

        while (2.0 <= r)
        {
            r /= 2.0;
            l += 1;
        }

        while (r < 1.0)
        {
            r *= 2.0;
            l -= 1;
        }
    }

    public static double r8_fractional(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_FRACTIONAL returns the fractional part of an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_FRACTIONAL, the fractional part of X.
        //
    {
        double value = Math.Abs(x) - (int) Math.Abs(x);

        return value;
    }

    public static double r8_haversine(double a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_HAVERSINE computes the haversine of an angle.
        //
        //  Discussion:
        //
        //    haversine(A) = ( 1 - cos ( A ) ) / 2
        //
        //    The haversine is useful in spherical trigonometry.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 November 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the angle.
        //
        //    Output, double R8_HAVERSINE, the haversine of the angle.
        //
    {
        double value = (1.0 - Math.Cos(a)) / 2.0;

        return value;
    }

    public static double r8_heaviside(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_HEAVISIDE evaluates the Heaviside function.
        //
        //  Discussion:
        //
        //    The Heaviside function is 0 for x < 0, 1 for x > 0, and 1/2 for x = 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 November 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_HEAVISIDE, the value.
        //
    {
        double value = x switch
        {
            < 0.0 => 0.0,
            0.0 => 0.5,
            _ => 1.0
        };

        return value;
    }

    public static double r8_hypot(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_HYPOT returns the value of sqrt ( X^2 + Y^2 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the arguments.
        //
        //    Output, double R8_HYPOT, the value of sqrt ( X^2 + Y^2 ).
        //
    {
        double a;
        double b;

        if (Math.Abs(x) < Math.Abs(y))
        {
            a = Math.Abs(y);
            b = Math.Abs(x);
        }
        else
        {
            a = Math.Abs(x);
            b = Math.Abs(y);
        }

        double value = a switch
        {
            //
            //  A contains the larger value.
            //
            0.0 => 0.0,
            _ => a * Math.Sqrt(1.0 + b / a * (b / a))
        };

        return value;
    }


    public static double r8_walsh_1d ( double x, int digit )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_WALSH_1D evaluates the Walsh function of a real scalar argument.
        //
        //  Discussion:
        //
        //    Consider the binary representation of X, and number the digits
        //    in descending order, from leading to lowest, with the units digit
        //    being numbered 0.
        //
        //    The Walsh function W(J)(X) is equal to the J-th binary digit of X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the Walsh function.
        //
        //    Input, int DIGIT, the index of the Walsh function.
        //
        //    Output, double R8_WALSH_1D, the value of the Walsh function.
        //
    {
        //
        //  Hide the effect of the sign of X.
        //
        x = Math.Abs ( x );
        //
        //  If DIGIT is positive, divide by 2 DIGIT times.
        //  If DIGIT is negative, multiply by 2 (-DIGIT) times.
        //
        x /= Math.Pow ( 2.0, digit );
        //
        //  Make it an integer.
        //  Because it's positive, and we're using INT, we don't change the
        //  units digit.
        //
        int n = ( int ) x;
        double value = (n % 2) switch
        {
            //
            //  Is the units digit odd or even?
            //
            0 => 0.0,
            _ => 1.0
        };

        return value;
    }

    public static double r8_wrap(double r, double rlo, double rhi)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_WRAP forces an R8 to lie between given limits by wrapping.
        //
        //  Discussion:
        //
        //    An R8 is a double value.
        //
        //  Example:
        //
        //    RLO = 4.0, RHI = 8.0
        //
        //     R  Value
        //
        //    -2     8
        //    -1     4
        //     0     5
        //     1     6
        //     2     7
        //     3     8
        //     4     4
        //     5     5
        //     6     6
        //     7     7
        //     8     8
        //     9     4
        //    10     5
        //    11     6
        //    12     7
        //    13     8
        //    14     4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 December 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, a value.
        //
        //    Input, double RLO, RHI, the desired bounds.
        //
        //    Output, double R8_WRAP, a "wrapped" version of the value.
        //
    {
        double rhi2;
        double rlo2;
        double value;
        //
        //  Guarantee RLO2 < RHI2.
        //
        if (rlo <= rhi)
        {
            rlo2 = rlo;
            rhi2 = rhi;
        }
        else
        {
            rlo2 = rhi;
            rhi2 = rlo;
        }

        //
        //  Find the width.
        //
        double rwide = rhi2 - rlo2;
        switch (rwide)
        {
            //
            //  Add enough copies of (RHI2-RLO2) to R so that the
            //  result ends up in the interval RLO2 - RHI2.
            //
            case 0.0:
                value = rlo;
                break;
            default:
            {
                int n;
                if (r < rlo2)
                {
                    n = (int) ((rlo2 - r) / rwide) + 1;
                    value = r + n * rwide;
                    if (Math.Abs(value - rhi) <= typeMethods.r8_epsilon())
                    {
                        value = rlo;
                    }
                }
                else
                {
                    n = (int) ((r - rlo2) / rwide);
                    value = r - n * rwide;
                    if (Math.Abs(value - rlo) <= typeMethods.r8_epsilon())
                    {
                        value = rhi;
                    }
                }

                break;
            }
        }

        return value;
    }
        
    public static double r8_atan(double y, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ATAN computes the inverse tangent of the ratio Y / X.
        //
        //  Discussion:
        //
        //    R8_ATAN returns an angle whose tangent is ( Y / X ), a job which
        //    the built in functions ATAN and ATAN2 already do.
        //
        //    However:
        //
        //    * R8_ATAN always returns a positive angle, between 0 and 2 PI,
        //      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
        //      and [-PI,+PI] respectively;
        //
        //    * R8_ATAN accounts for the signs of X and Y, (as does ATAN2).  The ATAN
        //     function by contrast always returns an angle in the first or fourth
        //     quadrants.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double Y, X, two quantities which represent the tangent of
        //    an angle.  If Y is not zero, then the tangent is (Y/X).
        //
        //    Output, double R8_ATAN, an angle between 0 and 2 * PI, whose tangent is
        //    (Y/X), and which lies in the appropriate quadrant so that the signs
        //    of its cosine and sine match those of X and Y.
        //
    {
        double theta = 0;
        switch (x)
        {
            //
            //  Special cases:
            //
            case 0.0:
                theta = y switch
                {
                    > 0.0 => Math.PI / 2.0,
                    < 0.0 => 3.0 * Math.PI / 2.0,
                    0.0 => 0.0,
                    _ => theta
                };

                break;
            default:
            {
                switch (y)
                {
                    case 0.0:
                        theta = x switch
                        {
                            > 0.0 => 0.0,
                            < 0.0 => Math.PI,
                            _ => theta
                        };

                        break;
                    //
                    default:
                        double abs_y = Math.Abs(y);
                        double abs_x = Math.Abs(x);

                        double theta_0 = Math.Atan2(abs_y, abs_x);

                        theta = x switch
                        {
                            > 0.0 when 0.0 < y => theta_0,
                            < 0.0 when 0.0 < y => Math.PI - theta_0,
                            < 0.0 when y < 0.0 => Math.PI + theta_0,
                            > 0.0 when y < 0.0 => 2.0 * Math.PI - theta_0,
                            _ => theta
                        };

                        break;
                }

                break;
            }
        }

        return theta;
    }


        
    public static double r8_nth_root ( double x, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_NTH_ROOT returns the nth-root of an R8.
        //
        //  Discussion:
        //
        //    The nth root of X is x^(1/n)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, real X, the number whose nth root is desired.
        //
        //    Input, integer N, the index of the root.
        //
        //    Output, real VALUE, the Nth root of X.
        //
    {
        double value = 0;
        switch (x)
        {
            //
            //  Potential Error 1: 0^0
            //  But we will use it as 1.
            //
            case 0.0 when n == 0:
                value = 1.0;
                return value;
            //
            //  Error 2: 0^(negative power)
            //
            case 0.0 when n < 0:
                value = double.NaN;
                return value;
            //
            //  Error 3: (negative)^(even strictly positive root)
            //
            case < 0.0 when n % 2 == 0 && 0 < n:
                value = double.NaN;
                return value;
        }

        switch (n)
        {
            //
            //  X^0 = 1
            //
            case 0:
                value = 1.0;
                break;
            //
            //  X^1 = X
            //
            case 1:
                value = x;
                break;
            //
            //  X^(-1) = 1/X
            //
            case -1:
                value = 1.0 / x;
                break;
            default:
            {
                double e = 1.0 / Math.Abs ( n );

                value = n switch
                {
                    < 0 => 1.0 / value,
                    _ => x switch
                    {
                        > 0.0 => Math.Pow(x, e),
                        0.0 => 0.0,
                        _ => -Math.Pow(-x, e)
                    }
                };

                break;
            }
        }

        return value;
    }

    public static double r8_pi ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_PI returns the value of PI as an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double R8_PI, the value of PI.
        //
    {
        return Math.PI;
    }

    public static double r8_pi_sqrt ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_PI_SQRT returns the square root of PI as an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double R8_PI_SQRT, the square root of PI.
        //
    {
        return Math.Sqrt(Math.PI);
    }


    public static void r8_print(double r, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_PRINT prints an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the value to print.
        //
        //    Input, string TITLE, a title.
        //
    {
        Console.WriteLine(title + "  " + r + "");
    }

    public static double r8_random(int[] iseed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_RANDOM returns a uniformly distributed random number between 0 and 1.
        //
        //  Discussion:
        //
        //    This routine uses a multiplicative congruential method with modulus
        //    2**48 and multiplier 33952834046453 (see G.S.Fishman,
        //    'Multiplicative congruential random number generators with modulus
        //    2**b: an exhaustive analysis for b = 32 and a partial analysis for
        //    b = 48', Math. Comp. 189, pp 331-344, 1990).
        //
        //    48-bit integers are stored in 4 integer array elements with 12 bits
        //    per element. Hence the routine is portable across machines with
        //    integers of 32 bits or more.
        //
        //  Parameters:
        //
        //    Input/output, integer ISEED(4).
        //    On entry, the seed of the random number generator; the array
        //    elements must be between 0 and 4095, and ISEED(4) must be odd.
        //    On exit, the seed is updated.
        //
        //    Output, double R8_RANDOM, the next pseudorandom number.
        //
    {
        const int ipw2 = 4096;
        const int m1 = 494;
        const int m2 = 322;
        const int m3 = 2508;
        const int m4 = 2549;
        const double r = 1.0 / 4096.0;
        //
        //  Multiply the seed by the multiplier modulo 2^48.
        //
        int it4 = iseed[3] * m4;
        int it3 = it4 / ipw2;
        it4 -= ipw2 * it3;
        it3 = it3 + iseed[2] * m4 + iseed[3] * m3;
        int it2 = it3 / ipw2;
        it3 -= ipw2 * it2;
        it2 = it2 + iseed[1] * m4 + iseed[2] * m3 + iseed[3] * m2;
        int it1 = it2 / ipw2;
        it2 -= ipw2 * it1;
        it1 = it1 + iseed[0] * m4 + iseed[1] * m3 + iseed[2] * m2 + iseed[3] * m1;
        it1 %= ipw2;
        //
        //  Return updated seed
        //
        iseed[0] = it1;
        iseed[1] = it2;
        iseed[2] = it3;
        iseed[3] = it4;
        //
        //  Convert 48-bit integer to a real number in the interval (0,1)
        //
        double value = r * (it1
                            + r * (it2
                                   + r * (it3
                                          + r * it4)));

        return value;
    }
}