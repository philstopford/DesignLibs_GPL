using System;

namespace Burkardt.FEM;

public static class FEM_Test_Methods
{
    public static double a00(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A00 evaluates A function #0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A00, the value of A(X).
        //
    {
        double value = 1.0;

        return value;
    }

    public static double c00(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C00 evaluates C function #0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C00, the value of C(X).
        //
    {
        double value = 1.0;

        return value;
    }

    public static double exact00(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT00 evaluates exact solution #00.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT00, the value of U(X).
        //
    {
        double value = x - Math.Sinh(x) / Math.Sinh(1.0);

        return value;
    }

    public static double exact_ux00(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX00 evaluates the derivative of exact solution #00.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX00, the value of dUdX(X).
        //
    {
        double value = 1.0 - Math.Cosh(x) / Math.Sinh(1.0);

        return value;
    }

    public static double f00(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F00 evaluates right hand side function #00.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F00, the value of F(X).
        //
    {
        double value = x;

        return value;
    }

    public static double a1(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A1 evaluates A function #1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A1, the value of A(X).
        //
    {
        double value = 1.0;

        return value;
    }

    public static double c1(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C1 evaluates C function #1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C1, the value of C(X).
        //
    {
        double value = 0.0;

        return value;
    }

    public static double exact1(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT1 evaluates exact solution #1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT1, the value of U(X).
        //
    {
        double value = x * (1.0 - x) * Math.Exp(x);

        return value;
    }

    public static double exact_ux1(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX1 evaluates the derivative of exact solution #1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX1, the value of dUdX(X).
        //
    {
        double value = (1.0 - x - x * x) * Math.Exp(x);

        return value;
    }

    public static double f1(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F1 evaluates right hand side function #1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F1, the value of F(X).
        //
    {
        double value = x * (x + 3.0) * Math.Exp(x);

        return value;
    }


    public static double a2(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A2 evaluates A function #2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A2, the value of A(X).
        //
    {
        double value = 1.0;

        return value;
    }

    public static double c2(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C2 evaluates C function #2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C2, the value of C(X).
        //
    {
        double value = 2.0;

        return value;
    }

    public static double exact2(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT2 evaluates exact solution #2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT2, the value of U(X).
        //
    {
        double value = x * (1.0 - x) * Math.Exp(x);

        return value;
    }

    public static double exact_ux2(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX2 evaluates the derivative of exact solution #2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX2, the value of dUdX(X).
        //
    {
        double value = (1.0 - x - x * x) * Math.Exp(x);

        return value;
    }

    public static double f2(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F2 evaluates right hand side function #2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F2, the value of F(X).
        //
    {
        double value = x * (5.0 - x) * Math.Exp(x);

        return value;
    }


    public static double a3(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A3 evaluates A function #3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A3, the value of A(X).
        //
    {
        double value = 1.0;

        return value;
    }

    public static double c3(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C3 evaluates C function #3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C3, the value of C(X).
        //
    {
        double value = 2.0 * x;

        return value;
    }

    public static double exact3(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT3 evaluates exact solution #3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT3, the value of U(X).
        //
    {
        double value = x * (1.0 - x) * Math.Exp(x);

        return value;
    }

    public static double exact_ux3(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX3 evaluates the derivative of exact solution #3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX3, the value of dUdX(X).
        //
    {
        double value = (1.0 - x - x * x) * Math.Exp(x);

        return value;
    }

    public static double f3(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F3 evaluates right hand side function #3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F3, the value of F(X).
        //
    {
        double value = -x * (2.0 * x * x - 3.0 * x - 3.0) * Math.Exp(x);

        return value;
    }

    public static double a4(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A4 evaluates A function #4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A4, the value of A(X).
        //
    {
        double value = 1.0 + x * x;

        return value;
    }

    public static double c4(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C4 evaluates C function #4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C4, the value of C(X).
        //
    {
        double value = 0.0;

        return value;
    }

    public static double exact4(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT4 evaluates exact solution #4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT4, the value of U(X).
        //
    {
        double value = x * (1.0 - x) * Math.Exp(x);

        return value;
    }

    public static double exact_ux4(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX4 evaluates the derivative of exact solution #4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX4, the value of dUdX(X).
        //
    {
        double value = (1.0 - x - x * x) * Math.Exp(x);

        return value;
    }

    public static double f4(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F4 evaluates right hand side function #4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F4, the value of F(X).
        //
    {
        double value = (x + 3.0 * x * x + 5.0 * x * x * x + x * x * x * x) * Math.Exp(x);

        return value;
    }

    public static double a5(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A5 evaluates A function #5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A5, the value of A(X).
        //
    {
        double value = x switch
        {
            <= 1.0 / 3.0 => 1.0 + x * x,
            _ => x + 7.0 / 9.0
        };

        return value;
    }

    public static double c5(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C5 evaluates C function #5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C5, the value of C(X).
        //
    {
        double value = 0.0;

        return value;
    }

    public static double exact5(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT5 evaluates exact solution #5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT5, the value of U(X).
        //
    {
        double value = x * (1.0 - x) * Math.Exp(x);

        return value;
    }

    public static double exact_ux5(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX5 evaluates the derivative of exact solution #5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX5, the value of dUdX(X).
        //
    {
        double value = (1.0 - x - x * x) * Math.Exp(x);

        return value;
    }

    public static double f5(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F5 evaluates right hand side function #5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F5, the value of F(X).
        //
    {
        double value = x switch
        {
            <= 1.0 / 3.0 => (x + 3.0 * x * x + 5.0 * x * x * x + x * x * x * x) * Math.Exp(x),
            _ => (-1.0 + 10.0 / 3.0 * x + 43.0 / 9.0 * x * x + x * x * x) * Math.Exp(x)
        };

        return value;
    }

    public static double a6(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A6 evaluates A function #6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A6, the value of A(X).
        //
    {
        double value = 1.0;

        return value;
    }

    public static double c6(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C6 evaluates C function #6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C6, the value of C(X).
        //
    {
        double value = 0.0;

        return value;
    }

    public static double exact6(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT6 returns exact solution #6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT6, the value of U(X).
        //
    {
        double value = Math.Sin(Math.PI * x);

        return value;
    }

    public static double exact_ux6(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX6 returns the derivative of exact solution #6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX6, the value of U(X).
        //
    {
        double value = Math.PI * Math.Cos(Math.PI * x);

        return value;
    }

    public static double f6(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F6 evaluates right hand side function #6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F6, the value of F(X).
        //
    {
        double value = Math.PI * Math.PI * Math.Sin(Math.PI * x);

        return value;
    }

    public static double a7(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A7 evaluates A function #7.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A7, the value of A(X).
        //
    {
        double alpha = 30.0;
        double x0 = 1.0 / 3.0;
        double value = 1.0 / alpha + alpha * Math.Pow(x - x0, 2);

        return value;
    }

    public static double c7(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C7 evaluates C function #7.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C7, the value of C(X).
        //
    {
        double value = 0.0;

        return value;
    }

    public static double exact7(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT7 returns exact solution #7.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT7, the value of U(X).
        //
    {
        double alpha = 30.0;
        double x0 = 1.0 / 3.0;
        double value = (1.0 - x)
                       * (Math.Atan(alpha * (x - x0)) + Math.Atan(alpha * x0));

        return value;
    }

    public static double exact_ux7(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX7 returns the derivative of exact solution #7.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX7, the value of U(X).
        //
    {
        double alpha = 30.0;
        double x0 = 1.0 / 3.0;
        double value = -Math.Atan(alpha * (x - x0)) - Math.Atan(alpha * x0)
                       + (1.0 - x) * alpha / (1.0 + alpha * alpha * Math.Pow(x - x0, 2));


        return value;
    }

    public static double f7(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F7 evaluates right hand side function #7.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F7, the value of F(X).
        //
    {
        double alpha = 30.0;
        double x0 = 1.0 / 3.0;
        double value = 2.0 * (1.0 + alpha * (x - x0) *
            (Math.Atan(alpha * (x - x0)) + Math.Atan(alpha * x0)));

        return value;
    }

    public static double a8(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A8 evaluates A function #8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A8, the value of A(X).
        //
    {
        double value = 1.0;

        return value;
    }

    public static double c8(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8 evaluates C function #8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C8, the value of C(X).
        //
    {
        double value = 0.0;

        return value;
    }

    public static double exact8(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT8 evaluates exact solution #8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT8, the value of U(X).
        //
    {
        double value = x switch
        {
            <= 2.0 / 3.0 => x * (1.0 - x) * Math.Exp(x),
            _ => x * (1.0 - x) * Math.Exp(2.0 / 3.0)
        };

        return value;
    }

    public static double exact_ux8(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX8 evaluates the derivative of exact solution #8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX8, the value of dUdX(X).
        //
    {
        double value = x switch
        {
            <= 2.0 / 3.0 => (1.0 - x - x * x) * Math.Exp(x),
            _ => (1.0 - 2.0 * x) * Math.Exp(2.0 / 3.0)
        };

        return value;
    }

    public static double f8(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F8 evaluates right hand side function #8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F8, the value of F(X).
        //
    {
        double value = x switch
        {
            <= 2.0 / 3.0 => x * (x + 3.0) * Math.Exp(x),
            _ => 2.0 * Math.Exp(2.0 / 3.0)
        };

        return value;
    }

    public static double a9(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A9 evaluates A function #9.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A8, the value of A(X).
        //
    {
        double value = 1.0;

        return value;
    }

    public static double c9(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C9 evaluates C function #9.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C9, the value of C(X).
        //
    {
        double value = 0.0;

        return value;
    }

    public static double exact9(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT9 evaluates exact solution #9.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT9, the value of U(X).
        //
    {
        double value = x switch
        {
            <= 2.0 / 3.0 => x * (1.0 - x) * Math.Exp(x),
            _ => x * (1.0 - x)
        };

        return value;
    }

    public static double exact_ux9(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX9 evaluates the derivative of exact solution #9.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX9, the value of dUdX(X).
        //
    {
        double value = x switch
        {
            <= 2.0 / 3.0 => (1.0 - x - x * x) * Math.Exp(x),
            _ => 1.0 - 2.0 * x
        };

        return value;
    }

    public static double f9(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F9 evaluates right hand side function #9.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F9, the value of F(X).
        //
    {
        double value = x switch
        {
            <= 2.0 / 3.0 => x * (x + 3.0) * Math.Exp(x),
            _ => 2.0
        };

        return value;
    }

    public static double a10(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A10 evaluates A function #10.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A10, the value of A(X).
        //
    {
        double value = 1.0;

        return value;
    }

    public static double c10(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C10 evaluates C function #10.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C10, the value of C(X).
        //
    {
        double value = 1.0;

        return value;
    }

    public static double exact10(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT10 evaluates exact solution #10.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT10, the value of U(X).
        //
    {
        double value = x - Math.Sinh(x) / Math.Sinh(1.0);

        return value;
    }

    public static double exact_ux10(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX10 evaluates the derivative of exact solution #10.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX10, the value of dUdX(X).
        //
    {
        double value = 1.0 - Math.Cosh(x) / Math.Sinh(1.0);

        return value;
    }

    public static double f10(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F10 evaluates right hand side function #10.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F10, the value of F(X).
        //
    {
        double value = x;

        return value;
    }
        
    public static double k1 ( double x )
//****************************************************************************80
//
//  Purpose:
//
//    K1 evaluates K function #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double K1, the value of K(X).
//
    {
        double value = 1.0;

        return value;
    }

    public static double f ( double x )
//****************************************************************************80
//
//  Purpose:
//
//    F evaluates the right hand side function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 November 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double F, the value of the right hand side at X.
//
    {
        double value = x;

        return value;
    }

    public static double[] exact ( int x_num, double[] x )
//****************************************************************************80
//
//  Purpose:
//
//    EXACT returns the exact solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 November 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X_NUM, the number of nodes.
//
//    Input, double X[X_NUM], the nodes.
//
//    Output, double UE[X_NUM], the exact solution at the nodes.
//
    {
        int x_i;

        double[] ue = new double[x_num];

        for ( x_i = 0; x_i < x_num; x_i++ )
        {
            ue[x_i] = x[x_i] - Math.Sinh ( x[x_i] ) / Math.Sinh ( 1.0 );
        }
        return ue;
    }        
        
    public static double u_exact ( double x, int problem )
//****************************************************************************80
//
//  Purpose:
//
//    U_EXACT returns the value of the exact solution at a point X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, int PROBLEM, indicates which problem to be solved.
//    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
//    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
//    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
//    u(0) = 0, u'(1)=1.
//
//    Output, double U_EXACT, the value of the exact solution at X.
//
    {
        double value = problem switch
        {
            //
            //  Test problem 1
            //
            1 => x,
            //
            //  Test problem 2
            //
            2 => 2.0 * (1.0 - Math.Cos(0.5 * Math.PI * x)) / Math.PI,
            _ => 0
        };

        return value;
    }
        
        
}