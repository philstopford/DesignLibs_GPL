using System;
using Burkardt.Types;

namespace Burkardt.IntegralNS;

public static class Filon
{
    public static double filon_fun_cos(int n, Func< int, double[], double[] > f, double a,
            double b, double t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FILON_FUN_COS uses Filon's method on integrals with a coMath.Sine factor.
        //
        //  Discussion:
        //
        //    The integral to be approximated has the form:
        //
        //      Integral ( A <= X <= B ) F(X) * COS(T*X) dX
        //
        //    where T is user specified.
        //
        //    The function is interpolated over each subinterval by
        //    a parabolic arc.
        //
        //  LicenMath.Sing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Stephen Chase, Lloyd Fosdick,
        //    An Algorithm for Filon Quadrature,
        //    Communications of the Association for Computing Machinery,
        //    Volume 12, Number 8, August 1969, pages 453-457.
        //
        //    Stephen Chase, Lloyd Fosdick,
        //    Algorithm 353:
        //    Filon Quadrature,
        //    Communications of the Association for Computing Machinery,
        //    Volume 12, Number 8, August 1969, pages 457-458.
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, int N, the number of data points.
        //    N must be odd, and greater than 1.
        //
        //    Input, double *F ( int n, double x[] ), the function which evaluates the 
        //    integrand.
        //
        //    Input, double A, B, the limits of integration.
        //
        //    Input, double T, the multiplier of the X argument of the coMath.Sine.
        //
        //    Output, double FILON_FUN_COS, the approximate value of the integral.
        //
    {
        double alpha;
        double beta;
        double gamma;
        int i;
        double value;

        if (Math.Abs(a - b) <= typeMethods.r8_epsilon())
        {
            value = 0.0;
            return value;
        }

        switch (n)
        {
            case <= 1:
                Console.WriteLine("");
                Console.WriteLine("FILON_FUN_COS - Fatal error!");
                Console.WriteLine("  N < 2");
                Console.WriteLine("  N = " + n + "");
                return 1;
        }

        if (n % 2 != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("FILON_FUN_COS - Fatal error!");
            Console.WriteLine("  N must be odd.");
            Console.WriteLine("  N = " + n + "");
            return 1;
        }

        //
        //  Set the X values.
        //
        double[] x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = ((n - i - 1) * a
                    + i * b)
                   / (n - 1);
        }

        double h = (b - a) / (n - 1);
        double theta = t * h;
        double sint = Math.Sin(theta);
        double cost = Math.Cos(theta);

        switch (6.0 * Math.Abs(theta))
        {
            case <= 1.0:
                alpha = 2.0 * Math.Pow(theta, 3) / 45.0
                        - 2.0 * Math.Pow(theta, 5) / 315.0
                        + 2.0 * Math.Pow(theta, 7) / 4725.0;

                beta = 2.0 / 3.0
                       + 2.0 * Math.Pow(theta, 2) / 15.0
                       - 4.0 * Math.Pow(theta, 4) / 105.0
                       + 2.0 * Math.Pow(theta, 6) / 567.0
                       - 4.0 * Math.Pow(theta, 8) / 22275.0;

                gamma = 4.0 / 3.0
                        - 2.0 * Math.Pow(theta, 2) / 15.0
                        + Math.Pow(theta, 4) / 210.0
                        - Math.Pow(theta, 6) / 11340.0;
                break;
            default:
                alpha = (Math.Pow(theta, 2) + theta * sint * cost - 2.0 * sint * sint)
                        / Math.Pow(theta, 3);

                beta = (2.0 * theta + 2.0 * theta * cost * cost
                        - 4.0 * sint * cost) / Math.Pow(theta, 3);

                gamma = 4.0 * (sint - theta * cost) / Math.Pow(theta, 3);
                break;
        }

        //
        //  Tabulate the function.
        //
        double[] ftab = f(n, x);

        double c2n = 0.5 * ftab[0] * Math.Cos(t * x[0]);
        for (i = 2; i < n - 1; i += 2)
        {
            c2n += ftab[i] * Math.Cos(t * x[i]);
        }

        c2n += 0.5 * ftab[n - 1] * Math.Cos(t * x[n - 1]);

        double c2nm1 = 0.0;
        for (i = 1; i <= n - 2; i += 2)
        {
            c2nm1 += ftab[i] * Math.Cos(t * x[i]);
        }

        value = h * (
            alpha * (ftab[n - 1] * Math.Sin(t * x[n - 1])
                     - ftab[0] * Math.Sin(t * x[0]))
            + beta * c2n
            + gamma * c2nm1);

        return value;
    }

    public static double filon_tab_cos(int n, double[] ftab, double a, double b, double t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FILON_TAB_COS uses Filon's method on integrals with a coMath.Sine factor.
        //
        //  Discussion:
        //
        //    The integral to be approximated has the form:
        //
        //      Integral ( A <= X <= B ) F(X) * COS(T*X) dX
        //
        //    where T is user specified.
        //
        //    The function is interpolated over each subinterval by
        //    a parabolic arc.
        //
        //  LicenMath.Sing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Stephen Chase, Lloyd Fosdick,
        //    An Algorithm for Filon Quadrature,
        //    Communications of the Association for Computing Machinery,
        //    Volume 12, Number 8, August 1969, pages 453-457.
        //
        //    Stephen Chase, Lloyd Fosdick,
        //    Algorithm 353:
        //    Filon Quadrature,
        //    Communications of the Association for Computing Machinery,
        //    Volume 12, Number 8, August 1969, pages 457-458.
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, int N, the number of data points.
        //    N must be odd, and greater than 1.
        //
        //    Input, double FTAB[N], contains the value of the function
        //    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(N-1).
        //
        //    Input, double A, B, the limits of integration.
        //
        //    Input, double T, the multiplier of the X argument of the coMath.Sine.
        //
        //    Output, double FILON_TAB_COS, the approximate value of the integral.
        //
    {
        double alpha;
        double beta;
        double gamma;
        int i;
        double value;

        if (Math.Abs(a - b) <= typeMethods.r8_epsilon())
        {
            value = 0.0;
            return value;
        }

        switch (n)
        {
            case <= 1:
                Console.WriteLine("");
                Console.WriteLine("FILON_TAB_COS - Fatal error!");
                Console.WriteLine("  N < 2");
                Console.WriteLine("  N = " + n + "");
                return 1;
        }

        if (n % 2 != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("FILON_TAB_COS - Fatal error!");
            Console.WriteLine("  N must be odd.");
            Console.WriteLine("  N = " + n + "");
            return 1;
        }

        //
        //  Set the X values.
        //
        double[] x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = ((n - i - 1) * a
                    + i * b)
                   / (n - 1);
        }

        double h = (b - a) / (n - 1);
        double theta = t * h;
        double sint = Math.Sin(theta);
        double cost = Math.Cos(theta);

        switch (6.0 * Math.Abs(theta))
        {
            case <= 1.0:
                alpha = 2.0 * Math.Pow(theta, 3) / 45.0
                        - 2.0 * Math.Pow(theta, 5) / 315.0
                        + 2.0 * Math.Pow(theta, 7) / 4725.0;

                beta = 2.0 / 3.0
                       + 2.0 * Math.Pow(theta, 2) / 15.0
                       - 4.0 * Math.Pow(theta, 4) / 105.0
                       + 2.0 * Math.Pow(theta, 6) / 567.0
                       - 4.0 * Math.Pow(theta, 8) / 22275.0;

                gamma = 4.0 / 3.0
                        - 2.0 * Math.Pow(theta, 2) / 15.0
                        + Math.Pow(theta, 4) / 210.0
                        - Math.Pow(theta, 6) / 11340.0;
                break;
            default:
                alpha = (Math.Pow(theta, 2) + theta * sint * cost - 2.0 * sint * sint)
                        / Math.Pow(theta, 3);

                beta = (2.0 * theta + 2.0 * theta * cost * cost
                        - 4.0 * sint * cost) / Math.Pow(theta, 3);

                gamma = 4.0 * (sint - theta * cost) / Math.Pow(theta, 3);
                break;
        }

        double c2n = +0.5 * ftab[0] * Math.Cos(t * x[0]);
        for (i = 2; i < n - 1; i += 2)
        {
            c2n += ftab[i] * Math.Cos(t * x[i]);
        }

        c2n += 0.5 * ftab[n - 1] * Math.Cos(t * x[n - 1]);

        double c2nm1 = 0.0;
        for (i = 1; i <= n - 2; i += 2)
        {
            c2nm1 += ftab[i] * Math.Cos(t * x[i]);
        }

        value = h * (
            alpha * (ftab[n - 1] * Math.Sin(t * x[n - 1])
                     - ftab[0] * Math.Sin(t * x[0]))
            + beta * c2n
            + gamma * c2nm1);
            
        return value;
    }

    public static double filon_fun_sin(int n, Func < int, double[], double[] > f, double a,
            double b, double t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FILON_FUN_SIN uses Filon's method on integrals with a Math.Sine factor.
        //
        //  Discussion:
        //
        //    The integral to be approximated has the form
        //
        //      Integral ( A <= X <= B ) F(X) * SIN(T*X) dX
        //
        //    where T is user specified.
        //
        //    The function is interpolated over each subinterval by
        //    a parabolic arc.
        //
        //  LicenMath.Sing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Stephen Chase, Lloyd Fosdick,
        //    An Algorithm for Filon Quadrature,
        //    Communications of the Association for Computing Machinery,
        //    Volume 12, Number 8, August 1969, pages 453-457.
        //
        //    Stephen Chase, Lloyd Fosdick,
        //    Algorithm 353:
        //    Filon Quadrature,
        //    Communications of the Association for Computing Machinery,
        //    Volume 12, Number 8, August 1969, pages 457-458.
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, int N, the number of data points, 
        //    including the endpoints.  N must be odd, and greater than 1.
        //
        //    Input, external F, the subroutine which evaluates the integrand,
        //    of the form subroutine F ( N, X, FX ).
        //
        //    Input, double A, B, the limits of integration.
        //
        //    Input, double T, multiplier of the X argument of the Math.Sine.
        //
        //    Output, double FILON_FUN_SIN, the approximate value of the integral.
        //
    {
        double alpha;
        double beta;
        double gamma;
        int i;
        double value;

        if (Math.Abs(a - b) <= typeMethods.r8_epsilon())
        {
            value = 0.0;
            return value;
        }

        switch (n)
        {
            case <= 1:
                Console.WriteLine("");
                Console.WriteLine("FILON_FUN_SIN - Fatal error!");
                Console.WriteLine("  N < 2");
                Console.WriteLine("  N = " + n + "");
                return 1;
        }

        if (n % 2 != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("FILON_FUN_SIN - Fatal error!");
            Console.WriteLine("  N must be odd.");
            Console.WriteLine("  N = " + n + "");
            return 1;
        }

        //
        //  Set the X values.
        //
        double[] x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = ((n - i - 1) * a
                    + i * b)
                   / (n - 1);
        }

        double h = (b - a) / (n - 1);
        double theta = t * h;
        double sint = Math.Sin(theta);
        double cost = Math.Cos(theta);

        switch (6.0 * Math.Abs(theta))
        {
            case <= 1.0:
                alpha = 2.0 * Math.Pow(theta, 3) / 45.0
                        - 2.0 * Math.Pow(theta, 5) / 315.0
                        + 2.0 * Math.Pow(theta, 7) / 4725.0;

                beta = 2.0 / 3.0
                       + 2.0 * Math.Pow(theta, 2) / 15.0
                       - 4.0 * Math.Pow(theta, 4) / 105.0
                       + 2.0 * Math.Pow(theta, 6) / 567.0
                       - 4.0 * Math.Pow(theta, 8) / 22275.0;

                gamma = 4.0 / 3.0
                        - 2.0 * Math.Pow(theta, 2) / 15.0
                        + Math.Pow(theta, 4) / 210.0
                        - Math.Pow(theta, 6) / 11340.0;
                break;
            default:
                alpha = (Math.Pow(theta, 2) + theta * sint * cost
                         - 2.0 * sint * sint) / Math.Pow(theta, 3);

                beta = (2.0 * theta + 2.0 * theta * cost * cost
                        - 4.0 * sint * cost) / Math.Pow(theta, 3);

                gamma = 4.0 * (sint - theta * cost) / Math.Pow(theta, 3);
                break;
        }

        //
        //  Tabulate the function.
        //
        double[] ftab = f(n, x);

        double s2n = +0.5 * ftab[0] * Math.Sin(t * x[0]);
        for (i = 2; i < n - 1; i += 2)
        {
            s2n += ftab[i] * Math.Sin(t * x[i]);
        }

        s2n += 0.5 * ftab[n - 1] * Math.Sin(t * x[n - 1]);

        double s2nm1 = 0.0;
        for (i = 1; i <= n - 2; i += 2)
        {
            s2nm1 += ftab[i] * Math.Sin(t * x[i]);
        }

        value = h * (
            alpha * (ftab[0] * Math.Cos(t * x[0])
                     - ftab[n - 1] * Math.Cos(t * x[n - 1]))
            + beta * s2n
            + gamma * s2nm1);

        return value;
    }

    public static double filon_tab_sin(int n, double[] ftab, double a, double b, double t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FILON_TAB_SIN uses Filon's method on integrals with a Math.Sine factor.
        //
        //  Discussion:
        //
        //    The integral to be approximated has the form
        //
        //      Integral ( A <= X <= B ) F(X) * SIN(T*X) dX
        //
        //    where T is user specified.
        //
        //    The function is interpolated over each subinterval by
        //    a parabolic arc.
        //
        //  LicenMath.Sing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Stephen Chase, Lloyd Fosdick,
        //    An Algorithm for Filon Quadrature,
        //    Communications of the Association for Computing Machinery,
        //    Volume 12, Number 8, August 1969, pages 453-457.
        //
        //    Stephen Chase, Lloyd Fosdick,
        //    Algorithm 353:
        //    Filon Quadrature,
        //    Communications of the Association for Computing Machinery,
        //    Volume 12, Number 8, August 1969, pages 457-458.
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, int N, the number of data points, 
        //    including the endpoints.  N must be odd, and greater than 1.
        //
        //    Input, double FTAB[N], contains the value of the function
        //    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(N-1).
        //
        //    Input, double A, B, the limits of integration.
        //
        //    Input, double T, multiplier of the X argument of the Math.Sine.
        //
        //    Output, double FILON_TAB_SIN, the approximate value of the integral.
        //
    {
        double alpha;
        double beta;
        double gamma;
        int i;
        double value;

        if (Math.Abs(a - b) <= typeMethods.r8_epsilon())
        {
            value = 0.0;
            return value;
        }

        switch (n)
        {
            case <= 1:
                Console.WriteLine("");
                Console.WriteLine("FILON_TAB_SIN - Fatal error!");
                Console.WriteLine("  N < 2");
                Console.WriteLine("  N = " + n + "");
                return 1;
        }

        if (n % 2 != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("FILON_TAB_SIN - Fatal error!");
            Console.WriteLine("  N must be odd.");
            Console.WriteLine("  N = " + n + "");
            return 1;
        }

        //
        //  Set the X values.
        //
        double[] x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = ((n - i - 1) * a
                    + i * b)
                   / (n - 1);
        }

        double h = (b - a) / (n - 1);
        double theta = t * h;
        double sint = Math.Sin(theta);
        double cost = Math.Cos(theta);

        switch (6.0 * Math.Abs(theta))
        {
            case <= 1.0:
                alpha = 2.0 * Math.Pow(theta, 3) / 45.0
                        - 2.0 * Math.Pow(theta, 5) / 315.0
                        + 2.0 * Math.Pow(theta, 7) / 4725.0;

                beta = 2.0 / 3.0
                       + 2.0 * Math.Pow(theta, 2) / 15.0
                       - 4.0 * Math.Pow(theta, 4) / 105.0
                       + 2.0 * Math.Pow(theta, 6) / 567.0
                       - 4.0 * Math.Pow(theta, 8) / 22275.0;

                gamma = 4.0 / 3.0
                        - 2.0 * Math.Pow(theta, 2) / 15.0
                        + Math.Pow(theta, 4) / 210.0
                        - Math.Pow(theta, 6) / 11340.0;
                break;
            default:
                alpha = (Math.Pow(theta, 2) + theta * sint * cost
                         - 2.0 * sint * sint) / Math.Pow(theta, 3);

                beta = (2.0 * theta + 2.0 * theta * cost * cost
                        - 4.0 * sint * cost) / Math.Pow(theta, 3);

                gamma = 4.0 * (sint - theta * cost) / Math.Pow(theta, 3);
                break;
        }

        double s2n = +0.5 * ftab[0] * Math.Sin(t * x[0]);
        for (i = 2; i < n - 1; i += 2)
        {
            s2n += ftab[i] * Math.Sin(t * x[i]);
        }

        s2n += 0.5 * ftab[n - 1] * Math.Sin(t * x[n - 1]);

        double s2nm1 = 0.0;
        for (i = 1; i <= n - 2; i += 2)
        {
            s2nm1 += ftab[i] * Math.Sin(t * x[i]);
        }

        value = h * (
            alpha * (ftab[0] * Math.Cos(t * x[0])
                     - ftab[n - 1] * Math.Cos(t * x[n - 1]))
            + beta * s2n
            + gamma * s2nm1);

        return value;
    }
}