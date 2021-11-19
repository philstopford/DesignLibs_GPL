using System;
using Burkardt.IntegralNS;

namespace CauchyPrincipalValueTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CAUCHY_PRINCIPAL_VALUE_TEST tests the CAUCHY_PRINCIPAL_VALUE library.
        //
        //  Location:
        //
        //    http://people.sc.fsu.edu/~jburkardt/cpp_src/cauchy_principal_value/cauchy_principal_value_test.cpp
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CAUCHY_PRINCIPAL_VALUE_TEST");
        Console.WriteLine("  Test the CAUCHY_PRINCIPAL_VALUE library.");

        cpv_test01();
        cpv_test02();

        Console.WriteLine("");
        Console.WriteLine("CAUCHY_PRINCIPAL_VALUE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void cpv_test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CPV_TEST01 seeks the CPV of Integral ( -1 <= t <= 1 ) exp(t) / t dt
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        double exact;
        int n;
        double value = 0;

        Console.WriteLine("");
        Console.WriteLine("CPV_TEST01:");
        Console.WriteLine("  CPV of Integral ( -1 <= t <= 1 ) exp(t) / t dt");

        Console.WriteLine("");
        Console.WriteLine("   N           Estimate             Error");
        Console.WriteLine("");

        exact = 2.11450175075;
        a = -1.0;
        b = +1.0;
        for (n = 2; n <= 8; n += 2)
        {
            value = CauchyPrincipalValue.cpv(f01, a, b, n);
            Console.WriteLine("  " + n.ToString().PadLeft(2)
                                   + "  " + value.ToString().PadLeft(24)
                                   + "  " + Math.Abs(value - exact).ToString().PadLeft(14) + "");
        }

    }

    private static double f01(double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F01 evaluates the integrand of Integral ( -1 <= t <= 1 ) exp(t) / t dt
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, real T, the argument.
        //
        //    Output, real VALUE, the value of the integrand.
        //
    {
        double value = Math.Exp(t);

        return value;
    }

    private static void cpv_test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CPV_TEST02 is another test.
        //
        // Discussion:
        //
        //    We seek
        //      CPV ( Integral ( 1-delta <= t <= 1+delta ) 1/(1-t)^3 dt )
        //    which we must rewrite as
        //      CPV ( Integral ( 1-delta <= t <= 1+delta ) 1/(1+t+t^2) 1/(1-t) dt )
        //    so that our "integrand" is 1/(1+t+t^2).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        double delta;
        double exact;
        int k;
        int n;
        double r1;
        double r2;
        double r3;
        double r4;
        double value = 0;

        Console.WriteLine("");
        Console.WriteLine("CPV_TEST02:");
        Console.WriteLine("  Compute CPV ( Integral ( 1-delta <= t <= 1+delta ) 1/(1-t)^3 dt )");
        Console.WriteLine("  Try this for delta = 1, 1/2, 1/4.");
        Console.WriteLine("");
        Console.WriteLine("   N          Estimate                  Exact                  Error");
        delta = 1.0;
        for (k = 1; k <= 3; k++)
        {
            Console.WriteLine("");
            r1 = Math.Pow(delta + 1.5, 2) + 0.75;
            r2 = Math.Pow(-delta + 1.5, 2) + 0.75;
            r3 = Math.Atan(Math.Sqrt(0.75) / (delta + 1.5));
            r4 = Math.Atan(Math.Sqrt(0.75) / (-delta + 1.5));
            exact = -Math.Log(r1 / r2) / 6.0 + (r3 - r4) / Math.Sqrt(3.0);
            for (n = 2; n <= 8; n += 2)
            {
                a = 1.0 - delta;
                b = 1.0 + delta;
                value = CauchyPrincipalValue.cpv(f02, a, b, n);
                Console.WriteLine("  " + n.ToString().PadLeft(2)
                                       + "  " + value.ToString().PadLeft(24)
                                       + "  " + exact.ToString().PadLeft(24)
                                       + "  " + Math.Abs(value - exact).ToString().PadLeft(14) + "");
            }

            delta /= 2.0;
        }
    }

    private static double f02(double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F02: integrand of Integral ( 1-delta <= t <= 1+delta ) 1/(1-t^3) dt
        //
        //  Discussion:
        //
        //    1/(1-t^3) = 1/(1+t+t^2) * 1/(1-t)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T, the evaluation point.
        //
        //    Output, double F02, the value of the integrand at T.
        //
    {
        double value = 1.0 / (1.0 + t + t * t);

        return value;
    }
}