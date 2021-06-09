using System;
using Burkardt.ZeroRC;

namespace Burkardt.ZeroRCTest
{
  class Program
  {
    static void Main(string[] args)
      //****************************************************************************80
      //
      //  Purpose:
      //
      //    zero_rc_test() tests zero_rc().
      //
      //  Licensing:
      //
      //    This code is distributed under the GNU LGPL license.
      //
      //  Modified:
      //
      //    29 May 2021
      //
      //  Author:
      //
      //    John Burkardt
      //
    {
      double a;
      double b;

      Console.WriteLine("");
      Console.WriteLine("zero_rc_test()");
      Console.WriteLine("  C++ version");
      Console.WriteLine("  zero_rc() seeks a root of a function F(X)");
      Console.WriteLine("  in an interval [A,B] using reverse communication.");

      a = 1.0;
      b = 2.0;
      example_test(a, b, f_01, "f_01(x) = sin ( x ) - x / 2");

      a = 0.0;
      b = 1.0;
      example_test(a, b, f_02, "f_02(x) = 2 * x - exp ( - x )");

      a = -1.0;
      b = 0.5;
      example_test(a, b, f_03, "f_03(x) = x * exp ( - x )");

      a = 0.0001;
      b = 20.0;
      example_test(a, b, f_04, "f_04(x) = exp ( x ) - 1 / ( 100 * x * x )");

      a = -5.0;
      b = 2.0;
      example_test(a, b, f_05, "f_05(x) = (x+3) * (x-1) * (x-1)");
      //
      //  Terminate.
      //
      Console.WriteLine("");
      Console.WriteLine("zero_rc_test()");
      Console.WriteLine("  Normal end of execution.");
      Console.WriteLine("");

    }

    static void example_test(double a, double b, Func<double, double> f, string title )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    example_test() tests zero_rc() on one test function.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    29 May 2021
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Input:
    //
    //    double A, B, the endpoints of the change of sign interval.
    //
    //    double F ( double x ), the function whose zero is being sought.
    //
    //    string TITLE, a title for the problem.
    //
    {
      double arg = 0;
      double value = 0;

      Console.WriteLine("");
      Console.WriteLine("  " + title + "");
      Console.WriteLine("");
      Console.WriteLine("    STATUS      X               F(X)");
      Console.WriteLine("");

      double t = Double.Epsilon;

      int status = 0;

      zerorc_data data = new zerorc_data();
      
      for (;;)
      {
        zerorc.zero_rc(ref data, a, b, t, ref arg, ref status, ref value);

        if (status < 0)
        {
          Console.WriteLine("");
          Console.WriteLine("  zero_rc() returned an error flag!");
          break;
        }

        value = f(arg);

        Console.WriteLine("  " + status.ToString().PadLeft(8)
          + "  " + arg.ToString("0.################").PadLeft(24)
          + "  " + value.ToString("0.######").PadLeft(14) + "");

        if (status == 0)
        {
          break;
        }
      }
    }

    static double f_01(double x)

      //****************************************************************************80
      //
      //  Purpose:
      //
      //    f_01() evaluates sin ( x ) - x / 2.
      //
      //  Licensing:
      //
      //    This code is distributed under the GNU LGPL license.
      //
      //  Modified:
      //
      //    13 April 2008
      //
      //  Author:
      //
      //    John Burkardt
      //
      //  Input:
      //
      //    double X, the evaluation point.
      //
      //  Output:
      //
      //    double F_01, the value of the function at X.
      //
    {
      double value = Math.Sin(x) - 0.5 * x;

      return value;
    }

    static double f_02(double x)

      //****************************************************************************80
      //
      //  Purpose:
      //
      //    f_02() evaluates 2*x-exp(-x).
      //
      //  Licensing:
      //
      //    This code is distributed under the GNU LGPL license.
      //
      //  Modified:
      //
      //    13 April 2008
      //
      //  Author:
      //
      //    John Burkardt
      //
      //  Input:
      //
      //    double X, the evaluation point.
      //
      //  Output:
      //
      //    double F_02, the value of the function at X.
      //
    {
      double value = 2.0 * x - Math.Exp(-x);

      return value;
    }

    static double f_03(double x)

      //****************************************************************************80
      //
      //  Purpose:
      //
      //    f_03() evaluates x*exp(-x).
      //
      //  Licensing:
      //
      //    This code is distributed under the GNU LGPL license.
      //
      //  Modified:
      //
      //    13 April 2008
      //
      //  Author:
      //
      //    John Burkardt
      //
      //  Input:
      //
      //    double X, the evaluation point.
      //
      //  Output:
      //
      //    double F_03, the value of the function at X.
      //
    {
      double value = x * Math.Exp(-x);

      return value;
    }

    static double f_04(double x)

      //****************************************************************************80
      //
      //  Purpose:
      //
      //    f_04() evaluates exp(x) - 1 / (100*x*x).
      //
      //  Licensing:
      //
      //    This code is distributed under the GNU LGPL license.
      //
      //  Modified:
      //
      //    13 April 2008
      //
      //  Author:
      //
      //    John Burkardt
      //
      //  Input:
      //
      //    double X, the evaluation point.
      //
      //  Output:
      //
      //    double F_04, the value of the function at X.
      //
    {
      double value = Math.Exp(x) - 1.0 / 100.0 / x / x;

      return value;
    }

    static double f_05(double x)

      //****************************************************************************80
      //
      //  Purpose:
      //
      //    f_05() evaluates (x+3)*(x-1)*(x-1).
      //
      //  Licensing:
      //
      //    This code is distributed under the GNU LGPL license.
      //
      //  Modified:
      //
      //    13 April 2008
      //
      //  Author:
      //
      //    John Burkardt
      //
      //  Input:
      //
      //    double X, the evaluation point.
      //
      //  Output:
      //
      //    double F_05, the value of the function at X.
      //
    {
      double value = (x + 3.0) * (x - 1.0) * (x - 1.0);

      return value;
    }
  }
}