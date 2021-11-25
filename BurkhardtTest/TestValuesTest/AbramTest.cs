using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class AbramsTest
{
  public static void abram0_values_test()
    //****************************************************************************80
    //
    //  Purpose:
    //
    //    ABRAM0_VALUES_TEST tests ABRAM0_VALUES.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //
    //    30 June 2006
    //
    //  Author:
    //
    //    John Burkardt
    //
  {
    double fx = 0;
    double x = 0;
    Console.WriteLine("");
    Console.WriteLine("ABRAM0_VALUES_TEST:");
    Console.WriteLine("  ABRAM0_VALUES stores values of ");
    Console.WriteLine("  the Abramowitz function of order 0.");
    Console.WriteLine("");
    Console.WriteLine("                X                   ABRAM0(X)");
    Console.WriteLine("");
    int n_data = 0;
    for (;;)
    {
      Abrams.abram0_values(ref n_data, ref x, ref fx);
      if (n_data == 0)
      {
        break;
      }

      Console.WriteLine("  "
                        + x.ToString("0.################").PadLeft(27) + "  "
                        + fx.ToString("0.################").PadLeft(27) + "");
    }
  }

  public static void abram1_values_test()
    //****************************************************************************80
    //
    //  Purpose:
    //
    //    ABRAM1_VALUES_TEST tests ABRAM1_VALUES.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //
    //    20 January 2007
    //
    //  Author:
    //
    //    John Burkardt
    //
  {
    double fx = 0;
    double x = 0;
    Console.WriteLine("");
    Console.WriteLine("ABRAM1_VALUES_TEST:");
    Console.WriteLine("  ABRAM1_VALUES stores values of ");
    Console.WriteLine("  the Abramowitz function of order 1.");
    Console.WriteLine("");
    Console.WriteLine("                X                   ABRAM1(X)");
    Console.WriteLine("");
    int n_data = 0;
    for (;;)
    {
      Abrams.abram1_values(ref n_data, ref x, ref fx);
      if (n_data == 0)
      {
        break;
      }

      Console.WriteLine("  "
                        + x.ToString("0.################").PadLeft(27) + "  "
                        + fx.ToString("0.################").PadLeft(27) + "");
    }
  }

  public static void abram2_values_test()
    //****************************************************************************80
    //
    //  Purpose:
    //
    //    ABRAM2_VALUES_TEST tests ABRAM2_VALUES.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //
    //    20 January 2007
    //
    //  Author:
    //
    //    John Burkardt
    //
  {
    double fx = 0;
    double x = 0;
    Console.WriteLine("");
    Console.WriteLine("ABRAM2_VALUES_TEST:");
    Console.WriteLine("  ABRAM2_VALUES stores values of ");
    Console.WriteLine("  the Abramowitz function of order 2.");
    Console.WriteLine("");
    Console.WriteLine("                X                   ABRAM3(X)");
    Console.WriteLine("");
    int n_data = 0;
    for (;;)
    {
      Abrams.abram2_values(ref n_data, ref x, ref fx);
      if (n_data == 0)
      {
        break;
      }

      Console.WriteLine("  "
                        + x.ToString("0.################").PadLeft(27) + "  "
                        + fx.ToString("0.################").PadLeft(27) + "");
    }
  }
}