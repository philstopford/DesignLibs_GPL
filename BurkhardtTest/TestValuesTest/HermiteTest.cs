using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class HermiteTest
{
    public static void hermite_function_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_FUNCTION_VALUES_TEST tests HERMITE_FUNCTION_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("HERMITE_FUNCTION_VALUES_TEST");
        Console.WriteLine("  HERMITE_FUNCTION_VALUES stores values of");
        Console.WriteLine("  the Hermite function.");
        Console.WriteLine("");
        Console.WriteLine("     N      X            Hf(N,X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Hermite.hermite_function_values(ref n_data, ref n, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString().PadLeft(6) + "  "
                              + x.ToString().PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void hermite_poly_phys_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_POLY_PHYS_VALUES_TEST tests HERMITE_POLY_PHYS_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("HERMITE_POLY_PHYS_VALUES_TEST");
        Console.WriteLine("  HERMITE_POLY_PHYS_VALUES stores values of");
        Console.WriteLine("  the physicist's Hermite polynomials.");
        Console.WriteLine("");
        Console.WriteLine("     N      X            H(N,X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Hermite.hermite_poly_phys_values(ref n_data, ref n, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString().PadLeft(6) + "  "
                              + x.ToString().PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void hermite_poly_prob_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_POLY_PROB_VALUES_TEST tests HERMITE_POLY_PROB_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("HERMITE_POLY_PROB_VALUES_TEST");
        Console.WriteLine("  HERMITE_POLY_PROB_VALUES stores values of");
        Console.WriteLine("  the probabilist's Hermite polynomials.");
        Console.WriteLine("");
        Console.WriteLine("     N      X            He(N,X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Hermite.hermite_poly_prob_values(ref n_data, ref n, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString().PadLeft(6) + "  "
                              + x.ToString().PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

}