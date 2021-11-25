using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class CosineTest
{
    public static void ci_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CI_VALUES_TEST tests CI_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("CI_VALUES_TEST:");
        Console.WriteLine("  CI_VALUES stores values of");
        Console.WriteLine("  the Cosine Integral function CI(X).");
        Console.WriteLine("");
        Console.WriteLine("      X            CI(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Cosine.ci_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void cin_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIN_VALUES_TEST tests CIN_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("CIN_VALUES_TEST:");
        Console.WriteLine("  CIN_VALUES stores values of");
        Console.WriteLine("  the Cosine Integral function CIN(X).");
        Console.WriteLine("");
        Console.WriteLine("      X            CIN(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Cosine.cin_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void cinh_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CINH_VALUES_TEST tests CINH_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("CINH_VALUES_TEST:");
        Console.WriteLine("  CINH_VALUES stores values of");
        Console.WriteLine("  the Hyperbolic Cosine Integral function CINH(X).");
        Console.WriteLine("");
        Console.WriteLine("      X            CINH(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Cosine.cinh_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void cos_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COS_VALUES_TEST tests COS_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("COS_VALUES_TEST:");
        Console.WriteLine("   COS_VALUES stores values of the cosine function.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Cosine.cos_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void cos_degree_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COS_DEGREE_VALUES_TEST tests COS_DEGREE_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("COS_DEGREE_VALUES_TEST:");
        Console.WriteLine("   COS_DEGREE_VALUES stores values of the cosine function.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Cosine.cos_degree_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + x + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void cos_power_int_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COS_POWER_INT_VALUES_TEST tests COS_POWER_INT_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double fx = 0;
        int n = 0;
        Console.WriteLine("");
        Console.WriteLine("COS_POWER_INT_VALUES_TEST:");
        Console.WriteLine("  COS_POWER_INT_VALUES returns values of");
        Console.WriteLine("  the integral of the N-th power of the cosine function.");
        Console.WriteLine("");
        Console.WriteLine("         A         B       N        FX");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Cosine.cos_power_int_values(ref n_data, ref a, ref b, ref n, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + a.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + b.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void cosh_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COSH_VALUES_TEST tests COSH_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("COSH_VALUES_TEST:");
        Console.WriteLine("   COSH_VALUES stores values of the hyperbolic cosine function.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Cosine.cosh_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

}