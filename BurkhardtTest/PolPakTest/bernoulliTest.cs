using System;
using Burkardt.PolynomialNS;

namespace PolPakTest;

public static class bernoulliTest
{
    public static void bernoulli_number_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNOULLI_NUMBER_TEST tests BERNOULLI_NUMBER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double c0 = 0;
        double[] c1 = new double[31];
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("BERNOULLI_NUMBER_TEST");
        Console.WriteLine("  BERNOULLI_NUMBER computes Bernoulli numbers;");

        Console.WriteLine("");
        Console.WriteLine("   I      Exact     BERNOULLI_NUMBER");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Bernoulli.bernoulli_number_values(ref n_data, ref n, ref c0);

            if (n_data == 0)
            {
                break;
            }

            Burkardt.Sequence.Bernoulli.bernoulli_number(n, ref c1);

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                              + c0.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + c1[n].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }

    public static void bernoulli_number2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNOULLI_NUMBER2_TEST tests BERNOULLI_NUMBER2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double c0 = 0;
        double[] c1 = new double[31];
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("BERNOULLI_NUMBER2_TEST");
        Console.WriteLine("  BERNOULLI_NUMBER2 computes Bernoulli numbers;");
        Console.WriteLine("");
        Console.WriteLine("   I      Exact     BERNOULLI_NUMBER2");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Bernoulli.bernoulli_number_values(ref n_data, ref n, ref c0);

            if (n_data == 0)
            {
                break;
            }

            Burkardt.Sequence.Bernoulli.bernoulli_number2(n, ref c1);

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                              + c0.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + c1[n].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }

    public static void bernoulli_number3_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNOULLI_NUMBER3_TEST tests BERNOULLI_NUMBER3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double c0 = 0;
        double c1 = 0;
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("BERNOULLI_NUMBER3_TEST");
        Console.WriteLine("  BERNOULLI_NUMBER3 computes Bernoulli numbers.");
        Console.WriteLine("");
        Console.WriteLine("   I      Exact     BERNOULLI_NUMBER3");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Bernoulli.bernoulli_number_values(ref n_data, ref n, ref c0);

            if (n_data == 0)
            {
                break;
            }

            c1 = Burkardt.Sequence.Bernoulli.bernoulli_number3(n);

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                              + c0.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                              + c1.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        }

    }

    public static void bernoulli_poly_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNOULLI_POLY_TEST tests BERNOULLI_POLY;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double bx;
        int i;
        int n = 15;
        double x;

        x = 0.2;

        Console.WriteLine("");
        Console.WriteLine("BERNOULLI_POLY_TEST");
        Console.WriteLine("  BERNOULLI_POLY evaluates Bernoulli polynomials;");
        Console.WriteLine("");
        Console.WriteLine("  X = " + x + "");
        Console.WriteLine("");
        Console.WriteLine("  I          BX");
        Console.WriteLine("");

        for (i = 1; i <= n; i++)
        {
            bx = Bernoulli.bernoulli_poly(i, x);

            Console.WriteLine("  "
                              + i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + bx.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }

    public static void bernoulli_poly2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNOULLI_POLY2_TEST tests BERNOULLI_POLY2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double bx;
        int i;
        int n = 15;
        double x;

        x = 0.2;

        Console.WriteLine("");
        Console.WriteLine("BERNOULLI_POLY2_TEST");
        Console.WriteLine("  BERNOULLI_POLY2 evaluates Bernoulli polynomials.");
        Console.WriteLine("");
        Console.WriteLine("  X = " + x + "");
        Console.WriteLine("");
        Console.WriteLine("  I          BX");
        Console.WriteLine("");

        for (i = 1; i <= n; i++)
        {
            bx = Bernoulli.bernoulli_poly2(i, x);

            Console.WriteLine("  "
                              + i.ToString(CultureInfo.InvariantCulture).PadLeft(2) + "  "
                              + bx.ToString(CultureInfo.InvariantCulture).PadLeft(16) + "");
        }

    }
        
    public static void poly_bernoulli_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLY_BERNOULLI_TEST tests POLY_BERNOULLI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int b;
        int k;
        int n;

        Console.WriteLine("");
        Console.WriteLine("POLY_BERNOULLI_TEST");
        Console.WriteLine("  POLY_BERNOULLI computes the poly-Bernoulli numbers");
        Console.WriteLine("  of negative index, B_n^(-k)");
        Console.WriteLine("");
        Console.WriteLine("   N   K    B_N^(-K)");
        Console.WriteLine("");

        for ( k = 0; k <= 6; k++ )
        {
            Console.WriteLine("");
            for ( n = 0; n <= 6; n++ )
            {
                b = Bernoulli.poly_bernoulli ( n, k );

                Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }

    }

}