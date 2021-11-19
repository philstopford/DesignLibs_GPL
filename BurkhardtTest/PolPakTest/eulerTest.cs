using System;

namespace PolPakTest;

using Sequence = Burkardt.Sequence.Euler;
using Polynomial = Burkardt.PolynomialNS.Euler;

public static class eulerTest
{
    public static void euler_number_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EULER_NUMBER_TEST tests EULER_NUMBER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int c1 = 0;
        int[] c2 = new int[13];
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("EULER_NUMBER_TEST");
        Console.WriteLine("  EULER_NUMBER computes Euler numbers.");
        Console.WriteLine("");
        Console.WriteLine("  N  exact   EULER_NUMBER");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Euler.euler_number_values(ref n_data, ref n, ref c1);

            if (n_data == 0)
            {
                break;
            }

            Sequence.euler_number(n, ref c2);

            Console.WriteLine("  "
                              + n.ToString().PadLeft(4) + "  "
                              + c1.ToString().PadLeft(12) + "  "
                              + c2[n].ToString().PadLeft(12) + "");

        }

    }

    public static void euler_number2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EULER_NUMBER2_TEST tests EULER_NUMBER2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int c1 = 0;
        double c2 = 0;
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("EULER_NUMBER2_TEST");
        Console.WriteLine("  EULER_NUMBER2 computes Euler numbers.");
        Console.WriteLine("");
        Console.WriteLine("  N  exact   EULER_NUMBER2");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Euler.euler_number_values(ref n_data, ref n, ref c1);

            if (n_data == 0)
            {
                break;
            }

            c2 = Sequence.euler_number2(n);

            Console.WriteLine("  "
                              + n.ToString().PadLeft(4) + "  "
                              + c1.ToString().PadLeft(12) + "  "
                              + c2.ToString().PadLeft(14) + "");

        }

    }

    public static void euler_poly_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EULER_POLY_TEST tests EULER_POLY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double f;
        int i;
        int n = 15;
        double x;

        x = 0.5;

        Console.WriteLine("");
        Console.WriteLine("EULER_POLY_TEST");
        Console.WriteLine("  EULER_POLY evaluates Euler polynomials.");
        Console.WriteLine("");
        Console.WriteLine("  N         X              F(X)");
        Console.WriteLine("");

        for (i = 0; i <= n; i++)
        {
            f = Polynomial.euler_poly(i, x);

            Console.WriteLine("  "
                              + i.ToString().PadLeft(2) + "  "
                              + x.ToString().PadLeft(14) + "  "
                              + f.ToString().PadLeft(14) + "");
        }

    }

    public static void eulerian_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EULERIAN_TEST tests EULERIAN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 7;

        int[] e = new int[N * N];
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("EULERIAN_TEST");
        Console.WriteLine("  EULERIAN evaluates Eulerian numbers.");
        Console.WriteLine("");

        Sequence.eulerian(N, ref e);

        for (i = 0; i < N; i++)
        {
            string cout = "";
            for (j = 0; j < N; j++)
            {
                cout += e[i + j * N].ToString().PadLeft(6) + "  ";
            }

            Console.WriteLine(cout);
        }

    }

}