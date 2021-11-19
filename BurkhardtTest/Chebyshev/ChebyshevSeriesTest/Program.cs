using System;
using Burkardt;
using Burkardt.Chebyshev;

namespace ChebyshevSeriestTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CHEBYSHEV_SERIES_TEST.
        //
        //  Discussion:
        //
        //    CHEBYSHEV_SERIES_TEST tests the CHEBYSHEV_SERIES library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 January 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV_SERIES_TEST:");
        Console.WriteLine("  Test the CHEBYSHEV_SERIES libary.");

        test01();
        test02();
        test03();

        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV_SERIES_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 considers an even Chebyshev series for EXP(X).
        //
        //  Discussion:
        //
        //    Table 5 is from Clenshaw, and contains 18 terms of the Chebyshev
        //    series for exp(x) over [-1,+1].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 January 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
        //  Reference:
        //
        //    Charles Clenshaw,
        //    Mathematical Tables, Volume 5,
        //    Chebyshev series for mathematical functions,
        //    London, 1962.
        //
    {
        double[] table5 =  {
                2.53213175550401667120,
                1.13031820798497005442,
                0.27149533953407656237,
                0.04433684984866380495,
                0.00547424044209373265,
                0.00054292631191394375,
                0.00004497732295429515,
                0.00000319843646240199,
                0.00000019921248066728,
                0.00000001103677172552,
                0.00000000055058960797,
                0.00000000002497956617,
                0.00000000000103915223,
                0.00000000000003991263,
                0.00000000000000142376,
                0.00000000000000004741,
                0.00000000000000000148,
                0.00000000000000000004
            }
            ;
        double s1 = 0;
        double s2 = 0;
        double s3 = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  ECHEBSER3 computes a Chebyshev series approximation");
        Console.WriteLine("  and the first three derivatives.");
        Console.WriteLine("");
        Console.WriteLine("  Errors of a Chebyshev series for exp(x)");
        Console.WriteLine("");
        Console.WriteLine("    x        err(y)       err(y')       err(y\")      err(y\"')");
        Console.WriteLine("");

        for (int i = -10; i <= 10; i++)
        {
            double x = i / 10.0;
            double s = ChebyshevSeries.echebser3(x, table5, 18, ref s1, ref s2, ref s3);
            double y = Math.Exp(x);
            s -= y;
            s1 -= y;
            s2 -= y;
            s3 -= y;

            Console.WriteLine(x.ToString().PadLeft(5)
                              + s.ToString().PadLeft(14)
                              + s1.ToString().PadLeft(14)
                              + s2.ToString().PadLeft(14)
                              + s3.ToString().PadLeft(14) + "");
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 considers an even Chebyshev series for COS(PI*X/2).
        //
        //  Discussion:
        //
        //    TABLE1 contains the even Chebyshev series coefficients for
        //    cos(pi*x/2) over [-1,+1].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 January 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
        //  Reference:
        //
        //    Charles Clenshaw,
        //    Mathematical Tables, Volume 5,
        //    Chebyshev series for mathematical functions,
        //    London, 1962.
        //
    {
        double s;
        double s1 = 0;
        double s2 = 0;
        double[] table1 =  {
                +0.94400243153646953490,
                -0.49940325827040708740,
                +0.02799207961754761751,
                -0.00059669519654884650,
                +0.00000670439486991684,
                -0.00000004653229589732,
                +0.00000000021934576590,
                -0.00000000000074816487,
                +0.00000000000000193230,
                -0.00000000000000000391,
                +0.00000000000000000001
            }
            ;
        double y = 0;
        double y1 = 0;
        double y2 = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  EVENCHEBSER2 computes an even Chebyshev series");
        Console.WriteLine("  and its first two derivatives.");
        Console.WriteLine("");
        Console.WriteLine("  Errors of an even Chebyshev series for cos(pi*x/2):");
        Console.WriteLine("");
        Console.WriteLine("    x        err(y)       err(y')       err(y\")");
        Console.WriteLine("");

        for (int i = 0; i <= 10; i++)
        {
            double x = i / 10.0;
            s = ChebyshevSeries.evenchebser2(x, table1, 11, ref s1, ref s2);
            Helpers.sincos(Math.PI/2 * x, ref y1, ref y);
            y1 = -y1 * Math.PI/2;
            y2 = -y * (Math.PI/2* Math.PI/2);
            s -= y;
            s1 -= y1;
            s2 -= y2;

            Console.WriteLine(x.ToString().PadLeft(5)
                              + s.ToString().PadLeft(14)
                              + s1.ToString().PadLeft(14)
                              + s2.ToString().PadLeft(14) + "");
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 considers an odd Chebyshev series for SINH(X).
        //
        //  Discussion:
        //
        //    TABLE5ODD contains the odd Chebyshev series coefficients for
        //    sinh(x) over -1 <= x <= 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 January 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
        //  Reference:
        //
        //    Charles Clenshaw,
        //    Mathematical Tables, Volume 5,
        //    Chebyshev series for mathematical functions,
        //    London, 1962.
        //
    {
        double s1 = 0;
        double s2 = 0;
        double[] table5odd =  {
                1.13031820798497005442,
                0.04433684984866380495,
                0.00054292631191394375,
                0.00000319843646240199,
                0.00000001103677172552,
                0.00000000002497956617,
                0.00000000000003991263,
                0.00000000000000004741,
                0.00000000000000000004
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST03:");
        Console.WriteLine("  ODDCHEBSER2 computes an odd Chebyshev series approximation.");
        Console.WriteLine("  and its first two derivatives.");
        Console.WriteLine("");
        Console.WriteLine("  Errors of an odd Chebyshev series y(x) approximating sinh(x):");
        Console.WriteLine("");
        Console.WriteLine("    x        err(y)       err(y')       err(y\")");
        Console.WriteLine("");

        for (int i = 0; i <= 10; i++)
        {
            double x = i / 10.0;
            double s = ChebyshevSeries.oddchebser2(x, table5odd, 9, ref s1, ref s2);
            double y = Math.Sinh(x);
            double y1 = Math.Cosh(x);
            s -= y;
            s1 -= y1;
            s2 -= y;
            Console.WriteLine(x.ToString().PadLeft(5)
                              + s.ToString().PadLeft(14)
                              + s1.ToString().PadLeft(14)
                              + s2.ToString().PadLeft(14) + "");
        }
    }
}