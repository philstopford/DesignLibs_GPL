using System;
using Burkardt.Cycle;

namespace CycleBrentTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CYCLE_BRENT_TEST.
        //
        //  Discussion:
        //
        //    CYCLE_BRENT_TEST tests the CYCLE_BRENT library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CYCLE_BRENT_TEST");
            
        Console.WriteLine("  Test the CYCLE_BRENT library.");

        test01();
        test02();
        test03();
        test04();
        test05();

        Console.WriteLine("");
        Console.WriteLine("CYCLE_BRENT_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests CYCLE_BRENT for a tiny example.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lam = 0;
        int mu = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Test CYCLE_BRENT on F1().");
        Console.WriteLine("  f1(0) = 6.");
        Console.WriteLine("  f1(1) = 6.");
        Console.WriteLine("  f1(2) = 0.");
        Console.WriteLine("  f1(3) = 1.");
        Console.WriteLine("  f1(4) = 4.");
        Console.WriteLine("  f1(5) = 3.");
        Console.WriteLine("  f1(6) = 3.");
        Console.WriteLine("  f1(7) = 4.");
        Console.WriteLine("  f1(8) = 0.");

        const int x0 = 2;
        Console.WriteLine("");
        Console.WriteLine("  Starting argument X0 = " + x0 + "");

        Brent.cycle_brent(f1, x0, ref lam, ref mu);

        Console.WriteLine("");
        Console.WriteLine("  Reported cycle length is " + lam + "");
        Console.WriteLine("  Expected value is 3");
        Console.WriteLine("");
        Console.WriteLine("  Reported distance to first cycle element is " + mu + "");
        Console.WriteLine("  Expected value is 2");
    }

    private static int f1(int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F1 is the iteration function for example 1.
        //
        //  Discussion:
        //
        //    This function has two cycles:
        //
        //    6, 3, 1, of length 3
        //    4, of length 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the argument of the function.
        //
        //    Output, int F1, the value of the function.
        //
    {
        int[] f_table = {
                6, 6, 0, 1, 4, 3, 3, 4, 0
            }
            ;

        int value = f_table[i];

        return value;
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests CYCLE_BRENT for F2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lam = 0;
        int mu = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Test CYCLE_BRENT for F2().");
        Console.WriteLine("  f2(i) = mod ( 22 * i + 1, 72 ).");

        const int x0 = 0;
        Console.WriteLine("");
        Console.WriteLine("  Starting argument X0 = " + x0 + "");

        Brent.cycle_brent(f2, x0, ref lam, ref mu);

        Console.WriteLine("");
        Console.WriteLine("  Reported cycle length is " + lam + "");
        Console.WriteLine("  Expected value is 9");
        Console.WriteLine("");
        Console.WriteLine("  Reported distance to first cycle element is " + mu + "");
        Console.WriteLine("  Expected value is 3");
    }

    private static int f2(int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F2 is the iteration function for example 2.
        //
        //  Discussion:
        //
        //    This function has a cycle
        //
        //    3, 67, 35, 51, 43, 11, 27, 19, 59, of length 9
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the argument of the function.
        //
        //    Output, int F2, the value of the function.
        //
    {
        int value = (22 * i + 1) % 72;

        return value;
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests CYCLE_BRENT for F3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lam = 0;
        int mu = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Test CYCLE_BRENT for F3().");
        Console.WriteLine("  f3(i) = mod ( 123 * i + 456, 100000 ).");

        const int x0 = 789;
        Console.WriteLine("");
        Console.WriteLine("  Starting argument X0 = " + x0 + "");

        Brent.cycle_brent(f3, x0, ref lam, ref mu);

        Console.WriteLine("");
        Console.WriteLine("  Reported cycle length is " + lam + "");
        Console.WriteLine("  Expected value is 50000");
        Console.WriteLine("");
        Console.WriteLine("  Reported distance to first cycle element is " + mu + "");
        Console.WriteLine("  Expected value is 0");
    }

    private static int f3(int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F3 is the iteration function for example 3.
        //
        //  Discussion:
        //
        //    This function has a cycle of length 50000
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the argument of the function.
        //
        //    Output, int F3, the value of the function.
        //
    {
        int value = (123 * i + 456) % 1000000;

        return value;
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests CYCLE_BRENT for F4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lam = 0;
        int mu = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Test CYCLE_BRENT for F4().");
        Console.WriteLine("  f4(i) = mod ( 31421 * i + 6927, 65536 ).");

        const int x0 = 1;
        Console.WriteLine("");
        Console.WriteLine("  Starting argument X0 = " + x0 + "");

        Brent.cycle_brent(f4, x0, ref lam, ref mu);

        Console.WriteLine("");
        Console.WriteLine("  Reported cycle length is " + lam + "");
        Console.WriteLine("  Expected value is 65536");
        Console.WriteLine("");
        Console.WriteLine("  Reported distance to first cycle element is " + mu + "");
        Console.WriteLine("  Expected value is 0");
    }

    private static int f4(int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F4 is the iteration function for example 4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the argument of the function.
        //
        //    Output, int F4, the value of the function.
        //
    {
        int value = (31421 * i + 6927) % 65536;

        return value;
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests CYCLE_BRENT for F5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lam = 0;
        int mu = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  Test CYCLE_BRENT for F5().");
        Console.WriteLine("  f5(i) = mod ( 16383 * i + 1, 65536 ).");

        int x0 = 1;
        Console.WriteLine("");
        Console.WriteLine("  Starting argument X0 = " + x0 + "");

        Brent.cycle_brent(f5, x0, ref lam, ref mu);

        Console.WriteLine("");
        Console.WriteLine("  Reported cycle length is " + lam + "");
        Console.WriteLine("  Expected value is 8");
        Console.WriteLine("");
        Console.WriteLine("  Reported distance to first cycle element is " + mu + "");
        Console.WriteLine("  Expected value is 0");

        int i = 0;
        x0 = 1;
        Console.WriteLine("  " + i + "  " + x0 + "");
        for (i = 1; i <= 10; i++)
        {
            x0 = f5(x0);
            Console.WriteLine("  " + i + "  " + x0 + "");
        }
    }

    private static int f5(int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F5 is the iteration function for example 5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the argument of the function.
        //
        //    Output, int F5, the value of the function.
        //
    {
        int value = (16383 * i + 1) % 65536;

        return value;
    }
}