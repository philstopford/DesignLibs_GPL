using System;
using Burkardt.FullertonFnLib;
using Burkardt.SubsetNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace TOMS515Test;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TOMS515_TEST.
        //
        //  Discussion:
        //
        //    TOMS515_TEST tests the TOMS515 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TOMS515_TEST");
        Console.WriteLine("  Test the TOMS515 library.");

        test01();
        test02();
        test03();
        test04();
        test05();

        Console.WriteLine("");
        Console.WriteLine("TOMS515_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests COMB by generating all 3-subsets of a 5 set.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int k = 3;
        int l;
        const int n = 5;

        int lmax = FullertonLib.i4_binom(n, k);

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Generate all K-subsets of an N set.");
        Console.WriteLine("  K = " + k + "");
        Console.WriteLine("  N = " + n + "");
        Console.WriteLine("  LMAX =" + lmax + "");

        if (!typeMethods.i4_choose_check(n, k))
        {
            Console.WriteLine("");
            Console.WriteLine("TEST01 - Warning!");
            Console.WriteLine("  The binomial coefficient cannot be");
            Console.WriteLine("  computed in integer arithmetic for");
            Console.WriteLine("  this choice of parameters.");
            return;
        }

        Console.WriteLine("");

        for (l = 1; l <= lmax; l++)
        {
            int[] c = Comb.comb(n, k, l);
            string cout = "  " + l.ToString().PadLeft(4) + ":  ";
            int i;
            for (i = 0; i < k; i++)
            {
                cout += "  " + c[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests COMB by generating 10 random 3-subsets of a 10 set.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        const int k = 3;
        const int n = 5;

        int lmax = FullertonLib.i4_binom(n, k);

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Generate 10 random K-subsets of an N set.");
        Console.WriteLine("  K = " + k + "");
        Console.WriteLine("  N = " + n + "");
        Console.WriteLine("  LMAX =" + lmax + "");

        if (!typeMethods.i4_choose_check(n, k))
        {
            Console.WriteLine("");
            Console.WriteLine("TEST02 - Warning!");
            Console.WriteLine("  The binomial coefficient cannot be");
            Console.WriteLine("  computed in integer arithmetic for");
            Console.WriteLine("  this choice of parameters.");
            return;
        }

        Console.WriteLine("");

        int seed = 123456789;

        for (j = 1; j <= 10; j++)
        {
            int l = UniformRNG.i4_uniform_ab(1, lmax, ref seed);
            int[] c = Comb.comb(n, k, l);
            string cout = "  " + l.ToString().PadLeft(4) + ":  ";
            int i;
            for (i = 0; i < k; i++)
            {
                cout += "  " + c[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests COMB by generating 10 random 3-subsets of a 25 set.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        const int k = 3;
        const int n = 25;

        int lmax = FullertonLib.i4_binom(n, k);

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Generate 10 random K-subsets of an N set.");
        Console.WriteLine("  K = " + k + "");
        Console.WriteLine("  N = " + n + "");
        Console.WriteLine("  LMAX =" + lmax + "");

        if (!typeMethods.i4_choose_check(n, k))
        {
            Console.WriteLine("");
            Console.WriteLine("TEST03 - Warning!");
            Console.WriteLine("  The binomial coefficient cannot be");
            Console.WriteLine("  computed in integer arithmetic for");
            Console.WriteLine("  this choice of parameters.");
            return;
        }

        Console.WriteLine("");

        int seed = 123456789;

        for (j = 1; j <= 10; j++)
        {
            int l = UniformRNG.i4_uniform_ab(1, lmax, ref seed);
            int[] c = Comb.comb(n, k, l);
            string cout = "  " + l.ToString().PadLeft(4) + ":  ";
            int i;
            for (i = 0; i < k; i++)
            {
                cout += "  " + c[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests COMB by generating 10 random 3-subsets of a 100 set.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        const int k = 3;
        const int n = 100;

        int lmax = FullertonLib.i4_binom(n, k);

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Generate 10 random K-subsets of an N set.");
        Console.WriteLine("  K = " + k + "");
        Console.WriteLine("  N = " + n + "");
        Console.WriteLine("  LMAX =" + lmax + "");

        if (!typeMethods.i4_choose_check(n, k))
        {
            Console.WriteLine("");
            Console.WriteLine("TEST04 - Warning!");
            Console.WriteLine("  The binomial coefficient cannot be");
            Console.WriteLine("  computed in integer arithmetic for");
            Console.WriteLine("  this choice of parameters.");
            return;
        }

        Console.WriteLine("");

        int seed = 123456789;

        for (j = 1; j <= 10; j++)
        {
            int l = UniformRNG.i4_uniform_ab(1, lmax, ref seed);
            int[] c = Comb.comb(n, k, l);
            string cout = "  " + l.ToString().PadLeft(4) + ":  ";
            int i;
            for (i = 0; i < k; i++)
            {
                cout += "  " + c[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests COMB by generating 10 random 10-subsets of a 100 set.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        const int k = 10;
        const int n = 100;

        int lmax = FullertonLib.i4_binom(n, k);

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  Generate 10 random K-subsets of an N set.");
        Console.WriteLine("  K = " + k + "");
        Console.WriteLine("  N = " + n + "");
        Console.WriteLine("  LMAX =" + lmax + "");
        Console.WriteLine("");
        Console.WriteLine("  Note that this function is already");
        Console.WriteLine("  failing because LMAX is negative.");
        Console.WriteLine("  The combinatorial coefficient C(100,10)");
        Console.WriteLine("  is too large to store in an integer.");
        Console.WriteLine("");
        Console.WriteLine("  Although the program continues to give");
        Console.WriteLine("  results, they cannot be relied on.");

        if (!typeMethods.i4_choose_check(n, k))
        {
            Console.WriteLine("");
            Console.WriteLine("TEST05 - Warning!");
            Console.WriteLine("  The binomial coefficient cannot be");
            Console.WriteLine("  computed in integer arithmetic for");
            Console.WriteLine("  this choice of parameters.");
            return;
        }

        Console.WriteLine("");

        int seed = 123456789;

        for (j = 1; j <= 10; j++)
        {
            int l = UniformRNG.i4_uniform_ab(1, lmax, ref seed);
            int[] c = Comb.comb(n, k, l);
            string cout = "  " + l.ToString().PadLeft(4) + ":  ";
            int i;
            for (i = 0; i < k; i++)
            {
                cout += "  " + c[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }
}