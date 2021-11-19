using System;
using Burkardt.Sequence;

namespace VanderCorputAdvancedTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for VAN_DER_CORPUT_ADVANCED_TEST.
        //
        //  Discussion:
        //
        //    VAN_DER_CORPUT_ADVANCED_TEST tests the VAN_DER_CORPUT_ADVANCED library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("VAN_DER_CORPUT_ADVANCED_TEST");
        Console.WriteLine("  Test the VAN_DER_CORPUT_ADVANCED library.");

        test01();
        test02();
        test03();
        test04();
        test045();
        test05();
        test06();
        test09();

        Console.WriteLine("");
        Console.WriteLine("VAN_DER_CORPUT_ADVANCED_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests VAN_DER_CORPUT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        double r;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  VAN_DER_CORPUT computes the elements of a");
        Console.WriteLine("  van der Corput sequence.");
        Console.WriteLine("  Each call produces the next value.  By default,");
        Console.WriteLine("  the base is 2, and the sequence starts at element 1.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we call VAN_DER_CORPUT several times.");
        Console.WriteLine("");
        Console.WriteLine("  Seed   van der Corput");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            seed = VanDerCorput.van_der_corput_seed_get();

            r = VanDerCorput.van_der_corput();

            Console.WriteLine(seed.ToString().PadLeft(6) + "  "
                                                         + r.ToString().PadLeft(10) + "");
        }

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests VAN_DER_CORPUT_SEQUENCE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 0;

        int base_;
        int i;
        double[] r = new double[N];
        int seed;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  VAN_DER_CORPUT_SEQUENCE computes several elements of");
        Console.WriteLine("  a van der Corput sequence on a single call.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we call VAN_DER_CORPUT_SEQUENCE once.");

        base_ = 2;
        VanDerCorput.van_der_corput_base_set(base_);

        seed = 0;
        VanDerCorput.van_der_corput_seed_set(seed);

        VanDerCorput.van_der_corput_sequence(N, ref r);

        Console.WriteLine("");
        Console.WriteLine("  Element   van der Corput");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            Console.WriteLine(i.ToString().PadLeft(6) + "  "
                                                      + r[i].ToString().PadLeft(10) + "");
        }

    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests VAN_DER_CORPUT_SEQUENCE, VAN_DER_CORPUT_SEED_GET, VAN_DER_CORPUT_SEED_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NMAX = 10;

        int i;
        int n;
        double[] r = new double[NMAX];
        int seed;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  VAN_DER_CORPUT_SEED_SET specifies the next element of");
        Console.WriteLine("    the van der Corput sequence to compute.");
        Console.WriteLine("  VAN_DER_CORPUT_SEED_GET reports the next element of the");
        Console.WriteLine("    van der Corput sequence that will be computed.");
        Console.WriteLine("");
        Console.WriteLine("  By default, the sequence starts at element 1.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we demonstrate computing elements");
        Console.WriteLine("  affects the seed, and how resetting the seed determines");
        Console.WriteLine("  the next element computed.");
        Console.WriteLine("");
        Console.WriteLine("  We start at element 0 and compute 10 elements.");
        Console.WriteLine("");

        seed = 0;
        VanDerCorput.van_der_corput_seed_set(seed);

        n = NMAX;
        VanDerCorput.van_der_corput_sequence(n, ref r);

        for (i = 0; i < n; i++)
        {
            Console.WriteLine((seed + i).ToString().PadLeft(6) + " "
                                                               + r[i].ToString().PadLeft(10) + "");
        }

        seed = VanDerCorput.van_der_corput_seed_get();

        Console.WriteLine("");
        Console.WriteLine("  The current seed is " + seed + "");

        Console.WriteLine("");
        Console.WriteLine("  We jump back to element 6 and compute 10 elements.");
        Console.WriteLine("");

        seed = 6;
        VanDerCorput.van_der_corput_seed_set(seed);

        n = NMAX;
        VanDerCorput.van_der_corput_sequence(n, ref r);

        for (i = 0; i < n; i++)
        {
            Console.WriteLine((seed + i).ToString().PadLeft(6) + " "
                                                               + r[i].ToString().PadLeft(10) + "");
        }

        seed = VanDerCorput.van_der_corput_seed_get();

        Console.WriteLine("");
        Console.WriteLine("  The current seed is " + seed + "");

        Console.WriteLine("");
        Console.WriteLine("  We restart at element 0 and compute 6 elements.");
        Console.WriteLine("");

        seed = 0;
        VanDerCorput.van_der_corput_seed_set(seed);

        n = 6;
        VanDerCorput.van_der_corput_sequence(n, ref r);

        for (i = 0; i < n; i++)
        {
            Console.WriteLine((seed + i).ToString().PadLeft(6) + " "
                                                               + r[i].ToString().PadLeft(10) + "");
        }

        seed = VanDerCorput.van_der_corput_seed_get();

        Console.WriteLine("");
        Console.WriteLine("  The current seed is " + seed + "");

        Console.WriteLine("");
        Console.WriteLine("  We jump to element 100 and compute 5 elements.");
        Console.WriteLine("");

        seed = 100;
        VanDerCorput.van_der_corput_seed_set(seed);

        n = 5;
        VanDerCorput.van_der_corput_sequence(n, ref r);

        for (i = 0; i < n; i++)
        {
            Console.WriteLine((seed + i).ToString().PadLeft(6) + " "
                                                               + r[i].ToString().PadLeft(10) + "");
        }

        seed = VanDerCorput.van_der_corput_seed_get();

        Console.WriteLine("");
        Console.WriteLine("  The current seed is " + seed + "");

    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests VAN_DER_CORPUT, VAN_DER_CORPUT_BASE_GET, VAN_DER_CORPUT_BASE_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int base_;
        int i;
        double r;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  VAN_DER_CORPUT_BASE_GET gets the current base.");
        Console.WriteLine("  VAN_DER_CORPUT_BASE_SET sets the current base.");
        Console.WriteLine("");
        Console.WriteLine("  The van der Corput base is usually a prime, but this is");
        Console.WriteLine("  not required.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we compute a van der Corput sequence");
        Console.WriteLine("  with the default base, then change the base,");
        Console.WriteLine("  reset the seed, and recompute the sequence.");

        seed = 0;
        VanDerCorput.van_der_corput_seed_set(seed);

        base_ = VanDerCorput.van_der_corput_base_get();

        Console.WriteLine("");
        Console.WriteLine("  VAN_DER_CORPUT_BASE_GET: Current base is " + base_ + "");
        Console.WriteLine("");
        Console.WriteLine("  Seed   van der Corput");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            seed = VanDerCorput.van_der_corput_seed_get();
            r = VanDerCorput.van_der_corput();

            Console.WriteLine(seed.ToString().PadLeft(6) + "  "
                                                         + r.ToString().PadLeft(10) + "");
        }

        base_ = 3;
        VanDerCorput.van_der_corput_base_set(base_);

        seed = 0;
        VanDerCorput.van_der_corput_seed_set(seed);

        Console.WriteLine("");
        Console.WriteLine("  Reset base to " + base_ + "");
        Console.WriteLine("  Reset seed to " + seed + "");

        Console.WriteLine("");
        Console.WriteLine("  Seed   van der Corput");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            seed = VanDerCorput.van_der_corput_seed_get();
            r = VanDerCorput.van_der_corput();

            Console.WriteLine(seed.ToString().PadLeft(6) + "  "
                                                         + r.ToString().PadLeft(10) + "");
        }

        base_ = 4;
        VanDerCorput.van_der_corput_base_set(base_);

        seed = 0;
        VanDerCorput.van_der_corput_seed_set(seed);

        Console.WriteLine("");
        Console.WriteLine("  Set BASE = " + base_ + "");
        Console.WriteLine("  Set SEED = " + seed + "");

        Console.WriteLine("");
        Console.WriteLine("  Seed   van der Corput");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            seed = VanDerCorput.van_der_corput_seed_get();
            r = VanDerCorput.van_der_corput();

            Console.WriteLine(seed.ToString().PadLeft(6) + "  "
                                                         + r.ToString().PadLeft(10) + "");
        }

    }

    private static void test045()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST045 tests VAN_DER_CORPUT, VAN_DER_CORPUT_SEED_GET, VAN_DER_CORPUT_SEED_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int base_;
        int i;
        double r;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("TEST045");
        Console.WriteLine("  VAN_DER_CORPUT_SEED_GET gets the current seed.");
        Console.WriteLine("  VAN_DER_CORPUT_SEED_SET sets the current seed.");
        Console.WriteLine("");
        Console.WriteLine("  The van der Corput base is usually a prime, but this is");
        Console.WriteLine("  not required.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we compute a van der Corput sequence");
        Console.WriteLine("  starting with the default seed, then check the seed,");
        Console.WriteLine("  reset the seed, and recompute the sequence.");

        base_ = 2;
        VanDerCorput.van_der_corput_base_set(base_);

        seed = 0;
        VanDerCorput.van_der_corput_seed_set(seed);

        Console.WriteLine("");
        Console.WriteLine("  All computations will use base " + base_ + ".");
        Console.WriteLine("");

        seed = VanDerCorput.van_der_corput_seed_get();

        Console.WriteLine("");
        Console.WriteLine("  Set SEED = " + seed + "");
        Console.WriteLine("");
        Console.WriteLine("  Seed   van der Corput");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            seed = VanDerCorput.van_der_corput_seed_get();

            r = VanDerCorput.van_der_corput();

            Console.WriteLine(seed.ToString().PadLeft(6) + "  "
                                                         + r.ToString().PadLeft(10) + "");
        }

        seed = VanDerCorput.van_der_corput_seed_get();

        Console.WriteLine("");
        Console.WriteLine("  Current seed is " + seed + "");
        Console.WriteLine("");

        seed = 100;
        VanDerCorput.van_der_corput_seed_set(seed);

        Console.WriteLine("");
        Console.WriteLine("  Set SEED = " + seed + "");
        Console.WriteLine("");
        Console.WriteLine("  Seed   van der Corput");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            seed = VanDerCorput.van_der_corput_seed_get();

            r = VanDerCorput.van_der_corput();

            Console.WriteLine(seed.ToString().PadLeft(6) + "  "
                                                         + r.ToString().PadLeft(10) + "");
        }

        seed = VanDerCorput.van_der_corput_seed_get();

        Console.WriteLine("");
        Console.WriteLine("  Current seed is " + seed + "");
        Console.WriteLine("");

        seed = 3;
        VanDerCorput.van_der_corput_seed_set(seed);

        Console.WriteLine("");
        Console.WriteLine("  Reset seed to " + seed + "");
        Console.WriteLine("");
        Console.WriteLine("  Seed   van der Corput");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            seed = VanDerCorput.van_der_corput_seed_get();

            r = VanDerCorput.van_der_corput();

            Console.WriteLine(seed.ToString().PadLeft(6) + "  "
                                                         + r.ToString().PadLeft(10) + "");
        }

        ;
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests I4_TO_VAN_DER_CORPUT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int base_;
        double r;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  I4_TO_VAN_DER_CORPUT returns the I-th element");
        Console.WriteLine("  of a van der Corput sequence to a given base.");
        Console.WriteLine("");
        Console.WriteLine("");
        Console.WriteLine("  Base    Seed   R");
        Console.WriteLine("");

        for (base_ = 2; base_ <= 5; base_++)
        {

            Console.WriteLine("");
            Console.WriteLine(base_.ToString().PadLeft(6) + "");

            for (seed = 0; seed <= 10; seed++)
            {
                r = VanDerCorput.i4_to_van_der_corput(seed, base_);

                Console.WriteLine("        "
                                  + seed.ToString().PadLeft(6) + "  "
                                  + r.ToString().PadLeft(10) + "");
            }

        }
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests I4_TO_VAN_DER_CORPUT_SEQUENCE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;

        int base_;
        int i;
        double[] r = new double[N];
        int seed;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  I4_TO_VAN_DER_CORPUT_SEQUENCE returns N elements");
        Console.WriteLine("  of a van der Corput sequence to a given base.");
        Console.WriteLine("");
        Console.WriteLine("");
        Console.WriteLine("  Base    Seed   R");
        Console.WriteLine("");

        for (base_ = 2; base_ <= 5; base_++)
        {

            Console.WriteLine("");
            Console.WriteLine(base_.ToString().PadLeft(6) + "");

            seed = 0;

            VanDerCorput.i4_to_van_der_corput_sequence(seed, base_, N, ref r);

            for (i = 0; i < N; i++)
            {

                Console.WriteLine("        "
                                  + (seed + i).ToString().PadLeft(6) + "  "
                                  + r[i].ToString().PadLeft(10) + "");
            }

        }
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests VDC_NUMERATOR_SEQUENCE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n;
        int[] r;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  VDC_NUMERATOR_SEQUENCE returns N elements");
        Console.WriteLine("  of a van der Corput numerator sequence in base 2.");
        Console.WriteLine("");
        Console.WriteLine("   N:  Sequence");
        Console.WriteLine("");

        for (n = 1; n <= 20; n++)
        {
            r = VanDerCorput.vdc_numerator_sequence(n);
            string cout = "  " + n.ToString().PadLeft(2) + ":";
            for (i = 0; i < n; i++)
            {
                cout += "  " + r[i].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);

        }
    }
}