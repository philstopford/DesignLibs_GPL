using System;
using System.Globalization;
using Burkardt.Ziggurat;

namespace ZigguratTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ZIGGURAT_TEST.
        //
        //  Discussion:
        //
        //    ZIGGURAT_TEST tests the ZIGGURAT library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int sample_num = 1000000;

        Console.WriteLine("");
        Console.WriteLine("ZIGGURAT_TEST");
        Console.WriteLine("  Test the ZIGGURAT library.");
        //
        //  Make sure that SEED controls the sequence, and can restart it.
        //
        shr3_seeded_test1();
        test02();
        test03();
        test04();
        //
        //  Measure the time it takes to generate 10,000 variables.
        //
        shr3_seeded_test2(sample_num);
        test06(sample_num);
        test07(sample_num);
        test08(sample_num);
        //
        //  Sample 10 values of the unsigned integer 32 bit generators.
        //
        test09();
        test10();
        test11();
        shr3_seeded_test3();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("ZIGGURAT_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void shr3_seeded_test1()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SHR3_SEEDED_TEST1 tests SHR3_SEEDED.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int seed;
        int value;

        Console.WriteLine("");
        Console.WriteLine("SHR3_SEEDED_TEST1");
        Console.WriteLine("  SHR3_SEEDED returns pseudorandom uniformly distributed");
        Console.WriteLine("  unsigned 32 bit integers.");
        Console.WriteLine("");
        Console.WriteLine("    Step        Signed      Unsigned   SHR3_SEEDED");
        Console.WriteLine("                  Seed          Seed");
        Console.WriteLine("");

        for (j = 0; j < 3; j++)
        {
            seed = (j % 2) switch
            {
                0 => 123456789,
                _ => 987654321
            };

            Console.WriteLine("");
            Console.WriteLine("  " + 0.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + seed.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                value = SHR3.shr3_seeded(ref seed);

                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                       + "  " + seed.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                       + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests R4_UNI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int seed;
        float value;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  R4_UNI returns pseudorandom uniformly distributed");
        Console.WriteLine("  floats (single precision reals) between 0 and 1.");

        for (j = 0; j < 3; j++)
        {
            seed = (j % 2) switch
            {
                0 => 123456789,
                _ => 987654321
            };

            Console.WriteLine("");
            Console.WriteLine("  " + 0.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + seed.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                value = r4.r4_uni(ref seed);

                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                       + "  " + seed.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                       + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests R4_NOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        float[] fn = new float[128];
        int i;
        int j;
        int[] kn = new int[128];
        int seed;
        float value;
        float[] wn = new float[128];

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  R4_NOR returns pseudorandom normally distributed");
        Console.WriteLine("  floats (single precision reals) between 0 and 1.");

        r4.r4_nor_setup(ref kn, ref fn, ref wn);

        for (j = 0; j < 3; j++)
        {
            seed = (j % 2) switch
            {
                0 => 123456789,
                _ => 987654321
            };

            Console.WriteLine("");
            Console.WriteLine("  " + 0.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + seed.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                value = r4.r4_nor(ref seed, kn, fn, wn);

                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                       + "  " + seed.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                       + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests R4_EXP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        float[] fe = new float[256];
        int i;
        int j;
        int[] ke = new int[256];
        int seed;
        float value;
        float[] we = new float[256];

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  R4_EXP returns pseudorandom exponentially distributed");
        Console.WriteLine("  floats (single precision reals) between 0 and 1.");

        r4.r4_exp_setup(ref ke, ref fe, ref we);

        for (j = 0; j < 3; j++)
        {
            seed = (j % 2) switch
            {
                0 => 123456789,
                _ => 987654321
            };

            Console.WriteLine("");
            Console.WriteLine("  " + 0.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + seed.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                value = r4.r4_exp(ref seed, ke, fe, we);

                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                       + "  " + seed.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                       + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    private static void shr3_seeded_test2(int sample_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SHR3_SEEDED_TEST2 times SHR3_SEEDED.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        DateTime ctime;
        int sample;
        int seed;
        int value = 0;

        Console.WriteLine("");
        Console.WriteLine("SHR3_SEEDED_TEST2");
        Console.WriteLine("  Measure the time it takes SHR3_SEEDED to generate");
        Console.WriteLine("  " + sample_num + " unsigned 32 bit integers.");

        seed = 123456789;

        ctime = DateTime.Now;

        for (sample = 0; sample < sample_num; sample++)
        {
            value = SHR3.shr3_seeded(ref seed);
        }

        TimeSpan ctime_elapsed = DateTime.Now - ctime;

        Console.WriteLine("");
        Console.WriteLine("  Final value = " + value + "");
        Console.WriteLine("  " + ctime_elapsed.TotalSeconds + " seconds.");
    }

    private static void test06(int sample_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 times R4_UNI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        DateTime ctime;
        int sample;
        int seed;
        float value = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  Measure the time it takes R4_UNI to generate");
        Console.WriteLine("  " + sample_num + " uniformly random floats.");

        seed = 123456789;

        ctime = DateTime.Now;

        for (sample = 0; sample < sample_num; sample++)
        {
            value = r4.r4_uni(ref seed);
        }

        TimeSpan ctime_elapsed = DateTime.Now - ctime;

        Console.WriteLine("");
        Console.WriteLine("  Final value = " + value + "");
        Console.WriteLine("  " + ctime_elapsed.TotalSeconds + " seconds.");
    }

    private static void test07(int sample_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 times R8_NOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        DateTime ctime;
        float[] fn = new float[128];
        int[] kn = new int[128];
        int sample;
        int seed;
        float value = 0;
        float[] wn = new float[129];

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  Measure the time it takes R8_NOR to generate");
        Console.WriteLine("  " + sample_num + " normal random floats.");

        r4.r4_nor_setup(ref kn, ref fn, ref wn);

        seed = 123456789;

        ctime = DateTime.Now;

        for (sample = 0; sample < sample_num; sample++)
        {
            value = r4.r4_nor(ref seed, kn, fn, wn);
        }

        TimeSpan ctime_elapsed = DateTime.Now - ctime;

        Console.WriteLine("");
        Console.WriteLine("  Final value = " + value + "");
        Console.WriteLine("  " + ctime_elapsed.TotalSeconds + " seconds.");

    }

    private static void test08(int sample_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 times R4_EXP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        DateTime ctime;
        float[] fe = new float[256];
        int[] ke = new int[256];
        int sample;
        int seed;
        float value = 0;
        float[] we = new float[256];

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  Measure the time it takes R4_EXP to generate");
        Console.WriteLine("  " + sample_num + " exponential random floats.");

        r4.r4_exp_setup(ref ke, ref fe, ref we);

        seed = 123456789;

        ctime = DateTime.Now;

        for (sample = 0; sample < sample_num; sample++)
        {
            value = r4.r4_exp(ref seed, ke, fe, we);
        }

        TimeSpan ctime_elapsed = DateTime.Now - ctime;

        Console.WriteLine("");
        Console.WriteLine("  Final value = " + value + "");
        Console.WriteLine("  " + ctime_elapsed.TotalSeconds + " seconds.");

    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests CONG_SEEDED.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        int jcong_new;
        int jcong_in;
        int jcong_old;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  CONG_SEEDED is a generator of pseudorandom uniformly");
        Console.WriteLine("  distributed unsigned 32 bit integers.");
        Console.WriteLine("");
        Console.WriteLine("    Input Seed   Output Seed  Output Value");
        Console.WriteLine("");

        jcong_new = 234567891;

        for (j = 1; j <= 10; j++)
        {
            jcong_old = jcong_new;
            jcong_in = jcong_new;
            jcong_new = Congruential.cong_seeded(ref jcong_in);
            Console.WriteLine("  " + jcong_old.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + jcong_in.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + jcong_new.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests KISS_SEEDED.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        int jcong_in;
        int jcong_old;
        int jsr_in;
        int jsr_old;
        int w_in;
        int w_old;
        int value;
        int z_in;
        int z_old;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  KISS_SEEDED is a generator of pseudorandom uniformly");
        Console.WriteLine("  distributed unsigned 32 bit integers.");
        Console.WriteLine("");
        Console.WriteLine("              JCONG           JSR             W             Z         Value");
        Console.WriteLine("");

        jcong_in = 234567891;
        jsr_in = 123456789;
        w_in = 345678912;
        z_in = 456789123;

        for (j = 1; j <= 10; j++)
        {
            jcong_old = jcong_in;
            jsr_old = jsr_in;
            w_old = w_in;
            z_old = z_in;
            value = KISS.kiss_seeded(ref jcong_in, ref jsr_in, ref w_in, ref z_in);
            Console.WriteLine("  In "
                              + "  " + jcong_old.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                              + "  " + jsr_old.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                              + "  " + w_old.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                              + "  " + z_old.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            Console.WriteLine("  Out"
                              + "  " + jcong_in.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                              + "  " + jsr_in.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                              + "  " + w_in.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                              + "  " + z_in.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                              + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests MWC_SEEDED.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        int w_in;
        int w_old;
        int value;
        int z_in;
        int z_old;

        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  MWC_SEEDED is a generator of pseudorandom uniformly");
        Console.WriteLine("  distributed unsigned 32 bit integers.");
        Console.WriteLine("");
        Console.WriteLine("       Input W       Input Z      Output W      Output Z  Output Value");
        Console.WriteLine("");

        w_in = 345678912;
        z_in = 456789123;

        for (j = 1; j <= 10; j++)
        {
            w_old = w_in;
            z_old = z_in;
            value = MultiplyWithCarry.mwc_seeded(ref w_in, ref z_in);
            Console.WriteLine("  " + w_old.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + z_old.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + w_in.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + z_in.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void shr3_seeded_test3()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SHR3_SEEDED_TEST3 tests SHR3_SEEDED.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        int jsr_new;
        int jsr_in;
        int jsr_old;

        Console.WriteLine("");
        Console.WriteLine("SHR3_SEEDED_TEST3");
        Console.WriteLine("  SHR3_SEEDED is a generator of pseudorandom uniformly");
        Console.WriteLine("  distributed unsigned 32 bit integers.");
        Console.WriteLine("");
        Console.WriteLine("    Input Seed   Output Seed  Output Value");
        Console.WriteLine("");

        jsr_new = 123456789;

        for (j = 1; j <= 10; j++)
        {
            jsr_old = jsr_new;
            jsr_in = jsr_new;
            jsr_new = SHR3.shr3_seeded(ref jsr_in);
            Console.WriteLine("  " + jsr_old.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + jsr_in.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + jsr_new.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }
}