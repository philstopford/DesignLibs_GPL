﻿using System;
using System.Globalization;
using Burkardt.Ziggurat;

namespace ZigguratTest;

internal static class Program
{
    private static void Main()
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
        int j;

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
            int seed = (j % 2) switch
            {
                0 => 123456789,
                _ => 987654321
            };

            Console.WriteLine("");
            Console.WriteLine("  " + 0.ToString().PadLeft(6)
                                   + "  " + seed.ToString().PadLeft(12) + "");
            Console.WriteLine("");

            int i;
            for (i = 1; i <= 10; i++)
            {
                int value = SHR3.shr3_seeded(ref seed);

                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + seed.ToString().PadLeft(12)
                                       + "  " + value.ToString().PadLeft(12) + "");
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
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  R4_UNI returns pseudorandom uniformly distributed");
        Console.WriteLine("  floats (single precision reals) between 0 and 1.");

        for (j = 0; j < 3; j++)
        {
            int seed = (j % 2) switch
            {
                0 => 123456789,
                _ => 987654321
            };

            Console.WriteLine("");
            Console.WriteLine("  " + 0.ToString().PadLeft(6)
                                   + "  " + seed.ToString().PadLeft(12) + "");
            Console.WriteLine("");

            int i;
            for (i = 1; i <= 10; i++)
            {
                float value = r4.r4_uni(ref seed);

                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + seed.ToString().PadLeft(12)
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
        int j;
        int[] kn = new int[128];
        float[] wn = new float[128];

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  R4_NOR returns pseudorandom normally distributed");
        Console.WriteLine("  floats (single precision reals) between 0 and 1.");

        r4.r4_nor_setup(ref kn, ref fn, ref wn);

        for (j = 0; j < 3; j++)
        {
            int seed = (j % 2) switch
            {
                0 => 123456789,
                _ => 987654321
            };

            Console.WriteLine("");
            Console.WriteLine("  " + 0.ToString().PadLeft(6)
                                   + "  " + seed.ToString().PadLeft(12) + "");
            Console.WriteLine("");

            int i;
            for (i = 1; i <= 10; i++)
            {
                float value = r4.r4_nor(ref seed, kn, fn, wn);

                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + seed.ToString().PadLeft(12)
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
        int j;
        int[] ke = new int[256];
        float[] we = new float[256];

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  R4_EXP returns pseudorandom exponentially distributed");
        Console.WriteLine("  floats (single precision reals) between 0 and 1.");

        r4.r4_exp_setup(ref ke, ref fe, ref we);

        for (j = 0; j < 3; j++)
        {
            int seed = (j % 2) switch
            {
                0 => 123456789,
                _ => 987654321
            };

            Console.WriteLine("");
            Console.WriteLine("  " + 0.ToString().PadLeft(6)
                                   + "  " + seed.ToString().PadLeft(12) + "");
            Console.WriteLine("");

            int i;
            for (i = 1; i <= 10; i++)
            {
                float value = r4.r4_exp(ref seed, ke, fe, we);

                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + seed.ToString().PadLeft(12)
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
        int sample;
        int value = 0;

        Console.WriteLine("");
        Console.WriteLine("SHR3_SEEDED_TEST2");
        Console.WriteLine("  Measure the time it takes SHR3_SEEDED to generate");
        Console.WriteLine("  " + sample_num + " unsigned 32 bit integers.");

        int seed = 123456789;

        DateTime ctime = DateTime.Now;

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
        int sample;
        float value = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  Measure the time it takes R4_UNI to generate");
        Console.WriteLine("  " + sample_num + " uniformly random floats.");

        int seed = 123456789;

        DateTime ctime = DateTime.Now;

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
        float[] fn = new float[128];
        int[] kn = new int[128];
        int sample;
        float value = 0;
        float[] wn = new float[129];

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  Measure the time it takes R8_NOR to generate");
        Console.WriteLine("  " + sample_num + " normal random floats.");

        r4.r4_nor_setup(ref kn, ref fn, ref wn);

        int seed = 123456789;

        DateTime ctime = DateTime.Now;

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
        float[] fe = new float[256];
        int[] ke = new int[256];
        int sample;
        float value = 0;
        float[] we = new float[256];

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  Measure the time it takes R4_EXP to generate");
        Console.WriteLine("  " + sample_num + " exponential random floats.");

        r4.r4_exp_setup(ref ke, ref fe, ref we);

        int seed = 123456789;

        DateTime ctime = DateTime.Now;

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

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  CONG_SEEDED is a generator of pseudorandom uniformly");
        Console.WriteLine("  distributed unsigned 32 bit integers.");
        Console.WriteLine("");
        Console.WriteLine("    Input Seed   Output Seed  Output Value");
        Console.WriteLine("");

        int jcong_new = 234567891;

        for (j = 1; j <= 10; j++)
        {
            int jcong_old = jcong_new;
            int jcong_in = jcong_new;
            jcong_new = Congruential.cong_seeded(ref jcong_in);
            Console.WriteLine("  " + jcong_old.ToString().PadLeft(12)
                                   + "  " + jcong_in.ToString().PadLeft(12)
                                   + "  " + jcong_new.ToString().PadLeft(12) + "");
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

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  KISS_SEEDED is a generator of pseudorandom uniformly");
        Console.WriteLine("  distributed unsigned 32 bit integers.");
        Console.WriteLine("");
        Console.WriteLine("              JCONG           JSR             W             Z         Value");
        Console.WriteLine("");

        int jcong_in = 234567891;
        int jsr_in = 123456789;
        int w_in = 345678912;
        int z_in = 456789123;

        for (j = 1; j <= 10; j++)
        {
            int jcong_old = jcong_in;
            int jsr_old = jsr_in;
            int w_old = w_in;
            int z_old = z_in;
            int value = KISS.kiss_seeded(ref jcong_in, ref jsr_in, ref w_in, ref z_in);
            Console.WriteLine("  In "
                              + "  " + jcong_old.ToString().PadLeft(12)
                              + "  " + jsr_old.ToString().PadLeft(12)
                              + "  " + w_old.ToString().PadLeft(12)
                              + "  " + z_old.ToString().PadLeft(12) + "");
            Console.WriteLine("  Out"
                              + "  " + jcong_in.ToString().PadLeft(12)
                              + "  " + jsr_in.ToString().PadLeft(12)
                              + "  " + w_in.ToString().PadLeft(12)
                              + "  " + z_in.ToString().PadLeft(12)
                              + "  " + value.ToString().PadLeft(12) + "");
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

        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  MWC_SEEDED is a generator of pseudorandom uniformly");
        Console.WriteLine("  distributed unsigned 32 bit integers.");
        Console.WriteLine("");
        Console.WriteLine("       Input W       Input Z      Output W      Output Z  Output Value");
        Console.WriteLine("");

        int w_in = 345678912;
        int z_in = 456789123;

        for (j = 1; j <= 10; j++)
        {
            int w_old = w_in;
            int z_old = z_in;
            int value = MultiplyWithCarry.mwc_seeded(ref w_in, ref z_in);
            Console.WriteLine("  " + w_old.ToString().PadLeft(12)
                                   + "  " + z_old.ToString().PadLeft(12)
                                   + "  " + w_in.ToString().PadLeft(12)
                                   + "  " + z_in.ToString().PadLeft(12)
                                   + "  " + value.ToString().PadLeft(12) + "");
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

        Console.WriteLine("");
        Console.WriteLine("SHR3_SEEDED_TEST3");
        Console.WriteLine("  SHR3_SEEDED is a generator of pseudorandom uniformly");
        Console.WriteLine("  distributed unsigned 32 bit integers.");
        Console.WriteLine("");
        Console.WriteLine("    Input Seed   Output Seed  Output Value");
        Console.WriteLine("");

        int jsr_new = 123456789;

        for (j = 1; j <= 10; j++)
        {
            int jsr_old = jsr_new;
            int jsr_in = jsr_new;
            jsr_new = SHR3.shr3_seeded(ref jsr_in);
            Console.WriteLine("  " + jsr_old.ToString().PadLeft(12)
                                   + "  " + jsr_in.ToString().PadLeft(12)
                                   + "  " + jsr_new.ToString().PadLeft(12) + "");
        }
    }
}