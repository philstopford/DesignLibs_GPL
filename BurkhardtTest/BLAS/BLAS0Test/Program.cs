using System;
using BLASTestData;
using Burkardt.Types;
using Burkardt.Uniform;

namespace BLAS0Test;

internal class Program
{
    private static void Main(string[] args)
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BLAS0_TEST.
//
//  Discussion:
//
//    BLAS0_TEST tests the BLAS0 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 March 2014
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("BLAS0_TEST");
        Console.WriteLine("  Test the BLAS0 library.");

        dmach_test();
        test01();
        test015();
        test02();
        test03();

        Console.WriteLine("");
        Console.WriteLine("BLAS0_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void dmach_test()
//****************************************************************************80
//
//  Purpose:
//
//    DMACH_TEST demonstrates DMACH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
    {
        int job;

        Console.WriteLine("");
        Console.WriteLine("DMACH_TEST");
        Console.WriteLine("  DMACH returns some approximate machine numbers.");
        Console.WriteLine("");
        job = 1;
        Console.WriteLine("  DMACH(1) = EPS =  " + BLASData.dmach(job) + "");
        job = 2;
        Console.WriteLine("  DMACH(2) = TINY = " + BLASData.dmach(job) + "");
        job = 3;
        Console.WriteLine("  DMACH(3) = HUGE = " + BLASData.dmach(job) + "");
    }

    private static void test01()

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests R4_ABS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 March 2014
//
//  Author:
//
//    John Burkardt
//
    {
        float r4_hi = 5.0f;
        float r4_lo = -3.0f;
        int test_num = 10;

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  R4_ABS returns the absolute value of an R4.");
        Console.WriteLine("");

        for (int test = 1; test <= test_num; test++)
        {
            float r4 = UniformRNG.r4_uniform_ab(r4_lo, r4_hi, ref seed);
            float r4_absolute = Math.Abs(r4);
            Console.WriteLine("  " + r4.ToString().PadLeft(10)
                                   + "  " + r4_absolute.ToString().PadLeft(10) + "");
        }
    }

    private static void test015()

//****************************************************************************80
//
//  Purpose:
//
//    TEST015 tests R4_SIGN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 March 2014
//
//  Author:
//
//    John Burkardt
//
    {
        int TEST_NUM = 5;

        float[] x_test =  {
                -1.25f, -0.25f, 0.0f, +0.5f, +9.0f
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST015");
        Console.WriteLine("  R4_SIGN returns the sign of a number.");
        Console.WriteLine("");

        for (int test = 0; test < TEST_NUM; test++)
        {
            float x = x_test[test];
            Console.WriteLine("  " + x.ToString().PadLeft(8)
                                   + "  " + typeMethods.r4_sign(x).ToString().PadLeft(8) + "");
        }
    }

    private static void test02()
//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests R8_ABS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 March 2014
//
//  Author:
//
//    John Burkardt
//
    {
        double r8_hi = 5.0;
        double r8_lo = -3.0;
        int test_num = 10;

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  R8_ABS returns the absolute value of an R8.");
        Console.WriteLine("");
        Console.WriteLine("      X         R8_ABS(X)");
        Console.WriteLine("");

        for (int test = 1; test <= test_num; test++)
        {
            double r8 = UniformRNG.r8_uniform_ab(r8_lo, r8_hi, ref seed);
            double r8_absolute = Math.Abs(r8);
            Console.WriteLine("  " + r8.ToString().PadLeft(10)
                                   + "  " + r8_absolute.ToString().PadLeft(10) + "");
        }
    }

    private static void test03()

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests R8_SIGN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 March 2014
//
//  Author:
//
//    John Burkardt
//
    {
        int TEST_NUM = 5;

        double[] x_test =  {
                -1.25, -0.25, 0.0, +0.5, +9.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  R8_SIGN returns the sign of a number.");
        Console.WriteLine("");

        for (int test = 0; test < TEST_NUM; test++)
        {
            double x = x_test[test];
            Console.WriteLine("  " + x.ToString().PadLeft(8)
                                   + "  " + typeMethods.r8_sign(x).ToString().PadLeft(8) + "");
        }
    }
}