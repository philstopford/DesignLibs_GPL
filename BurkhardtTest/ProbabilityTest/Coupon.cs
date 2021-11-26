﻿using System;
using System.Globalization;
using Burkardt.Probability;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void coupon_complete_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    COUPON_COMPLETE_PDF_TEST tests COUPON_COMPLETE_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int type_num;

        Console.WriteLine("");
        Console.WriteLine("COUPON_COMPLETE_PDF_TEST");
        Console.WriteLine("  COUPON_COMPLETE_PDF evaluates the coupon collector's");
        Console.WriteLine("  complete collection pdf.");
        Console.WriteLine("");

        for (type_num = 2; type_num <= 4; type_num++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Number of coupon types is " + type_num + "");
            Console.WriteLine("");
            Console.WriteLine("   BOX_NUM      PDF             CDF");
            Console.WriteLine("");
            double cdf = 0.0;
            int box_num;
            for (box_num = 1; box_num <= 20; box_num++)
            {
                double pdf = Coupon.coupon_complete_pdf(type_num, box_num);
                cdf += pdf;
                Console.WriteLine("  " + box_num.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    private static void coupon_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    COUPON_SAMPLE_TEST tests COUPON_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        const int N_TRIAL = 10;
        const int MAX_TYPE = 25;

        int[] coupon = new int[MAX_TYPE];
        int n_coupon = 0;
        int n_type;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("COUPON_SAMPLE_TEST:");
        Console.WriteLine("  COUPON_SAMPLE samples the coupon PDF.");
        Console.WriteLine("");

        for (n_type = 5; n_type <= MAX_TYPE; n_type += 5)
        {
            Console.WriteLine("");
            Console.WriteLine("  Number of coupon types is " + n_type + "");
            double expect = n_type * Math.Log(n_type);
            Console.WriteLine("  Expected wait is about " + expect + "");
            Console.WriteLine("");

            double average = 0.0;
            int i;
            for (i = 1; i <= N_TRIAL; i++)
            {
                Coupon.coupon_sample(n_type, ref seed, ref coupon, ref n_coupon);

                Console.WriteLine("  "
                                  + i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                  + n_coupon.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");

                average += n_coupon;
            }

            average /= N_TRIAL;
            Console.WriteLine("");
            Console.WriteLine("  Average wait was " + average + "");
        }
    }

}