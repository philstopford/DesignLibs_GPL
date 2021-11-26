using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void student_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_CDF_TEST tests STUDENT_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("STUDENT_CDF_TEST");
        Console.WriteLine("  STUDENT_CDF evaluates the Student CDF;");
        Console.WriteLine("  STUDENT_PDF evaluates the Student PDF;");
        Console.WriteLine("  STUDENT_SAMPLE samples the Student PDF;");

        const double a = 0.5;
        const double b = 2.0;
        const double c = 6.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A = " + a + "");
        Console.WriteLine("  PDF parameter B = " + b + "");
        Console.WriteLine("  PDF parameter C = " + c + "");

        if (!Student.student_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("STUDENT_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameter values are illegal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Student.student_sample(a, b, c, ref seed);
            double pdf = Student.student_pdf(x, a, b, c);
            double cdf = Student.student_cdf(x, a, b, c);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void student_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_SAMPLE_TEST tests STUDENT_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 March 2016
//
//  Author:
//
//    John Burkardt
//
    {
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        double[] x = new double [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("STUDENT_SAMPLE_TEST");
        Console.WriteLine("  STUDENT_MEAN evaluates the Student mean;");
        Console.WriteLine("  STUDENT_SAMPLE samples the Student PDF;");
        Console.WriteLine("  STUDENT_VARIANCE computes the Student variance;");

        double a = 0.5;
        double b = 2.0;
        double c = 6.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A = " + a + "");
        Console.WriteLine("  PDF parameter B = " + b + "");
        Console.WriteLine("  PDF parameter C = " + c + "");

        if (!Student.student_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("STUDENT_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameter values are illegal.");
            return;
        }

        double mean = Student.student_mean(a, b, c);
        double variance = Student.student_variance(a, b, c);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Student.student_sample(a, b, c, ref seed);
        }

        mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
        variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
        double xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
        double xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");
    }

    private static void student_noncentral_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_NONCENTRAL_CDF_TEST tests STUDENT_NONCENTRAL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("STUDENT_NONCENTRAL_CDF_TEST");
        Console.WriteLine("  STUDENT_NONCENTRAL_CDF evaluates the Student Noncentral CDF;");

        const double x = 0.50;
        const int idf = 10;
        const double b = 1.0;

        double cdf = Student.student_noncentral_cdf(x, idf, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF argument X =              " + x + "");
        Console.WriteLine("  PDF parameter IDF =           " + idf + "");
        Console.WriteLine("  PDF parameter B =             " + b + "");
        Console.WriteLine("  CDF value =                   " + cdf + "");

    }

}