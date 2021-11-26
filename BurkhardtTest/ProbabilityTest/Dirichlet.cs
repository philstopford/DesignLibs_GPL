using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void dirichlet_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_SAMPLE_TEST tests DIRICHLET_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        const int N = 3;
        const int SAMPLE_NUM = 1000;

        double[] a =  {
                0.250, 0.500, 1.250
            }
            ;
        int i;
        int j;
        int seed = 123456789;
        double[] x = new double[N * SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("DIRICHLET_SAMPLE");
        Console.WriteLine("  DIRICHLET_MEAN computes the Dirichlet mean;");
        Console.WriteLine("  DIRICHLET_SAMPLE samples the Dirichlet distribution;");
        Console.WriteLine("  DIRICHLET_VARIANCE computes the Dirichlet variance;");

        Console.WriteLine("");
        Console.WriteLine("  Number of components N = " + N + "");

        typeMethods.r8vec_print(N, a, "  PDF parameter A:");

        if (!Dirichlet.dirichlet_check(N, a))
        {
            Console.WriteLine("");
            Console.WriteLine("DIRICHLET_SAMPLE - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double[] mean = Dirichlet.dirichlet_mean(N, a);
        double[] variance = Dirichlet.dirichlet_variance(N, a);
        typeMethods.r8vec_print(N, mean, "  PDF mean:");
        typeMethods.r8vec_print(N, variance, "  PDF variance:");

        double[] m2 = Dirichlet.dirichlet_moment2(N, a);

        typeMethods.r8mat_print(N, N, m2, "  Second moment matrix:");

        for (j = 0; j < SAMPLE_NUM; j++)
        {
            double[] y = Dirichlet.dirichlet_sample(N, a, ref seed);
            for (i = 0; i < N; i++)
            {
                x[i + j * N] = y[i];
            }
        }

        mean = typeMethods.r8row_mean(N, SAMPLE_NUM, x);
        variance = typeMethods.r8row_variance(N, SAMPLE_NUM, x);
        double[] xmax = typeMethods.r8row_max(N, SAMPLE_NUM, x);
        double[] xmin = typeMethods.r8row_min(N, SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("");
        Console.WriteLine("  Component Mean, Variance, Min, Max:");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            Console.WriteLine("  "
                              + i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + mean[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + variance[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + xmax[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + xmin[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void dirichlet_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_PDF_TEST tests DIRICHLET_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        const int N = 3;

        double[] a =  {
                0.250, 0.500, 1.250
            }
            ;
        double[] x =  {
                0.500, 0.125, 0.375
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("DIRICHLET_PDF_TEST");
        Console.WriteLine("  DIRICHLET_PDF evaluates the Dirichlet PDF;");

        Console.WriteLine("");
        Console.WriteLine("  Number of components N = " + N + "");

        typeMethods.r8vec_print(N, a, "  PDF parameter A:");

        if (!Dirichlet.dirichlet_check(N, a))
        {
            Console.WriteLine("");
            Console.WriteLine("DIRICHLET_PDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        typeMethods.r8vec_print(N, x, "  PDF argument X:");

        double pdf = Dirichlet.dirichlet_pdf(x, N, a);

        Console.WriteLine("");
        Console.WriteLine("  PDF value = " + pdf + "");
    }

    private static void dirichlet_mix_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_MIX_SAMPLE_TEST tests DIRICHLET_MIX_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        const int COMP_NUM = 2;
        const int ELEM_NUM = 3;
        const int SAMPLE_NUM = 1000;

        double[] a =  {
                0.250, 0.500, 1.250,
                1.500, 0.500, 2.000
            }
            ;
        int comp = 0;
        double[] comp_weight =  {
                1.0, 2.0
            }
            ;
        int i;
        int j;
        int seed = 123456789;
        double[] x = new double [ELEM_NUM * SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("DIRICHLET_MIX_SAMPLE_TEST");
        Console.WriteLine("  DIRICHLET_MIX_SAMPLE samples the Dirichlet Mix distribution;");
        Console.WriteLine("  DIRICHLET_MIX_MEAN computes the Dirichlet Mix mean;");

        Console.WriteLine("");
        Console.WriteLine("  Number of elements ELEM_NUM =   " + ELEM_NUM + "");
        Console.WriteLine("  Number of components COMP_NUM = " + COMP_NUM + "");
        typeMethods.r8mat_print(ELEM_NUM, COMP_NUM, a, "  PDF parameters A(ELEM,COMP):");
        typeMethods.r8vec_print(COMP_NUM, comp_weight, "  Component weights");

        if (!Dirichlet.dirichlet_mix_check(COMP_NUM, ELEM_NUM, a, comp_weight))
        {
            Console.WriteLine("");
            Console.WriteLine("DIRICHLET_MIX_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double[] mean = Dirichlet.dirichlet_mix_mean(COMP_NUM, ELEM_NUM, a, comp_weight);

        typeMethods.r8vec_print(ELEM_NUM, mean, "  PDF mean:");

        for (j = 0; j < SAMPLE_NUM; j++)
        {
            double[] y = Dirichlet.dirichlet_mix_sample(COMP_NUM, ELEM_NUM, a, comp_weight, ref seed,
                ref comp);

            for (i = 0; i < ELEM_NUM; i++)
            {
                x[i + j * ELEM_NUM] = y[i];
            }
        }

        mean = typeMethods.r8row_mean(ELEM_NUM, SAMPLE_NUM, x);
        double[] variance = typeMethods.r8row_variance(ELEM_NUM, SAMPLE_NUM, x);
        double[] xmax = typeMethods.r8row_max(ELEM_NUM, SAMPLE_NUM, x);
        double[] xmin = typeMethods.r8row_min(ELEM_NUM, SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("");
        Console.WriteLine("  Component Mean, Variance, Max, Min:");
        Console.WriteLine("");

        for (i = 0; i < ELEM_NUM; i++)
        {
            Console.WriteLine("  "
                              + i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + mean[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + variance[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + xmax[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + xmin[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void dirichlet_mix_pdf_test()
//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_MIX_PDF_TEST tests DIRICHLET_MIX_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        const int COMP_NUM = 2;
        const int ELEM_NUM = 3;

        double[] a =  {
                0.250, 0.500, 1.250,
                1.500, 0.500, 2.000
            }
            ;
        double[] comp_weight =  {
                1.0, 2.0
            }
            ;
        double[] x =  {
                0.500, 0.125, 0.375
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("DIRICHLET_MIX_PDF_TEST");
        Console.WriteLine("  DIRICHLET_MIX_PDF evaluates the Dirichlet Mix PDF.");

        Console.WriteLine();

        Console.WriteLine("");
        Console.WriteLine("  Number of elements ELEM_NUM =   " + ELEM_NUM + "");
        Console.WriteLine("  Number of components COMP_NUM = " + COMP_NUM + "");
        typeMethods.r8mat_print(ELEM_NUM, COMP_NUM, a, "  PDF parameters A(ELEM,COMP):");
        typeMethods.r8vec_print(COMP_NUM, comp_weight, "  Component weights");

        Console.WriteLine();

        if (!Dirichlet.dirichlet_mix_check(COMP_NUM, ELEM_NUM, a, comp_weight))
        {
            Console.WriteLine("");
            Console.WriteLine("DIRICHLET_MIX_PDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        typeMethods.r8vec_print(ELEM_NUM, x, "  PDF argument X:");

        Console.WriteLine();

        double pdf = Dirichlet.dirichlet_mix_pdf(x, COMP_NUM, ELEM_NUM, a, comp_weight);

        Console.WriteLine("");
        Console.WriteLine("  PDF value =           " + pdf + "");

    }

}