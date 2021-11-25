using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class MultinomialTest
{
    public static void multinomial_pdf_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTINOMIAL_PDF_VALUES_TEST tests MULTINOMIAL_PDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m = 0;
        int n = 0;
        double pdf = 0;
        Console.WriteLine("");
        Console.WriteLine("MULTINOMIAL_PDF_VALUES_TEST:");
        Console.WriteLine("  MULTINOMIAL_PDF_VALUES stores values of the Multinomial PDF.");
        Console.WriteLine("  Given M possible outcomes on a single trial,");
        Console.WriteLine("  with each outcome having probability P,");
        Console.WriteLine("  PDF is the probability that after N trials,");
        Console.WriteLine("  outcome I occurred X(I) times.");
        Console.WriteLine("");
        Console.WriteLine("     N     M     I      P        X        PDF()");
        int n_data1 = 0;
        int n_data2 = 0;
        for (;;)
        {
            Multinomial.multinomial_pdf_sizes(ref n_data1, ref m);
            if (n_data1 == 0)
            {
                break;
            }

            double[] p = new double[m];
            int[] x = new int[m];
            Multinomial.multinomial_pdf_values(ref n_data2, m, ref n, ref p, ref x, ref pdf);
            Console.WriteLine("");
            int i;
            for (i = 0; i < m; i++)
            {
                Console.WriteLine("              " + i.ToString().PadLeft(4)
                                                   + "  " + p[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                                   + "  " + x[i].ToString().PadLeft(4) + "");
            }

            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + m.ToString().PadLeft(4)
                                   + "                        " + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }
}