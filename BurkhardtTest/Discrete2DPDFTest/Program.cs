using System;
using Burkardt.PDFLib;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Discrete2DPDFTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for DISCRETE_PDF_SAMPLE_2D.
            //
            //  Discussion:
            //
            //    This program is an example of how discrete sample or density data
            //    can be used to define a PDF (probability density function). 
            //
            //    In this function and the functions it calls, we assume that we have
            //    data for an array of 20 by 20 square subcells of the unit square.
            //    We wish to derive a PDF that will allow us to sample an arbitrary
            //    number of points from this region.
            //
            //    In particular, we intend to use the discrete data to generate a PDF
            //    which we will then use to generate sample points.
            //
            //    Roughly speaking, we have kept track of how many fish we caught in
            //    each part of a lake, and now we want to simulate catching N fish
            //    under the same conditions.
            //
            //    The statistics for each simulation should be governed by the discrete
            //    PDF, but with random variation.  In other words, the actual number
            //    of points taken from each subregion is random, and the actual location of
            //    each point in a subregion is random, but over many simulations, the
            //    statistics of the sample points should reproduce the statistics of
            //    the original discrete sample that defined the PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 December 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;

            Console.WriteLine("");
            Console.WriteLine("DISCRETE_PDF_SAMPLE_2D:");

            Console.WriteLine("  Generate sample data using a discrete PDF.");

            n = 1000;
            test01(n);
            test02(n);

            Console.WriteLine("");
            Console.WriteLine("DISCRETE_PDF_SAMPLE_2D:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 looks at a 20x20 region.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 December 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of sample points to be generated.
            //
        {
            double[] cdf;
            string filename;
            int n1 = 0;
            int n2 = 0;
            double[] pdf;
            int seed;
            double[] u;
            double[] xy;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  Consider data skewed toward the upper left corner of the unit square.");
            Console.WriteLine("  Generate " + n + " samples");
            //
            //  Get the dimensions of the PDF data.
            //
            Discrete2D.get_discrete_pdf_size1(ref n1, ref n2);
            Console.WriteLine("  PDF data is on a " + n1 + " by " + n2 + " grid.");
            //
            //  Construct a PDF from the data.
            //
            pdf = Discrete2D.get_discrete_pdf_data1(n1, n2);
            //
            //  "Integrate" the data over rows and columns of the region to get the CDF.
            //
            cdf = Discrete2D.set_discrete_cdf(n1, n2, pdf);
            //
            //  Choose N CDF values at random.
            //
            seed = 123456789;

            u = UniformRNG.r8vec_uniform_01_new(n, ref seed);
            //
            //  Find the cell corresponding to each CDF value,
            //  and choose a random point in that cell.
            //
            xy = Discrete2D.discrete_cdf_to_xy(n1, n2, cdf, n, u, ref seed);
            //
            //  Write data to a file for examination, plotting, or analysis.
            //
            filename = "test01.txt";
            typeMethods.r8mat_write(filename, 2, n, xy);

            Console.WriteLine("");
            Console.WriteLine("  Wrote sample data to file \"" + filename + "\".");
        }

        static void test02(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 looks at a 12x8 region.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 December 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of sample points to be generated.
            //
        {
            double[] cdf;
            string filename;
            int n1 = 0;
            int n2 = 0;
            double[] pdf;
            int seed;
            double[] u;
            double[] xy;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  Consider data suggested by the shape and density of Iowa.");
            Console.WriteLine("  Generate " + n + " samples");
            //
            //  Get the dimensions of the PDF data.
            //
            Discrete2D.get_discrete_pdf_size2(ref n1, ref n2);
            Console.WriteLine("  PDF data is on a " + n1 + " by " + n2 + " grid.");
            //
            //  Construct a PDF from the data.
            //
            pdf = Discrete2D.get_discrete_pdf_data2(n1, n2);
            //
            //  "Integrate" the data over rows and columns of the region to get the CDF.
            //
            cdf = Discrete2D.set_discrete_cdf(n1, n2, pdf);
            //
            //  Choose N CDF values at random.
            //
            seed = 123456789;

            u = UniformRNG.r8vec_uniform_01_new(n, ref seed);
            //
            //  Find the cell corresponding to each CDF value,
            //  and choose a random point in that cell.
            //
            xy = Discrete2D.discrete_cdf_to_xy(n1, n2, cdf, n, u, ref seed);
            //
            //  Write data to a file for examination, plotting, or analysis.
            //
            filename = "test02.txt";
            typeMethods.r8mat_write(filename, 2, n, xy);

            Console.WriteLine("");
            Console.WriteLine("  Wrote sample data to file \"" + filename + "\".");
        }
    }
}