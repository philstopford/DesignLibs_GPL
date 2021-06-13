using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.PDFLib
{
    public static class Discrete2D
    {
        public static double[] discrete_cdf_to_xy(int n1, int n2, double[] cdf, int n, double[] u,
        ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISCRETE_CDF_TO_XY finds XY points corresponding to discrete CDF values.
        //
        //  Discussion:
        //
        //    This program is given a discrete CDF function and a set of N random
        //    values U.  Each value of U corresponds to a particular (I,J) subregion
        //    whose CDF value just exceeds the value of U.  Inside that subregion,
        //    we pick a point at random - this is equivalent to assuming the PDF
        //    is constant over the subregion.
        //
        //    This function is part of an example program, for which various
        //    assumptions have been made.  In particular, the region is the unit
        //    square, and the subregions are formed by a 20 by 20 grid of equal
        //    subsquares.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the number of rows and columns of data.
        //
        //    Input, double CDF[N1*N2], the CDF values associated with each 
        //    subcell.  A particular ordering has been given to the subcells so that the
        //    CDF is a monotonoe function when the subcells are listed in that order.
        //
        //    Input, int N, the number of sample points.
        //
        //    Input, double U[N], N random values.
        //
        //    Input/output, int *SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double DISCRETE_CDF_TO_XY[2*N], the sample points.
        //
        {
            double high;
            int i;
            int j;
            int k;
            double low;
            double[] r;
            double[] xy;

            xy = new double[2 * n];

            low = 0.0;
            for (j = 0; j < n2; j++)
            {
                for (i = 0; i < n1; i++)
                {
                    high = cdf[i + j * n1];
                    for (k = 0; k < n; k++)
                    {
                        if (low <= u[k] && u[k] <= high)
                        {
                            r = UniformRNG.r8vec_uniform_01_new(2, ref seed);
                            xy[0 + k * 2] = ((double) (i) + r[0]) / (double) n1;
                            xy[1 + k * 2] = ((double) (j) + r[1]) / (double) n2;
                        }
                    }

                    low = high;
                }
            }

            return xy;
        }

        public static double[] get_discrete_pdf_data1(int n1, int n2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GET_DISCRETE_PDF_DATA1 returns the value of the discrete PDF function in each cell.
            //
            //  Discussion:
            //
            //    Cell (I,J) extends from 
            //
            //      (I-1) * H < Y < I * H
            //      (J-1) * H < X < J * H
            //
            //    We have data for each cell, representing the integral of some PDF
            //    over that cell.  The function pdf(x,y) must be nonnegative.  However,
            //    we don't impose any other conditions on it.
            //
            //    The array PDF(:,:) contains the integral of pdf(x,y) over each cell,
            //    or, almost as good, simply a sample or average value.
            //
            //    We load the array PDF, and then we normalize it so that the sum of
            //    all the entries is 1.  
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
            //    Input, int N1, N2, the number of rows and columns of data.
            //
            //    Output, double GET_DISCRETE_PDF[N1*N2].  PDF(I,J) is the discrete PDF 
            //    for the cell (I,J), normalized so that the sum over all cells is 1.
            //
        {
            double[] pdf;
            double[] pdf_save =  {
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, 0.0000, 0.0001, 0.0001, 0.0002,
                0.0002, 0.0002, 0.0003, 0.0003, 0.0003,
                0.0003, 0.0003, 0.0002, 0.0002, 0.0002,
                0.0002, 0.0001, 0.0001, 0.0001, 0.0000,
                0.0000, 0.0001, 0.0002, 0.0003, 0.0004,
                0.0004, 0.0005, 0.0006, 0.0006, 0.0006,
                0.0006, 0.0006, 0.0005, 0.0005, 0.0004,
                0.0003, 0.0003, 0.0002, 0.0001, 0.0000,
                0.0000, 0.0002, 0.0003, 0.0005, 0.0006,
                0.0008, 0.0009, 0.0009, 0.0010, 0.0010,
                0.0010, 0.0009, 0.0008, 0.0008, 0.0007,
                0.0006, 0.0004, 0.0003, 0.0002, 0.0000,
                0.0000, 0.0003, 0.0005, 0.0008, 0.0010,
                0.0012, 0.0014, 0.0015, 0.0015, 0.0015,
                0.0015, 0.0014, 0.0013, 0.0011, 0.0010,
                0.0008, 0.0006, 0.0005, 0.0003, 0.0000,
                0.0000, 0.0004, 0.0009, 0.0013, 0.0016,
                0.0019, 0.0021, 0.0023, 0.0023, 0.0023,
                0.0021, 0.0020, 0.0018, 0.0016, 0.0013,
                0.0011, 0.0009, 0.0007, 0.0004, 0.0000,
                0.0000, 0.0007, 0.0014, 0.0020, 0.0025,
                0.0030, 0.0033, 0.0034, 0.0034, 0.0033,
                0.0031, 0.0028, 0.0025, 0.0022, 0.0018,
                0.0015, 0.0012, 0.0009, 0.0006, 0.0000,
                0.0000, 0.0011, 0.0021, 0.0031, 0.0039,
                0.0045, 0.0049, 0.0051, 0.0050, 0.0047,
                0.0043, 0.0039, 0.0034, 0.0029, 0.0024,
                0.0019, 0.0015, 0.0011, 0.0007, 0.0000,
                0.0000, 0.0017, 0.0033, 0.0048, 0.0060,
                0.0069, 0.0074, 0.0074, 0.0072, 0.0066,
                0.0059, 0.0052, 0.0045, 0.0037, 0.0031,
                0.0025, 0.0019, 0.0014, 0.0009, 0.0000,
                0.0000, 0.0025, 0.0050, 0.0073, 0.0091,
                0.0104, 0.0109, 0.0107, 0.0101, 0.0091,
                0.0080, 0.0068, 0.0057, 0.0047, 0.0038,
                0.0030, 0.0023, 0.0017, 0.0011, 0.0000,
                0.0000, 0.0038, 0.0075, 0.0110, 0.0136,
                0.0153, 0.0157, 0.0151, 0.0138, 0.0121,
                0.0104, 0.0087, 0.0071, 0.0058, 0.0046,
                0.0036, 0.0027, 0.0019, 0.0012, 0.0000,
                0.0000, 0.0055, 0.0110, 0.0160, 0.0198,
                0.0218, 0.0219, 0.0205, 0.0182, 0.0155,
                0.0129, 0.0106, 0.0085, 0.0068, 0.0053,
                0.0041, 0.0031, 0.0022, 0.0014, 0.0000,
                0.0000, 0.0077, 0.0154, 0.0224, 0.0276,
                0.0299, 0.0293, 0.0266, 0.0229, 0.0190,
                0.0154, 0.0123, 0.0098, 0.0077, 0.0059,
                0.0045, 0.0034, 0.0024, 0.0015, 0.0000,
                0.0000, 0.0100, 0.0202, 0.0295, 0.0362,
                0.0385, 0.0368, 0.0324, 0.0271, 0.0219,
                0.0174, 0.0137, 0.0107, 0.0082, 0.0063,
                0.0048, 0.0035, 0.0025, 0.0016, 0.0000,
                0.0000, 0.0120, 0.0244, 0.0356, 0.0432,
                0.0455, 0.0426, 0.0366, 0.0298, 0.0236,
                0.0184, 0.0143, 0.0110, 0.0084, 0.0064,
                0.0048, 0.0035, 0.0025, 0.0016, 0.0000,
                0.0000, 0.0134, 0.0266, 0.0382, 0.0461,
                0.0480, 0.0445, 0.0376, 0.0301, 0.0235,
                0.0181, 0.0139, 0.0106, 0.0081, 0.0061,
                0.0046, 0.0033, 0.0023, 0.0015, 0.0000,
                0.0000, 0.0151, 0.0261, 0.0362, 0.0436,
                0.0447, 0.0412, 0.0347, 0.0276, 0.0214,
                0.0164, 0.0125, 0.0095, 0.0072, 0.0054,
                0.0041, 0.0029, 0.0021, 0.0013, 0.0000,
                0.0000, 0.0174, 0.0220, 0.0295, 0.0349,
                0.0361, 0.0333, 0.0281, 0.0225, 0.0175,
                0.0134, 0.0102, 0.0078, 0.0059, 0.0044,
                0.0033, 0.0024, 0.0017, 0.0010, 0.0000,
                0.0000, 0.0097, 0.0152, 0.0200, 0.0235,
                0.0244, 0.0227, 0.0193, 0.0156, 0.0122,
                0.0094, 0.0072, 0.0055, 0.0041, 0.0031,
                0.0023, 0.0017, 0.0012, 0.0007, 0.0000,
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000
            }
            ;
            double scale;
            double total;

            pdf = typeMethods.r8mat_copy_new(n1, n2, pdf_save);
            //
            //  Normalize to get an integral of 1.
            //
            total = typeMethods.r8mat_sum(n1, n2, pdf);

            Console.WriteLine("");
            Console.WriteLine("  PDF data sums to " + total + "");

            scale = 1.0 / total;

            typeMethods.r8mat_scale(n1, n2, scale, ref pdf);

            return pdf;
        }

        public static double[] get_discrete_pdf_data2(int n1, int n2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GET_DISCRETE_PDF_DATA2 returns the discrete PDF data array.
            //
            //  Discussion:
            //
            //    Cell (I,J) extends from 
            //
            //      (I-1) * H < Y < I * H
            //      (J-1) * H < X < J * H
            //
            //    We have data for each cell, representing the integral of some PDF
            //    over that cell.  The function pdf(x,y) must be nonnegative.  However,
            //    we don't impose any other conditions on it.
            //
            //    The array PDF(:,:) contains the integral of pdf(x,y) over each cell,
            //    or, almost as good, simply a sample or average value.
            //
            //    We load the array PDF, and then we normalize it so that the sum of
            //    all the entries is 1.  
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
            //    Input, integer ( kind = 4 ) N1, N2, the number of rows and columns
            //    of PDF data.
            //
            //    Output, real ( kind = 8 ) PDF(N1,N2).  PDF(I,J) is the discrete PDF 
            //    for the cell (I,J), normalized so that the sum over all cells is 1.
            //
        {
            double[] pdf;
            double[] pdf_save = {
                10.0, 20.0, 10.0, 10.0, 20.0, 10.0, 30.0, 10.0, 10.0, 10.0, 10.0, 50.0,
                25.0, 30.0, 10.0, 25.0, 30.0, 40.0, 30.0, 20.0, 10.0, 20.0, 30.0, 40.0,
                25.0, 30.0, 20.0, 10.0, 40.0, 200.0, 50.0, 40.0, 10.0, 30.0, 60.0, 40.0,
                20.0, 30.0, 40.0, 10.0, 75.0, 100.0, 100.0, 30.0, 25.0, 25.0, 90.0, 30.0,
                75.0, 15.0, 20.0, 10.0, 75.0, 50.0, 40.0, 10.0, 100.0, 25.0, 25.0, 80.0,
                25.0, 50.0, 50.0, 10.0, 25.0, 25.0, 15.0, 10.0, 25.0, 25.0, 10.0, 10.0,
                100.0, 50.0, 50.0, 10.0, 10.0, 10.0, 5.0, 100.0, 50.0, 50.0, 10.0, 10.0,
                10.0, 10.0, 25.0, 50.0, 10.0, 50.0, 10.0, 50.0, 25.0, 25.0, 25.0, 10.0
            }
            ;
            double scale;
            double total;

            pdf = typeMethods.r8mat_copy_new(n1, n2, pdf_save);
            //
            //  Normalize to get an integral of 1.
            //
            total = typeMethods.r8mat_sum(n1, n2, pdf);

            Console.WriteLine("");
            Console.WriteLine("  PDF data sums to " + total + "");

            scale = 1.0 / total;

            typeMethods.r8mat_scale(n1, n2, scale, ref pdf);

            return pdf;
        }

        public static void get_discrete_pdf_size1(ref int n1, ref int n2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GET_DISCRETE_PDF_SIZE1 returns the dimension of the discrete PDF data
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 December 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int &N1, &N2, the number of rows and columns
        //    of data.
        //
        {
            n1 = 20;
            n2 = 20;
        }

        public static void get_discrete_pdf_size2(ref int n1, ref int n2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GET_DISCRETE_PDF_SIZE2 returns the dimension of the discrete PDF data
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 December 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int &N1, &N2, the number of rows and columns
        //    of data.
        //
        {
            n1 = 12;
            n2 = 8;
        }
        
        public static double[] set_discrete_cdf ( int n1, int n2, double[] pdf )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SET_DISCRETE_CDF sets a CDF from a discrete PDF.
            //
            //  Discussion:
            //
            //    Here, we proceed from cell (1,1) to (2,1) to (20,1), (1,2), (2,2)...(20,20).
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
            //    Input, int N1, N2, the number of rows and columns of data.
            //
            //    Input, double PDF[N1,N2], the discrete PDF for the cell (I,J),
            //    normalized so that the sum over all cells is 1.
            //
            //    Output, double CDF[N1,N2], the discrete CDF for the cell (I,J).
            //    The last entry of CDF should be 1.
            //
        {
            double[] cdf;
            int i;
            int j;
            double total;

            cdf = new double[n1*n2];

            total = 0.0;
            for ( j = 0; j < n2; j++ )
            {
                for ( i = 0; i < n1; i++ )
                {
                    total = total + pdf[i+j*n1];
                    cdf[i+j*n1] = total;
                }
            }
            return cdf;
        }
        
    }
}