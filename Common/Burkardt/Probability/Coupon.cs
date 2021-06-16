using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Coupon
    {
        public static double coupon_complete_pdf(int type_num, int box_num)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COUPON_COMPLETE_PDF evaluates the Complete Coupon Collection PDF.
        //
        //  Discussion:
        //
        //    PDF(TYPE_NUM;BOX_NUM) is the probability that, given an inexhaustible
        //    supply of boxes, inside each of which there is one of TYPE_NUM distinct
        //    coupons, which are uniformly distributed among the boxes, that it will
        //    require opening exactly BOX_NUM boxes to achieve at least one of each
        //    kind of coupon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Herbert Wilf,
        //    Some New Aspects of the Coupon Collector's Problem,
        //    SIAM Review,
        //    Volume 48, Number 3, September 2006, pages 549-565.
        //
        //  Parameters:
        //
        //    Input, int BOX_NUM, the number of boxes that had to be opened
        //    in order to just get at least one of each coupon.
        //    0 <= BOX_NUM.  If BOX_NUM < TYPE_NUM, then PDF is surely 0.
        //
        //    Input, int TYPE_NUM, the number of distinct coupons.
        //    1 <= TYPE_NUM.
        //
        //    Output, double COUPON_COMPLETE_PDF, the value of the PDF.
        //
        {
            double pdf;
            //
            //  Nonsense cases.
            //
            if (box_num < 0)
            {
                pdf = 0.0;
            }
            else if (type_num < 1)
            {
                pdf = 0.0;
            }
            //
            //  Degenerate but meaningful case.
            //
            else if (type_num == 1)
            {
                if (box_num == 1)
                {
                    pdf = 1.0;
                }
                else
                {
                    pdf = 0.0;
                }
            }
            //
            //  Easy cases.
            //
            else if (box_num < type_num)
            {
                pdf = 0.0;
            }
            //
            //  General case.
            //
            else
            {
                double factor = 1.0;
                for (int i = 1; i <= type_num; i++)
                {
                    factor = factor * (double) (i) / (double) (type_num);
                }

                for (int i = type_num + 1; i <= box_num; i++)
                {
                    factor = factor / (double) (type_num);
                }

                pdf = factor * (double) (Misc.stirling2_value(box_num - 1, type_num - 1));
            }

            return pdf;
        }

        public static double coupon_mean(int j, int type_num)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COUPON_MEAN returns the mean of the Coupon PDF.
        //
        //  Discussion:
        //
        //    In this version of the coupon collector's problem, we assume
        //    that each box contains 1 coupon, that there are TYPE_NUM distinct types
        //    of coupon, uniformly distributed among an inexhaustible supply
        //    of boxes, and that the collector's goal is to get J distinct
        //    types of coupons by opening one box after another.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int J, the number of distinct coupons to be collected.
        //    J must be between 1 and N.
        //
        //    Input, int TYPE_NUM, the number of distinct coupons.
        //
        //    Output, double COUPON_MEAN, the mean number of boxes that
        //    must be opened in order to just get J distinct kinds.
        //
        {
            if (type_num < j)
            {
                Console.WriteLine(" ");
                Console.WriteLine("COUPON_MEAN - Fatal error!");
                Console.WriteLine("  Number of distinct coupons desired must be no more");
                Console.WriteLine("  than the total number of boxes opened.");
                return (1);
            }

            double mean = 0.0;

            for (int i = 1; i <= j; i++)
            {
                mean = mean + 1.0 / (double) (type_num - i + 1);
            }

            mean = mean * (double) (type_num);

            return mean;
        }

        public static void coupon_sample(int type_num, ref int seed, ref int[] coupon, ref int box_num )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COUPON_SAMPLE simulates the coupon collector's problem.
        //
        //  Discussion:
        //
        //    The coupon collector needs to collect one of each of TYPE_NUM
        //    coupons.  The collector may draw one coupon (or, in some settings,
        //    open one box) on each trial, and takes as many trials as necessary
        //    to complete the task.  On each trial, the probability of picking
        //    any particular type of coupon is always 1 / TYPE_NUM.
        //
        //    Interesting questions include;
        //
        //    * what is the expected number of drawings necessary to complete
        //      the collection?
        //
        //    * How does the expected number of drawings necessary to complete
        //      the collection vary as TYPE_NUM increases?
        //
        //    * What is the distribution of the numbers of each type of coupon
        //      in a typical collection when it is just completed?
        //
        //    As TYPE_NUM increases, the number of coupons necessary to be
        //    collected in order to get a complete set in any simulation
        //    strongly tends to the value TYPE_NUM * LOG ( TYPE_NUM ).
        //
        //    If TYPE_NUM is 1, the simulation ends with a single drawing.
        //
        //    If TYPE_NUM is 2, then we may call the coupon taken on the first drawing
        //    a "Head", say, and the process then is similar to the question of the
        //    length, plus one, of a run of Heads or Tails in coin flipping.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int TYPE_NUM, the number of types of coupons.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int COUPON[TYPE_NUM], the number of coupons of each type
        //    that were collected during the simulation.
        //
        //    Output, int *BOX_NUM, the total number of boxes opened.
        //
        {
            int i;
            int box_max = 2000;

            for (i = 0; i < type_num; i++)
            {
                coupon[i] = 0;
            }

            int straight = 0;
            box_num = 0;
            //
            //  Draw another coupon.
            //
            while (box_num < box_max)
            {
                i = UniformRNG.i4_uniform_ab(1, type_num, ref seed);
                //
                //  Increment the number of I coupons.
                //
                coupon[i - 1] = coupon[i - 1] + 1;
                box_num = box_num + 1;
                //
                //  If I is the next one we needed, increase STRAIGHT by 1.
                //
                if (i == straight + 1)
                {
                    for (;;)
                    {
                        straight = straight + 1;
                        //
                        //  If STRAIGHT = TYPE_NUM, we have all of them.
                        //
                        if (type_num <= straight)
                        {
                            return;
                        }

                        //
                        //  If the next coupon has not been collected, our straight is over.
                        //
                        if (coupon[straight] <= 0)
                        {
                            break;
                        }
                    }
                }
            }

            Console.WriteLine(" ");
            Console.WriteLine("COUPON_SAMPLE - Fatal error!");
            Console.WriteLine("  Maximum number of coupons drawn without success.");
        }

        public static double coupon_variance(int j, int type_num)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COUPON_VARIANCE returns the variance of the Coupon PDF.
        //
        //  Discussion:
        //
        //    In this version of the coupon collector's problem, we assume
        //    that each box contains 1 coupon, that there are TYPE_NUM distinct types
        //    of coupon, uniformly distributed among an inexhaustible supply
        //    of boxes, and that the collector's goal is to get J distinct
        //    types of coupons by opening one box after another.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int J, the number of distinct coupons to be collected.
        //    J must be between 1 and N.
        //
        //    Input, int TYPE_NUM, the number of distinct coupons.
        //
        //    Output, double COUPON_VARIANCE, the variance of the number of
        //    boxes that must be opened in order to just get J distinct kinds.
        //
        {
            if (type_num < j)
            {
                Console.WriteLine(" ");
                Console.WriteLine("COUPON_VARIANCE - Fatal error!");
                Console.WriteLine("  Number of distinct coupons desired must be no more");
                Console.WriteLine("  than the total number of distinct coupons.");
                return (1);
            }

            double variance = 0.0;
            for (int i = 1; i <= j; i++)
            {
                variance = variance + (double) (i - 1) /
                    Math.Pow((double) (type_num - i + 1), 2);
            }

            variance = variance * (double) (type_num);

            return variance;
        }

    }
}