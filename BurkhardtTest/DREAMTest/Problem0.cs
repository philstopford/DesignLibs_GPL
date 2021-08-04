using System;
using Burkardt.DREAM;
using Burkardt.PDFLib;

namespace DREAMTest
{
    public static class Problem0
    {
        public static void test()
        {
            ProblemSize ps = new ProblemSize()
            {
                chain_num = 10,
                cr_num = 3,
                gen_num = 10,
                pair_num = 3,
                par_num = 10
            };

            ProblemValue pv = new ProblemValue()
            {
                chain_filename = "problem0_chain00.txt",
                gr_filename = "problem0_gr.txt",
                gr_threshold = 1.2,
                jumpstep = 5,
                printstep = 10,
                restart_read_filename = "",
                restart_write_filename = "problem0_restart.txt"
            };

            pv.limits = new double[ps.par_num * 2 ];
            for (int j = 0; j < ps.par_num; j++)
            {
                pv.limits[0 + j * 2] = -10.0;
                pv.limits[1 + j * 2] = +10.0;
            }
            
            Dream.dream(ref ps, ref pv, prior_sample, prior_density, sample_likelihood );

        }
        
        static Dream.DensityResult prior_density(int par_num, double[] zp, int zpIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PRIOR_DENSITY evaluates the prior density function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 May 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int PAR_NUM, the total number of parameters.
            //    1 <= PAR_NUM.
            //
            //    Input, double ZP[PAR_NUM], the argument of the density
            //    function.
            //
            //    Output, real PRIOR_DENSITY, the value of the prior density function.
            //
        {
            double[] a =  {
                -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0
            }
            ;
            double[] b =  {
                +10.0, +10.0, +10.0, +10.0, +10.0, +10.0, +10.0, +10.0, +10.0, +10.0
            }
            ;
            int i;
            double value;

            value = 1.0;

            for (i = 0; i < par_num; i++)
            {
                value = value * PDF.r8_uniform_pdf(a[i], b[i], zp[zpIndex + i]);
            }

            return new Dream.DensityResult() { result = value };
        }

        static Dream.SampleResult prior_sample(int par_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PRIOR_SAMPLE samples from the prior distribution.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 May 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int PAR_NUM, the total number of parameters.
            //    1 <= PAR_NUM.
            //
            //    Output, double PRIOR_SAMPLE[PAR_NUM], the sample from the distribution.
            //
        {
            double[] a =  {
                -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0
            }
            ;
            double[] b =  {
                +10.0, +10.0, +10.0, +10.0, +10.0, +10.0, +10.0, +10.0, +10.0, +10.0
            }
            ;
            int i;
            double[] zp;

            zp = new double[par_num];

            for (i = 0; i < par_num; i++)
            {
                zp[i] = PDF.r8_uniform_sample(a[i], b[i]);
            }

            return new Dream.SampleResult() {result = zp};
        }

        static double sample_likelihood(int par_num, double[] zp)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SAMPLE_LIKELIHOOD computes the log likelihood function.
            //
            //  Discussion:
            //
            //    This is a one mode Gaussian.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 May 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int PAR_NUM, the total number of parameters.
            //    1 <= PAR_NUM.
            //
            //    Input, double ZP[PAR_NUM], a sample.
            //
            //    Output, double SAMPLE_LIKELIHOOD, the log likelihood function 
            //    for the sample.
            //
        {
            int i;
            const double pi = 3.141592653589793;
            double t1;
            double t2;
            double value;

            t1 = 0.0;
            for (i = 0; i < par_num; i++)
            {
                t1 = t1 + Math.Pow(zp[i] + 5.0, 2);
            }

            t1 = Math.Log(1.0 / 3.0 / (2.0 * pi)) - 0.5 * t1;

            t2 = 0.0;
            for (i = 0; i < par_num; i++)
            {
                t2 = t2 + Math.Pow(zp[i] - 5.0, 2);
            }

            t2 = Math.Log(2.0 / 3.0 / (2.0 * pi)) - 0.5 * t2;

            if (t1 < t2)
            {
                value = t2;
            }
            else
            {
                value = t1;
            }

            return value;
        }        
    }
}