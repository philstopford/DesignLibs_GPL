using System;
using Burkardt.DREAM;
using Burkardt.PDFLib;
using Burkardt.Types;

namespace DREAMTest;

public static partial class Problem1c
{
    private static Covariance c;

    public static void test()
    {
        ProblemSize ps = new()
        {
            chain_num = 10,
            cr_num = 3,
            gen_num = 10,
            pair_num = 3,
            //par_num = 100,
            par_num = 5,
        };

        c = new Covariance(ps.par_num);
            
        ProblemValue pv = new()
        {
            chain_filename = "problem1_chain00.txt",
            gr_filename = "problem1_gr.txt",
            gr_threshold = 1.2,
            jumpstep = 5,
            printstep = 10,
            restart_read_filename = "",
            restart_write_filename = "problem1_restart.txt",
            limits = new double[ps.par_num * 2]
        };

        for (int j = 0; j < ps.par_num; j++)
        {
            pv.limits[0 + j * 2] = 9.9;
            pv.limits[1 + j * 2] = +10.0;
        }

        Dream.dream(ref ps, ref pv, prior_sample, prior_density, sample_likelihood );

    }

    private static Dream.DensityResult prior_density(int par_num, double[] zp, int zpIndex = 0)

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
        //    26 June 2013
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
        double det;
        double[] factor;
        double[] mean;
        double value = 0;

        mean = c.mean_get();
        factor = c.factor_get();
        det = c.det_get();

        value = PDF.r8vec_multinormal_pdf(par_num, mean, factor, det, zp, zpIndex);

        return new Dream.DensityResult {result = value};
    }

    private static Dream.SampleResult prior_sample(int par_num)

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
        //    07 August 2013
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
        double[] factor;
        int i;
        double[] mean;
        double[] x;
        double[] zp;

        x = new double[par_num];

        for (i = 0; i < par_num; i++)
        {
            x[i] = PDF.r8_normal_01_sample();
        }

        factor = c.factor_get();
        zp = typeMethods.r8mat_mtv_new(par_num, par_num, factor, x);

        mean = c.mean_get();
        for (i = 0; i < par_num; i++)
        {
            zp[i] += mean[i];
        }

        return new Dream.SampleResult {result = zp};
    }

    private static double sample_likelihood(int par_num, double[] zp)

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
        //    26 June 2013
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
        double det;
        double[] factor;
        int i;
        double[] mean;
        const double pi = 3.141592653589793;
        double value = 0;
        double[] x;
        double xcx;
        double[] y;

        mean = c.mean_get();
        x = new double[par_num];

        for (i = 0; i < par_num; i++)
        {
            x[i] = zp[i] - mean[i];
        }

        factor = c.factor_get();
        y = typeMethods.r8mat_utsol(par_num, factor, x);
        //
        //  Compute:
        //    (x-mu)' * inv(C)          * (x-mu)
        //  = (x-mu)' * inv(R'*R)       * (x-mu)
        //  = (x-mu)' * inv(R) * inv(R) * (x-mu)
        //  = y' * y.
        //
        xcx = typeMethods.r8vec_dot_product(par_num, y, y);

        det = c.det_get();
        value = -0.5 * par_num * Math.Log(2.0 * pi)
                - 0.5 * Math.Log(det)
                - 0.5 * xcx;

        return value;
    }
}