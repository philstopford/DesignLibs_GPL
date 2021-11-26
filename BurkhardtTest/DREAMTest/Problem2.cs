﻿using System;
using Burkardt.DREAM;
using Burkardt.PDFLib;

namespace DREAMTest;

public static class Problem2
{
    public static void test()
    {
        ProblemSize ps = new()
        {
            chain_num = 10,
            cr_num = 3,
            gen_num = 10,
            pair_num = 3,
            par_num = 10
        };

        ProblemValue pv = new()
        {
            chain_filename = "problem2_chain00.txt",
            gr_filename = "problem2_gr.txt",
            gr_threshold = 1.2,
            jumpstep = 5,
            printstep = 10,
            restart_read_filename = "",
            restart_write_filename = "problem2_restart.txt",
            limits = new double[ps.par_num * 2]
        };

        for (int j = 0; j < ps.par_num; j++)
        {
            pv.limits[0 + j * 2] = -10.0;
            pv.limits[1 + j * 2] = +10.0;
        }

        Dream.dream(ref ps, ref pv, prior_sample, prior_density, sample_likelihood);

    }

    private static Dream.DensityResult prior_density(int par_num, double[] zp, int zpIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRIOR_DENSITY evaluates the prior density function.
        //
        //  Discussion:
        //
        //    "The initial sample was generated from a normal distribution with
        //    variance-covariance matrix 5 * Id."
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 July 2013
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
        int i;
        Dream.DensityData data = new() { zp_mean_1d = 0.0, zp_sd_1d = Math.Sqrt(5.0) };

        double value = 1.0;
        for (i = 0; i < par_num; i++)
        {
            value *= PDF.r8_normal_pdf(data.zp_mean_1d, data.zp_sd_1d, zp[i]);
        }

        return new Dream.DensityResult { result = value, data = data };
    }

    private static Dream.SampleResult prior_sample(int par_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRIOR_SAMPLE samples from the prior distribution.
        //
        //  Discussion:
        //
        //    "The initial sample was generated from a normal distribution with
        //    variance-covariance matrix 5 * Id."
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2018
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
        int i;

        Dream.SampleData data = new()
        {
            zp_mean_1d = 0.0,
            zp_sd_1d = Math.Sqrt(5.0)
        };
            
        double[] zp = new double[par_num];

        for (i = 0; i < par_num; i++)
        {
            zp[i] = PDF.r8_normal_sample(data.zp_mean_1d, data.zp_sd_1d);
        }

        return new Dream.SampleResult {result = zp, data = data};
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
        //    "A 10-dimensional twisted Gaussian density function, given by the
        //    unnormalized density
        //      pi_b(x) proportional to pi(phi_b(x)),
        //    with
        //      phi_b(x) = (x1, x2+bx1^2-100b,x3,...,x10).
        //    Here, pi signifies the density of a multivariate normal distribution,
        //    Nd(0,Sigma), with Sigma=diag(100,1,...,1), and phi_b is a function that
        //    is used to transform pi to a twisted distribution."
        //
        //    The value b = 0.01 corresponds to a mildly nonlinear problem,
        //    and b = 0.1 to a highly nonlinear problem.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 July 2013
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
        const double b = 0.01;
        int i;

        double y2 = zp[1] + b * zp[0] * zp[0] - 100.0 * b;

        double xcx = zp[0] * zp[0] / 100.0 + y2 * y2;
        for (i = 2; i < par_num; i++)
        {
            xcx += zp[i] * zp[i];
        }

        double value = -0.5 * par_num * Math.Log(2.0 * Math.PI)
                       - 0.5 * Math.Log(100.0)
                       - 0.5 * xcx;

        return value;
    }
}