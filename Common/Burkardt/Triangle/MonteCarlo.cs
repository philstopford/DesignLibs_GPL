using System;
using Burkardt.FEM;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.TriangleNS;

public static class MonteCarlo
{
    public class tusData
    {
        public int seed { get; set; }
        public double[] result { get; set; }
    }

    public static double[] triangle_monte_carlo(double[] t, int p_num, int f_num,
            Func<int, int, tusData> triangle_unit_sample,
            Func<int, double[], int, double[]> triangle_integrand, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_MONTE_CARLO applies the Monte Carlo rule to integrate a function.
        //
        //  Discussion:
        //
        //    The function f(x,y) is to be integrated over a triangle T.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Input, int P_NUM, the number of sample points.
        //
        //    Input, int F_NUM, the number of functions to integrate.
        //
        //    Input, external TRIANGLE_UNIT_SAMPLE, the sampling routine.
        //
        //    Input, external TRIANGLE_INTEGRAND, the integrand routine.
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, dobule TRIANGLE_MONTE_CARLO[F_NUM], the approximate integrals.
        //
    {
        int i;

        double area = Integrals.triangle_area(t);

        tusData data = triangle_unit_sample(p_num, seed);
        seed = data.seed;
        double[] p = data.result;

        double[] p2 = new double[2 * p_num];

        Reference.reference_to_physical_t3(t, p_num, p, ref p2);

        double[] fp = triangle_integrand(p_num, p2, f_num);

        double[] result = new double[f_num];

        for (i = 0; i < f_num; i++)
        {
            double fp_sum = 0.0;
            int j;
            for (j = 0; j < p_num; j++)
            {
                fp_sum += fp[i + j * f_num];
            }

            result[i] = area * fp_sum / p_num;
        }

        return result;
    }

    public static double[] triangle_integrand_01(int p_num, double[] p, int f_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_INTEGRAND_01 evaluates 1 integrand function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input, double P[2*P_NUM], the evaluation points.
        //
        //    Input, int F_NUM, the number of integrands.
        //
        //    Output, double FP[F_NUM*P_NUM], the integrand values.
        //
    {
        int j;

        double[] fp = new double[f_num * p_num];

        for (j = 0; j < p_num; j++)
        {
            fp[0 + j * f_num] = 1.0;
        }

        return fp;
    }

    public static double[] triangle_integrand_02(int p_num, double[] p, int f_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_INTEGRAND_02 evaluates 2 integrand functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input, double P[2*P_NUM], the evaluation points.
        //
        //    Input, int F_NUM, the number of integrands.
        //
        //    Output, double FP[F_NUM*P_NUM], the integrand values.
        //
    {
        int j;

        double[] fp = new double[f_num * p_num];

        for (j = 0; j < p_num; j++)
        {
            fp[0 + j * f_num] = p[0 + j * 2];
            fp[1 + j * f_num] = p[1 + j * 2];
        }

        return fp;
    }

    public static double[] triangle_integrand_03(int p_num, double[] p, int f_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_INTEGRAND_03 evaluates 3 integrand functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input, double P[2*P_NUM], the evaluation points.
        //
        //    Input, int F_NUM, the number of integrands.
        //
        //    Output, double FP[F_NUM*P_NUM], the integrand values.
        //
    {
        int j;

        double[] fp = new double[f_num * p_num];

        for (j = 0; j < p_num; j++)
        {
            fp[0 + j * f_num] = p[0 + j * 2] * p[0 + j * 2];
            fp[1 + j * f_num] = p[0 + j * 2] * p[1 + j * 2];
            fp[2 + j * f_num] = p[1 + j * 2] * p[1 + j * 2];
        }

        return fp;
    }

    public static double[] triangle_integrand_04(int p_num, double[] p, int f_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_INTEGRAND_04 evaluates 4 integrand functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input, double P[2*P_NUM], the evaluation points.
        //
        //    Input, int F_NUM, the number of integrands.
        //
        //    Output, double FP[F_NUM*P_NUM], the integrand values.
        //
    {
        int j;

        double[] fp = new double[f_num * p_num];

        for (j = 0; j < p_num; j++)
        {
            fp[0 + j * f_num] = p[0 + j * 2] * p[0 + j * 2] * p[0 + j * 2];
            fp[1 + j * f_num] = p[0 + j * 2] * p[0 + j * 2] * p[1 + j * 2];
            fp[2 + j * f_num] = p[0 + j * 2] * p[1 + j * 2] * p[1 + j * 2];
            fp[3 + j * f_num] = p[1 + j * 2] * p[1 + j * 2] * p[1 + j * 2];
        }

        return fp;
    }

    public static double[] triangle_integrand_05(int p_num, double[] p, int f_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_INTEGRAND_05 evaluates 5 integrand functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input, double P[2*P_NUM], the evaluation points.
        //
        //    Input, int F_NUM, the number of integrands.
        //
        //    Output, double FP[F_NUM*P_NUM], the integrand values.
        //
    {
        int j;

        double[] fp = new double[f_num * p_num];

        for (j = 0; j < p_num; j++)
        {
            fp[0 + j * f_num] = p[0 + j * 2] * p[0 + j * 2] * p[0 + j * 2] * p[0 + j * 2];
            fp[1 + j * f_num] = p[0 + j * 2] * p[0 + j * 2] * p[0 + j * 2] * p[1 + j * 2];
            fp[2 + j * f_num] = p[0 + j * 2] * p[0 + j * 2] * p[1 + j * 2] * p[1 + j * 2];
            fp[3 + j * f_num] = p[0 + j * 2] * p[1 + j * 2] * p[1 + j * 2] * p[1 + j * 2];
            fp[4 + j * f_num] = p[1 + j * 2] * p[1 + j * 2] * p[1 + j * 2] * p[1 + j * 2];
        }

        return fp;
    }

    public static tusData triangle_unit_sample_01(int p_num, int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_SAMPLE_01 selects points from the unit triangle.
        //
        //  Discussion:
        //
        //    The unit triangle has vertices (1,0), (0,1), (0,0).
        //
        //    Any point in the unit simplex CAN be chosen by this algorithm.
        //
        //    However, the points that are chosen tend to be clustered near
        //    the centroid.
        //
        //    This routine is supplied as an example of "bad" sampling.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double TRIANGLE_UNIT_SAMPLE_01[2*P_NUM], the points.
        //
    {
        int j;

        double[] x = new double[2 * p_num];

        for (j = 0; j < p_num; j++)
        {
            double[] e = UniformRNG.r8vec_uniform_01_new(3, ref seed);

            double e_sum = typeMethods.r8vec_sum(3, e);

            int i;
            for (i = 0; i < 2; i++)
            {
                x[i + j * 2] = e[i] / e_sum;
            }

        }

        tusData ret = new()
        {
            seed = seed,
            result = x
        };

        return ret;
    }

    public static tusData triangle_unit_sample_02(int p_num, int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_SAMPLE_02 selects points from the unit triangle.
        //
        //  Discussion:
        //
        //    The unit triangle has vertices (1,0), (0,1), (0,0).
        //
        //    The sampling is uniform.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Greg Turk,
        //    Generating Random Points in a Triangle,
        //    in Graphics Gems,
        //    edited by Andrew Glassner,
        //    AP Professional, 1990, pages 24-28.
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double TRIANGLE_UNIT_SAMPLE_02[2*P_NUM], the points.
        //
    {
        int j;

        double[] x = new double[2 * p_num];

        for (j = 0; j < p_num; j++)
        {
            double[] r = UniformRNG.r8vec_uniform_01_new(2, ref seed);

            switch (r[0] + r[1])
            {
                case > 1.0:
                    r[0] = 1.0 - r[0];
                    r[1] = 1.0 - r[1];
                    break;
            }

            x[0 + j * 2] = r[0];
            x[1 + j * 2] = r[1];

        }

        tusData ret = new()
        {
            seed = seed,
            result = x
        };

        return ret;
    }

    public static tusData triangle_unit_sample_03(int p_num, int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_SAMPLE_03 selects points from the unit triangle.
        //
        //  Discussion:
        //
        //    The unit triangle has vertices (1,0), (0,1), (0,0).
        //
        //    This routine uses Turk's rule 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Greg Turk,
        //    Generating Random Points in a Triangle,
        //    in Graphics Gems,
        //    edited by Andrew Glassner,
        //    AP Professional, 1990, pages 24-28.
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double TRIANGLE_UNIT_SAMPLE_03[2*P_NUM], the points.
        //
    {
        int j;

        double[] x = new double[2 * p_num];

        for (j = 0; j < p_num; j++)
        {
            double[] r = UniformRNG.r8vec_uniform_01_new(2, ref seed);

            r[1] = Math.Sqrt(r[1]);

            double a = 1.0 - r[1];
            double b = (1.0 - r[0]) * r[1];

            x[0 + j * 2] = a;
            x[1 + j * 2] = b;

        }

        tusData ret = new()
        {
            seed = seed,
            result = x
        };

        return ret;
    }

    public static tusData triangle_unit_sample_04(int p_num, int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_SAMPLE_04 selects points from the unit triangle.
        //
        //  Discussion:
        //
        //    The unit triangle has vertices (1,0), (0,1), (0,0).
        //
        //    The sampling is uniform.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Reuven Rubinstein,
        //    Monte Carlo Optimization, Simulation, and Sensitivity
        //    of Queueing Networks,
        //    Krieger, 1992,
        //    ISBN: 0894647644,
        //    LC: QA298.R79.
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double TRIANGLE_UNIT_SAMPLE_04[2*P_NUM], the points.
        //
    {
        int j;
        //
        //  The construction begins by sampling DIM_NUM+1 points from the
        //  exponential distribution with parameter 1.
        //
        double[] x = new double[2 * p_num];

        for (j = 0; j < p_num; j++)
        {
            double[] e = UniformRNG.r8vec_uniform_01_new(3, ref seed);

            int i;
            for (i = 0; i <= 2; i++)
            {
                e[i] = -Math.Log(e[i]);
            }

            double e_sum = typeMethods.r8vec_sum(3, e);

            for (i = 0; i < 2; i++)
            {
                x[i + 2 * j] = e[i] / e_sum;
            }
        }

        tusData ret = new()
        {
            seed = seed,
            result = x
        };

        return ret;
    }

}