using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.TetrahedronNS;

public static class MonteCarlo
{
    public class TetrahedronSampleResult
    {
        public int seed;
        public double[] result;
    }

    public static TetrahedronSampleResult tetrahedron_monte_carlo(double[] t, int p_num, int f_num,
            Func<int, int, TetrahedronSampleResult> tetrahedron_unit_sample,
            Func<int, double[], int, double[]> tetrahedron_integrand, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_MONTE_CARLO applies the Monte Carlo rule to integrate a function.
        //
        //  Discussion:
        //
        //    The function f(x,y,z) is to be integrated over a tetrahedron T.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[3*4], the tetrahedron vertices.
        //
        //    Input, int P_NUM, the number of sample points.
        //
        //    Input, int F_NUM, the number of functions to integrate.
        //
        //    Input, external TETRAHEDRON_UNIT_SAMPLE, the sampling routine.
        //
        //    Input, external TETRAHEDRON_INTEGRAND, the integrand routine.
        //
        //    Input/output, int *SEED, a seed for the random 
        //    number generator.
        //
        //    Output, dobule TETRAHEDRON_MONTE_CARLO[F_NUM], the approximate integrals.
        //
    {
        int i;
        double[] p2 = new double[1];

        double volume = Tetrahedron.tetrahedron_volume(t);

        TetrahedronSampleResult tmp = tetrahedron_unit_sample(p_num, seed);
        double[] p = tmp.result;
        seed = tmp.seed;
        p2 = Tetrahedron.reference_to_physical_tet4(t, p_num, p);

        double[] fp = tetrahedron_integrand(p_num, p2, f_num);

        double[] result = new double[f_num];

        for (i = 0; i < f_num; i++)
        {
            double fp_sum = 0.0;
            int j;
            for (j = 0; j < p_num; j++)
            {
                fp_sum += fp[i + j * f_num];
            }

            result[i] = volume * fp_sum / p_num;
        }

        return new TetrahedronSampleResult { result = result, seed = seed };
    }

    public static TetrahedronSampleResult tetrahedron_unit_sample_01(int p_num, int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_UNIT_SAMPLE_01 selects points from the unit tetrahedron.
        //
        //  Discussion:
        //
        //    The unit tetrahedron has vertices (1,0,0), (0,1,0), (0,0,1), (0,0,0).
        //
        //    Any point in the unit tetrahedron CAN be chosen by this algorithm.
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
        //    16 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input/output, int *SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double TETRAHEDRON_UNIT_SAMPLE_01[3*P_NUM], the points.
        //
    {
        int j;

        double[] p = new double[3 * p_num];

        for (j = 0; j < p_num; j++)
        {
            double[] e = UniformRNG.r8vec_uniform_01_new(4, ref seed);

            double e_sum = typeMethods.r8vec_sum(4, e);

            int i;
            for (i = 0; i < 3; i++)
            {
                p[i + j * 3] = e[i] / e_sum;
            }
        }

        return new TetrahedronSampleResult { result = p, seed = seed };
    }

    public static TetrahedronSampleResult tetrahedron_unit_sample_02(int p_num, int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_UNIT_SAMPLE_02 selects points from the unit tetrahedron.
        //
        //  Discussion:
        //
        //    The unit tetrahedron has vertices (1,0,0), (0,1,0), (0,0,1), (0,0,0).
        //
        //    The sampling is uniform.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Claudio Rocchini, Paolo Cignoni,
        //    Generating Random Points in a Tetrahedron,
        //    Journal of Graphics Tools,
        //    Volume 5, Number 5, 2000, pages 9-12.
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input/output, int *SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double TETRAHEDRON_UNIT_SAMPLE_02[3*P_NUM], the points.
        //
    {
        int j;

        double[] p = new double[3 * p_num];

        for (j = 0; j < p_num; j++)
        {
            double[] c = UniformRNG.r8vec_uniform_01_new(3, ref seed);

            switch (c[0] + c[1])
            {
                case > 1.0:
                    c[0] = 1.0 - c[0];
                    c[1] = 1.0 - c[1];
                    break;
            }

            double t;
            switch (c[1] + c[2])
            {
                case > 1.0:
                    t = c[2];
                    c[2] = 1.0 - c[0] - c[1];
                    c[1] = 1.0 - t;
                    break;
                default:
                {
                    switch (c[0] + c[1] + c[2])
                    {
                        case > 1.0:
                            t = c[2];
                            c[2] = c[0] + c[1] + c[2] - 1.0;
                            c[0] = 1.0 - c[1] - t;
                            break;
                    }

                    break;
                }
            }

            p[0 + j * 3] = c[0];
            p[1 + j * 3] = c[1];
            p[2 + j * 3] = c[2];

        }

        return new TetrahedronSampleResult { result = p, seed = seed };
    }

    public static TetrahedronSampleResult tetrahedron_unit_sample_03(int p_num, int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_UNIT_SAMPLE_03 selects points from the unit tetrahedron.
        //
        //  Discussion:
        //
        //    The unit tetrahedron has vertices (1,0,0), (0,1,0), (0,0,1), (0,0,0).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 August 2009
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
        //    Input/output, int *SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double TETRAHEDRON_UNIT_SAMPLE_01[3*P_NUM], the points.
        //
    {
        int j;

        double[] p = new double[3 * p_num];

        for (j = 0; j < p_num; j++)
        {
            double[] r = UniformRNG.r8vec_uniform_01_new(3, ref seed);

            double e = Math.Pow(r[0], 1.0 / 3.0);
            double f = Math.Sqrt(r[1]);
            double g = r[2];

            double a = 1.0 - e;
            double b = (1.0 - f) * e;
            double c = (1.0 - g) * f * e;

            p[0 + j * 3] = a;
            p[1 + j * 3] = b;
            p[2 + j * 3] = c;

        }

        return new TetrahedronSampleResult { result = p, seed = seed };
    }

    public static TetrahedronSampleResult tetrahedron_unit_sample_04(int p_num, int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_UNIT_SAMPLE_04 selects points from the unit tetrahedron.
        //
        //  Discussion:
        //
        //    The unit tetrahedron has vertices (1,0,0), (0,1,0), (0,0,1), (0,0,0).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 August 2009
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
        //    Input/output, int *SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double TETRAHEDRON_UNIT_SAMPLE_01[3*P_NUM], the points.
        //
    {
        int j;

        double[] p = new double[3 * p_num];
        //
        //  The construction begins by sampling DIM_NUM+1 points from the
        //  exponential distribution with parameter 1.
        //
        for (j = 0; j < p_num; j++)
        {
            double[] e = UniformRNG.r8vec_uniform_01_new(4, ref seed);

            int i;
            for (i = 0; i < 4; i++)
            {
                e[i] = -Math.Log(e[i]);
            }

            double e_sum = typeMethods.r8vec_sum(4, e);

            for (i = 0; i < 3; i++)
            {
                p[i + j * 3] = e[i] / e_sum;
            }
        }

        return new TetrahedronSampleResult { result = p, seed = seed };
    }
}