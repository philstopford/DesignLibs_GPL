using System;
using System.Globalization;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace SandiaRules2Test;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SANDIA_RULES2_TEST.
        //
        //  Discussion:
        //
        //    SANDIA_RULES2_TEST tests the SANDIA_RULES2 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SANDIA_RULES2_TEST");
        Console.WriteLine("  Test the SANDIA_RULES2 library.");

        test185();

        Console.WriteLine("");
        Console.WriteLine("SANDIA_RULES2_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static JacobiQuadrature.ParameterData parameter(JacobiQuadrature.ParameterData data, int dim, int offset)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARAMETER is a user-supplied routine to retrieve parameters.
        //
        //  Discussion:
        //
        //    The declaration for this function is in SANDIA_RULES.H
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM, the spatial dimension.
        //
        //    Input, int OFFSET, the offset of the parameter within the 
        //    spatial dimension.
        //
        //    Output, double PARAMETER, the value of the OFFSET-th parameter
        //    associated with the DIM-th dimension.
        //
    {
        int i;

        int j = 0;
        for (i = 0; i < dim; i++)
        {
            j += data.NP[i];
        }

        double value = data.P[j + offset];

        data.value = value;
        return data;
    }

    private static void test185()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST185 tests JACOBI_POINTS and JACOBI_WEIGHTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int TEST_NUM = 3;

        double[] alpha_test = {0.5, 1.0, 2.5};
        double[] beta_test = {0.5, 1.0, 2.5};
        const int order_max = 10;
        int test1;
        JacobiQuadrature.ParameterData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST185");
        Console.WriteLine("  JACOBI_POINTS and JACOBI_WEIGHTS compute a Gauss-Jacobi rule");
        Console.WriteLine("  which is appropriate for integrands of the form");
        Console.WriteLine("");
        Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.");
        Console.WriteLine("");
        Console.WriteLine("  For technical reasons, the parameters ALPHA and BETA are to");
        Console.WriteLine("  supplied by a function called PARAMETER.");
        Console.WriteLine("");
        Console.WriteLine("   Order       ALPHA        BETA   ||X1-X2||   ||W1-W2||");
        Console.WriteLine("");

        for (test1 = 0; test1 < TEST_NUM; test1++)
        {
            double alpha = alpha_test[test1];

            int test2;
            for (test2 = 0; test2 < TEST_NUM; test2++)
            {
                double beta = beta_test[test2];

                int dim = 0;
                data.NP = new int[1];
                data.NP[0] = 2;
                data.P = new double[2];
                data.P[0] = alpha;
                data.P[1] = beta;

                int order;
                for (order = 1; order <= order_max; order++)
                {
                    double[] w1 = new double[order];
                    double[] w2 = new double[order];
                    double[] x1 = new double[order];
                    double[] x2 = new double[order];

                    JacobiQuadrature.jacobi_compute(order, alpha, beta, ref x1, ref w1);

                    data = JacobiQuadrature.jacobi_points(parameter, data, order, dim, ref x2);
                    data = JacobiQuadrature.jacobi_weights(parameter, data, order, dim, ref w2);

                    double x_diff = typeMethods.r8vec_diff_norm_li(order, x1, x2);
                    double w_diff = typeMethods.r8vec_diff_norm_li(order, w1, w2);

                    Console.WriteLine("  " + order.ToString().PadLeft(6)
                                           + "  " + alpha.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                           + "  " + beta.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                           + "  " + x_diff.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                           + "  " + w_diff.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
                }
            }
        }
    }
}