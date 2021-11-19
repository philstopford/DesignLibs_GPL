using System;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace SandiaRules2Test;

internal class Program
{
    private static void Main(string[] args)
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
        int j;
        double value = 0;

        j = 0;
        for (i = 0; i < dim; i++)
        {
            j += data.NP[i];
        }

        value = data.P[j + offset];

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
        int TEST_NUM = 3;

        double alpha;
        double[] alpha_test = {0.5, 1.0, 2.5};
        double beta;
        double[] beta_test = {0.5, 1.0, 2.5};
        int dim;
        int order;
        int order_max = 10;
        int test1;
        int test2;
        double w_diff;
        double[] w1;
        double[] w2;
        double x_diff;
        double[] x1;
        double[] x2;
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
            alpha = alpha_test[test1];

            for (test2 = 0; test2 < TEST_NUM; test2++)
            {
                beta = beta_test[test2];

                dim = 0;
                data.NP = new int[1];
                data.NP[0] = 2;
                data.P = new double[2];
                data.P[0] = alpha;
                data.P[1] = beta;

                for (order = 1; order <= order_max; order++)
                {
                    w1 = new double[order];
                    w2 = new double[order];
                    x1 = new double[order];
                    x2 = new double[order];

                    JacobiQuadrature.jacobi_compute(order, alpha, beta, ref x1, ref w1);

                    data = JacobiQuadrature.jacobi_points(parameter, data, order, dim, ref x2);
                    data = JacobiQuadrature.jacobi_weights(parameter, data, order, dim, ref w2);

                    x_diff = typeMethods.r8vec_diff_norm_li(order, x1, x2);
                    w_diff = typeMethods.r8vec_diff_norm_li(order, w1, w2);

                    Console.WriteLine("  " + order.ToString().PadLeft(6)
                                           + "  " + alpha.ToString().PadLeft(10)
                                           + "  " + beta.ToString().PadLeft(10)
                                           + "  " + x_diff.ToString().PadLeft(10)
                                           + "  " + w_diff.ToString().PadLeft(10) + "");
                }
            }
        }
    }
}