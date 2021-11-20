using System;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace Burkardt.PyramidNS;

public static class QuadratureRule
{
    public static void pyramid_handle(int legendre_order, int jacobi_order, string filename)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_HANDLE computes the requested pyramid rule and outputs it.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int LEGENDRE_ORDER, JACOBI_ORDER, the orders 
        //    of the component Legendre and Jacobi rules.
        //
        //    Input, string FILENAME, the rootname for the files,  
        //    write files 'file_w.txt' and 'file_x.txt', and 'file_r.txt', weights,
        //    abscissas, and region.
        //
    {
        const int DIM_NUM = 3;

        int k;
        double[] pyramid_r =
        {
            -1.0, -1.0, 0.0,
            +1.0, -1.0, 0.0,
            -1.0, +1.0, 0.0,
            +1.0, +1.0, 0.0,
            0.0, 0.0, 1.0
        };
        //
        //  Compute the factor rules.
        //
        double[] legendre_w = new double[legendre_order];
        double[] legendre_x = new double[legendre_order];

        Legendre.QuadratureRule.legendre_compute(legendre_order, ref legendre_x, ref legendre_w);

        double[] jacobi_w = new double[jacobi_order];
        double[] jacobi_x = new double[jacobi_order];

        const double jacobi_alpha = 2.0;
        const double jacobi_beta = 0.0;

        JacobiQuadrature.jacobi_compute(jacobi_order, jacobi_alpha, jacobi_beta, ref jacobi_x, ref jacobi_w);
        //
        //  Compute the pyramid rule.
        //
        int pyramid_order = legendre_order * legendre_order * jacobi_order;

        double[] pyramid_w = new double[pyramid_order];
        double[] pyramid_x = new double[DIM_NUM * pyramid_order];

        const double volume = 4.0 / 3.0;

        int l = 0;
        for (k = 0; k < jacobi_order; k++)
        {
            double xk = (jacobi_x[k] + 1.0) / 2.0;
            double wk = jacobi_w[k] / 2.0;
            int j;
            for (j = 0; j < legendre_order; j++)
            {
                double xj = legendre_x[j];
                double wj = legendre_w[j];
                int i;
                for (i = 0; i < legendre_order; i++)
                {
                    double xi = legendre_x[i];
                    double wi = legendre_w[i];
                    pyramid_w[l] = wi * wj * wk / 4.0 / volume;
                    pyramid_x[0 + l * 3] = xi * (1.0 - xk);
                    pyramid_x[1 + l * 3] = xj * (1.0 - xk);
                    pyramid_x[2 + l * 3] = xk;
                    l += 1;
                }
            }
        }

        //
        //  Write the rule to files.
        //
        string filename_w = filename + "_w.txt";
        string filename_x = filename + "_x.txt";
        string filename_r = filename + "_r.txt";

        Console.WriteLine("");
        Console.WriteLine("  Creating quadrature files.");
        Console.WriteLine("");
        Console.WriteLine("  \"Root\" file name is   \"" + filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Weight file will be   \"" + filename_w + "\".");
        Console.WriteLine("  Abscissa file will be \"" + filename_x + "\".");
        Console.WriteLine("  Region file will be   \"" + filename_r + "\".");

        typeMethods.r8mat_write(filename_w, 1, pyramid_order, pyramid_w);
        typeMethods.r8mat_write(filename_x, DIM_NUM, pyramid_order, pyramid_x);
        typeMethods.r8mat_write(filename_r, DIM_NUM, 5, pyramid_r);

    }
}