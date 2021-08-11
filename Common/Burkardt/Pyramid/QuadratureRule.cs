using System;
using Burkardt.Types;

namespace Burkardt.Pyramid
{
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
            int DIM_NUM = 3;

            string filename_r;
            string filename_w;
            string filename_x;
            int i;
            int j;
            double jacobi_alpha;
            double jacobi_beta;
            double[] jacobi_w;
            double[] jacobi_x;
            int k;
            int l;
            double[] legendre_w;
            double[] legendre_x;
            int pyramid_order;
            double[] pyramid_r =
            {
                -1.0, -1.0, 0.0,
                +1.0, -1.0, 0.0,
                -1.0, +1.0, 0.0,
                +1.0, +1.0, 0.0,
                0.0, 0.0, 1.0
            };
            double[] pyramid_w;
            double[] pyramid_x;
            double volume;
            double wi;
            double wj;
            double wk;
            double xi;
            double xj;
            double xk;
            //
            //  Compute the factor rules.
            //
            legendre_w = new double[legendre_order];
            legendre_x = new double[legendre_order];

            Legendre.QuadratureRule.legendre_compute(legendre_order, ref legendre_x, ref legendre_w);

            jacobi_w = new double[jacobi_order];
            jacobi_x = new double[jacobi_order];

            jacobi_alpha = 2.0;
            jacobi_beta = 0.0;

            JacobiQuadrature.jacobi_compute(jacobi_order, jacobi_alpha, jacobi_beta, ref jacobi_x, ref jacobi_w);
            //
            //  Compute the pyramid rule.
            //
            pyramid_order = legendre_order * legendre_order * jacobi_order;

            pyramid_w = new double[pyramid_order];
            pyramid_x = new double[DIM_NUM * pyramid_order];

            volume = 4.0 / 3.0;

            l = 0;
            for (k = 0; k < jacobi_order; k++)
            {
                xk = (jacobi_x[k] + 1.0) / 2.0;
                wk = jacobi_w[k] / 2.0;
                for (j = 0; j < legendre_order; j++)
                {
                    xj = legendre_x[j];
                    wj = legendre_w[j];
                    for (i = 0; i < legendre_order; i++)
                    {
                        xi = legendre_x[i];
                        wi = legendre_w[i];
                        pyramid_w[l] = wi * wj * wk / 4.0 / volume;
                        pyramid_x[0 + l * 3] = xi * (1.0 - xk);
                        pyramid_x[1 + l * 3] = xj * (1.0 - xk);
                        pyramid_x[2 + l * 3] = xk;
                        l = l + 1;
                    }
                }
            }

            //
            //  Write the rule to files.
            //
            filename_w = filename + "_w.txt";
            filename_x = filename + "_x.txt";
            filename_r = filename + "_r.txt";

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
}