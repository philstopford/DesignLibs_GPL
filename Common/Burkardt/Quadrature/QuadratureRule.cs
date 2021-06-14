using System;
using Burkardt.Types;

namespace Burkardt
{
    public static class QuadratureRule
    {
        public static void rule_write ( int order, string filename, double[] x, double[] w, 
        double[] r )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RULE_WRITE writes a quadrature rule to three files.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Input, double A, the left endpoint.
        //
        //    Input, double B, the right endpoint.
        //
        //    Input, string FILENAME, specifies the output filenames.
        //    "filename_w.txt", "filename_x.txt", "filename_r.txt" 
        //    defining weights, abscissas, and region.
        // 
        {
            string filename_r;
            string filename_w;
            string filename_x;

            filename_w = filename + "_w.txt";
            filename_x = filename + "_x.txt";
            filename_r = filename + "_r.txt";

            Console.WriteLine("");
            Console.WriteLine("  Creating quadrature files.");
            Console.WriteLine("");
            Console.WriteLine("  Root file name is     \"" + filename   + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Weight file will be   \"" + filename_w + "\".");
            Console.WriteLine("  Abscissa file will be \"" + filename_x + "\".");
            Console.WriteLine("  Region file will be   \"" + filename_r + "\".");
            
            typeMethods.r8mat_write ( filename_w, 1, order, w );
            typeMethods.r8mat_write ( filename_x, 1, order, x );
            typeMethods.r8mat_write ( filename_r, 1, 2,     r );
  
            return;
        }
    }
}