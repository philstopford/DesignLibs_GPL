using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void compare ( int node_num, double[] node_xy, int[] indx, int nunk,
            double[] f, Func<double, double, ExactResult> exact )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMPARE compares the exact and computed solution at the nodes.
        //
        //  Discussion:
        //
        //    This is a rough comparison, done only at the nodes.  Such a pointwise
        //    comparison is easy, because the value of the finite element
        //    solution is exactly the value of the finite element coefficient
        //    associated with that node.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the nodes.
        //
        //    Input, int INDX[NODE_NUM], the index of the unknown in the finite
        //    element linear system.
        //
        //    Input, int NUNK, the number of unknowns in the finite element system.
        //
        //    Input, double F[NUNK], the solution vector of the finite
        //    element system.
        //
    {
        int node;

        Console.WriteLine("");
        Console.WriteLine("COMPARE:");
        Console.WriteLine("  Compare computed and exact solutions at the nodes.");
        Console.WriteLine("");
        Console.WriteLine("         X           Y          U           U");
        Console.WriteLine("                             computed     exact");
        Console.WriteLine("");

        for ( node = 0; node < node_num; node++ )
        {
            double x = node_xy[0+node*2];
            double y = node_xy[1+node*2];

            ExactResult res = exact ( x, y );
            double u = res.u;

            int i = indx[node];
            double uh = f[i-1];

            Console.WriteLine(x.ToString().PadLeft(12)  + "  "
                                                        + y.ToString().PadLeft(12)  + "  "
                                                        + uh.ToString().PadLeft(12) + "  "
                                                        + u.ToString().PadLeft(12)  + "");
        }
    }
    public static void solution_write ( double[] f, int[] indx, int node_num, int nunk,
            string output_filename, double[] node_xy, Func<double, double, ExactResult> exact )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SOLUTION_WRITE writes the solution to a file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double F[NUNK], the coefficients of the solution.
        //
        //    Input, int INDX[NODE_NUM], gives the index of the unknown quantity
        //    associated with the given node.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int NUNK, the number of unknowns.
        //
        //    Input, string OUTPUT_FILENAME, the name of the file
        //    in which the data should be stored.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
        //
    {
        int node;
        List<string> output = new();

        for ( node = 0; node < node_num; node++ )
        {
            double x = node_xy[0+node*2];
            double y = node_xy[1+node*2];

            double u;
            switch (indx[node])
            {
                case > 0:
                    u = f[indx[node]-1];
                    break;
                default:
                {
                    ExactResult res = exact ( x, y );
                    u = res.u;
                    break;
                }
            }

            output.Add(u.ToString().PadLeft(14));
        }

        try
        {
            File.WriteAllLines(output_filename, output);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("SOLUTION_WRITE - Warning!");
            Console.WriteLine("  Could not write the solution file.");
        }
    }
}