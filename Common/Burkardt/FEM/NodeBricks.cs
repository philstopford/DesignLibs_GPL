using Burkardt.Types;

namespace Burkardt.FEM
{
    public class NodeBricks
    {
        public static double[] nodes_brick8()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NODES_BRICK8 returns the natural coordinates of the BRICK8 element.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 March 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double NODES_BRICK8[3*8], the coordinates.
            //
        {
            double[] p;
            double[] p_save = {
                -1.0, -1.0, -1.0,
                +1.0, -1.0, -1.0,
                +1.0, +1.0, -1.0,
                -1.0, +1.0, -1.0,
                -1.0, -1.0, +1.0,
                +1.0, -1.0, +1.0,
                +1.0, +1.0, +1.0,
                -1.0, +1.0, +1.0
            }
            ;

            p = typeMethods.r8mat_copy_new(3, 8, p_save);

            return p;
        }

        public static double[] nodes_brick20()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NODES_BRICK20 returns the natural coordinates of the BRICK20 element.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 March 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double NODES_BRICK20[3*20], the coordinates.
            //
        {
            double[] p;
            double[] p_save = {
                -1.0, -1.0, -1.0,
                +1.0, -1.0, -1.0,
                +1.0, +1.0, -1.0,
                -1.0, +1.0, -1.0,
                -1.0, -1.0, +1.0,
                +1.0, -1.0, +1.0,
                +1.0, +1.0, +1.0,
                -1.0, +1.0, +1.0,
                0.0, -1.0, -1.0,
                +1.0, 0.0, -1.0,
                0.0, +1.0, -1.0,
                -1.0, 0.0, -1.0,
                -1.0, -1.0, 0.0,
                +1.0, -1.0, 0.0,
                +1.0, +1.0, 0.0,
                -1.0, +1.0, 0.0,
                0.0, -1.0, +1.0,
                +1.0, 0.0, +1.0,
                0.0, +1.0, +1.0,
                -1.0, 0.0, +1.0
            }
            ;

            p = typeMethods.r8mat_copy_new(3, 20, p_save);

            return p;
        }

        public static double[] nodes_brick27()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NODES_BRICK27 returns the natural coordinates of the BRICK27 element.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 March 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double NODES_BRICK27[3*27], the coordinates.
            //
        {
            double[] p;
            double[] p_save = {
                -1.0, -1.0, -1.0,
                +1.0, -1.0, -1.0,
                +1.0, +1.0, -1.0,
                -1.0, +1.0, -1.0,
                -1.0, -1.0, +1.0,
                +1.0, -1.0, +1.0,
                +1.0, +1.0, +1.0,
                -1.0, +1.0, +1.0,
                0.0, -1.0, -1.0,
                +1.0, 0.0, -1.0,
                0.0, +1.0, -1.0,
                -1.0, 0.0, -1.0,
                -1.0, -1.0, 0.0,
                +1.0, -1.0, 0.0,
                +1.0, +1.0, 0.0,
                -1.0, +1.0, 0.0,
                0.0, -1.0, +1.0,
                +1.0, 0.0, +1.0,
                0.0, +1.0, +1.0,
                -1.0, 0.0, +1.0,
                0.0, 0.0, -1.0,
                0.0, -1.0, 0.0,
                +1.0, 0.0, 0.0,
                0.0, +1.0, 0.0,
                -1.0, 0.0, 0.0,
                0.0, 0.0, +1.0,
                0.0, 0.0, 0.0
            }
            ;

            p = typeMethods.r8mat_copy_new(3, 27, p_save);

            return p;
        }
    }
}