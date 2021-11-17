using System;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace GeometryTest;

public static class BasisTest
{
    public static void basis_map_3d_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_MAP_3D_TEST tests BASIS_MAP_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a = null;
        double[] c = new double[3*3];
        int i;
        int j;
        int k;
        double[] u = {
            1.0, 2.0, 3.0,
            0.0, 0.0, 1.0,
            1.0, 0.0, 2.0 };
        double[] v = {
            14.0, 4.0, 4.0,
            3.0, 1.0, 0.0,
            7.0, 3.0, 2.0 };

        Console.WriteLine("");
        Console.WriteLine("BASIS_MAP_3D_TEST");
        Console.WriteLine("  BASIS_MAP_3D computes the linear transform A");
        Console.WriteLine("  which maps vectors U1, U2 and U3 to vectors");
        Console.WriteLine("  V1, V2 and V3.");

        typeMethods.r8mat_print ( 3, 3, u, "  The matrix U" );

        typeMethods.r8mat_print ( 3, 3, v, "  The matrix V" );

        a = Matrix.basis_map_3d ( u, v );

        switch (a)
        {
            case null:
                Console.WriteLine("");
                Console.WriteLine("  The matrix [ U1 | U2 | U3 ] was singular.");
                Console.WriteLine("  No transformation was computed.");
                return;
        }

        typeMethods.r8mat_print ( 3, 3, a, "  The transformation matrix" );

        for ( i = 0; i < 3; i++ )
        {
            for ( j = 0; j < 3; j++ )
            {
                c[i+j*3] = 0.0;
                for ( k = 0; k < 3; k++ )
                {
                    c[i+j*3] += a[i+k*3] * u[k+j*3];
                }
            }
        }

        typeMethods.r8mat_print ( 3, 3, c, "  The product matrix A * [ U1 | U2 | U3 ]" );
    }
}