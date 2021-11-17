using System;
using Burkardt.Types;
using Burkardt.Vector;

namespace GeometryTest;

public static class ProvecTest
{
    public static void test170 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST170 tests PROVEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 4;
        int N = 2;

        double[] base_ = {
            4.0, 3.0, 2.0, 1.0,
            1.0, 2.0, 3.0, 4.0 };
        double[] vecm = { 1.0, 1.0, 1.0, 2.0 };
        double[] vecn = new double[N];
        double[] vecnm = new double[M];

        Console.WriteLine("");
        Console.WriteLine("TEST170");
        Console.WriteLine("  PROVEC projects a vector onto a subspace.");

        typeMethods.r8mat_transpose_print ( M, N, base_, "  Base vectors" );

        typeMethods.r8vec_print ( M, vecm, "  Vector to be projected:" );

        Geometry.provec ( M, N, base_, vecm, ref vecn, ref vecnm );

        typeMethods.r8vec_print ( N, vecn, "  Projected vector in BASE coordinates:" );

        typeMethods.r8vec_print ( M, vecnm, "  Projected vector in original coordinates:" );

    }

}