using System;
using System.Globalization;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.FEM;

public static class BasisBrick
{
    public static double[] basis_brick8(int n, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_BRICK8: BRICK8 basis functions at natural coordinates.
        //
        //  Discussion:
        //
        //      8------7        t  s
        //     /|     /|        | /
        //    5------6 |        |/
        //    | |    | |        0-------r
        //    | 4----|-3        
        //    |/     |/        
        //    1------2        
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
        //    Input, int N, the number of evaluation points.
        //
        //    Input, double P[3*N], natural coordinates of evaluation
        //    points.
        //
        //    Output, double BASIS_BRICK8[8*N], the basis function values.
        //
    {
        int j;

        double[] phi = new double[8 * n];

        for (j = 0; j < n; j++)
        {
            phi[0 + j * 8] =
                (1.0 - p[0 + j * 3]) * (1.0 - p[1 + j * 3]) * (1.0 - p[2 + j * 3]) / 8.0;
            phi[1 + j * 8] =
                (1.0 + p[0 + j * 3]) * (1.0 - p[1 + j * 3]) * (1.0 - p[2 + j * 3]) / 8.0;
            phi[2 + j * 8] =
                (1.0 + p[0 + j * 3]) * (1.0 + p[1 + j * 3]) * (1.0 - p[2 + j * 3]) / 8.0;
            phi[3 + j * 8] =
                (1.0 - p[0 + j * 3]) * (1.0 + p[1 + j * 3]) * (1.0 - p[2 + j * 3]) / 8.0;
            phi[4 + j * 8] =
                (1.0 - p[0 + j * 3]) * (1.0 - p[1 + j * 3]) * (1.0 + p[2 + j * 3]) / 8.0;
            phi[5 + j * 8] =
                (1.0 + p[0 + j * 3]) * (1.0 - p[1 + j * 3]) * (1.0 + p[2 + j * 3]) / 8.0;
            phi[6 + j * 8] =
                (1.0 + p[0 + j * 3]) * (1.0 + p[1 + j * 3]) * (1.0 + p[2 + j * 3]) / 8.0;
            phi[7 + j * 8] =
                (1.0 - p[0 + j * 3]) * (1.0 + p[1 + j * 3]) * (1.0 + p[2 + j * 3]) / 8.0;
        }

        return phi;
    }

    public static void basis_brick8_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_BRICK8_TEST verifies BASIS_BRICK8.
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
        //    None
        //
    {
        int j;
        const int node_num = 8;
        const int test_num = 5;

        Console.WriteLine("");
        Console.WriteLine("BASIS_BRICK8_TEST:");
        Console.WriteLine("  Verify basis functions for element BRICK8.");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + node_num + "");
        Console.WriteLine("");
        Console.WriteLine("  The basis function values at basis nodes");
        Console.WriteLine("  should form the identity matrix.");
        Console.WriteLine("");

        int n = node_num;

        double[] p = NodeBricks.nodes_brick8();

        double[] phi = basis_brick8(n, p);

        for (j = 0; j < n; j++)
        {
            string cout = "  ";
            int i;
            for (i = 0; i < node_num; i++)
            {
                cout += phi[i + j * node_num].ToString(CultureInfo.InvariantCulture).PadLeft(7);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  The basis function values at ANY point P");
        Console.WriteLine("  should sum to 1:");
        Console.WriteLine("");
        Console.WriteLine("    ------------P-------------     PHI_SUM");
        Console.WriteLine("");

        n = test_num;
        int seed = 123456789;

        p = UniformRNG.r8mat_uniform_01_new(3, n, ref seed);

        phi = basis_brick8(n, p);

        for (j = 0; j < n; j++)
        {
            double phi_sum = typeMethods.r8vec_sum(node_num, phi, aIndex: + j * node_num);
            Console.WriteLine("  " + p[0 + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + p[1 + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + p[2 + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + phi_sum.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }
    }

    public static double[] basis_brick20(int n, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_BRICK20: BRICK20 basis functions at natural coordinates.
        //
        //  Discussion:
        //
        //        8----19---7
        //       /|        /|
        //     20 |      18 |        t   s
        //     /  16     /  15       |  /
        //    5----17---6   |        | /
        //    |   |     |   |        |/
        //    |   4--11-|---3        0---------r
        //   13  /     14  /        
        //    | 12      | 10       
        //    |/        |/        
        //    1----9----2
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
        //    Input, int N, the number of evaluation points.
        //
        //    Input, double P[3*N], natural coordinates of evaluation
        //    points.
        //
        //    Output, double BASIS_BRICK20[20*N], the basis function values.
        //
    {
        int j;

        double[] phi = new double[20 * n];

        for (j = 0; j < n; j++)
        {
            phi[0 + j * 20] =
                (1.0 - p[0 + j * 3]) * (1.0 - p[1 + j * 3]) * (1.0 - p[2 + j * 3])
                * (-p[0 + j * 3] - p[1 + j * 3] - p[2 + j * 3] - 2.0) / 8.0;
            phi[1 + j * 20] =
                (1.0 + p[0 + j * 3]) * (1.0 - p[1 + j * 3]) * (1.0 - p[2 + j * 3])
                * (+p[0 + j * 3] - p[1 + j * 3] - p[2 + j * 3] - 2.0) / 8.0;
            phi[2 + j * 20] =
                (1.0 + p[0 + j * 3]) * (1.0 + p[1 + j * 3]) * (1.0 - p[2 + j * 3])
                * (+p[0 + j * 3] + p[1 + j * 3] - p[2 + j * 3] - 2.0) / 8.0;
            phi[3 + j * 20] =
                (1.0 - p[0 + j * 3]) * (1.0 + p[1 + j * 3]) * (1.0 - p[2 + j * 3])
                * (-p[0 + j * 3] + p[1 + j * 3] - p[2 + j * 3] - 2.0) / 8.0;
            phi[4 + j * 20] =
                (1.0 - p[0 + j * 3]) * (1.0 - p[1 + j * 3]) * (1.0 + p[2 + j * 3])
                * (-p[0 + j * 3] - p[1 + j * 3] + p[2 + j * 3] - 2.0) / 8.0;
            phi[5 + j * 20] =
                (1.0 + p[0 + j * 3]) * (1.0 - p[1 + j * 3]) * (1.0 + p[2 + j * 3])
                * (+p[0 + j * 3] - p[1 + j * 3] + p[2 + j * 3] - 2.0) / 8.0;
            phi[6 + j * 20] =
                (1.0 + p[0 + j * 3]) * (1.0 + p[1 + j * 3]) * (1.0 + p[2 + j * 3])
                * (+p[0 + j * 3] + p[1 + j * 3] + p[2 + j * 3] - 2.0) / 8.0;
            phi[7 + j * 20] =
                (1.0 - p[0 + j * 3]) * (1.0 + p[1 + j * 3]) * (1.0 + p[2 + j * 3])
                * (-p[0 + j * 3] + p[1 + j * 3] + p[2 + j * 3] - 2.0) / 8.0;

            phi[8 + j * 20] = (1.0 + p[0 + j * 3]) * (1.0 - p[0 + j * 3])
                                                   * (1.0 - p[1 + j * 3]) * (1.0 - p[2 + j * 3]) / 4.0;
            phi[9 + j * 20] = (1.0 + p[0 + j * 3]) * (1.0 + p[1 + j * 3])
                                                   * (1.0 - p[1 + j * 3]) * (1.0 - p[2 + j * 3]) / 4.0;
            phi[10 + j * 20] = (1.0 + p[0 + j * 3]) * (1.0 - p[0 + j * 3])
                                                    * (1.0 + p[1 + j * 3]) * (1.0 - p[2 + j * 3]) / 4.0;
            phi[11 + j * 20] = (1.0 - p[0 + j * 3]) * (1.0 + p[1 + j * 3])
                                                    * (1.0 - p[1 + j * 3]) * (1.0 - p[2 + j * 3]) / 4.0;
            phi[12 + j * 20] = (1.0 - p[0 + j * 3]) * (1.0 - p[1 + j * 3])
                                                    * (1.0 + p[2 + j * 3]) * (1.0 - p[2 + j * 3]) / 4.0;
            phi[13 + j * 20] = (1.0 + p[0 + j * 3]) * (1.0 - p[1 + j * 3])
                                                    * (1.0 + p[2 + j * 3]) * (1.0 - p[2 + j * 3]) / 4.0;
            phi[14 + j * 20] = (1.0 + p[0 + j * 3]) * (1.0 + p[1 + j * 3])
                                                    * (1.0 + p[2 + j * 3]) * (1.0 - p[2 + j * 3]) / 4.0;
            phi[15 + j * 20] = (1.0 - p[0 + j * 3]) * (1.0 + p[1 + j * 3])
                                                    * (1.0 + p[2 + j * 3]) * (1.0 - p[2 + j * 3]) / 4.0;
            phi[16 + j * 20] = (1.0 + p[0 + j * 3]) * (1.0 - p[0 + j * 3])
                                                    * (1.0 - p[1 + j * 3]) * (1.0 + p[2 + j * 3]) / 4.0;
            phi[17 + j * 20] = (1.0 + p[0 + j * 3]) * (1.0 + p[1 + j * 3])
                                                    * (1.0 - p[1 + j * 3]) * (1.0 + p[2 + j * 3]) / 4.0;
            phi[18 + j * 20] = (1.0 + p[0 + j * 3]) * (1.0 - p[0 + j * 3])
                                                    * (1.0 + p[1 + j * 3]) * (1.0 + p[2 + j * 3]) / 4.0;
            phi[19 + j * 20] = (1.0 - p[0 + j * 3]) * (1.0 + p[1 + j * 3])
                                                    * (1.0 - p[1 + j * 3]) * (1.0 + p[2 + j * 3]) / 4.0;
        }

        return phi;
    }
    //****************************************************************************80

    public static void basis_brick20_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_BRICK20_TEST verifies BASIS_BRICK20.
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
        //    None
        //
    {
        int j;
        const int node_num = 20;
        const int test_num = 5;

        Console.WriteLine("");
        Console.WriteLine("BASIS_BRICK20_TEST:");
        Console.WriteLine("  Verify basis functions for element BRICK20.");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + node_num + "");
        Console.WriteLine("");
        Console.WriteLine("  The basis function values at basis nodes");
        Console.WriteLine("  should form the identity matrix.");
        Console.WriteLine("");

        int n = node_num;

        double[] p = NodeBricks.nodes_brick20();

        double[] phi = basis_brick20(n, p);

        for (j = 0; j < n; j++)
        {
            string cout = "  ";
            int i;
            for (i = 0; i < node_num; i++)
            {
                cout += phi[i + j * node_num].ToString(CultureInfo.InvariantCulture).PadLeft(7);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  The basis function values at ANY point P");
        Console.WriteLine("  should sum to 1:");
        Console.WriteLine("");
        Console.WriteLine("    ------------P-------------     PHI_SUM");
        Console.WriteLine("");

        n = test_num;
        int seed = 123456789;

        p = UniformRNG.r8mat_uniform_01_new(3, n, ref seed);

        phi = basis_brick20(n, p);

        for (j = 0; j < n; j++)
        {
            double phi_sum = typeMethods.r8vec_sum(node_num, phi, aIndex: + j * node_num);
            Console.WriteLine("  " + p[0 + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + p[1 + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + p[2 + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + phi_sum.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }
    }

    public static double[] basis_brick27(int n, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_BRICK20: BRICK20 basis functions at natural coordinates.
        //
        //  Discussion:
        //
        //        8----19---7
        //       /|         /
        //     20 |   26   /|
        //     /          / |
        //    5----17----6  |
        //    |   |      |  |
        //    |  16---24-|-15
        //    |  /|      | /|
        //    |25 |  27  |23|        t
        //    |/         |/ |        |   s
        //   13----22---14  |        |  /
        //    |   |      |  |        | /
        //    |   |      |  |        |/
        //    |   4--11--|--3        0---------r
        //    |  /       | /        
        //    | 12   21  |10       
        //    |/         |/        
        //    1----9-----2
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
        //    Input, int N, the number of evaluation points.
        //
        //    Input, double P[3*N], natural coordinates of evaluation
        //    points.
        //
        //    Output, double BASIS_BRICK27[27*N)], the basis function values.
        //
    {
        int j;

        double[] phi = new double[27 * n];

        for (j = 0; j < n; j++)
        {
            double rm = p[0 + j * 3] + 1.0;
            double rz = p[0 + j * 3];
            double rp = p[0 + j * 3] - 1.0;

            double sm = p[1 + j * 3] + 1.0;
            double sz = p[1 + j * 3];
            double sp = p[1 + j * 3] - 1.0;

            double tm = p[2 + j * 3] + 1.0;
            double tz = p[2 + j * 3];
            double tp = p[2 + j * 3] - 1.0;

            phi[0 + j * 27] = rz * rp * sz * sp * tz * tp / 8.0;
            phi[1 + j * 27] = rm * rz * sz * sp * tz * tp / 8.0;
            phi[2 + j * 27] = rm * rz * sm * sz * tz * tp / 8.0;
            phi[3 + j * 27] = rz * rp * sm * sz * tz * tp / 8.0;
            phi[4 + j * 27] = rz * rp * sz * sp * tm * tz / 8.0;
            phi[5 + j * 27] = rm * rz * sz * sp * tm * tz / 8.0;
            phi[6 + j * 27] = rm * rz * sm * sz * tm * tz / 8.0;
            phi[7 + j * 27] = rz * rp * sm * sz * tm * tz / 8.0;

            phi[8 + j * 27] = -rm * rp * sz * sp * tz * tp / 4.0;
            phi[9 + j * 27] = -rm * rz * sm * sp * tz * tp / 4.0;
            phi[10 + j * 27] = -rm * rp * sm * sz * tz * tp / 4.0;
            phi[11 + j * 27] = -rz * rp * sm * sp * tz * tp / 4.0;
            phi[12 + j * 27] = -rz * rp * sz * sp * tm * tp / 4.0;
            phi[13 + j * 27] = -rm * rz * sz * sp * tm * tp / 4.0;
            phi[14 + j * 27] = -rm * rz * sm * sz * tm * tp / 4.0;
            phi[15 + j * 27] = -rz * rp * sm * sz * tm * tp / 4.0;
            phi[16 + j * 27] = -rm * rp * sz * sp * tm * tz / 4.0;
            phi[17 + j * 27] = -rm * rz * sm * sp * tm * tz / 4.0;
            phi[18 + j * 27] = -rm * rp * sm * sz * tm * tz / 4.0;
            phi[19 + j * 27] = -rz * rp * sm * sp * tm * tz / 4.0;

            phi[20 + j * 27] = rm * rp * sm * sp * tz * tp / 2.0;
            phi[21 + j * 27] = rm * rp * sz * sp * tm * tp / 2.0;
            phi[22 + j * 27] = rm * rz * sm * sp * tm * tp / 2.0;
            phi[23 + j * 27] = rm * rp * sm * sz * tm * tp / 2.0;
            phi[24 + j * 27] = rz * rp * sm * sp * tm * tp / 2.0;
            phi[25 + j * 27] = rm * rp * sm * sp * tm * tz / 2.0;

            phi[26 + j * 27] = -rm * rp * sm * sp * tm * tp;
        }

        return phi;
    }

    public static void basis_brick27_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_BRICK27_TEST verifies BASIS_BRICK27.
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
        //    None
        //
    {
        int j;
        const int node_num = 27;
        const int test_num = 5;

        Console.WriteLine("");
        Console.WriteLine("BASIS_BRICK27_TEST:");
        Console.WriteLine("  Verify basis functions for element BRICK27.");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + node_num + "");
        Console.WriteLine("");
        Console.WriteLine("  The basis function values at basis nodes");
        Console.WriteLine("  should form the identity matrix.");
        Console.WriteLine("");

        int n = node_num;

        double[] p = NodeBricks.nodes_brick27();

        double[] phi = basis_brick27(n, p);

        for (j = 0; j < n; j++)
        {
            string cout = "  ";
            int i;
            for (i = 0; i < node_num; i++)
            {
                cout += phi[i + j * node_num].ToString(CultureInfo.InvariantCulture).PadLeft(7);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  The basis function values at ANY point P");
        Console.WriteLine("  should sum to 1:");
        Console.WriteLine("");
        Console.WriteLine("    ------------P-------------     PHI_SUM");
        Console.WriteLine("");

        n = test_num;
        int seed = 123456789;

        p = UniformRNG.r8mat_uniform_01_new(3, n, ref seed);

        phi = basis_brick27(n, p);

        for (j = 0; j < n; j++)
        {
            double phi_sum = typeMethods.r8vec_sum(node_num, phi, aIndex: + j * node_num);
            Console.WriteLine("  " + p[0 + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + p[1 + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + p[2 + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + phi_sum.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }
    }
}