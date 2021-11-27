using System;
using System.Globalization;
using Burkardt.FEM;
using Burkardt.Types;
using Burkardt.Uniform;

namespace FEM1DPackTest;

internal static class Program
{
    private static void Main()
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM1D_PACK_TEST.
//
//  Discussion:
//
//    FEM1D_PACK_TEST tests the FEM1D_PACK library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 March 2011
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("FEM1D_PACK_TEST");
        Console.WriteLine("  Test the FEM1D_PACK library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("FEM1D_PACK_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
            
    }

    private static void test01()
//****************************************************************************80
//
//  Purpose:
//
//    TEST01 verifies LOCAL_BASIS_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2011
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
        const int NODE_NUM = 4;

        double[] node_x =  {
                1.0, 2.0, 4.0, 4.5
            }
            ;
        double[] phi;
        double[] phi_matrix = new double[NODE_NUM * NODE_NUM];
        double x;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  LOCAL_BASIS_1D evaluates the local basis functions");
        Console.WriteLine("  for a 1D element.");
        Console.WriteLine("");
        Console.WriteLine("  Test that the basis functions, evaluated at the nodes,");
        Console.WriteLine("  form the identity matrix.");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + NODE_NUM + "");

        Console.WriteLine("");
        Console.WriteLine("  Node coordinates:");
        Console.WriteLine("");
        for (int j = 0; j < NODE_NUM; j++)
        {
            Console.WriteLine("  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + node_x[j].ToString(CultureInfo.InvariantCulture).PadLeft(7) + "");
        }

        for (int j = 0; j < NODE_NUM; j++)
        {
            x = node_x[j];
            phi = LocalBasis.local_basis_1d(NODE_NUM, node_x, x);
            for (int i = 0; i < NODE_NUM; i++)
            {
                phi_matrix[i + j * NODE_NUM] = phi[i];
            }
        }

        typeMethods.r8mat_print(NODE_NUM, NODE_NUM, phi_matrix, "  A(I,J) = PHI(I) at node (J):");

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("  The PHI functions should sum to 1 at random X values:");
        Console.WriteLine("");
        Console.WriteLine("       X        Sum ( PHI(:)(X) )");
        Console.WriteLine("");

        for (int j = 1; j <= 5; j++)
        {
            x = UniformRNG.r8_uniform(1.0, 4.5, ref seed);
            phi = LocalBasis.local_basis_1d(NODE_NUM, node_x, x);
            double s = typeMethods.r8vec_sum(NODE_NUM, phi);
            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + s.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }
}