using System;
using Burkardt.FEM;

namespace FEM3DPackTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FEM3D_PACK_TEST.
        //
        //  Discussion:
        //
        //    FEM3D_PACK_TEST tests the FEM3D_PACK library.
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
    {
        Console.WriteLine("");
        Console.WriteLine("FEM3D_PACK_TEST:");
        Console.WriteLine("  Test the FEM3D_PACK library.");

        Basis_mn.basis_mn_tet4_test();
        Basis_mn.basis_mn_tet10_test();
        BasisBrick.basis_brick8_test();
        BasisBrick.basis_brick20_test();
        BasisBrick.basis_brick27_test();

        test03();
        Console.WriteLine("");
        Console.WriteLine("FEM3D_PACK_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests PHYSICAL_TO_REFERENCE_TET4 and REFERENCE_TO_PHYSICAL_TET4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        int n = 10;
        double[] phy;
        double[] ref_;
        double[] ref2;
        int seed;
        double[] t =  {
                1.0, 2.0, 3.0,
                4.0, 1.0, 2.0,
                2.0, 4.0, 4.0,
                3.0, 2.0, 5.0
            }
            ;

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  For an order 4 tetrahedron,");
        Console.WriteLine("  PHYSICAL_TO_REFERENCE_TET4 maps a physical point to");
        Console.WriteLine("    a reference point.");
        Console.WriteLine("  REFERENCE_TO_PHYSICAL_TET4 maps a reference point to");
        Console.WriteLine("    a physical point.");
        Console.WriteLine("");
        Console.WriteLine("    ( R       S       T ) ==> ( X       Y       Z ) " + 
                          "==> ( R2      S2      T2 )");
        Console.WriteLine("");

        ref_ = Reference.reference_tet4_uniform(n, ref seed);

        phy = Reference.reference_to_physical_tet4(t, n, ref ref_);
        ref2 = PhysicalToRef.physical_to_reference_tet4(t, n, phy);

        for (j = 0; j < n; j++)
        {
            Console.WriteLine("  " + ref_[0 + j * 3].ToString().PadLeft(8)
                                   + "  " + ref_[1 + j * 3].ToString().PadLeft(8)
                                   + "  " + ref_[2 + j * 3].ToString().PadLeft(8)
                                   + "  " + phy[0 + j * 3].ToString().PadLeft(8)  
                                   + "  " + phy[1 + j * 3].ToString().PadLeft(8) 
                                   + "  " + phy[2 + j * 3].ToString().PadLeft(8) 
                                   + "  " + ref2[0 + j * 3].ToString().PadLeft(8) 
                                   + "  " + ref2[1 + j * 3].ToString().PadLeft(8)
                                   + "  " + ref2[2 + j * 3].ToString().PadLeft(8) + "");
        }
    }
}