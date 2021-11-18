using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SubsetTestNS;

public static class UnsignTest
{
    public static void nim_sum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NIM_SUM_TEST tests NIM_SUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 32;

        int i;
        uint[] i1vec = new uint[N];
        uint[] i2vec = new uint[N];
        uint[] i3vec = new uint[N];
        int ihi = 1000;
        int ilo = 0;
        int ntest = 5;
        int seed;
        uint ui1;
        uint ui2;
        uint ui3;

        Console.WriteLine("");
        Console.WriteLine("NIM_SUM_TEST");
        Console.WriteLine("  NIM_SUM computes the Nim sum of two integers.");
        Console.WriteLine("");
        Console.WriteLine("    I    J    Nim(I+J)");
        Console.WriteLine("");

        seed = 123456789;

        for (i = 1; i <= ntest; i++)
        {
            ui1 = (uint)UniformRNG.i4_uniform_ab(ilo, ihi, ref seed);
            typeMethods.ui4_to_ubvec(ui1, N, ref i1vec);

            ui2 = (uint)UniformRNG.i4_uniform_ab(ilo, ihi, ref seed);
            typeMethods.ui4_to_ubvec(ui2, N, ref i2vec);

            ui3 = typeMethods.nim_sum(ui1, ui2);
            typeMethods.ui4_to_ubvec(ui3, N, ref i3vec);

            Console.WriteLine("");
            Console.WriteLine("  I1, I2, I3 in decimal:");
            Console.WriteLine("");

            Console.WriteLine("  "
                              + ui1.ToString(CultureInfo.InvariantCulture).PadLeft(5) + "");
            Console.WriteLine("  "
                              + ui2.ToString(CultureInfo.InvariantCulture).PadLeft(5) + "");
            Console.WriteLine("  "
                              + ui3.ToString(CultureInfo.InvariantCulture).PadLeft(5) + "");

            Console.WriteLine("");
            Console.WriteLine("  I1, I2, I3 in binary:");
            Console.WriteLine("");

            typeMethods.ubvec_print(N, i1vec, " ");
            typeMethods.ubvec_print(N, i2vec, " ");
            typeMethods.ubvec_print(N, i3vec, " ");
        }

    }

    public static void ubvec_add_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UBVEC_ADD_TEST tests UBVEC_ADD;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;

        uint[] bvec1 = new uint[N];
        uint[] bvec2 = new uint[N];
        uint[] bvec3 = new uint[N];
        uint ui;
        uint uj;
        uint uk;
        int seed = 123456789;
        int test;
        int test_num = 10;

        Console.WriteLine("");
        Console.WriteLine("UBVEC_ADD_TEST");
        Console.WriteLine("  UBVEC_ADD adds unsigned binary vectors ");
        Console.WriteLine("  representing uintegers;");
        Console.WriteLine("");
        Console.WriteLine("        I        J        K = I + J");
        Console.WriteLine("");

        for (test = 1; test <= test_num; test++)
        {
            ui = (uint ) UniformRNG.i4_uniform_ab(0, 100, ref seed);
            uj = (uint ) UniformRNG.i4_uniform_ab(0, 100, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  " + ui.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + uj.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

            uk = ui + uj;

            Console.WriteLine("  Directly:         "
                              + "  " + uk.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

            typeMethods.ui4_to_ubvec(ui, N, ref bvec1);
            typeMethods.ui4_to_ubvec(uj, N, ref bvec2);

            typeMethods.ubvec_add(N, bvec1, bvec2, ref bvec3);
            uk = typeMethods.ubvec_to_ui4(N, bvec3);

            Console.WriteLine(" UBVEC_ADD          "
                              + "  " + uk.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }
    }

    public static void ubvec_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UBVEC_PRINT_TEST tests UBVEC_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 10;
        uint[] ubvec =  {
                1, 0, 0, 1, 0, 1, 1, 1, 0, 0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("UBVEC_PRINT_TEST");
        Console.WriteLine("  UBVEC_PRINT prints an unsigned binary vector.");

        typeMethods.ubvec_print(n, ubvec, "  UBVEC:");
    }

    public static void ubvec_to_ui4_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UBVEC_TO_UI4_TEST tests UBVEC_TO_UI4;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;

        uint[] bvec = new uint[N];
        uint ui1;
        uint ui2;
        int j;

        Console.WriteLine("");
        Console.WriteLine("UBVEC_TO_UI4_TEST");
        Console.WriteLine("  UBVEC_TO_UI4 converts an unsigned binary vector");
        Console.WriteLine("  to an uinteger;");
        Console.WriteLine("");
        Console.WriteLine("  I --> BVEC  -->  I");
        Console.WriteLine("");

        for (ui1 = 0; ui1 <= 10; ui1++)
        {
            typeMethods.ui4_to_ubvec(ui1, N, ref bvec);

            ui2 = typeMethods.ubvec_to_ui4(N, bvec);

            string cout = ui1.ToString(CultureInfo.InvariantCulture).PadLeft(3) + "  ";
            for (j = 0; j < N; j++)
            {
                cout += bvec[j].ToString(CultureInfo.InvariantCulture).PadLeft(1);
            }

            cout += "  ";
            Console.WriteLine(cout + ui2.ToString(CultureInfo.InvariantCulture).PadLeft(3) + "");
        }
    }

    public static void ubvec_xor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UBVEC_XOR_TEST tests UBVEC_XOR;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        uint[] bvec1 = new uint[10];
        uint[] bvec2 = new uint[10];
        uint[] bvec3 = new uint[10];
        int n = 10;
        uint ui;
        uint uj;
        uint uk;
        int seed = 123456789;
        int test;
        int test_num = 10;

        Console.WriteLine("");
        Console.WriteLine("UBVEC_XOR_TEST");
        Console.WriteLine("  UBVEC_XOR exclusive-ors unsigned binary vectors ");
        Console.WriteLine("  representing uintegers;");
        Console.WriteLine("");
        Console.WriteLine("        I        J        K = I XOR J");
        Console.WriteLine("");

        for (test = 1; test <= test_num; test++)
        {
            ui = (uint ) UniformRNG.i4_uniform_ab(0, 100, ref seed);
            uj = (uint ) UniformRNG.i4_uniform_ab(0, 100, ref seed);
            typeMethods.ui4_to_ubvec(ui, n, ref bvec1);
            typeMethods.ui4_to_ubvec(uj, n, ref bvec2);
            typeMethods.ubvec_xor(n, bvec1, bvec2, bvec3);
            uk = typeMethods.ubvec_to_ui4(n, bvec3);

            Console.WriteLine("  " + ui.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + uj.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + uk.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }
    }

    public static void ui4_to_ubvec_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UI4_TO_UBVEC_TEST tests UI4_TO_UBVEC;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;

        uint[] bvec = new uint[N];
        uint ui1;
        uint ui2;
        int j;

        Console.WriteLine("");
        Console.WriteLine("UI4_TO_UBVEC_TEST");
        Console.WriteLine("  UI4_TO_UBVEC converts an uinteger to an ");
        Console.WriteLine("  unsigned binary vector;");
        Console.WriteLine("");
        Console.WriteLine("  I --> BVEC  -->  I");
        Console.WriteLine("");

        for (ui1 = 0; ui1 <= 10; ui1++)
        {
            typeMethods.ui4_to_ubvec(ui1, N, ref bvec);

            ui2 = typeMethods.ubvec_to_ui4(N, bvec);

            string cout = ui1.ToString(CultureInfo.InvariantCulture).PadLeft(3) + "  ";
            for (j = 0; j < N; j++)
            {
                cout += bvec[j].ToString(CultureInfo.InvariantCulture).PadLeft(1);
            }

            cout += "  ";
            Console.WriteLine(cout + ui2.ToString(CultureInfo.InvariantCulture).PadLeft(3) + "");
        }
    }

}