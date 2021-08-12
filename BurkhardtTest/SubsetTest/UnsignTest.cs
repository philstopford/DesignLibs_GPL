using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SubsetTest
{
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
                                  + ui1.ToString().PadLeft(5) + "");
                Console.WriteLine("  "
                                  + ui2.ToString().PadLeft(5) + "");
                Console.WriteLine("  "
                                  + ui3.ToString().PadLeft(5) + "");

                Console.WriteLine("");
                Console.WriteLine("  I1, I2, I3 in binary:");
                Console.WriteLine("");

                typeMethods.ubvec_print(N, i1vec, " ");
                typeMethods.ubvec_print(N, i2vec, " ");
                typeMethods.ubvec_print(N, i3vec, " ");
            }

        }
    }
}