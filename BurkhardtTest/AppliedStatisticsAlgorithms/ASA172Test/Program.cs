using System;
using Burkardt.AppliedStatistics;

namespace ASA172Test;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA172_TEST.
        //
        //  Discussion:
        //
        //    ASA172_TEST tests the ASA172 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ASA172_TEST:");
        Console.WriteLine("  Test the ASA172 library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("ASA172_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }


    private static void test01()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 compares indices computed by a triple loop.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int KDIM = 3;

        int ifault;
        int[] iprod = new int[KDIM];
        int[] ivec= new int[KDIM];
        int jsub;
        int kdim = KDIM;
        int n;
        int[] nr =  {
                3, 2, 4
            }
            ;
        bool qfor;
        bool qind;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  SIMDO can convert between compressed and");
        Console.WriteLine("  vector indices representing a nested loop.");
        Console.WriteLine("");
        Console.WriteLine("  Here, we set QFOR = FALSE, meaning we do");
        Console.WriteLine("  NOT want to convert from FORTRAN ordering");
        Console.WriteLine("  to lexical ordering.");
        Console.WriteLine("");
        Console.WriteLine("  Here, we actually carry out a triple loop");
        Console.WriteLine("  list the indices, and then compare.");

        qfor = false;
        //
        //  If QFOR is FALSE, then the definition of IPROD is reversed...
        //
        iprod[0] = nr[kdim - 1];
        for (int i = 1; i < kdim; i++)
        {
            iprod[i] = iprod[i - 1] * nr[kdim - 1 - i];
        }

        n = iprod[kdim - 1];
        //
        //  Carry out the nested loops, and use JSUB to count each iteration.
        //  In the inmost loop, print JSUB and the corresponding (I1,I2,I3) vector.
        //
        jsub = 0;

        Console.WriteLine("");
        Console.WriteLine("  #1: Generate JSUB by counting as we DO the loops:");
        Console.WriteLine("");
        Console.WriteLine("  DO I1 = 1, N1");
        Console.WriteLine("    DO I2 = 1, N2");
        Console.WriteLine("      DO I3 = 1, N3");
        Console.WriteLine("");
        Console.WriteLine("      JSUB            I1        I2        I3");
        Console.WriteLine("");
        for (int i = 1; i <= nr[0]; i++)
        {
            ivec[0] = i;
            for (int j = 1; j <= nr[1]; j++)
            {
                ivec[1] = j;
                for (int k = 1; k <= nr[2]; k++)
                {
                    ivec[2] = k;
                    jsub += 1;
                    Console.WriteLine("  " + jsub.ToString().PadLeft(8) + "    "
                                      + "  " + i.ToString().PadLeft(8)
                                      + "  " + j.ToString().PadLeft(8)
                                      + "  " + k.ToString().PadLeft(8) + "");
                }
            }
        }

        //
        //  Now for each value of JSUB, retrieve the corresponding index subscript.
        //  In order to use the QFOR = .FALSE. switch, I apparently have to reverse
        //  the sense of the NR vector//
        //
        qind = true;

        Console.WriteLine("");
        Console.WriteLine("  #2: Loop on JSUB, retrieve loop indices");
        Console.WriteLine("      QIND = TRUE J ->I(J)");
        Console.WriteLine("      QFOR = FALSE");
        Console.WriteLine("");
        Console.WriteLine("      JSUB            I1        I2        I3");
        Console.WriteLine("");

        for (int j = 1; j <= n; j++)
        {
            jsub = j;
            ifault = Algorithms.simdo(qind, qfor, iprod, kdim, ref jsub, ref ivec);
            if (ifault != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("Fatal error!");
                Console.WriteLine("  simdo returns nonzero ifault.");
            }

            Console.WriteLine("  " + jsub.ToString().PadLeft(8) + "    "
                              + "  " + ivec[0].ToString().PadLeft(8)
                              + "  " + ivec[1].ToString().PadLeft(8)
                              + "  " + ivec[2].ToString().PadLeft(8) + "");
        }

        //
        //  Carry out the nested loops, and DO NOT compute JSUB.
        //  Have SIMDO determine JSUB.
        //
        qind = false;

        Console.WriteLine("");
        Console.WriteLine("  #3: For any set of loop indices, retrieve JSUB");
        Console.WriteLine("      QIND = FALSE I(J) -> J");
        Console.WriteLine("      QFOR = FALSE");
        Console.WriteLine("");
        Console.WriteLine("      JSUB            I1        I2        I3");
        Console.WriteLine("");
        for (int i = 1; i <= nr[0]; i++)
        {
            ivec[0] = i;
            for (int j = 1; j <= nr[1]; j++)
            {
                ivec[1] = j;
                for (int k = 1; k <= nr[2]; k++)
                {
                    ivec[2] = k;
                    ifault = Algorithms.simdo(qind, qfor, iprod, kdim, ref jsub, ref ivec);
                    Console.WriteLine("  " + jsub.ToString().PadLeft(8) + "    "
                                      + "  " + i.ToString().PadLeft(8)
                                      + "  " + j.ToString().PadLeft(8)
                                      + "  " + k.ToString().PadLeft(8) + "");
                }
            }
        }
    }

    private static void test02()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 compares indices computed by a triple loop.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int KDIM = 3;

        int ifault;
        int[] iprod = new int[KDIM];
        int[] ivec = new int[KDIM];
        int jsub;
        int kdim = KDIM;
        int n;
        int[] nr =  {
            3, 2, 4
        };
        bool qfor;
        bool qind;

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  SIMDO can convert between compressed and");
        Console.WriteLine("  vector indices representing a nested loop.");
        Console.WriteLine("");
        Console.WriteLine("  Here, we set QFOR = TRUE, meaning we DO");
        Console.WriteLine("  want to convert from the FORTRAN ");
        Console.WriteLine("  ordering to lexical convention.");
        Console.WriteLine("");
        Console.WriteLine("  Here, we actually carry out a triple loop");
        Console.WriteLine("  list the indices, and then compare.");

        qfor = true;

        iprod[0] = nr[0];
        for (int i = 1; i < kdim; i++)
        {
            iprod[i] = iprod[i - 1] * nr[i];
        }

        n = iprod[kdim - 1];
        //
        //  Carry out the nested loops, and use JSUB to count each iteration.
        //  In the inmost loop, print JSUB and the corresponding (I1,I2,I3) vector.
        //
        jsub = 0;

        Console.WriteLine("");
        Console.WriteLine("  #1: Generate JSUB by counting as we do the loops.");
        Console.WriteLine("");
        Console.WriteLine("  DO I3 = 1, N3");
        Console.WriteLine("    DO I2 = 1, N2");
        Console.WriteLine("      DO I1 = 1, N1");
        Console.WriteLine("");
        Console.WriteLine("      JSUB            I1        I2        I3");
        Console.WriteLine("");
        for (int i = 1; i <= nr[2]; i++)
        {
            ivec[2] = i;
            for (int j = 1; j <= nr[1]; j++)
            {
                ivec[1] = j;
                for (int k = 1; k <= nr[0]; k++)
                {
                    ivec[0] = k;
                    jsub += 1;
                    Console.WriteLine("  " + jsub.ToString().PadLeft(8) + "    "
                                      + "  " + k.ToString().PadLeft(8)
                                      + "  " + j.ToString().PadLeft(8)
                                      + "  " + i.ToString().PadLeft(8) + "");
                }
            }
        }

        //
        //  Reverse the order, so that the loop indices are generated in lexical order.
        //
        qind = true;

        Console.WriteLine("");
        Console.WriteLine("  #2: Setting QFOR false means loop indices");
        Console.WriteLine("  are generated in lexical order.");
        Console.WriteLine("      QIND = TRUE J -> I(J)");
        Console.WriteLine("      QFOR = TRUE");
        Console.WriteLine("");
        Console.WriteLine("      JSUB            I1        I2        I3");
        Console.WriteLine("");

        for (int j = 1; j <= n; j++)
        {
            jsub = j;
            ifault = Algorithms.simdo(qind, qfor, iprod, kdim, ref jsub, ref ivec);

            if (ifault != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("Fatal error!");
                Console.WriteLine("  simdo returns nonzero ifault.");
            }

            Console.WriteLine("  " + jsub.ToString().PadLeft(8) + "    "
                              + "  " + ivec[0].ToString().PadLeft(8)
                              + "  " + ivec[1].ToString().PadLeft(8)
                              + "  " + ivec[2].ToString().PadLeft(8) + "");
        }

        //
        //  Carry out the nested loops, and DO NOT compute JSUB.
        //  Have SIMDO determine JSUB.
        //
        qind = false;

        Console.WriteLine("");
        Console.WriteLine("  #3: For any set of loop indices, retrieve JSUB");
        Console.WriteLine("      QIND = FALSE I(J) -> J");
        Console.WriteLine("      QFOR = TRUE");
        Console.WriteLine("");
        Console.WriteLine("      JSUB            I1        I2        I3");
        Console.WriteLine("");
        for (int i = 1; i <= nr[2]; i++)
        {
            ivec[2] = i;
            for (int j = 1; j <= nr[1]; j++)
            {
                ivec[1] = j;
                for (int k = 1; k <= nr[0]; k++)
                {
                    ivec[0] = k;
                    ifault = Algorithms.simdo(qind, qfor, iprod, kdim, ref jsub, ref ivec);
                    Console.WriteLine("  " + jsub.ToString().PadLeft(8) + "    "
                                      + "  " + k.ToString().PadLeft(8)
                                      + "  " + j.ToString().PadLeft(8)
                                      + "  " + i.ToString().PadLeft(8) + "");
                }
            }
        }
    }

}