using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.IO;
using Burkardt.Types;

namespace HBIOTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for HB_IO_TEST.
        //
        //  Discussion:
        //
        //    HB_IO_TEST tests the HB_IO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("HB_IO_TEST");
        Console.WriteLine("  Test the HB_IO library.");

        test01("rua_32.txt");
        test01("rse_5.txt");
        test02();
        test03("rua_32.txt");
        test03("rse_5.txt");
        test04();
        test05("rua_32.txt");
        test05("rse_5.txt");
        test06();
        test07("rua_32.txt");
        test08();
        test09("rua_32.txt");
        test10();
        test11();
        test12();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("HB_IO_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(string input_file)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests HB_HEADER_READ;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int indcrd = 0;
        char[] indfmt = null;
        string[] input;
        char[] key = null;
        char[] mxtype = null;
        int ncol = 0;
        int neltvl = 0;
        int nnzero = 0;
        int nrhs = 0;
        int nrhsix = 0;
        int nrow = 0;
        int ptrcrd = 0;
        char[] ptrfmt = null;
        int rhscrd = 0;
        char[] rhsfmt = null;
        char[] rhstyp = null;
        char[] title = null;
        int totcrd = 0;
        int valcrd = 0;
        char[] valfmt = null;
        int inputIndex = 0;

        Console.WriteLine(" ");
        Console.WriteLine("TEST01");
        Console.WriteLine("  HB_HEADER_READ reads the header of an HB file.");
        Console.WriteLine("");
        Console.WriteLine("  Reading the file '" + input_file + "'.");

        try
        {
            input = File.ReadAllLines(input_file);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("TEST01 - Warning!");
            Console.WriteLine("  Error opening the file.");
            return;
        }

        HB.hb_header_read(input, ref inputIndex, ref title, ref key, ref totcrd, ref ptrcrd, ref indcrd,
            ref valcrd, ref rhscrd, ref mxtype, ref nrow, ref ncol, ref nnzero, ref neltvl, ref ptrfmt,
            ref indfmt, ref valfmt, ref rhsfmt, ref rhstyp, ref nrhs, ref nrhsix);

        //
        //  Print out the  header information.
        //
        HB.hb_header_print(title, key, totcrd, ptrcrd, indcrd, valcrd,
            rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, valfmt,
            rhsfmt, rhstyp, nrhs, nrhsix);

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests HB_HEADER_WRITE;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int indcrd = 8;
        char[] indfmt = "(16I5)".ToCharArray();
        Array.Resize(ref indfmt, 17);

        char[] key = "RUA_32".ToCharArray();
        Array.Resize(ref key, 9);

        char[] mxtype = "PUA".ToCharArray();
        Array.Resize(ref mxtype, 4);

        const int ncol = 32;
        const int neltvl = 0;
        const int nnzero = 126;
        const int nrhs = 0;
        const int nrhsix = 0;
        const int nrow = 32;
        List<string> output = new();
        const string output_file = "rua_32_header.txt";
        const int ptrcrd = 3;
        char[] ptrfmt = new char[17]; // "(16I5)"
        ptrfmt[0] = '(';
        ptrfmt[1] = '1';
        ptrfmt[2] = '6';
        ptrfmt[3] = 'I';
        ptrfmt[4] = '5';
        ptrfmt[5] = ')';
        const int rhscrd = 0;
        char[] rhsfmt = " ".ToCharArray();
        Array.Resize(ref rhsfmt, 21);
        char[] rhstyp = "   ".ToCharArray();
        Array.Resize(ref rhstyp, 4);
        char[] title = "1Real unsymmetric assembled matrix based on IBM32".ToCharArray();
        Array.Resize(ref title, 73);
        const int totcrd = 11;
        const int valcrd = 0;
        char[] valfmt = " ".ToCharArray();
        Array.Resize(ref valfmt, 21);

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  HB_HEADER_WRITE writes the header of an HB file.");
        Console.WriteLine("");
        Console.WriteLine("  Writing the file '" + output_file + "'.");

        HB.hb_header_write(ref output, title, key, totcrd, ptrcrd, indcrd,
            valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt,
            valfmt, rhsfmt, rhstyp, nrhs, nrhsix);

        try
        {
            File.WriteAllLines(output_file, output);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("TEST02 - Warning!");
            Console.WriteLine("  Error opening the file.");
        }
    }

    private static void test03(string input_file)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests HB_STRUCTURE_READ;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int indcrd = 0;
        char[] indfmt = null;
        string[] input;
        char[] key = null;
        char[] mxtype = null;
        int ncol = 0;
        int neltvl = 0;
        int nnzero = 0;
        int nrhs = 0;
        int nrhsix = 0;
        int nrow = 0;
        int ptrcrd = 0;
        char[] ptrfmt = null;
        int rhscrd = 0;
        char[] rhsfmt = null;
        char[] rhstyp = null;
        int[] rowind;
        char[] title = null;
        int totcrd = 0;
        int valcrd = 0;
        char[] valfmt = null;
        int inputIndex = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  HB_STRUCTURE_READ reads the structure of an HB file.");
        Console.WriteLine("");
        Console.WriteLine("  Reading the file '" + input_file + "'.");

        try
        {
            input = File.ReadAllLines(input_file);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("TEST03 - Warning!");
            Console.WriteLine("  Error opening the file.");
            return;
        }

        Console.WriteLine("  Reading the header.");

        HB.hb_header_read(input, ref inputIndex, ref title, ref key, ref totcrd, ref ptrcrd, ref indcrd,
            ref valcrd, ref rhscrd, ref mxtype, ref nrow, ref ncol, ref nnzero, ref neltvl, ref ptrfmt,
            ref indfmt, ref valfmt, ref rhsfmt, ref rhstyp, ref nrhs, ref nrhsix);

        int[] colptr = new int[ncol + 1];

        switch (mxtype[2])
        {
            case 'A':
                rowind = new int[nnzero];
                break;
            case 'E':
                rowind = new int[neltvl];
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("TEST03 - Warning!");
                Console.WriteLine("  Illegal value of MXTYPE character 3 = " + mxtype[2] + "");
                return;
        }

        Console.WriteLine("  Reading the structure.");

        HB.hb_structure_read(input, ref inputIndex, ncol, mxtype, nnzero, neltvl,
            ptrcrd, ptrfmt, indcrd, indfmt, ref colptr, ref rowind);

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("  KEY =    '" + key + "'.");
        Console.WriteLine("");
        Console.WriteLine("  NROW =   " + nrow + "");
        Console.WriteLine("  NCOL =   " + ncol + "");
        Console.WriteLine("  NNZERO = " + nnzero + "");
        Console.WriteLine("  NELTVL = " + neltvl + "");

        HB.hb_structure_print(ncol, mxtype, nnzero, neltvl, colptr, rowind);

    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests HB_STRUCTURE_WRITE;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int NCOL = 32;
        const int NELTVL = 0;
        const int NNZERO = 126;

        int[] colptr =
        {
            1, 7, 12, 18, 22, 26, 29, 34, 39, 46,
            53, 58, 61, 63, 65, 68, 71, 74, 79, 82,
            85, 88, 90, 94, 97, 102, 106, 110, 112, 117,
            121, 124, 127
        };
        Array.Resize(ref colptr, NCOL + 1);
        char[] indfmt = "(16I5)".ToCharArray();
        Array.Resize(ref indfmt, 17);

        char[] mxtype = "RUA".ToCharArray();
        Array.Resize(ref mxtype, 4);

        List<string> output = new();
        string output_file = "rua_32_structure.txt";

        char[] ptrfmt = "(16I5)".ToCharArray();
        Array.Resize(ref ptrfmt, 17);

        int[] rowind =
        {
            1, 2, 3, 4, 7, 26, 1, 2, 9, 21,
            28, 2, 3, 6, 8, 9, 29, 3, 4, 5,
            12, 3, 5, 23, 27, 1, 6, 16, 3, 7,
            14, 21, 31, 1, 8, 12, 17, 27, 7, 9,
            10, 13, 19, 23, 27, 1, 10, 11, 21, 23,
            25, 27, 2, 11, 15, 18, 29, 6, 12, 24,
            11, 13, 3, 14, 2, 15, 20, 4, 16, 22,
            4, 16, 17, 6, 10, 18, 20, 30, 1, 19,
            26, 8, 16, 20, 3, 21, 32, 11, 22, 2,
            17, 21, 23, 12, 24, 26, 6, 15, 18, 24,
            25, 13, 18, 22, 26, 5, 24, 26, 27, 9,
            28, 3, 5, 27, 29, 32, 12, 17, 23, 30,
            13, 14, 31, 24, 28, 32
        };
        Array.Resize(ref rowind, NNZERO);

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  HB_STRUCTURE_WRITE writes the structure of an HB file.");
        Console.WriteLine("");
        Console.WriteLine("  Writing the file '" + output_file + "'.");

        HB.hb_structure_write(ref output, NCOL, mxtype, NNZERO, NELTVL,
            ptrfmt, indfmt, colptr, rowind);

        try
        {
            File.WriteAllLines(output_file, output);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("TEST04 - Warning!");
            Console.WriteLine("  Error opening the file.");
        }

    }

    private static void test05(string input_file)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests HB_VALUES_READ;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int indcrd = 0;
        char[] indfmt = null;
        string[] input;
        char[] key = null;
        char[] mxtype = null;
        int ncol = 0;
        int neltvl = 0;
        int nnzero = 0;
        int nrhs = 0;
        int nrhsix = 0;
        int nrow = 0;
        int ptrcrd = 0;
        char[] ptrfmt = null;
        int rhscrd = 0;
        char[] rhsfmt = null;
        char[] rhstyp = null;
        int[] rowind;
        char[] title = null;
        int totcrd = 0;
        int valcrd = 0;
        char[] valfmt = null;
        double[] values;
        int inputIndex = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  HB_VALUES_READ reads the values of an HB file.");
        Console.WriteLine("");
        Console.WriteLine("  Reading the file '" + input_file + "'.");

        try
        {
            input = File.ReadAllLines(input_file);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("TEST05 - Warning!");
            Console.WriteLine("  Error opening the file.");
            return;
        }

        Console.WriteLine("  Reading the header.");

        HB.hb_header_read(input, ref inputIndex, ref title, ref key, ref totcrd, ref ptrcrd, ref indcrd,
            ref valcrd, ref rhscrd, ref mxtype, ref nrow, ref ncol, ref nnzero, ref neltvl, ref ptrfmt,
            ref indfmt, ref valfmt, ref rhsfmt, ref rhstyp, ref nrhs, ref nrhsix);

        int[] colptr = new int[ncol + 1];

        switch (mxtype[2])
        {
            case 'A':
                rowind = new int[nnzero];
                break;
            case 'E':
                rowind = new int[neltvl];
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("TEST05 - Warning!");
                Console.WriteLine("  Illegal value of MXTYPE character 3.");
                return;
        }

        Console.WriteLine("  Reading the structure.");

        HB.hb_structure_read(input, ref inputIndex, ncol, mxtype, nnzero, neltvl,
            ptrcrd, ptrfmt, indcrd, indfmt, ref colptr, ref rowind);

        switch (mxtype[2])
        {
            case 'A':
                values = new double[nnzero];
                break;
            case 'E':
                values = new double[neltvl];
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("TEST05 - Warning!");
                Console.WriteLine("  Illegal value of MXTYPE character 3 = " + mxtype[2] + "");
                return;
        }

        Console.WriteLine("  Reading the values.");

        HB.hb_values_read(input, ref inputIndex, valcrd, mxtype, nnzero, neltvl, valfmt, ref values);

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("  KEY =    '" + key + "'.");
        Console.WriteLine("");
        Console.WriteLine("  NROW =   " + nrow + "");
        Console.WriteLine("  NCOL =   " + ncol + "");
        Console.WriteLine("  NNZERO = " + nnzero + "");
        Console.WriteLine("  NELTVL = " + neltvl + "");

        HB.hb_values_print(ncol, colptr, mxtype, nnzero, neltvl, values);
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests HB_VALUES_WRITE;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int NELTVL = 0;
        const int NNZERO = 126;

        char[] mxtype = "RUA".ToCharArray();
        Array.Resize(ref mxtype, 4);
        List<string> output = new();
        const string output_file = "rua_32_values.txt";
        const int valcrd = 13;
        char[] valfmt = "(10F7.1)".ToCharArray();
        Array.Resize(ref valfmt, 21);
        double[] values =
        {
            101.0, 102.0, 103.0, 104.0, 107.0,
            126.0, 201.0, 202.0, 209.0, 221.0,
            228.0, 302.0, 303.0, 306.0, 308.0,
            309.0, 329.0, 403.0, 404.0, 405.0,
            412.0, 503.0, 505.0, 523.0, 527.0,
            601.0, 606.0, 616.0, 703.0, 707.0,
            714.0, 721.0, 731.0, 801.0, 808.0,
            812.0, 817.0, 827.0, 907.0, 909.0,
            910.0, 913.0, 919.0, 923.0, 927.0,
            1001.0, 1010.0, 1011.0, 1021.0, 1023.0,
            1025.0, 1027.0, 1102.0, 1111.0, 1115.0,
            1118.0, 1129.0, 1206.0, 1212.0, 1224.0,
            1311.0, 1313.0, 1403.0, 1414.0, 1502.0,
            1515.0, 1520.0, 1604.0, 1616.0, 1622.0,
            1704.0, 1716.0, 1717.0, 1806.0, 1810.0,
            1818.0, 1820.0, 1830.0, 1901.0, 1919.0,
            1926.0, 2008.0, 2016.0, 2020.0, 2103.0,
            2121.0, 2132.0, 2211.0, 2222.0, 2302.0,
            2317.0, 2321.0, 2323.0, 2412.0, 2424.0,
            2426.0, 2506.0, 2515.0, 2518.0, 2524.0,
            2525.0, 2613.0, 2618.0, 2622.0, 2626.0,
            2705.0, 2724.0, 2726.0, 2727.0, 2809.0,
            2828.0, 2903.0, 2905.0, 2927.0, 2929.0,
            2932.0, 3012.0, 3017.0, 3023.0, 3030.0,
            3113.0, 3114.0, 3131.0, 3224.0, 3228.0,
            3232.0
        };
        Array.Resize(ref values, NNZERO);

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  HB_VALUES_WRITE writes the values of an HB file.");
        Console.WriteLine("");
        Console.WriteLine("  Writing the file '" + output_file + "'.");


        HB.hb_values_write(ref output, valcrd, mxtype, NNZERO, NELTVL, valfmt, values);

        try
        {
            File.WriteAllLines(output_file, output);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("TEST06 - Warning!");
            Console.WriteLine("  Error opening the file.");
        }
    }

    private static void test07(string input_file)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests HB_RHS_READ, HB_GUESS_READ, HB_EXACT_READ;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int indcrd = 0;
        char[] indfmt = null;
        string[] input;
        char[] key = null;
        char[] mxtype = null;
        int ncol = 0;
        int neltvl = 0;
        int nnzero = 0;
        int nrhs = 0;
        int nrhsix = 0;
        int nrow = 0;
        int ptrcrd = 0;
        char[] ptrfmt = null;
        int rhscrd = 0;
        char[] rhsfmt = null;
        int[] rhsind = null;
        int[] rhsptr = null;
        char[] rhstyp = null;
        double[] rhsval = null;
        double[] rhsvec = null;
        int[] rowind;
        char[] title = null;
        int totcrd = 0;
        int valcrd = 0;
        char[] valfmt = null;
        double[] values;
        int inputIndex = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  HB_RHS_READ reads right hand sides from an HB file.");
        Console.WriteLine("  HB_GUESS_READ reads starting guesses from an HB file.");
        Console.WriteLine("  HB_EXACT_READ reads exact solutions from an HB file.");

        Console.WriteLine("");
        Console.WriteLine("  Reading the file '" + input_file + "'.");

        try
        {
            input = File.ReadAllLines(input_file);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("TEST07 - Warning!");
            Console.WriteLine("  Error opening the file.");
            return;
        }

        Console.WriteLine("  Reading the header.");

        HB.hb_header_read(input, ref inputIndex, ref title, ref key, ref totcrd, ref ptrcrd, ref indcrd,
            ref valcrd, ref rhscrd, ref mxtype, ref nrow, ref ncol, ref nnzero, ref neltvl, ref ptrfmt,
            ref indfmt, ref valfmt, ref rhsfmt, ref rhstyp, ref nrhs, ref nrhsix);

        int[] colptr = new int[ncol + 1];

        switch (mxtype[2])
        {
            case 'A':
                rowind = new int[nnzero];
                break;
            case 'E':
                rowind = new int[neltvl];
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("TEST07 - Warning!");
                Console.WriteLine("  Illegal value of MXTYPE character 3 = " + mxtype[2] + "");
                return;
        }

        Console.WriteLine("  Reading the structure.");

        HB.hb_structure_read(input, ref inputIndex, ncol, mxtype, nnzero, neltvl,
            ptrcrd, ptrfmt, indcrd, indfmt, ref colptr, ref rowind);

        switch (mxtype[2])
        {
            case 'A':
                values = new double[nnzero];
                break;
            case 'E':
                values = new double[neltvl];
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("TEST07 - Warning!");
                Console.WriteLine("  Illegal value of MXTYPE character 3 = " + mxtype[2] + "");
                return;
        }

        Console.WriteLine("  Reading the values.");

        HB.hb_values_read(input, ref inputIndex, valcrd, mxtype, nnzero, neltvl, valfmt, ref values);
        switch (rhscrd)
        {
            //
            //  Read the right hand sides.
            //
            case > 0:
            {
                Console.WriteLine("  Reading the right hand side.");

                switch (rhstyp[0])
                {
                    case 'F':
                        rhsval = new double[nrow * nrhs];
                        break;
                    case 'M' when mxtype[2] == 'A':
                        rhsptr = new int[nrhs + 1];
                        rhsind = new int [nrhsix];
                        rhsvec = new double[nrhsix];
                        break;
                    case 'M':
                    {
                        rhsval = mxtype[2] switch
                        {
                            'E' => new double [nnzero * nrhs],
                            _ => rhsval
                        };

                        break;
                    }
                }

                HB.hb_rhs_read(input, ref inputIndex, nrow, nnzero, nrhs, nrhsix,
                    rhscrd, ptrfmt, indfmt, rhsfmt, mxtype, rhstyp, ref rhsval,
                    ref rhsind, ref rhsptr, ref rhsvec);

                Console.WriteLine("  Done reading the right hand side.");

                switch (rhstyp[0])
                {
                    case 'F':
                        typeMethods.r8mat_print_some(nrow, nrhs, rhsval, 1, 1, 5, 5, "  Part of RHS");
                        break;
                    case 'M' when mxtype[2] == 'A':
                        typeMethods.i4vec_print_part(nrhs + 1, rhsptr, 10, "  Part of RHSPTR");
                        typeMethods.i4vec_print_part(nrhsix, rhsind, 10, "  Part of RHSIND");
                        typeMethods.r8vec_print_part(nrhsix, rhsvec, 10, "  Part of RHSVEC");
                        break;
                    case 'M' when mxtype[2] == 'E':
                        typeMethods.r8mat_print_some(nnzero, nrhs, rhsval, 1, 1, 5, 5, "  Part of RHS");
                        break;
                }

                switch (rhstyp[1])
                {
                    //
                    //  Read the starting guesses.
                    //
                    case 'G':
                        Console.WriteLine("  Reading the starting guesses.");

                        double[] guess = new double[nrow * nrhs];

                        HB.hb_guess_read(input, ref inputIndex, nrow, nrhs, rhscrd, rhsfmt, rhstyp, ref guess);

                        typeMethods.r8mat_print_some(nrow, nrhs, guess, 1, 1, 5, 5, "  Part of GUESS");
                        break;
                }

                switch (rhstyp[2])
                {
                    //
                    //  Read the exact solutions.
                    //
                    case 'X':
                        Console.WriteLine("  Reading the exact solutions.");

                        double[] exact = new double[nrow * nrhs];

                        HB.hb_exact_read(input, ref inputIndex, nrow, nrhs, rhscrd, rhsfmt, rhstyp, ref exact);

                        typeMethods.r8mat_print_some(nrow, nrhs, exact, 1, 1, 5, 5, "  Part of EXACT");
                        break;
                }

                break;
            }
        }

    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests HB_RHS_WRITE, HB_GUESS_WRITE, HB_EXACT_WRITE;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int NNZERO = 126;
        const int NRHS = 1;
        const int NRHSIX = 0;
        const int NROW = 32;

        double[] exact =
        {
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
            11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
            21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0,
            31.0, 32.0
        };
        Array.Resize(ref exact, NROW);

        double[] guess =
        {
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0
        };
        Array.Resize(ref guess, NROW * NRHS);

        char[] indfmt = "(16I5)".ToCharArray();
        Array.Resize(ref indfmt, 17);

        char[] mxtype = "RUA".ToCharArray();
        Array.Resize(ref mxtype, 4);

        List<string> output = new();
        const string output_file = "rua_32_rhs.txt";

        char[] ptrfmt = "(16I5)".ToCharArray();
        Array.Resize(ref ptrfmt, 17);

        const int rhscrd = 12;
        char[] rhsfmt = "(10F7.1)".ToCharArray();
        Array.Resize(ref rhsfmt, 21);

        int[] rhsind = null;
        int[] rhsptr = null;

        double[] rhsval =
        {
            101.0, 102.0, 103.0, 104.0, 107.0, 126.0, 201.0, 202.0, 209.0, 221.0,
            228.0, 302.0, 303.0, 306.0, 308.0, 309.0, 329.0, 403.0, 404.0, 405.0,
            412.0, 503.0, 505.0, 523.0, 527.0, 601.0, 606.0, 616.0, 703.0, 707.0,
            714.0, 721.0
        };
        Array.Resize(ref rhsval, NROW * NRHS);

        double[] rhsvec = null;

        char[] rhstyp = "FGX".ToCharArray();
        Array.Resize(ref rhstyp, 4);

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  HB_RHS_WRITE writes the right hand sides to an HB file.");
        Console.WriteLine("  HB_GUESS_WRITE writes starting guesses to an HB file.");
        Console.WriteLine("  HB_EXACT_WRITE writes exact solutions to an HB file.");
        Console.WriteLine("");
        Console.WriteLine("  Writing the file '" + output_file + "'.");

        //
        //  Write the right hand sides.
        //
        HB.hb_rhs_write(ref output, NROW, NNZERO, NRHS, NRHSIX,
            rhscrd, ptrfmt, indfmt, rhsfmt, mxtype, rhstyp, rhsval,
            rhsind, rhsptr, rhsvec);
        //
        //  Write the right hand sides.
        //
        HB.hb_guess_write(ref output, NROW, NRHS, rhscrd, rhsfmt, rhstyp, guess);
        //
        //  Write the right hand sides.
        //
        HB.hb_exact_write(ref output, NROW, NRHS, rhscrd, rhsfmt, rhstyp, exact);

        try
        {
            File.WriteAllLines(output_file, output);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("TEST08 - Warning!");
            Console.WriteLine("  Error opening the file.");
        }
    }

    private static void test09(string input_file)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests HB_FILE_READ;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] colptr = null;
        double[] exact = null;
        double[] guess = null;
        int indcrd = 0;
        char[] indfmt = null;
        char[] key = null;
        char[] mxtype = null;
        int ncol = 0;
        int neltvl = 0;
        int nnzero = 0;
        int nrhs = 0;
        int nrhsix = 0;
        int nrow = 0;
        int ptrcrd = 0;
        char[] ptrfmt = null;
        int rhscrd = 0;
        char[] rhsfmt = null;
        int[] rhsind = null;
        int[] rhsptr = null;
        char[] rhstyp = null;
        double[] rhsval = null;
        double[] rhsvec = null;
        int[] rowind = null;
        char[] title = null;
        int totcrd = 0;
        int valcrd = 0;
        char[] valfmt = null;
        double[] values = null;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  HB_FILE_READ reads all the data in an HB file.");
        Console.WriteLine("  HB_FILE_MODULE is the module that stores the data.");

        Console.WriteLine("");
        Console.WriteLine("  Reading the file '" + input_file + "'.");

        HB.hb_file_read(input_file, ref title, ref key, ref totcrd, ref ptrcrd, ref indcrd,
            ref valcrd, ref rhscrd, ref mxtype, ref nrow, ref ncol, ref nnzero, ref neltvl,
            ref ptrfmt, ref indfmt, ref valfmt, ref rhsfmt, ref rhstyp, ref nrhs, ref nrhsix,
            ref colptr, ref rowind, ref values, ref rhsval, ref rhsptr, ref rhsind, ref rhsvec,
            ref guess, ref exact);
        //
        //  Print out the header information.
        //
        HB.hb_header_print(title, key, totcrd, ptrcrd, indcrd, valcrd,
            rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, valfmt,
            rhsfmt, rhstyp, nrhs, nrhsix);
        //
        //  Print the structure information.
        //
        HB.hb_structure_print(ncol, mxtype, nnzero, neltvl, colptr, rowind);
        //
        //  Print the values.
        //
        HB.hb_values_print(ncol, colptr, mxtype, nnzero, neltvl, values);

        switch (rhscrd)
        {
            case > 0:
            {
                switch (rhstyp[0])
                {
                    //
                    //  Print a bit of the right hand sides.
                    //
                    case 'F':
                        typeMethods.r8mat_print_some(nrow, nrhs, rhsval, 1, 1, 5, 5, "  Part of RHS");
                        break;
                    case 'M' when mxtype[2] == 'A':
                        typeMethods.i4vec_print_part(nrhs + 1, rhsptr, 10, "  Part of RHSPTR");
                        typeMethods.i4vec_print_part(nrhsix, rhsind, 10, "  Part of RHSIND");
                        typeMethods.r8vec_print_part(nrhsix, rhsvec, 10, "  Part of RHSVEC");
                        break;
                    case 'M' when mxtype[2] == 'E':
                        typeMethods.r8mat_print_some(nnzero, nrhs, rhsval, 1, 1, 5, 5, "  Part of RHS");
                        break;
                }

                switch (rhstyp[1])
                {
                    //
                    //  Print a bit of the starting guesses.
                    //
                    case 'G':
                        typeMethods.r8mat_print_some(nrow, nrhs, guess, 1, 1, 5, 5, "  Part of GUESS");
                        break;
                }

                switch (rhstyp[2])
                {
                    //
                    //  Print a bit of the exact solutions.
                    //
                    case 'X':
                        typeMethods.r8mat_print_some(nrow, nrhs, exact, 1, 1, 5, 5, "  Part of EXACT");
                        break;
                }

                break;
            }
        }

    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests HB_FILE_WRITE;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int NCOL = 32;
        const int NELTVL = 0;
        const int NNZERO = 126;
        const int NRHS = 1;
        const int NRHSIX = 0;
        const int NROW = 32;

        int[] colptr =
        {
            1, 7, 12, 18, 22, 26, 29, 34, 39, 46,
            53, 58, 61, 63, 65, 68, 71, 74, 79, 82,
            85, 88, 90, 94, 97, 102, 106, 110, 112, 117,
            121, 124, 127
        };
        Array.Resize(ref colptr, NCOL + 1);

        double[] exact =
        {
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
            11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
            21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0,
            31.0, 32.0
        };
        Array.Resize(ref exact, NROW * NRHS);

        double[] guess =
        {
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0
        };
        Array.Resize(ref guess, NROW * NRHS);

        const int indcrd = 8;
        char[] indfmt = "(16I5)".ToCharArray();
        Array.Resize(ref indfmt, 17);

        char[] key = "RUA_32".ToCharArray();
        Array.Resize(ref key, 9);

        char[] mxtype = "RUA".ToCharArray();
        Array.Resize(ref mxtype, 4);
        const string output_file = "rua_32_file.txt";
        List<string> output = new();
        const int ptrcrd = 3;

        char[] ptrfmt = "(16I5)".ToCharArray();
        Array.Resize(ref ptrfmt, 17);

        const int rhscrd = 12;

        char[] rhsfmt = "(10F7.1)".ToCharArray();
        Array.Resize(ref rhsfmt, 21);

        int[] rhsind = null;
        int[] rhsptr = null;

        double[] rhsval =
        {
            101.0, 102.0, 103.0, 104.0, 107.0, 126.0, 201.0, 202.0, 209.0, 221.0,
            228.0, 302.0, 303.0, 306.0, 308.0, 309.0, 329.0, 403.0, 404.0, 405.0,
            412.0, 503.0, 505.0, 523.0, 527.0, 601.0, 606.0, 616.0, 703.0, 707.0,
            714.0, 721.0
        };
        Array.Resize(ref rhsval, NROW * NRHS);

        char[] rhstyp = "FGX".ToCharArray();
        Array.Resize(ref rhstyp, 4);

        double[] rhsvec = null;
        int[] rowind =
        {
            1, 2, 3, 4, 7, 26, 1, 2, 9, 21,
            28, 2, 3, 6, 8, 9, 29, 3, 4, 5,
            12, 3, 5, 23, 27, 1, 6, 16, 3, 7,
            14, 21, 31, 1, 8, 12, 17, 27, 7, 9,
            10, 13, 19, 23, 27, 1, 10, 11, 21, 23,
            25, 27, 2, 11, 15, 18, 29, 6, 12, 24,
            11, 13, 3, 14, 2, 15, 20, 4, 16, 22,
            4, 16, 17, 6, 10, 18, 20, 30, 1, 19,
            26, 8, 16, 20, 3, 21, 32, 11, 22, 2,
            17, 21, 23, 12, 24, 26, 6, 15, 18, 24,
            25, 13, 18, 22, 26, 5, 24, 26, 27, 9,
            28, 3, 5, 27, 29, 32, 12, 17, 23, 30,
            13, 14, 31, 24, 28, 32
        };
        Array.Resize(ref rowind, NNZERO);

        char[] title = "1Real unsymmetric assembled matrix based on IBM32".ToCharArray();
        Array.Resize(ref title, 73);

        const int totcrd = 36;
        const int valcrd = 13;

        char[] valfmt = "(10F7.1)".ToCharArray();
        Array.Resize(ref valfmt, 21);

        double[] values =
        {
            101.0, 102.0, 103.0, 104.0, 107.0,
            126.0, 201.0, 202.0, 209.0, 221.0,
            228.0, 302.0, 303.0, 306.0, 308.0,
            309.0, 329.0, 403.0, 404.0, 405.0,
            412.0, 503.0, 505.0, 523.0, 527.0,
            601.0, 606.0, 616.0, 703.0, 707.0,
            714.0, 721.0, 731.0, 801.0, 808.0,
            812.0, 817.0, 827.0, 907.0, 909.0,
            910.0, 913.0, 919.0, 923.0, 927.0,
            1001.0, 1010.0, 1011.0, 1021.0, 1023.0,
            1025.0, 1027.0, 1102.0, 1111.0, 1115.0,
            1118.0, 1129.0, 1206.0, 1212.0, 1224.0,
            1311.0, 1313.0, 1403.0, 1414.0, 1502.0,
            1515.0, 1520.0, 1604.0, 1616.0, 1622.0,
            1704.0, 1716.0, 1717.0, 1806.0, 1810.0,
            1818.0, 1820.0, 1830.0, 1901.0, 1919.0,
            1926.0, 2008.0, 2016.0, 2020.0, 2103.0,
            2121.0, 2132.0, 2211.0, 2222.0, 2302.0,
            2317.0, 2321.0, 2323.0, 2412.0, 2424.0,
            2426.0, 2506.0, 2515.0, 2518.0, 2524.0,
            2525.0, 2613.0, 2618.0, 2622.0, 2626.0,
            2705.0, 2724.0, 2726.0, 2727.0, 2809.0,
            2828.0, 2903.0, 2905.0, 2927.0, 2929.0,
            2932.0, 3012.0, 3017.0, 3023.0, 3030.0,
            3113.0, 3114.0, 3131.0, 3224.0, 3228.0,
            3232.0
        };
        Array.Resize(ref values, NNZERO);

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  HB_FILE_WRITE writes an HB file.");
        Console.WriteLine("");
        Console.WriteLine("  Writing the file '" + output_file + "'.");

        HB.hb_file_write(output_file, title, key, totcrd, ptrcrd, indcrd,
            valcrd, rhscrd, mxtype, NROW, NCOL, NNZERO, NELTVL, ptrfmt, indfmt,
            valfmt, rhsfmt, rhstyp, NRHS, NRHSIX, colptr, rowind, values,
            rhsval, rhsptr, rhsind, rhsvec, guess, exact);

    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests HB_FILE_WRITE;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int NCOL = 32;
        const int NELTVL = 0;
        const int NNZERO = 126;
        const int NRHS = 2;
        const int NRHSIX = 0;
        const int NROW = 32;

        int[] colptr =
        {
            1, 7, 12, 18, 22, 26, 29, 34, 39, 46,
            53, 58, 61, 63, 65, 68, 71, 74, 79, 82,
            85, 88, 90, 94, 97, 102, 106, 110, 112, 117,
            121, 124, 127
        };
        Array.Resize(ref colptr, NCOL + 1);

        double[] exact =
        {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0
        };
        Array.Resize(ref exact, NROW * NRHS);

        double[] guess =
        {
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0
        };
        Array.Resize(ref guess, NROW * NRHS);

        const int indcrd = 8;
        char[] indfmt = "(16I5)".ToCharArray();
        Array.Resize(ref indfmt, 17);

        char[] key = "RUA_32".ToCharArray();
        Array.Resize(ref key, 9);

        char[] mxtype = "RUA".ToCharArray();
        Array.Resize(ref mxtype, 4);

        const string output_file = "rua_32_ax.txt";
        List<string> output = new();

        const int ptrcrd = 3;
        char[] ptrfmt = "(16I5)".ToCharArray();
        Array.Resize(ref ptrfmt, 17);

        const int rhscrd = 12;
        char[] rhsfmt = "(10F7.1)".ToCharArray();
        Array.Resize(ref rhsfmt, 21);

        int[] rhsind = null;
        int[] rhsptr = null;

        char[] rhstyp = "FGX".ToCharArray();
        Array.Resize(ref rhstyp, 4);

        double[] rhsvec = null;

        int[] rowind =
        {
            1, 2, 3, 4, 7, 26, 1, 2, 9, 21,
            28, 2, 3, 6, 8, 9, 29, 3, 4, 5,
            12, 3, 5, 23, 27, 1, 6, 16, 3, 7,
            14, 21, 31, 1, 8, 12, 17, 27, 7, 9,
            10, 13, 19, 23, 27, 1, 10, 11, 21, 23,
            25, 27, 2, 11, 15, 18, 29, 6, 12, 24,
            11, 13, 3, 14, 2, 15, 20, 4, 16, 22,
            4, 16, 17, 6, 10, 18, 20, 30, 1, 19,
            26, 8, 16, 20, 3, 21, 32, 11, 22, 2,
            17, 21, 23, 12, 24, 26, 6, 15, 18, 24,
            25, 13, 18, 22, 26, 5, 24, 26, 27, 9,
            28, 3, 5, 27, 29, 32, 12, 17, 23, 30,
            13, 14, 31, 24, 28, 32
        };
        Array.Resize(ref rowind, NNZERO);

        char[] title = "1Real unsymmetric assembled matrix based on IBM32".ToCharArray();
        Array.Resize(ref title, 73);

        const int totcrd = 36;
        const int valcrd = 13;
        char[] valfmt = "(10F7.1)".ToCharArray();
        Array.Resize(ref valfmt, 21);

        double[] values =
        {
            101.0, 102.0, 103.0, 104.0, 107.0,
            126.0, 201.0, 202.0, 209.0, 221.0,
            228.0, 302.0, 303.0, 306.0, 308.0,
            309.0, 329.0, 403.0, 404.0, 405.0,
            412.0, 503.0, 505.0, 523.0, 527.0,
            601.0, 606.0, 616.0, 703.0, 707.0,
            714.0, 721.0, 731.0, 801.0, 808.0,
            812.0, 817.0, 827.0, 907.0, 909.0,
            910.0, 913.0, 919.0, 923.0, 927.0,
            1001.0, 1010.0, 1011.0, 1021.0, 1023.0,
            1025.0, 1027.0, 1102.0, 1111.0, 1115.0,
            1118.0, 1129.0, 1206.0, 1212.0, 1224.0,
            1311.0, 1313.0, 1403.0, 1414.0, 1502.0,
            1515.0, 1520.0, 1604.0, 1616.0, 1622.0,
            1704.0, 1716.0, 1717.0, 1806.0, 1810.0,
            1818.0, 1820.0, 1830.0, 1901.0, 1919.0,
            1926.0, 2008.0, 2016.0, 2020.0, 2103.0,
            2121.0, 2132.0, 2211.0, 2222.0, 2302.0,
            2317.0, 2321.0, 2323.0, 2412.0, 2424.0,
            2426.0, 2506.0, 2515.0, 2518.0, 2524.0,
            2525.0, 2613.0, 2618.0, 2622.0, 2626.0,
            2705.0, 2724.0, 2726.0, 2727.0, 2809.0,
            2828.0, 2903.0, 2905.0, 2927.0, 2929.0,
            2932.0, 3012.0, 3017.0, 3023.0, 3030.0,
            3113.0, 3114.0, 3131.0, 3224.0, 3228.0,
            3232.0
        };
        Array.Resize(ref values, NNZERO);

        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  HB_MATVEC_A_MEM multiplies a matrix times a vector.");
        Console.WriteLine("");
        Console.WriteLine("  This particular version assumes:");
        Console.WriteLine("  * the matrix is in \"A\" format (assembled),");
        Console.WriteLine("  * the matrix and vectors can fit in memory,");
        Console.WriteLine("  * the matrix and multiplicand have been read into");
        Console.WriteLine("    memory before the routine is called.");
        Console.WriteLine("");
        Console.WriteLine("  For this example, the first vector X is zero except");
        Console.WriteLine("  for a 1 in row 10.  This means A*X should return");
        Console.WriteLine("  column 10 of A.");
        Console.WriteLine("");
        Console.WriteLine("  The second vector X is all 1's.  A*X should be");
        Console.WriteLine("  the sum of the entries of each row.");

        double[] rhsval = HB.hb_matvec_a_mem(NROW, NCOL, NNZERO, NRHS, colptr, rowind, values,
            exact);

        typeMethods.r8mat_print(NROW, NRHS, rhsval, "  The product vectors A*X");

        Console.WriteLine("");
        Console.WriteLine("  Writing the file '" + output_file + "'.");

        HB.hb_file_write(output_file, title, key, totcrd, ptrcrd, indcrd,
            valcrd, rhscrd, mxtype, NROW, NCOL, NNZERO, NELTVL, ptrfmt, indfmt,
            valfmt, rhsfmt, rhstyp, NRHS, NRHSIX, colptr, rowind, values,
            rhsval, rhsptr, rhsind, rhsvec, guess, exact);

    }

    private static void test12()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tests HB_FILE_WRITE;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int NCOL = 32;
        const int NELTVL = 0;
        const int NNZERO = 126;
        const int NRHS = 2;
        const int NRHSIX = 0;
        const int NROW = 32;

        int[] colptr =
        {
            1, 7, 12, 18, 22, 26, 29, 34, 39, 46,
            53, 58, 61, 63, 65, 68, 71, 74, 79, 82,
            85, 88, 90, 94, 97, 102, 106, 110, 112, 117,
            121, 124, 127
        };
        Array.Resize(ref colptr, NCOL + 1);

        double[] exact =
        {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0
        };
        Array.Resize(ref exact, NROW * NRHS);

        double[] guess =
        {
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0
        };
        Array.Resize(ref guess, NROW * NRHS);

        const int indcrd = 8;
        char[] indfmt = "(16I5)".ToCharArray();
        Array.Resize(ref indfmt, 17);

        char[] key = "RUA_32".ToCharArray();
        Array.Resize(ref key, 9);

        char[] mxtype = "RUA".ToCharArray();
        Array.Resize(ref mxtype, 4);

        const string output_file = "rua_32_xa.txt";
        List<string> output = new();
        const int ptrcrd = 3;

        char[] ptrfmt = "(16I5)".ToCharArray();
        Array.Resize(ref ptrfmt, 17);

        const int rhscrd = 12;

        char[] rhsfmt = "(10F7.1)".ToCharArray();
        Array.Resize(ref rhsfmt, 21);

        int[] rhsind = null;
        int[] rhsptr = null;

        char[] rhstyp = "FGX".ToCharArray();
        Array.Resize(ref rhstyp, 4);

        double[] rhsvec = null;
        int[] rowind =
        {
            1, 2, 3, 4, 7, 26, 1, 2, 9, 21,
            28, 2, 3, 6, 8, 9, 29, 3, 4, 5,
            12, 3, 5, 23, 27, 1, 6, 16, 3, 7,
            14, 21, 31, 1, 8, 12, 17, 27, 7, 9,
            10, 13, 19, 23, 27, 1, 10, 11, 21, 23,
            25, 27, 2, 11, 15, 18, 29, 6, 12, 24,
            11, 13, 3, 14, 2, 15, 20, 4, 16, 22,
            4, 16, 17, 6, 10, 18, 20, 30, 1, 19,
            26, 8, 16, 20, 3, 21, 32, 11, 22, 2,
            17, 21, 23, 12, 24, 26, 6, 15, 18, 24,
            25, 13, 18, 22, 26, 5, 24, 26, 27, 9,
            28, 3, 5, 27, 29, 32, 12, 17, 23, 30,
            13, 14, 31, 24, 28, 32
        };
        Array.Resize(ref rowind, NNZERO);

        char[] title = "1Real unsymmetric assembled matrix based on IBM32".ToCharArray();
        Array.Resize(ref title, 73);
        const int totcrd = 36;
        const int valcrd = 13;
        char[] valfmt = "(10F7.1)".ToCharArray();
        Array.Resize(ref valfmt, 21);

        double[] values =
        {
            101.0, 102.0, 103.0, 104.0, 107.0,
            126.0, 201.0, 202.0, 209.0, 221.0,
            228.0, 302.0, 303.0, 306.0, 308.0,
            309.0, 329.0, 403.0, 404.0, 405.0,
            412.0, 503.0, 505.0, 523.0, 527.0,
            601.0, 606.0, 616.0, 703.0, 707.0,
            714.0, 721.0, 731.0, 801.0, 808.0,
            812.0, 817.0, 827.0, 907.0, 909.0,
            910.0, 913.0, 919.0, 923.0, 927.0,
            1001.0, 1010.0, 1011.0, 1021.0, 1023.0,
            1025.0, 1027.0, 1102.0, 1111.0, 1115.0,
            1118.0, 1129.0, 1206.0, 1212.0, 1224.0,
            1311.0, 1313.0, 1403.0, 1414.0, 1502.0,
            1515.0, 1520.0, 1604.0, 1616.0, 1622.0,
            1704.0, 1716.0, 1717.0, 1806.0, 1810.0,
            1818.0, 1820.0, 1830.0, 1901.0, 1919.0,
            1926.0, 2008.0, 2016.0, 2020.0, 2103.0,
            2121.0, 2132.0, 2211.0, 2222.0, 2302.0,
            2317.0, 2321.0, 2323.0, 2412.0, 2424.0,
            2426.0, 2506.0, 2515.0, 2518.0, 2524.0,
            2525.0, 2613.0, 2618.0, 2622.0, 2626.0,
            2705.0, 2724.0, 2726.0, 2727.0, 2809.0,
            2828.0, 2903.0, 2905.0, 2927.0, 2929.0,
            2932.0, 3012.0, 3017.0, 3023.0, 3030.0,
            3113.0, 3114.0, 3131.0, 3224.0, 3228.0,
            3232.0
        };
        Array.Resize(ref values, NNZERO);

        Console.WriteLine("");
        Console.WriteLine("TEST12");
        Console.WriteLine("  HB_VECMAT_A_MEM multiplies a vector times a matrix.");
        Console.WriteLine("");
        Console.WriteLine("  This particular version assumes:");
        Console.WriteLine("  * the matrix is in \"A\" format (assembled),");
        Console.WriteLine("  * the matrix and vectors can fit in memory,");
        Console.WriteLine("  * the matrix and multiplicand have been read into");
        Console.WriteLine("    memory before the routine is called.");
        Console.WriteLine("");
        Console.WriteLine("  For this example, the first vector X is zero except");
        Console.WriteLine("  for a 1 in row 10.  This means A'*X should return");
        Console.WriteLine("  row 10 of A.");
        Console.WriteLine("");
        Console.WriteLine("  The second vector X is all 1's.  A'*X should be");
        Console.WriteLine("  the sum of the entries of each column.");

        double[] rhsval = HB.hb_vecmat_a_mem(NROW, NCOL, NNZERO, NRHS, colptr, rowind, values,
            exact);

        typeMethods.r8mat_print(NCOL, NRHS, rhsval, "  The product vectors A'*X");

        Console.WriteLine("");
        Console.WriteLine("  Writing the file '" + output_file + "'.");

        HB.hb_file_write(output_file, title, key, totcrd, ptrcrd, indcrd,
            valcrd, rhscrd, mxtype, NROW, NCOL, NNZERO, NELTVL, ptrfmt, indfmt,
            valfmt, rhsfmt, rhstyp, NRHS, NRHSIX, colptr, rowind, values,
            rhsval, rhsptr, rhsind, rhsvec, guess, exact);
    }
}