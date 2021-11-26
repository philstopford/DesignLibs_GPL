﻿using System;
using Burkardt.IO;
using Burkardt.Types;

namespace CNFTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CNF_IO_TEST.
        //
        //  Discussion:
        //
        //    CNF_IO_TEST tests the CNF_IO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CNF_IO_TEST");

        Console.WriteLine("  Test the CNF_IO library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        test07();

        Console.WriteLine("");
        Console.WriteLine("CNF_IO_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 calls CNF_WRITE to write a small CNF example to a CNF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int C_NUM = 2;
        const int L_NUM = 5;

        const string cnf_file_name = "cnf_io_v3_c2.cnf";
        int[] l_c_num =  {
                2, 3
            }
            ;
        int[] l_val =  {
                1, -3, 2, 3, -1
            }
            ;
        const int v_num = 3;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  CNF_WRITE can write CNF data to a CNF file.");
        Console.WriteLine("");
        Console.WriteLine("  Here is the data:");

        CNF.cnf_print(v_num, C_NUM, L_NUM, l_c_num, l_val);

        Console.WriteLine("");
        Console.WriteLine("  Now we call CNF_WRITE to store this information");
        Console.WriteLine("  in the file \"" + cnf_file_name + "\".");

        CNF.cnf_write(v_num, C_NUM, L_NUM, l_c_num, l_val, cnf_file_name);

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 calls CNF_HEADER_READ to read the header of a small CNF example file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int c_num = 0;
        const string cnf_file_name = "cnf_io_v3_c2.cnf";
        int l_num = 0;
        int v_num = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  CNF_HEADER_READ reads the header of a CNF file.");

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + cnf_file_name + "\".");

        bool error = CNF.cnf_header_read(cnf_file_name, ref v_num, ref c_num, ref l_num);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("  The header information could not be read.");
                return;
        }

        Console.WriteLine("");
        Console.WriteLine("  The number of variables       V_NUM  = " + v_num + "");
        Console.WriteLine("  The number of clauses         C_NUM  = " + c_num + "");
        Console.WriteLine("  The number of signed literals L_NUM  = " + l_num + "");

    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 calls CNF_DATA_READ to read the data of a small CNF example file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int c_num = 0;
        const string cnf_file_name = "cnf_io_v3_c2.cnf";
        int l_num = 0;
        int v_num = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  CNF_DATA_READ reads the data of a CNF file.");

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + cnf_file_name + "\".");

        bool error = CNF.cnf_header_read(cnf_file_name, ref v_num, ref c_num, ref l_num);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("  The header information could not be read.");
                return;
        }

        Console.WriteLine("");
        Console.WriteLine("  The number of variables       V_NUM  = " + v_num + "");
        Console.WriteLine("  The number of clauses         C_NUM  = " + c_num + "");
        Console.WriteLine("  The number of signed literals L_NUM  = " + l_num + "");

        int[] l_c_num = new int[c_num];
        int[] l_val = new int[l_num];

        CNF.cnf_data_read(cnf_file_name, v_num, c_num, l_num, ref l_c_num, ref l_val);

        Console.WriteLine("");
        Console.WriteLine("  Here is the data as read from the file:");

        CNF.cnf_print(v_num, c_num, l_num, l_c_num, l_val);
    }

    private static void test04()

        //*****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 calls CNF_WRITE to write a CNF example to a CNF file.
        //
        //  Discussion:
        //
        //    This formula is used as an example in the Quinn reference.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int C_NUM = 18;
        const int L_NUM = 36;

        const string cnf_file_name = "cnf_io_v16_c18.cnf";
        int[] l_c_num =  {
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2
            }
            ;
        int[] l_val =  {
                1, 2,
                -2, -4,
                3, 4,
                -4, -5,
                5, -6,
                6, -7,
                6, 7,
                7, -16,
                8, -9,
                -8, -14,
                9, 10,
                9, -10,
                -10, -11,
                10, 12,
                11, 12,
                13, 14,
                14, -15,
                15, 16
            }
            ;
        const int v_num = 16;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  CNF_WRITE can write CNF data to a CNF file.");
        Console.WriteLine("");
        Console.WriteLine("  Here is the data to be written to the file:");

        CNF.cnf_print(v_num, C_NUM, L_NUM, l_c_num, l_val);

        Console.WriteLine("");
        Console.WriteLine("  Now we call CNF_WRITE to store this information");
        Console.WriteLine("  in the file \"" + cnf_file_name + "\".");

        CNF.cnf_write(v_num, C_NUM, L_NUM, l_c_num, l_val, cnf_file_name);

    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 calls CNF_HEADER_READ to read the header of a small CNF example file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int c_num = 0;
        const string cnf_file_name = "cnf_io_v16_c18.cnf";
        int l_num = 0;
        int v_num = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  CNF_HEADER_READ reads the header of a CNF file.");

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + cnf_file_name + "\".");

        bool error = CNF.cnf_header_read(cnf_file_name, ref v_num, ref c_num, ref l_num);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("  The header information could not be read.");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("  The number of variables       V_NUM  = " + v_num + "");
                Console.WriteLine("  The number of clauses         C_NUM  = " + c_num + "");
                Console.WriteLine("  The number of signed literals L_NUM  = " + l_num + "");
                break;
        }
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 calls CNF_DATA_READ to read the data of a small CNF example file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int c_num = 0;
        const string cnf_file_name = "cnf_io_v16_c18.cnf";
        int l_num = 0;
        int v_num = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  CNF_DATA_READ reads the data of a CNF file.");

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + cnf_file_name + "\".");

        bool error = CNF.cnf_header_read(cnf_file_name, ref v_num, ref c_num, ref l_num);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("  The header information could not be read.");
                return;
        }

        Console.WriteLine("");
        Console.WriteLine("  The number of variables       V_NUM  = " + v_num + "");
        Console.WriteLine("  The number of clauses         C_NUM  = " + c_num + "");
        Console.WriteLine("  The number of signed literals L_NUM  = " + l_num + "");

        int[] l_c_num = new int[c_num];
        int[] l_val = new int[l_num];

        CNF.cnf_data_read(cnf_file_name, v_num, c_num, l_num, ref l_c_num, ref l_val);

        Console.WriteLine("");
        Console.WriteLine("  Here is the data as read from the file:");

        CNF.cnf_print(v_num, c_num, l_num, l_c_num, l_val);

    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 calls CNF_EVALUATE to evaluate a formula.
        //
        //  Discussion:
        //
        //    This formula is used as an example in the Quinn reference.
        //    Here, we seek the logical inputs that make the formula true.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int C_NUM = 18;
        const int L_NUM = 36;
        const int V_NUM = 16;

        int i;
        int j;
        int[] l_c_num =  {
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2
            }
            ;
        int[] l_val =  {
                1, 2,
                -2, -4,
                3, 4,
                -4, -5,
                5, -6,
                6, -7,
                6, 7,
                7, -16,
                8, -9,
                -8, -14,
                9, 10,
                9, -10,
                -10, -11,
                10, 12,
                11, 12,
                13, 14,
                14, -15,
                15, 16
            }
            ;
        bool[] v_val = new bool[V_NUM];

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  Seek inputs to a circuit that produce a 1 (TRUE) output.");
        Console.WriteLine("");
        Console.WriteLine("  Here is the CNF data defining the formula:");

        CNF.cnf_print(V_NUM, C_NUM, L_NUM, l_c_num, l_val);
        //
        //  Initialize the logical vector.
        //
        for (j = 0; j < V_NUM; j++)
        {
            v_val[j] = false;
        }

        //
        //  Compute the number of binary vectors to check.
        //
        int ihi = (int)Math.Pow(2, V_NUM);

        Console.WriteLine("");
        Console.WriteLine("  Number of input vectors to check is " + ihi + "");
        Console.WriteLine("");
        Console.WriteLine("   #       Index    ---------Input Values----------");
        Console.WriteLine("");
        //
        //  Check every possible input vector.
        //
        int solution_num = 0;

        for (i = 0; i < ihi; i++)
        {
            bool f_val = CNF.cnf_evaluate(V_NUM, C_NUM, L_NUM, l_c_num, l_val, v_val);

            switch (f_val)
            {
                case true:
                {
                    solution_num += 1;
                    string cout = "  " + solution_num.ToString().PadLeft(2)
                                       + "  " + i.ToString().PadLeft(10);
                    for (j = 0; j < V_NUM; j++)
                    {

                        cout += v_val[j] ? 1 : 0; // v_val[j];
                    }

                    Console.WriteLine(cout);
                    break;
                }
            }

            typeMethods.lvec_next(V_NUM, ref v_val);
        }

        //
        //  Report.
        //
        Console.WriteLine("");
        Console.WriteLine("  Number of solutions found was " + solution_num + "");
    }
}