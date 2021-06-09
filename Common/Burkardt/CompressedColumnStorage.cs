using System;
using System.IO;
using System.Linq;
using Burkardt.Table;
using Burkardt.Types;

namespace Burkardt
{
    public static class CompressedColumnStorage
    {
        public static void ccs_data_read(string prefix, int ncc, int n, ref int[] icc, ref int[] ccc,
        ref double[] acc )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ccs_DATA_READ reads data about a sparse matrix in CCS format.
        //
        //  Discussion:
        //
        //    Three files are presumed to exist:
        //    * prefix_icc.txt contains NCC ICC values;
        //    * prefix_ccc.txt contains N+1 CCC values;
        //    * prefix_acc.txt contains NCC ACC values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string PREFIX, a common prefix for the filenames.
        //
        //    Input, int NCC, the number of CCS elements.
        //
        //    Input, int N, the number of columns in the matrix.
        //
        //    Output, int ICC[NCC], the CCS rows.
        //
        //    Output, int CCC[N+1], the compressed CCS columns.
        //
        //    Output, double ACC[NCC], the CCS values.
        //
        {
            string filename_icc = prefix + "_icc.txt";
            typeMethods.i4vec_data_read(filename_icc, ncc, ref icc);

            string filename_ccc = prefix + "_ccc.txt";
            typeMethods.i4vec_data_read(filename_ccc, n + 1, ref ccc);

            string filename_acc = prefix + "_acc.txt";
            typeMethods.r8vec_data_read(filename_acc, ncc, ref acc);
        }

        public static void ccs_header_read(string prefix, ref int ncc, ref int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ccs_HEADER_READ reads header information about a sparse matrix in CCS format.
        //
        //  Discussion:
        //
        //    Three files are presumed to exist:
        //    * prefix_icc.txt contains NCC ICC values;
        //    * prefix_ccc.txt contains N+1 CCC values;
        //    * prefix_acc.txt contains NCC ACC values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string PREFIX, a common prefix for the filenames.
        //
        //    Output, int &NCC, the number of CCS elements.
        //
        //    Output, int &N, the number of columns in the matrix.
        //
        {
            string filename_icc = prefix + "_icc.txt";
            ncc = TableMisc.file_row_count(filename_icc);

            string filename_ccc = prefix + "_ccc.txt";
            n = TableMisc.file_row_count(filename_ccc) - 1;
        }

        public static void ccs_print(int m, int n, int ncc, int[] icc, int[] ccc, double[] acc,
        string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ccs_PRINT prints a sparse matrix in CCS format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows in the matrix.
        //
        //    Input, int N, the number of columns in the matrix.
        //
        //    Input, int NCC, the number of CCS elements.
        //
        //    Input, int ICC[NCC], the CCS rows.
        //
        //    Input, int CCC[N+1], the compressed CCS columns.
        //
        //    Input, double ACC[NCC], the CCS values.
        //
        //    Input, string TITLE, a title.
        //
        {
            int i;
            int j;
            int k;

            Console.WriteLine("");
            Console.WriteLine(title + "");
            Console.WriteLine("     #     I     J         A");
            Console.WriteLine("  ----  ----  ----  ----------------");
            Console.WriteLine("");

            if (ccc[0] == 0)
            {
                j = 0;
                for (k = 0; k < ncc; k++)
                {
                    i = icc[k];
                    while (ccc[j + 1] <= k)
                    {
                        j = j + 1;
                    }

                    Console.WriteLine(k.ToString().PadLeft(4) + "  "
                        + i.ToString().PadLeft(4) + "  "
                        + j.ToString().PadLeft(4) + "  "
                        + acc[k].ToString().PadLeft(16) + "");
                }
            }
            else
            {
                j = 1;
                for (k = 0; k < ncc; k++)
                {
                    i = icc[k];
                    while (ccc[j] <= k + 1)
                    {
                        j = j + 1;
                    }

                    Console.WriteLine((k + 1).ToString().PadLeft(4) + "  "
                        + i.ToString().PadLeft(4) + "  "
                        + j.ToString().PadLeft(4) + "  "
                        + acc[k].ToString().PadLeft(16) + "");
                }
            }

            return;
        }

        public static void ccs_print_some(int i_min, int i_max, int j_min, int j_max, int ncc,
            int n, int[] icc, int[] ccc, double[] acc, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ccs_PRINT_SOME prints some of a sparse matrix in CCS format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I_MIN, IMAX, the first and last rows to print.
        //
        //    Input, int J_MIN, J_MAX, the first and last columns 
        //    to print.
        //
        //    Input, int NCC, the number of CCS elements.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, int ICC[NCC], the CCS rows.
        //
        //    Input, int CCC[N+1], the compressed CCS columns.
        //
        //    Input, double ACC[NCC], the CCS values.
        //
        //    Input, string TITLE, a title.
        //
        {
            int i;
            int j;
            int k;

            Console.WriteLine("");
            Console.WriteLine(title + "");
            Console.WriteLine("     #     I     J         A");
            Console.WriteLine("  ----  ----  ----  ----------------");
            Console.WriteLine("");

            if (ccc[0] == 0)
            {
                j = 0;
                for (k = 0; k < ncc; k++)
                {
                    i = icc[k];
                    while (ccc[j + 2] <= k)
                    {
                        j = j + 1;
                    }

                    if (i_min <= i && i <= i_max &&
                        j_min <= j && j <= j_max)
                    {
                        Console.WriteLine(k.ToString().PadLeft(4) + "  "
                            + i.ToString().PadLeft(4) + "  "
                            + j.ToString().PadLeft(4) + "  "
                            + acc[k].ToString().PadLeft(16) + "");
                    }
                }
            }
            else
            {
                j = 1;
                for (k = 0; k < ncc; k++)
                {
                    i = icc[k];
                    while (ccc[j + 1] <= k + 1)
                    {
                        j = j + 1;
                    }

                    if (i_min <= i && i <= i_max &&
                        j_min <= j && j <= j_max)
                    {
                        Console.WriteLine((k + 1).ToString().PadLeft(4) + "  "
                            + i.ToString().PadLeft(4) + "  "
                            + j.ToString().PadLeft(4) + "  "
                            + acc[k].ToString().PadLeft(16) + "");
                    }
                }
            }
        }

        public static void ccs_write(string prefix, int ncc, int n, int[] icc, int[] ccc,
        double[] acc )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ccs_WRITE writes a sparse matrix in CCS format to 3 files.
        //
        //  Discussion:
        //
        //    Three files will be created:
        //    * prefix_icc.txt contains NCC ICC values;
        //    * prefix_ccc.txt contains N+1 CCC values;
        //    * prefix_acc.txt contains NCC ACC values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string PREFIX, a common prefix for the filenames.
        //
        //    Input, int NCC, the number of CCS elements.
        //
        //    Input, int N, the number of columns in the matrix.
        //
        //    Input, int ICC[NCC], the CCS rows.
        //
        //    Input, int CCC[N+1], the compressed CCS columns.
        //
        //    Input, double ACC[NCC], the CCS values.
        //
        {
            string filename_icc = prefix + "_icc.txt";
            typeMethods.i4vec_write(filename_icc, ncc, icc);

            string filename_ccc = prefix + "_ccc.txt";
            typeMethods.i4vec_write(filename_ccc, n + 1, ccc);

            string filename_acc = prefix + "_acc.txt";
            typeMethods.r8vec_write(filename_acc, ncc, acc);
        }
    }
}