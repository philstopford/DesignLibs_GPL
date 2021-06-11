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
        
        public static void ccs_to_st ( int m, int n, int ncc, int[] icc, int[] ccc, double[] acc, 
        ref int nst, ref int[] ist, ref int[] jst, ref double[] ast )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ccs_TO_ST converts sparse matrix information from CCS to ST format.
        //
        //  Discussion:
        //
        //    Only JST actually needs to be computed.  The other three output 
        //    quantities are simply copies.  
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
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, int NCC, the number of CCS elements.
        //
        //    Input, int ICC[NCC], the CCS rows.
        //
        //    Input, int CCC[N+1], the CCS compressed columns.
        //
        //    Input, double ACC[NCC], the CCS values.
        //
        //    Output, int &NST, the number of ST elements.
        //
        //    Output, int IST[NST], JST[NST], the ST rows and columns.
        //
        //    Output, double AST[NST], the ST values.
        //
        {
            int j;
            int jhi;
            int jlo;
            int k;
            int khi;
            int klo;

            nst = 0;

            if ( ccc[0] == 0 )
            {
                jlo = 0;
                jhi = n - 1;
  
                for ( j = jlo; j <= jhi; j++ )
                {
                    klo = ccc[j];
                    khi = ccc[j+1] - 1;

                    for ( k = klo; k <= khi; k++ )
                    {
                        ist[nst] = icc[k];
                        jst[nst] = j;
                        ast[nst] = acc[k];
                        nst = nst + 1;
                    }
                }
            }
            else
            {
                jlo = 1;
                jhi = n;
  
                for ( j = jlo; j <= jhi; j++ )
                {
                    klo = ccc[j-1];
                    khi = ccc[j] - 1;

                    for ( k = klo; k <= khi; k++ )
                    {
                        ist[nst] = icc[k-1];
                        jst[nst] = j;
                        ast[nst] = acc[k-1];
                        nst = nst + 1;
                    }
                }
            }
        }
        
        public static double[] ccs_mv ( int m, int n, int ncc, int[] icc, int[] ccc, double[] acc, 
        double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ccs_MV multiplies a CCS matrix by a vector
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, int NCC, the number of CCS values.
        //
        //    Input, int ICC[NCC], the CCS rows.
        //
        //    Input, int CCC[N+1], the compressed CCS columns
        //
        //    Input, double ACC[NCC], the CCS values.
        //
        //    Input, double X[N], the vector to be multiplied.
        //
        //    Output, double ccs_MV[M], the product A*X.
        //
        {
            double[] b;
            int i;
            int j;
            int k;

            b = new double[m];

            for ( i = 0; i < m; i++ )
            {
                b[i] = 0.0;
            }

            for ( j = 0; j < n; j++ )
            {
                for ( k = ccc[j]; k < ccc[j+1]; k++ )
                {
                    i = icc[k];
                    b[i] = b[i] + acc[k] * x[j];
                }
            }

            return b;
        }
        
    }
}