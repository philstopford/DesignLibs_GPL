using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.SortNS;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8st_data_read(string input_filename, int m, int n, int nst,
            ref int[] ist, ref int[] jst, ref double[] ast )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ST_DATA_READ reads the data of an R8ST file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string INPUT_FILENAME, the name of the file.
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of columns .
        //
        //    Input, int NST, the number of nonzeros .
        //
        //    Output, int IST[NST], JST[NST], the row and column indices.
        //
        //    Output, double AST[NST], the nonzero values.
        //
        {
            int k = 0;
            try
            {
                string[] input = File.ReadAllLines(input_filename);
                
                for (k = 0; k < nst; k++)
                {
                    string[] tokens = Helpers.splitStringByWhitespace(input[k]);

                    i4 ti = s_to_i4(tokens[0]);
                    i4 tj = s_to_i4(tokens[1]);
                    r8 taij = s_to_r8(tokens[2]);

                    int i = ti.val;
                    int j = tj.val;
                    double aij = taij.val;
                    
                    ist[k] = i;
                    jst[k] = j;
                    ast[k] = aij;
                }

            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("R8ST_DATA_READ - Fatal error!");
                Console.WriteLine("  I/O error reading data index " + k + "");
            }

        }

        public static void r8st_header_print(int i_min, int i_max, int j_min, int j_max, int m,
                int n, int nst)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8ST_HEADER_PRINT prints the header of an R8ST file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int I_MIN, I_MAX, the minimum and maximum row indices.
            //
            //    Input, int J_MIN, J_MAX, the minimum and maximum column indices.
            //
            //    Input, int M, the number of rows.
            //
            //    Input, int N, the number of columns.
            //
            //    Input, int NST, the number of nonzeros.
            //
        {
            Console.WriteLine("");
            Console.WriteLine("  Header information:");
            Console.WriteLine("");
            Console.WriteLine("  Minimum row index I_MIN = " + i_min + "");
            Console.WriteLine("  Maximum row index I_MAX = " + i_max + "");
            Console.WriteLine("  Minimum col index J_MIN = " + j_min + "");
            Console.WriteLine("  Maximum col index J_MAX = " + j_max + "");
            Console.WriteLine("  Number of rows        M = " + m + "");
            Console.WriteLine("  Number of columns     N = " + n + "");
            Console.WriteLine("  Number of nonzeros  NST = " + nst + "");

        }

        public static void r8st_header_read(string input_filename, ref int i_min, ref int i_max,
        ref int j_min, ref int j_max, ref int m, ref int n, ref int nst )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ST_HEADER_READ reads the header of an R8ST file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string INPUT_FILENAME, the name of the file.
        //
        //    Output, int &I_MIN, &I_MAX, the minimum and maximum row indices.
        //
        //    Output, int &J_MIN, &J_MAX, the minimum and maximum column indices.
        //
        //    Output, int &M, the number of rows.
        //
        //    Output, int &N, the number of columns .
        //
        //    Output, int &NST, the number of nonzeros.
        //
        {
            try
            {
                string[] input = File.ReadAllLines(input_filename);

                nst = 0;
                i_min = +i4_huge();
                i_max = -i4_huge();
                j_min = +i4_huge();
                j_max = -i4_huge();

                for (;;)
                {
                    string[] tokens = Helpers.splitStringByWhitespace(input[nst]);

                    i4 ti = s_to_i4(tokens[0]);
                    i4 tj = s_to_i4(tokens[1]);
                    r8 taij = s_to_r8(tokens[2]);

                    int i = ti.val;
                    int j = tj.val;

                    nst = nst + 1;
                    i_min = Math.Min(i_min, i);
                    i_max = Math.Max(i_max, i);
                    j_min = Math.Min(j_min, j);
                    j_max = Math.Max(j_max, j);

                    if (nst == input.Length)
                    {
                        break;
                    }
                }
                
            }
            catch (Exception e)
            {
            }

            m = i_max - i_min + 1;
            n = j_max - j_min + 1;
        }

        public static void r8st_print(int m, int n, int nst, int[] ist, int[] jst, double[] ast,
        string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ST_PRINT prints a sparse matrix in R8ST format.
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
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, int NST, the number of elements.
        //
        //    Input, int IST[NST], JST[NST], the rows and columns.
        //
        //    Input, double AST[NST], the values.
        //
        //    Input, string TITLE, a title.
        //
        {
            Console.WriteLine("");
            Console.WriteLine(title + "");
            Console.WriteLine("  " + m + " rows and " + n + " columns.");
            Console.WriteLine("     #     I     J       A");
            Console.WriteLine("  ----  ----  ----  --------------");
            Console.WriteLine("");
            for (int k = 0; k < nst; k++)
            {
                Console.WriteLine(k.ToString().PadLeft(4) + "  "
                    + ist[k].ToString().PadLeft(4) + "  "
                    + jst[k].ToString().PadLeft(4) + "  "
                    + ast[k].ToString().PadLeft(16) + "");
            }
        }

        public static void r8st_print_some(int i_min, int i_max, int j_min, int j_max, int nst,
            int[] ist, int[] jst, double[] ast, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ST_PRINT_SOME prints some of a sparse matrix in R8ST format.
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
        //  Parameters:
        //
        //    Input, int I_MIN, IMAX, the first and last rows to print.
        //
        //    Input, int J_MIN, J_MAX, the first and last columns 
        //    to print.
        //
        //    Input, int NST, the number of elements.
        //
        //    Input, int IST[NST], JST[NST], the rows and columns.
        //
        //    Input, double AST[NST], the ST values.
        //
        //    Input, string TITLE, a title.
        //
        {
            Console.WriteLine("");
            Console.WriteLine(title + "");
            Console.WriteLine("     #     I     J       A");
            Console.WriteLine("  ----  ----  ----  --------------");
            Console.WriteLine("");
            for (int k = 0; k < nst; k++)
            {
                if (i_min <= ist[k] && ist[k] <= i_max &&
                    j_min <= jst[k] && jst[k] <= j_max)
                {
                    Console.WriteLine(k.ToString().PadLeft(4) + "  "
                        + ist[k].ToString().PadLeft(4) + "  "
                        + jst[k].ToString().PadLeft(4) + "  "
                        + ast[k].ToString().PadLeft(16) + "");
                }
            }
        }

        public static void r8st_sort_a(int m, int n, int nst, ref int[] ist, ref int[] jst, ref double[] ast )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ST_SORT_A sorts the entries of an R8ST matrix by column.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2014
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
        //    Input, int NST, the number of nonzeros.
        //
        //    Input/output, int IST[NST], JST[NST], the row and column indices.
        //
        //    Input/output, double AST[NST], the nonzero values.
        //
        {
            double aij;
            int cij;
            int i;
            int indx;
            int isgn;
            int j;
            int rij;
            //
            //  Initialize.
            //
            i = 0;
            indx = 0;
            isgn = 0;
            j = 0;
            //
            //  Call the external heap sorter.
            //
            SortHeapExternalData data = new SortHeapExternalData();
            for (;;)
            {
                Sort.sort_heap_external(ref data, nst, ref indx, ref i, ref j, isgn);
                //
                //  Interchange the I and J objects.
                //
                if (0 < indx)
                {
                    rij = ist[i - 1];
                    ist[i - 1] = ist[j - 1];
                    ist[j - 1] = rij;

                    cij = jst[i - 1];
                    jst[i - 1] = jst[j - 1];
                    jst[j - 1] = cij;

                    aij = ast[i - 1];
                    ast[i - 1] = ast[j - 1];
                    ast[j - 1] = aij;
                }
                //
                //  Compare the I and J objects.
                //
                else if (indx < 0)
                {
                    if (jst[i - 1] == jst[j - 1])
                    {
                        if (ist[i - 1] < ist[j - 1])
                        {
                            isgn = -1;
                        }
                        else if (ist[i - 1] == ist[j - 1])
                        {
                            isgn = 0;
                        }
                        else if (ist[j - 1] < ist[i - 1])
                        {
                            isgn = +1;
                        }
                    }
                    else if (jst[i - 1] < jst[j - 1])
                    {
                        isgn = -1;
                    }
                    else if (jst[j - 1] < jst[i - 1])
                    {
                        isgn = +1;
                    }
                }
                else if (indx == 0)
                {
                    break;
                }
            }

            return;
        }

        public static void r8st_transpose(ref int m, ref int n, int nst, ref int[] ist, ref int[] jst,
        double[] ast )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ST_TRANSPOSE transposes an R8ST matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &M, the number of rows.
        //
        //    Input/output, int &N, the number of columns.
        //
        //    Input, int NST, the number of nonzeros.
        //
        //    Input/output, int IST[NST], JST[NST], the row and column indices.
        //
        //    Input, double AST[NST], the nonzero values.
        //
        {
            int k;
            int t;

            t = m;
            m = n;
            n = t;

            for (k = 0; k < nst; k++)
            {
                t = ist[k];
                ist[k] = jst[k];
                jst[k] = t;
            }
        }

        public static void r8st_write(string output_filename, int m, int n, int nst, int[] ist,
        int[] jst, double[] ast )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ST_WRITE writes an R8ST file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string OUTPUT_FILENAME, the name of the file.
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, int NST, the number of nonzeros.
        //
        //    Input, int IST[NST], JST[NST], the row and column indices.
        //
        //    Input, double AST[NST], the nonzero values.
        //
        {
            try
            {
                List<string> output = new List<string>();
                for (int k = 0; k < nst; k++)
                {
                    output.Add("  " + ist[k]
                        + "  " + jst[k]
                        + "  " + ast[k] + "");
                }
                
                File.WriteAllLines(output_filename, output);

            }
            catch (Exception e)
            {
                Console.WriteLine(e);
                throw;
            }
        }
    }
}