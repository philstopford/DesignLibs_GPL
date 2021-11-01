using System;
using System.Collections.Generic;
using System.IO;
using Burkardt;
using Burkardt.IO;
using Burkardt.Table;
using Burkardt.Types;

namespace SVDBasisTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for SVD_BASIS.
            //
            //  Discussion:
            //
            //    SVD_BASIS forms a basis from the SVD of a set of data vectors.
            //
            //    This program uses the singular value decomposition (SVD) to analyze
            //    a set of data, and extract a number of dominant modes.
            //
            //    This program is intended as an intermediate application, in
            //    the following situation:
            //
            //    A) a "high fidelity" or "high resolution" PDE solver is used
            //       to determine many (say N = 500) solutions of a discretized
            //       PDE at various times, or parameter values.  Each solution
            //       may be regarded as an M vector.  Typically, each solution
            //       involves an M by M linear system, greatly reduced in
            //       complexity because of bandedness or sparsity.
            //
            //    B) This program is applied to extract L dominant modes from
            //       the N solutions.  This is done using the singular value
            //       decomposition of the M by N matrix, each of whose columns
            //       is one of the original solution vectors.
            //
            //    C) a "reduced order model" program may then attempt to solve
            //       a discretized version of the PDE, using the L dominant
            //       modes as basis vectors.  Typically, this means that a dense
            //       L by L linear system will be involved.
            //
            //    Thus, the program might read in 500 files, and write out
            //    5 or 10 files of the corresponding size and "shape", representing
            //    the dominant solution modes.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DATA_FILE_BASE_MAX = 20;

            string basis_file;
            int basis_num = 0;
            bool comment;
            char comment_char;
            int comp_num = 0;
            string data_file = "";
            int data_file_base_num = 0;
            string[] data_file_base = new string[DATA_FILE_BASE_MAX];
            int data_file_num = 0;
            int dim_num = 0;
            string file_name;
            int i;
            int ii;
            int j;
            int k;
            int l;
            int node_num = 0;
            double[] point;
            int point_num = 0;
            double[] sval;
            double[] table;

            Console.WriteLine("");
            Console.WriteLine("SVD_BASIS:");
            Console.WriteLine("");
            Console.WriteLine("  Given a PDE for which:");
            Console.WriteLine("    C is the number of components of the solution");
            Console.WriteLine("      at any single point,");
            Console.WriteLine("    P is the number of points where a solution is given,");
            Console.WriteLine("    N is the number of solution vectors,");
            Console.WriteLine("    L is the number of modes to be extracted.");
            Console.WriteLine("");
            Console.WriteLine("  Then we let M = C*P be the abstract spatial dimension.");
            Console.WriteLine("");
            Console.WriteLine("  Set up A, the M by N matrix of solution vectors,");
            Console.WriteLine("");
            Console.WriteLine("  Get A = U * S * V', the singular value decomposition.");
            Console.WriteLine("");
            Console.WriteLine("  The first L columns of U are our modes.");
            Console.WriteLine("");
            //
            //  What is the basis size?
            //
            Console.WriteLine("  How many basis vectors (L) are to be extracted?");
            basis_num = Convert.ToInt32(Console.ReadLine());

            Console.WriteLine("");
            Console.WriteLine("  L = " + basis_num + "");
            //
            //  Gather one or more "base" file names.
            //
            data_file_base_num = 0;

            for (;;)
            {
                if (DATA_FILE_BASE_MAX <= data_file_base_num)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  No more base file names can be entered.");
                    return;
                }

                //
                //  Get the next base file name.
                //
                Console.WriteLine("");
                Console.WriteLine("  You specify a consecutive sequence of file names");
                Console.WriteLine("  by giving the first \"base\" file name.");
                Console.WriteLine("");
                Console.WriteLine("  If there are no more sequences to enter,");
                Console.WriteLine("  just hit RETURN.");
                Console.WriteLine("");
                Console.WriteLine("  Enter a new base file name, or $ if done:");
                //
                //  CIN won't allow you to enter a blank line!
                //  I just don't have the energy today to replace CIN by GETLINE....
                //
                file_name = Console.ReadLine();

                if (file_name == "$")
                {
                    file_name = " ";
                }

                if (typeMethods.s_len_trim(file_name) <= 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  RETURN was entered.");
                    Console.WriteLine("  Presumably, there are no more file sequences.");
                    break;
                }

                data_file_base[data_file_base_num] = file_name;
                data_file_base_num = data_file_base_num + 1;

                Console.WriteLine("");
                Console.WriteLine(data_file_base_num + ":  \"" + file_name + "\"");
                //
                //  For the very first base file, get the data sizes.
                //
                if (data_file_base_num == 1)
                {
                    TableHeader h = typeMethods.r8mat_header_read(file_name);
                    comp_num = h.m;
                    node_num = h.n;

                    dim_num = comp_num * node_num;

                    Console.WriteLine("");
                    Console.WriteLine("  According to the first base file,");
                    Console.WriteLine("  The number of solution components C =   " + comp_num + "");
                    Console.WriteLine("  The number of solution points P =       " + node_num + "");
                    Console.WriteLine("  The \"size\" of each solution M = (C*P) = " + dim_num + "");
                    //
                    //  Idiocy check.  L must be less than or equal to M.
                    //
                    if (dim_num < basis_num)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("SVD_BASIS - Fatal error!");
                        Console.WriteLine("");
                        Console.WriteLine("  M < L.");
                        Console.WriteLine("");
                        Console.WriteLine("  That is, the number of modes requested (L) is greater");
                        Console.WriteLine("  than the spatial dimension (M).");
                        Console.WriteLine("  Technically, the program could pad out the answer");
                        Console.WriteLine("  with L-M zero vectors, but instead, we will stop");
                        Console.WriteLine("  assuming you made an error, or a misapprehension.");
                        Console.WriteLine("");
                        Console.WriteLine("SVD_BASIS:");
                        Console.WriteLine("  Abnormal end of execution.");
                        Console.WriteLine("");
                        return;
                    }

                }

            }

            //
            //  Count all the data files.
            //
            data_file_num = 0;

            for (i = 0; i < data_file_base_num; i++)
            {
                data_file = data_file_base[i];

                for (;;)
                {
                    if (!File.Exists(data_file))
                    {
                        break;
                    }

                    data_file_num = data_file_num + 1;

                    Files.file_name_inc_nowrap(ref data_file);
                }
            }

            if (data_file_num == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("SVD_BASIS - Fatal error!");
                Console.WriteLine("  There do not seem to be any solution files;");
                Console.WriteLine("  that is, files whose names are \"incremented\"");
                Console.WriteLine("  versions of the first file name.");
                Console.WriteLine("");
                Console.WriteLine("  The first file we looked for was \"" + data_file + "\"");
                Console.WriteLine("");
                Console.WriteLine("SVD_BASIS:");
                Console.WriteLine("  Abnormal end of execution.");
                Console.WriteLine("");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  The number of data files N = " + data_file_num + "");
            //
            //  Set up an array to hold all the data.
            //
            point_num = data_file_num;

            Console.WriteLine("");
            Console.WriteLine("  The data is stored in an M by N matrix A.");
            Console.WriteLine("");
            Console.WriteLine("  The \"spatial\" dimension M is " + dim_num + "");
            Console.WriteLine("  The number of data points N is " + point_num + "");
            //
            //  Allocate space for the POINT array.
            //
            point = new double[dim_num * point_num];
            //
            //  Read the data.
            //
            l = 0;

            for (ii = 0; ii < data_file_base_num; ii++)
            {
                data_file = data_file_base[ii];

                for (;;)
                {
                    if (!File.Exists(data_file))
                    {
                        break;
                    }

                    table = typeMethods.r8mat_data_read(data_file, comp_num, node_num);

                    k = 0;
                    for (j = 0; j < node_num; j++)
                    {
                        for (i = 0; i < comp_num; i++)
                        {
                            point[k + l * dim_num] = table[i + j * comp_num];
                            k = k + 1;
                        }
                    }

                    l = l + 1;
                    Files.file_name_inc_nowrap(ref data_file);

                }
            }

            Console.WriteLine("");
            Console.WriteLine("  The data has been read into the matrix A.");
            //
            //----------------------------------------------------------------------------
            //
            //  Compute the SVD of A.
            //
            //----------------------------------------------------------------------------
            //
            sval = new double[basis_num];

            SingleValueDecomposition.singular_vectors(dim_num, point_num, basis_num, ref point, ref sval);
            //
            //----------------------------------------------------------------------------
            //
            //  Write the first L left singular vectors (columns of U) to files.
            //
            //----------------------------------------------------------------------------
            //
            Console.WriteLine("");
            Console.WriteLine("SVD_BASIS:");
            Console.WriteLine("  Ready to write the left singular vectors to files.");
            Console.WriteLine("");
            Console.WriteLine("  Do you want comments in the header of the file?");
            Console.WriteLine("  (These begin with the \"#\" character.) (Y/N)");
            Console.WriteLine("");
            Console.WriteLine("  Enter \"Y\" or \"N\":");

            comment_char = Console.ReadLine().ToCharArray()[0];

            if (comment_char == 'Y' || comment_char == 'y')
            {
                comment = true;
            }
            else
            {
                comment = false;
            }

            basis_file = "svd_000.txt";

            for (j = 0; j < basis_num; j++)
            {
                Files.file_name_inc_nowrap(ref basis_file);

                if (j + 1 == 1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Writing first file " + basis_file + "");
                }

                if (j + 1 == basis_num)
                {
                    Console.WriteLine("  Writing last file  " + basis_file + "");
                }

                basis_write(basis_file, comp_num, node_num, sval[j],
                    point, comment, uIndex: +0 + j * dim_num);

                Console.WriteLine("");
                Console.WriteLine("SVD_BASIS:");
                Console.WriteLine("  Normal end of execution.");
                Console.WriteLine("");
            }

            static void basis_write(string file_out_name, int m, int n, double s, double[] u,
                    bool comment, int uIndex = 0)

                //****************************************************************************80
                //
                //  Purpose:
                //
                //    BASIS_WRITE writes a basis vector to a file.
                //
                //  Discussion:
                //
                //    The initial lines of the file are comments, which begin with a
                //    "#" character.
                //
                //  Licensing:
                //
                //    This code is distributed under the GNU LGPL license. 
                //
                //  Modified:
                //
                //    22 November 2011
                //
                //  Author:
                //
                //    John Burkardt
                //
                //  Parameters:
                //
                //    Input, string FILE_OUT_NAME, the name of the file to write.
                //
                //    Input, int M, the number of data components.
                //
                //    Input, int N, the number of data items.
                //
                //    Input, double S, the associated singular value.
                //
                //    Input, double U[M*N], the data values.
                //
                //    Input, bool COMMENT, is TRUE if comments are to be included.
                //
            {
                List<string> file_out = new List<string>();
                int i;
                int j;

                if (comment)
                {
                    file_out.Add("#  " + file_out_name + "");
                    file_out.Add("#  created by routine BASIS_WRITE.C" + "");
                    file_out.Add("#  part of SVD_BASIS.C." + "");
                    file_out.Add("#");
                    file_out.Add("#  Number of components M =  " + m.ToString().PadLeft(12) + "");
                    file_out.Add("#  Number of items N =       " + n.ToString().PadLeft(12) + "");
                    file_out.Add("#  Singular value S =        " + s.ToString().PadLeft(14) + "");
                    file_out.Add("#  EPSILON (unit roundoff) = " + typeMethods.r8_epsilon() + "");
                    file_out.Add("#");
                }

                for (j = 0; j < n; j++)
                {
                    string cout = "";
                    for (i = 0; i < m; i++)
                    {
                        cout += u[i + j * m].ToString().PadLeft(10) + "  ";
                    }

                    file_out.Add(cout);
                }

                try
                {
                    File.WriteAllLines(file_out_name, file_out);
                }
                catch (Exception e)
                {
                    Console.WriteLine("");
                    Console.WriteLine("BASIS_WRITE - Fatal error!");
                    Console.WriteLine("  Could not open the output file.");
                }
            }

        }
    }
}