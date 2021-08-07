using System;
using Burkardt.PentominoNS;
using Burkardt.PolyominoNS;
using Burkardt.Types;

namespace PolyominoTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for polyominoes_test.
            //
            //  Discussion:
            //
            //    polyominoes_test tests polyominoes.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 April 2020
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("polyominoes_test");
            Console.WriteLine("  Test polyominoes.");

            pentomino_matrix_test();
            pentomino_plot_test();
            polyomino_condense_test();
            polyomino_embed_number_test();
            polyomino_embed_list_test();
            polyomino_enumerate_test();
            polyomino_index_test();
            polyomino_lp_write_test();
            polyomino_transform_test();
            Console.WriteLine("");
            Console.WriteLine("polyominoes_test");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void pentomino_matrix_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PENTOMINO_MATRIX_TEST tests PENTOMINO_MATRIX.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 April 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int j;
            int k;
            int[] p = null;
            int p_m = 0;
            int p_n = 0;
            string[] pentominoes =
                { "F", "I", "L", "N", "P", "T", "U", "V", "W", "X", "Y", "Z" };

            Console.WriteLine("");
            Console.WriteLine("PENTOMINO_MATRIX_TEST");
            Console.WriteLine("  PENTOMINO_MATRIX returns a 0/1 matrix representing a pentomino.");

            for (k = 0; k < 12; k++)
            {
                Pentomino.pentomino_matrix(pentominoes[k], ref p_m, ref p_n, ref p);
                Console.WriteLine("");
                Console.WriteLine("  " + pentominoes[k] + " pentomino (" + p_m + "," + p_n + "):");
                Console.WriteLine("");
                for (i = 0; i < p_m; i++)
                {
                    string cout = "    ";
                    for (j = 0; j < p_n; j++)
                    {
                        cout += p[i * p_n + j];
                    }

                    Console.WriteLine(cout);
                }
            }

        }

        static void pentomino_plot_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PENTOMINO_PLOT_TEST tests PENTOMINO_PLOT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 April 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int k;
            string label;
            int[] p = null;
            int p_m = 0;
            int p_n = 0;
            string[] pentominoes =
                { "F", "I", "L", "N", "P", "T", "U", "V", "W", "X", "Y", "Z" };

            Console.WriteLine("");
            Console.WriteLine("PENTOMINO_PLOT_TEST");
            Console.WriteLine("  PENTOMINO_PLOT plots a pentomino.");

            for (k = 0; k < 12; k++)
            {
                Pentomino.pentomino_matrix(pentominoes[k], ref p_m, ref p_n, ref p);
                label = pentominoes[k];
                Pentomino.pentomino_plot(p_m, p_n, p, label);
            }

        }

        static void polyomino_condense_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYOMINO_CONDENSE_TEST tests POLYOMINO_CONDENSE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 April 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Local parameters:
            //
            //    Local, int MP, NP, the number of rows and columns in the representation
            //    of the polyomino P.
            //
            //    Local, int P[MP*NP], a matrix representing the polyomino.  
            //
        {
            int mp;
            int np;
            int[] p1 = { 0, 1, 1, 1, 1, 0, 0, 1, 0 };
            int[] p2 = { 0, 1, 2, 1, 3, 0, 0, -9, 0 };
            int[] p3 = { 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0 };
            int[] p4 = { 0, 0, 0, 0, 0, 0, 0, 0 };

            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_CONDENSE_TEST:");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  POLYOMINO_CONDENSE 'cleans up' a matrix that is supposed");
            Console.WriteLine("  to represent a polyomino:");
            Console.WriteLine("  * nonzero entries are set to 1;");
            Console.WriteLine("  * initial and final zero rows and columns are deleted.");
            //
            //  Nothing happens:
            //
            mp = 3;
            np = 3;
            polyomino_condense_demo(mp, np, p1);
            //
            //  Nonzero, but non-one entries are set to 1.
            //
            mp = 3;
            np = 3;
            polyomino_condense_demo(mp, np, p2);
            //
            //  Extraneous zero rows and columns are removed.
            //
            mp = 3;
            np = 4;
            polyomino_condense_demo(mp, np, p3);
            //
            //  Null matrices are detected.
            //
            mp = 2;
            np = 4;
            polyomino_condense_demo(mp, np, p4);
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_CONDENSE_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void polyomino_condense_demo(int mp, int np, int[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    polyomino_condense_demo demonstrates POLYOMINO_CONDENSE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 April 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int MP, NP, the number of rows and columns in the representation
            //    of the polyomino P.
            //
            //    Input, int P[MP*NP], a matrix representing the polyomino.  
            //
            //    Local, int MQ, NQ, the number of rows and columns in the representation
            //    of the condensed polyomino Q.
            //
            //    Local, int Q[MQ*NQ], a matrix representing the condensed polyomino.  
            //
        {
            string label;
            int mq = 0;
            int nq = 0;
            int[] q = null;

            label = "  The initial (" + (mp) + "," + (np) + ") polynomino P:";
            Polyomino.polyomino_print(mp, np, p, label);

            Polyomino.polyomino_condense(mp, np, p, ref mq, ref nq, ref q);

            label = "  The condensed (" + (mq) + "," + (nq) + ") polynomino Q:";
            Polyomino.polyomino_print(mq, nq, q, label);
        }

        static void polyomino_embed_list_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYOMINO_EMBED_LIST_TEST tests POLYOMINO_EMBED_LIST.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 May 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int j;
            int k;
            int[] list;
            int mk;
            int nk;
            int mp = 3;
            int mq;
            int mr = 4;
            int np = 2;
            int nq;
            int nr = 4;
            int number;
            int[] p =
            {
                0, 0, 1,
                1, 1, 1
            };
            int[] q = new int[4*4];
            int[] r =
            {
                0, 1, 1, 1,
                1, 1, 1, 0,
                1, 0, 1, 1,
                1, 1, 1, 1
            };

            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_EMBED_LIST_TEST:");
            Console.WriteLine("  POLYOMINO_EMBED_LIST lists the offsets used");
            Console.WriteLine("  to embed a fixed polyomino in a region.");

            Polyomino.polyomino_print(mr, nr, r, "  The given region R:");

            Polyomino.polyomino_print(mp, np, p, "  The given polyomino P:");
            //
            //  Get the number of embeddings.
            //
            number = Polyomino.polyomino_embed_number(mr, nr, r, mp, np, p);

            Console.WriteLine("");
            Console.WriteLine("  As a fixed polyomino, P can be embedded in R in " + number + " ways");
            /*
            Get the list of embeddings.
            */
            list = Polyomino.polyomino_embed_list(mr, nr, r, mp, np, p, number);

            for (k = 0; k < number; k++)
            {
                mk = list[k + 0 * number];
                nk = list[k + 1 * number];
                mq = mr;
                nq = nr;

                for (j = 0; j < nq; j++)
                {
                    for (i = 0; i < mq; i++)
                    {
                        q[i + j * mq] = r[i + j * mr];
                    }
                }

                for (j = 0; j < np; j++)
                {
                    for (i = 0; i < mp; i++)
                    {
                        q[i + mk + (j + nk) * mq] = q[i + mk + (j + nk) * mq] + p[i + j * mp];
                    }
                }

                Console.WriteLine("");
                Console.WriteLine("  Embedding number " + k + ":");
                Console.WriteLine("");
                for (i = 0; i < mq; i++)
                {
                    string cout = "";
                    for (j = 0; j < nq; j++)
                    {
                        cout += " " + q[i + j * mq];
                    }

                    Console.WriteLine(cout);
                }
            }
        }

        static void polyomino_embed_number_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYOMINO_EMBED_NUMBER_TEST tests POLYOMINO_EMBED_NUMBER.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 May 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int mp = 3;
            int mr = 4;
            int np = 2;
            int nr = 4;
            int number;
            int[] p =
            {
                0, 0, 1,
                1, 1, 1
            };
            int[] r =
            {
                0, 1, 1, 1,
                1, 1, 1, 0,
                1, 0, 1, 1,
                1, 1, 1, 1
            };

            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_EMBED_NUMBER_TEST:");
            Console.WriteLine("  POLYOMINO_EMBED_NUMBER reports the number of ways a");
            Console.WriteLine("  fixed polyomino can be embedded in a region.");

            Polyomino.polyomino_print(mr, nr, r, "  The given region R:");

            Polyomino.polyomino_print(mp, np, p, "  The given polyomino P:");

            number = Polyomino.polyomino_embed_number(mr, nr, r, mp, np, p);

            Console.WriteLine("");
            Console.WriteLine("  As a fixed polyomino, P can be embedded in R in " + number + " ways");

            return;
        }

        static void polyomino_enumerate_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYOMINO_ENUMERATE_TEST tests the POLYOMINO_ENUMERATE library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 May 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_ENUMERATE_TEST:");
            Console.WriteLine("  C++ version,");
            Console.WriteLine("  POLYOMINO_ENUMERATE enumerates various classes of polyomino.");

            polyomino_enumerate_chiral_test();
            polyomino_enumerate_fixed_test();
            polyomino_enumerate_free_test();
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_ENUMERATE_TEST:");
            Console.WriteLine("  Normal end of execution.");

            return;
        }

        static void polyomino_enumerate_chiral_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYOMINO_ENUMERATE_CHIRAL_TEST tests POLYOMINO_ENUMERATE_CHIRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 May 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n_data = 0;
            long number = 0;
            int order = 0;

            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_ENUMERATE_CHIRAL_TEST:");
            Console.WriteLine("  POLYOMINO_ENUMERATE_CHIRAL returns the number of chiral");
            Console.WriteLine("  polyominoes of a given order;");

            n_data = 0;
            Console.WriteLine("");
            Console.WriteLine("   Order      Number");
            Console.WriteLine("");

            for (;;)
            {
                Polyomino.polyomino_enumerate_chiral(ref n_data, ref order, ref number);

                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + order.ToString().PadLeft(4)
                                       + "  " + number.ToString().PadLeft(24) + "");
            }
        }

        static void polyomino_enumerate_fixed_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYOMINO_ENUMERATE_FIXED_TEST tests POLYOMINO_ENUMERATE_FIXED.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n_data = 0;
            long number = 0;
            int order = 0;

            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_ENUMERATE_FIXED_TEST:");
            Console.WriteLine("  POLYOMINO_ENUMERATE_FIXED returns the number of fixed");
            Console.WriteLine("  polyominoes of a given order;");

            n_data = 0;
            Console.WriteLine("");
            Console.WriteLine("   Order      Number");
            Console.WriteLine("");

            for (;;)
            {
                Polyomino.polyomino_enumerate_fixed(ref n_data, ref order, ref number);

                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + order.ToString().PadLeft(4)
                                       + "  " + number.ToString().PadLeft(24) + "");
            }
        }

        static void polyomino_enumerate_free_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYOMINO_ENUMERATE_FREE_TEST tests POLYOMINO_ENUMERATE_FREE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n_data = 0;
            long number = 0;
            int order = 0;

            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_ENUMERATE_FREE_TEST:");
            Console.WriteLine("  POLYOMINO_ENUMERATE_FREE returns the number of free");
            Console.WriteLine("  polyominoes of a given order;");

            n_data = 0;
            Console.WriteLine("");
            Console.WriteLine("   Order      Number");
            Console.WriteLine("");

            for (;;)
            {
                Polyomino.polyomino_enumerate_free(ref n_data, ref order, ref number);

                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + order.ToString().PadLeft(4)
                                       + "  " + number.ToString().PadLeft(24) + "");
            }
        }

        static void polyomino_index_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYOMINO_INDEX_TEST tests POLYOMINO_INDEX.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 April 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int j;
            int m = 3;
            int n = 4;
            //
            //  P is listed in column-major order;
            //
            int[] p =
            {
                1, 1, 0,
                0, 1, 1,
                1, 1, 1,
                1, 0, 0
            };
            int[] pin;

            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_INDEX_TEST");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  POLYOMINO_INDEX assigns an index to each nonzero entry");
            Console.WriteLine("  of a polyomino.");

            Polyomino.polyomino_print(m, n, p, "  The polyomino P:");

            pin = Polyomino.polyomino_index(m, n, p);

            Console.WriteLine("");
            Console.WriteLine("  PIN: Index vector for P:");
            Console.WriteLine("");
            for (i = 0; i < m; i++)
            {
                string cout = "";
                for (j = 0; j < n; j++)
                {
                    cout += " " + pin[i + j * m];
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_INDEX_TEST");
            Console.WriteLine("  Normal end of execution.");
        }

        static void polyomino_lp_write_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYOMINO_LP_WRITE_TEST tests POLYOMINO_LP_WRITE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 May 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] a;
            int[] b;
            string filename = "reid.lp";
            string label = "\\ LP file for the Reid example, created by POLYOMINO_LP_WRITE.";
            int m = 0;
            int n = 0;

            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_LP_WRITE_TEST:");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  POLYOMINO_LP_WRITE writes an LP file associated");
            Console.WriteLine("  with a binary programming problem for tiling a region");
            Console.WriteLine("  with copies of a single polyomino.");
            //
            //  Get the coefficients and right hand side for the Reid system.
            //
            polyomino_monohedral_example_reid_size(ref m, ref n);

            a = new int[m * n];
            b = new int[m];

            polyomino_monohedral_example_reid_system(m, n, a, b);
            //
            //  Create the LP file.
            //
            Polyomino.polyomino_lp_write(filename, label, m, n, a, b);

            Console.WriteLine("");
            Console.WriteLine("  Created the LP file '" + filename + "'");

            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_LP_WRITE_TEST:");
            Console.WriteLine("  Normal end of execution.");

            return;
        }

        static void polyomino_monohedral_example_reid_size(ref int m, ref int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYOMINO_MONOHEDRAL_EXAMPLE_REID_SIZE returns the size of the Reid system.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 May 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, int &M, &N, the number of equations and variables.
            //
        {
            m = 9;
            n = 10;

        }

        static void polyomino_monohedral_example_reid_system(int m, int n, int[] a, int[] b)

            //****************************************************************************80
            /*
            Purpose:
            
            POLYOMINO_MONOHEDRAL_EXAMPLE_REID_SYSTEM sets up the Reid linear system.
            
            Licensing:
            
            This code is distributed under the GNU LGPL license.
            
            Modified:
            
            17 May 2018
            
            Author:
            
            John Burkardt
            
            Parameters:
            
            Output, double A[9*10], the system matrix.
            
            Output, double B[9], the right hand side.
            */
        {
            /*
            Note that the matrix is specified in column major form.
            */
            int[] a_save =
            {
                1, 1, 0, 0, 0, 0, 0, 0, 2,
                0, 0, 1, 1, 0, 0, 0, 0, 2,
                0, 0, 0, 1, 1, 0, 0, 0, 2,
                0, 0, 0, 0, 0, 1, 1, 0, 2,
                0, 0, 0, 0, 0, 0, 1, 1, 2,
                1, 0, 1, 0, 0, 0, 0, 0, 2,
                0, 1, 0, 1, 0, 0, 0, 0, 2,
                0, 0, 1, 0, 0, 1, 0, 0, 2,
                0, 0, 0, 1, 0, 0, 1, 0, 2,
                0, 0, 0, 0, 1, 0, 0, 1, 2
            };

            int[] b_save = { 1, 1, 1, 1, 1, 1, 1, 1, 8 };

            typeMethods.i4mat_copy(9, 10, a_save, ref a);
            typeMethods.i4vec_copy(9, b_save, ref b);
        }

        static void polyomino_transform_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYOMINO_TRANSFORM_TEST tests POLYOMINO_TRANSFORM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 April 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Local parameters:
            //
            //    Local, int M, N, the number of rows and columns in the representation
            //    of the polyomino P.
            //
            //    Local, int P[M*N], a matrix of 0"s and 1"s representing the 
            //    polyomino.  The matrix should be "tight", that is, there should be a
            //    1 in row 1, and in column 1, and in row M, and in column N.
            //
        {
            string label;
            int m = 3;
            int mq = 0;
            int n = 3;
            int nq = 0;
            //
            //  P is given by columns, not rows.
            //
            int[] p =
            {
                0, 1, 0,
                1, 1, 1,
                1, 0, 0
            };
            int[] q = null;
            int reflect;
            int rotate;

            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_TRANSFORM_TEST:");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  POLYOMINO_TRANSFORM can transform a polyomino.");
            Console.WriteLine("  Generate all 8 combinations of rotation and reflection");
            Console.WriteLine("  applied to a polyomino represented by a binary matrix.");

            Polyomino.polyomino_print(m, n, p, "  The given polyomino P:");

            for (reflect = 0; reflect <= 1; reflect++)
            {
                for (rotate = 0; rotate <= 3; rotate++)
                {
                    Polyomino.polyomino_transform(m, n, p, rotate, reflect, ref mq, ref nq, ref q);

                    label = "  P after " + (reflect) + " reflections and " + (rotate) + " rotations:";

                    Polyomino.polyomino_print(mq, nq, q, label);
                }
            }

            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_TRANSFORM_TEST:");
            Console.WriteLine("  Normal end of execution.");

        }
    }
}