using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8vec2_print(int n, double[] a1, double[] a2, string title)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC2_PRINT prints an R8VEC2.
            //
            //  Discussion:
            //
            //    An R8VEC2 is a dataset consisting of N pairs of real values, stored
            //    as two separate vectors A1 and A2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 November 2002
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of components of the vector.
            //
            //    Input, double A1[N], double A2[N], the vectors to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            Console.WriteLine();
            Console.WriteLine(title);
            Console.WriteLine();
            for (int i = 0; i <= n - 1; i++)
            {
                Console.WriteLine(i.ToString().PadLeft(6)
                                  + ": " + a1[i].ToString().PadLeft(14)
                                  + "  " + a2[i].ToString().PadLeft(14));
            }
        }
        
        public static void r8vec3_print(int n, double[] a1, double[] a2, double[] a3, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC3_PRINT prints a triple of real vectors.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 September 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of components of the vector.
            //
            //    Input, double A1[N], double A2[N], double A3[N], the vectors
            //    to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            Console.WriteLine();
            Console.WriteLine(title);
            Console.WriteLine();
            for (int i = 0; i <= n - 1; i++)
            {
                Console.WriteLine(i.ToString().PadLeft(4) + ": "
                                                          + a1[i].ToString().PadLeft(10) + "  "
                                                          + a2[i].ToString().PadLeft(10) + "  "
                                                          + a3[i].ToString().PadLeft(10));
            }
        }
        
        public static void r8vec_print(int n, double[] a, string title, int aIndex = 0)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_PRINT prints an R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of components of the vector.
            //
            //    Input, double A[N], the vector to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {

            Console.WriteLine(title);
            Console.WriteLine();
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + ": " + a[aIndex + i].ToString().PadLeft(14) + "");
            }
        }
         public static void r8vec_print_part(int n, double[] a, int max_print, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_PRINT_PART prints "part" of an R8VEC.
        //
        //  Discussion:
        //
        //    The user specifies MAX_PRINT, the maximum number of lines to print.
        //
        //    If N, the size of the vector, is no more than MAX_PRINT, then
        //    the entire vector is printed, one entry per line.
        //
        //    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
        //    followed by a line of periods suggesting an omission,
        //    and the last entry.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries of the vector.
        //
        //    Input, double A[N], the vector to be printed.
        //
        //    Input, int MAX_PRINT, the maximum number of lines
        //    to print.
        //
        //    Input, string TITLE, a title.
        //
        {
            int i;

            if (max_print <= 0)
            {
                return;
            }

            if (n <= 0)
            {
                return;
            }

            Console.WriteLine("");
            Console.WriteLine(title + "");
            Console.WriteLine("");

            if (n <= max_print)
            {
                for (i = 0; i < n; i++)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                        + "  " + a[i].ToString().PadLeft(14) + "");
                }
            }
            else if (3 <= max_print)
            {
                for (i = 0; i < max_print - 2; i++)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + a[i].ToString().PadLeft(14) + "");
                }

                Console.WriteLine("  ........  ..............");
                i = n - 1;
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + "  " + a[i].ToString().PadLeft(14) + "");
            }
            else
            {
                for (i = 0; i < max_print - 1; i++)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + a[i].ToString().PadLeft(14) + "");
                }

                i = max_print - 1;
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                    + ": " + a[i].ToString().PadLeft(14)
                    + "  " + "...more entries...");
            }

            return;
        }
      
    }
}