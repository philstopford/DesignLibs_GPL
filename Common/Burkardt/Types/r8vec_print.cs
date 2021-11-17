using System;

namespace Burkardt.Types;

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

    public static void r8vec_print_part(int n, double[] a, int max_print, string title)

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

        switch (max_print)
        {
            case <= 0:
                return;
        }

        switch (n)
        {
            case <= 0:
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
        else
        {
            switch (max_print)
            {
                case >= 3:
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
                    break;
                }
                default:
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
                    break;
                }
            }
        }
    }

    public static void r8vec_print_16 ( int n, double[] a, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_PRINT_16 prints an R8VEC to 16 decimal places.
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
        //    29 May 2014
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
        int i;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");
        for ( i = 0; i < n; i++ )
        {
            Console.WriteLine("  " + i.ToString().PadLeft(8)
                                   + ": " + a[i].ToString("0.################").PadLeft(24)  + "");
        }
    }

    public static void r8vec_print_some ( int n, double[] a, int max_print, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_PRINT_SOME prints "some" of an R8VEC.
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

        switch (max_print)
        {
            case <= 0:
                return;
        }

        switch (n)
        {
            case <= 0:
                return;
        }

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");

        if ( n <= max_print )
        {
            for ( i = 0; i < n; i++ )
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + "  " + a[i].ToString().PadLeft(14) + "");
            }
        }
        else
        {
            switch (max_print)
            {
                case >= 3:
                {
                    for ( i = 0; i < max_print - 2; i++ )
                    {
                        Console.WriteLine("  " + i.ToString().PadLeft(8)
                                               + ": " + a[i].ToString().PadLeft(14) + "");
                    }
                    Console.WriteLine("  ........  ..............");
                    i = n - 1;
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + ": " + a[i].ToString().PadLeft(14) + "");
                    break;
                }
                default:
                {
                    for ( i = 0; i < max_print - 1; i++ )
                    {
                        Console.WriteLine("  " + i.ToString().PadLeft(8)
                                               + ": " + a[i].ToString().PadLeft(14) + "");
                    }
                    i = max_print - 1;
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + ": " + a[i].ToString().PadLeft(14)
                                           + "  " + "...more entries...");
                    break;
                }
            }
        }
    }
        
    public static void r8vec_print_some ( int n, double[] a, int i_lo, int i_hi, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_PRINT_SOME prints "some" of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8 values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2006
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
        //    Input, integer I_LO, I_HI, the first and last indices to print.
        //    The routine expects 1 <= I_LO <= I_HI <= N.
        //
        //    Input, string TITLE, a title.
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");
        for ( i = Math.Max ( 1, i_lo ); i <= Math.Min ( n, i_hi ); i++ )
        {
            Console.WriteLine("  " + i.ToString().PadLeft(8)       + "  "
                              + "  " + a[i-1].ToString().PadLeft(14)  + "");
        }
    }
        
    public static void r8vec_mask_print ( int n, double[] a, int mask_num, int[] mask,
            string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_MASK_PRINT prints a masked R8VEC.
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
        //    19 September 2005
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
        //    Input, int MASK_NUM, the number of masked elements.
        //
        //    Input, int MASK[MASK_NUM], the indices of the vector to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("  Masked vector printout:");

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");
        for ( i = 0; i < mask_num; i++ )
        {
            Console.WriteLine("  " + i.ToString().PadLeft(6)
                                   + ": " + mask[i].ToString().PadLeft(6)
                                   + "  " + a[mask[i]-1].ToString().PadLeft(12) + "");
        }
    }

    public static void r8vec2_print_some(int n, double[] x1, double[] x2, int max_print,
            string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC2_PRINT_SOME prints "some" of an R8VEC2.
        //
        //  Discussion:
        //
        //    An R8VEC2 is a dataset consisting of N pairs of real values, stored
        //    as two separate vectors A1 and A2.
        //
        //    The user specifies MAX_PRINT, the maximum number of lines to print.
        //
        //    If N, the size of the vectors, is no more than MAX_PRINT, then
        //    the entire vectors are printed, one entry of each per line.
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
        //    13 November 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries of the vectors.
        //
        //    Input, double X1[N], X2[N], the vector to be printed.
        //
        //    Input, int MAX_PRINT, the maximum number of lines to print.
        //
        //    Input, string TITLE, a title.
        //
    {
        int i;

        switch (max_print)
        {
            case <= 0:
                return;
        }

        switch (n)
        {
            case <= 0:
                return;
        }

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");

        if (n <= max_print)
        {
            for (i = 0; i < n; i++)
            {
                Console.WriteLine(i.ToString().PadLeft(6) + ": "
                                                          + x1[i].ToString().PadLeft(14) + "  "
                                                          + x2[i].ToString().PadLeft(14) + "");
            }
        }
        else
        {
            switch (max_print)
            {
                case >= 3:
                {
                    for (i = 0; i < max_print - 2; i++)
                    {
                        Console.WriteLine(i.ToString().PadLeft(6) + ": "
                                                                  + x1[i].ToString().PadLeft(14) + "  "
                                                                  + x2[i].ToString().PadLeft(14) + "");
                    }

                    Console.WriteLine("......  ..............  ..............");
                    i = n - 1;
                    Console.WriteLine(i.ToString().PadLeft(6) + ": "
                                                              + x1[i].ToString().PadLeft(14) + "  "
                                                              + x2[i].ToString().PadLeft(14) + "");
                    break;
                }
                default:
                {
                    for (i = 0; i < max_print - 1; i++)
                    {
                        Console.WriteLine(i.ToString().PadLeft(6) + ": "
                                                                  + x1[i].ToString().PadLeft(14) + "  "
                                                                  + x2[i].ToString().PadLeft(14) + "");
                    }

                    i = max_print - 1;
                    Console.WriteLine(i.ToString().PadLeft(6) + ": "
                                                              + x1[i].ToString().PadLeft(14) + "  "
                                                              + x2[i].ToString().PadLeft(14) + "...more entries...");
                    break;
                }
            }
        }
    }
}