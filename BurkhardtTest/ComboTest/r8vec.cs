using System;
using System.Globalization;
using Burkardt.Types;

namespace ComboTest;

internal static partial class Program
{
    private static void r8vec_backtrack_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_BACKTRACK_TEST tests R8VEC_BACKTRACK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int n = 8;
        const int maxstack = 100;
        int[] ncan = new int[8];
        double[] stacks = new double[100];
        double[] w =  {
                15.0, 22.0, 14.0, 26.0, 32.0, 9.0, 16.0, 8.0
            }
            ;
        double[] x = new double[8];

        Console.WriteLine("");
        Console.WriteLine("R8VEC_BACKTRACK_TEST");
        Console.WriteLine("  I4VEC_BACKTRACK uses backtracking, seeking a vector X of");
        Console.WriteLine("  N values which satisfies some condition.");

        Console.WriteLine("");
        Console.WriteLine("  In this demonstration, we have 8 values W(I).");
        Console.WriteLine("  We seek all subsets that sum to 53.0.");
        Console.WriteLine("  X(I) is 0.0 or 1.0 if the entry is skipped or used.");
        Console.WriteLine("");

        const double t = 53.0;

        for (i = 0; i < n; i++)
        {
            x[i] = 0.0;
        }

        int indx = 0;
        int k = 0;
        int nstack = 0;
        for (i = 0; i < n; i++)
        {
            ncan[i] = 0;
        }

        int found_num = 0;

        for (;;)
        {
            typeMethods.r8vec_backtrack(n, maxstack, stacks, ref x, ref indx, ref k, ref nstack, ref ncan);

            double total;
            if (indx == 1)
            {
                found_num += 1;
                string cout = "  " + found_num.ToString().PadLeft(2) + "   ";

                total = typeMethods.r8vec_dot_product(n, w, x);
                cout += "  " + total.ToString(CultureInfo.InvariantCulture).PadLeft(8) + ":  ";

                for (i = 0; i < n; i++)
                {
                    switch (x[i])
                    {
                        case 1.0:
                            cout += "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(8);
                            break;
                    }
                }

                Console.WriteLine(cout);
            }
            //
            //  Given that we've chose X(1:K-1), what are our choices for X(K)?
            //
            //     if T < TOTAL, 
            //       no choices
            //     if T = TOTAL, 
            //       X(K) = 0
            //     if T > TOTAL and K < N, 
            //       X(k) = 0
            //       If ( W(K)+TOTAL <= T ) X(K) = 1
            //     If T > TOTAL and K = N,
            //       If ( W(K) + TOTAL) = T ) X(K) = 1
            //
            else if (indx == 2)
            {
                total = typeMethods.r8vec_dot_product(k - 1, w, x);

                if (t < total)
                {
                    ncan[k - 1] = 0;
                }
                else if (Math.Abs(t - total) <= double.Epsilon)
                {
                    ncan[k - 1] += 1;
                    stacks[nstack] = 0.0;
                    nstack += 1;
                }
                else
                {
                    switch (total)
                    {
                        case < t when k < n:
                        {
                            ncan[k - 1] += 1;
                            stacks[nstack] = 0.0;
                            nstack += 1;

                            if (!(total + w[k - 1] <= t))
                            {
                                continue;
                            }

                            ncan[k - 1] += 1;
                            stacks[nstack] = 1.0;
                            nstack += 1;
                            break;
                        }
                        case < t when k == n:
                        {
                            if (Math.Abs(total + w[k - 1] - t) <= double.Epsilon)
                            {
                                ncan[k - 1] += 1;
                                stacks[nstack] = 1.0;
                                nstack += 1;
                            }

                            break;
                        }
                    }
                }
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  Done!");
                break;
            }
        }
    }
}