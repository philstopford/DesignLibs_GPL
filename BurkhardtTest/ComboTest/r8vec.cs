using System;
using Burkardt.Types;

namespace ComboTest;

internal partial class Program
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
        int found_num;
        int i;
        int indx;
        int k;
        int n = 8;
        int maxstack = 100;
        int[] ncan = new int[8];
        int nstack;
        double[] stacks = new double[100];
        double t;
        double total;
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

        t = 53.0;

        for (i = 0; i < n; i++)
        {
            x[i] = 0.0;
        }

        indx = 0;
        k = 0;
        nstack = 0;
        for (i = 0; i < n; i++)
        {
            ncan[i] = 0;
        }

        found_num = 0;

        for (;;)
        {
            typeMethods.r8vec_backtrack(n, maxstack, stacks, ref x, ref indx, ref k, ref nstack, ref ncan);

            if (indx == 1)
            {
                found_num += 1;
                string cout = "  " + found_num.ToString(CultureInfo.InvariantCulture).PadLeft(2) + "   ";

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
                else if (t == total)
                {
                    ncan[k - 1] += 1;
                    stacks[nstack] = 0.0;
                    nstack += 1;
                }
                else if (total < t && k < n)
                {
                    ncan[k - 1] += 1;
                    stacks[nstack] = 0.0;
                    nstack += 1;

                    if (total + w[k - 1] <= t)
                    {
                        ncan[k - 1] += 1;
                        stacks[nstack] = 1.0;
                        nstack += 1;
                    }
                }
                else if (total < t && k == n)
                {
                    if (total + w[k - 1] == t)
                    {
                        ncan[k - 1] += 1;
                        stacks[nstack] = 1.0;
                        nstack += 1;
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