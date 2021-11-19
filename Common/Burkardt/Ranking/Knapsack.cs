using System;
using Burkardt.Types;

namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static void knapsack_01(int n, double mass_limit, ref double[] p, ref double[] w, ref double[] x,
            ref double mass, ref double profit )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KNAPSACK_01 solves the 0/1 knapsack problem.
        // 
        //  Discussion:
        // 
        //    The 0/1 knapsack problem is as follows:
        // 
        //      Given:
        //        a set of N objects,
        //        a profit P(I) and weight W(I) associated with each object,
        //        and a weight limit MASS_LIMIT,
        //      Determine:
        //        a set of choices X(I) which are 0 or 1, that maximizes the profit
        //          P = Sum ( 1 <= I <= N ) P(I) * X(I)
        //        subject to the constraint
        //          Sum ( 1 <= I <= N ) W(I) * X(I) <= MASS_LIMIT.
        // 
        //    This routine assumes that the objects have already been sorted
        //    in order of decreasing "profit density", P(I)/W(I).
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    28 July 2011
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Reference:
        // 
        //    Donald Kreher, Douglas Simpson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998,
        //    ISBN: 0-8493-3988-X,
        //    LC: QA164.K73.
        // 
        //  Parameters:
        // 
        //    Input, int N, the number of objects.
        // 
        //    Input, double MASS_LIMIT, the weight limit of the
        //    chosen objects.
        // 
        //    Input/output, double P[N], the "profit" or value of each object.
        //    P is assumed to be nonnegative.
        // 
        //    Input/output, double W[N], the "weight" or cost of each object.
        //    W is assumed to be  nonnegative.
        // 
        //    Output, double X[N], the choice function for the objects.
        //    0, the object was not taken.
        //    1, the object was taken.
        // 
        //    Output, double &MASS, the total mass of the objects taken.
        // 
        //    Output, double &PROFIT, the total profit of the objects taken.
        // 
    {
        int i;
        int k = 0;
        double mass_2 = 0;
        const int maxstack = 100;
        double profit_2 = 0;

        int[] ncan = new int[n];
        double[] stack = new double[maxstack];
        double[] x_best = new double[n];

        int nstack = 0;
        // 
        //  Initialize the "best so far" data.
        // 
        for (i = 0; i < n; i++)
        {
            x_best[i] = 0.0;
        }

        double profit_best = 0.0;
        double mass_best = 0;
        // 
        //  Begin the backtracking solution.
        // 
        int indx = 0;

        for (;;)
        {
            typeMethods.r8vec_backtrack(n, maxstack, stack, ref x, ref indx, ref k, ref nstack, ref ncan);
            // 
            //  Got a new candidate.  Compare it to the best so far.
            // 
            if (indx == 1)
            {
                profit = typeMethods.r8vec_dot_product(n, p, x);
                mass = typeMethods.r8vec_dot_product(n, w, x);

                if (profit_best < profit || Math.Abs(profit - profit_best) <= double.Epsilon && mass < mass_best)
                {
                    profit_best = profit;
                    mass_best = mass;
                    for (i = 0; i < n; i++)
                    {
                        x_best[i] = x[i];
                    }
                }
            }
            // 
            //  Need candidates for X(K).
            // 
            //  X(K) = 1 is possible if:
            // 
            //    * adding W(K) to our mass doesn''t put us over our mass limit;
            //    * and adding P(K) to our current profit, and taking the best we
            //      could get using rational X for the remainder would put us over
            //      our current best.
            // 
            //  X(K) = 0 is always possible.
            // 
            else if (indx == 2)
            {
                ncan[k - 1] = 0;

                double mass_1 = w[k - 1];
                for (i = 0; i < k - 1; i++)
                {
                    mass_1 += w[i] * x[i];
                }

                if (mass_1 <= mass_limit)
                {
                    double mass_remaining = mass_limit - mass_1;

                    double profit_1 = p[k - 1];
                    for (i = 0; i < k - 1; i++)
                    {
                        profit_1 += p[i] * x[i];
                    }

                    if (k < n)
                    {
                        knapsack_rational(n - k, mass_remaining, p, w,
                            ref x, ref mass_2, ref profit_2, pIndex:k, wIndex:k, xIndex:k);
                    }
                    else
                    {
                        profit_2 = 0.0;
                    }

                    if (profit_best < profit_1 + profit_2)
                    {
                        if (maxstack <= nstack)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("KNAPSACK_01 - Fatal error!");
                            Console.WriteLine("  Exceeded stack space.");
                            return;
                        }

                        ncan[k - 1] += 1;
                        nstack += 1;
                        stack[nstack - 1] = 1.0;
                    }
                }

                if (maxstack <= nstack)
                {
                    Console.WriteLine("");
                    Console.WriteLine("KNAPSACK_01 - Fatal error!");
                    Console.WriteLine("  Exceeded stack space.");
                    return;
                }

                ncan[k - 1] += 1;
                nstack += 1;
                stack[nstack - 1] = 0.0;
            }
            // 
            //  Done.  Return the best solution.
            // 
            else
            {
                profit = profit_best;
                mass = mass_best;
                for (i = 0; i < n; i++)
                {
                    x[i] = x_best[i];
                }

                break;
            }
        }
    }

    public static void knapsack_rational(int n, double mass_limit, double[] p, double[] w,
            ref double[] x, ref double mass, ref double profit, int pIndex = 0, int wIndex = 0, int xIndex = 0 )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KNAPSACK_RATIONAL solves the rational knapsack problem.
        // 
        //  Discussion:
        // 
        //    The rational knapsack problem is a generalization of the 0/1 knapsack
        //    problem.  It is mainly used to derive a bounding function for the
        //    0/1 knapsack problem.
        // 
        //    The 0/1 knapsack problem is as follows:
        // 
        //      Given:
        //        a set of N objects,
        //        a profit P(I) and weight W(I) associated with each object,
        //        and a weight limit MASS_LIMIT,
        //      Determine:
        //        a set of choices X(I) which are 0 or 1, that maximizes the profit
        //          P = Sum ( 1 <= I <= N ) P(I) * X(I)
        //        subject to the constraint
        //          Sum ( 1 <= I <= N ) W(I) * X(I) <= MASS_LIMIT.
        // 
        //    By contrast, the rational knapsack problem allows the values X(I)
        //    to be any value between 0 and 1.  A solution for the rational knapsack
        //    problem is known.  Arrange the objects in order of their "profit density"
        //    ratios P(I)/W(I), and then take in order as many of these as you can.
        //    If you still have "room" in the weight constraint, then you should
        //    take the maximal fraction of the very next object, which will complete
        //    your weight limit, and maximize your profit.
        // 
        //    If should be obvious that, given the same data, a solution for
        //    the rational knapsack problem will always have a profit that is
        //    at least as high as for the 0/1 problem.  Since the rational knapsack
        //    maximum profit is easily computed, this makes it a useful bounding
        //    function.
        // 
        //    Note that this routine assumes that the objects have already been
        //    arranged in order of the "profit density".
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    25 July 2011
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Reference:
        // 
        //    Donald Kreher, Douglas Simpson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998,
        //    ISBN: 0-8493-3988-X,
        //    LC: QA164.K73.
        // 
        //  Parameters:
        // 
        //    Input, int N, the number of objects.
        // 
        //    Input, double MASS_LIMIT, the weight limit of the
        //    chosen objects.
        // 
        //    Input, double P[N], the "profit" or value of each object.
        //    The entries of P are assumed to be nonnegative.
        // 
        //    Input, double W[N], the "weight" or cost of each object.
        //    The entries of W are assumed to be nonnegative.
        // 
        //    Output, double X[N], the choice function for the objects.
        //    0.0, the object was not taken.
        //    1.0, the object was taken.
        //    R, where 0 < R < 1, a fractional amount of the object was taken.
        // 
        //    Output, double &MASS, the total mass of the objects taken.
        // 
        //    Output, double &PROFIT, the total profit of the objects taken.
        // 
    {
        int i;

        mass = 0.0;
        profit = 0.0;

        for (i = 0; i < n; i++)
        {
            if (mass_limit <= mass)
            {
                x[i + xIndex] = 0.0;
            }
            else if (mass + w[i] <= mass_limit)
            {
                x[i + xIndex] = 1.0;
                mass += w[i + wIndex];
                profit += p[i + pIndex];
            }
            else
            {
                x[i + xIndex] = (mass_limit - mass) / w[i + wIndex];
                mass = mass_limit;
                profit += p[i + pIndex] * x[i + xIndex];
            }
        }
    }

    public static void knapsack_reorder(int n, ref double[] p, ref double[] w )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KNAPSACK_REORDER reorders the knapsack data by "profit density".
        // 
        //  Discussion:
        // 
        //    This routine must be called to rearrange the data before calling
        //    routines that handle a knapsack problem.
        // 
        //    The "profit density" for object I is defined as P(I)/W(I).
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    26 July 2011
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Reference:
        // 
        //    Donald Kreher, Douglas Simpson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998,
        //    ISBN: 0-8493-3988-X,
        //    LC: QA164.K73.
        // 
        //  Parameters:
        // 
        //    Input, int N, the number of objects.
        // 
        //    Input/output, double P[N], the "profit" or value of each object.
        // 
        //    Input/output, double W[N], the "weight" or cost of each object.
        // 
    {
        int i;
        // 
        //  Rearrange the objects in order of "profit density".
        //
        for (i = 0; i < n; i++)
        {
            int j;
            for (j = i + 1; j < n; j++)
            {
                if (p[i] * w[j] < p[j] * w[i])
                {
                    double t = p[i];
                    p[i] = p[j];
                    p[j] = t;

                    t = w[i];
                    w[i] = w[j];
                    w[j] = t;
                }
            }
        }
    }
}