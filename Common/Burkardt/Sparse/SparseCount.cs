using System;
using Burkardt.Composition;
using Burkardt.Types;

namespace Burkardt.Sparse;

public class SparseCount
{
    public static int cc_se_size(int dim_num, int level_max)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CC_SE_SIZE: Clenshaw Curtis Slow Exponential Growth.
        //
        //  Discussion:
        //
        //    The grid is defined as the sum of the product rules whose LEVEL
        //    satisfies:
        //
        //      0 <= LEVEL <= LEVEL_MAX.
        //
        //    This calculation is much faster than a previous method.  It simply
        //    computes the number of new points that are added at each level in the
        //    1D rule, and then counts the new points at a given DIM_NUM dimensional
        //    level vector as the product of the new points added in each dimension.
        //
        //    This approach will work for nested families, and may be extensible
        //    to other families, and to mixed rules.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the maximum value of LEVEL.
        //
        //    Output, int CC_SE_SIZE, the number of points in the grid.
        //
    {
        int point_num;
        switch (level_max)
        {
            //
            //  Special case.
            //
            case < 0:
                point_num = 0;
                return point_num;
            case 0:
                point_num = 1;
                return point_num;
        }

        //
        //  Construct the vector that counts the new points in the 1D rule.
        //
        int[] new_1d = new int[level_max + 1];

        new_1d[0] = 1;
        new_1d[1] = 2;

        int o = 3;

        for (int l = 2; l <= level_max; l++)
        {
            int p = 2 * l + 1;
            if (o < p)
            {
                new_1d[l] = o - 1;
                o = 2 * o - 1;
            }
            else
            {
                new_1d[l] = 0;
            }
        }

        //
        //  Count the number of points by counting the number of new points 
        //  associated with each level vector.
        //
        int[] level_1d = new int[dim_num];

        point_num = 0;

        for (int level = 0; level <= level_max; level++)
        {
            bool more = false;
            int h = 0;
            int t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);

                int v = 1;
                for (int dim = 0; dim < dim_num; dim++)
                {
                    v *= new_1d[level_1d[dim]];
                }

                point_num += v;

                if (!more)
                {
                    break;
                }
            }
        }

        return point_num;
    }

    public static int cfn_e_size(int dim_num, int level_max)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CFN_E_SIZE: Closed Fully Nested Exponential Growth.
        //
        //  Discussion:
        //
        //    The grid is defined as the sum of the product rules whose LEVEL
        //    satisfies:
        //
        //      0 <= LEVEL <= LEVEL_MAX.
        //
        //    This calculation is much faster than a previous method.  It simply
        //    computes the number of new points that are added at each level in the
        //    1D rule, and then counts the new points at a given DIM_NUM dimensional
        //    level vector as the product of the new points added in each dimension.
        //
        //    This approach will work for nested families, and may be extensible
        //    to other families, and to mixed rules.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the maximum value of LEVEL.
        //
        //    Output, int CFN_E_SIZE, the number of points in the grid.
        //
    {
        int point_num;
        switch (level_max)
        {
            //
            //  Special case.
            //
            case < 0:
                point_num = 0;
                return point_num;
            case 0:
                point_num = 1;
                return point_num;
        }

        //
        //  Construct the vector that counts the new points in the 1D rule.
        //
        int[] new_1d = new int[level_max + 1];

        new_1d[0] = 1;
        new_1d[1] = 2;

        int j = 1;
        for (int l = 2; l <= level_max; l++)
        {
            j *= 2;
            new_1d[l] = j;
        }

        //
        //  Count the number of points by counting the number of new points 
        //  associated with each level vector.
        //
        int[] level_1d = new int[dim_num];

        point_num = 0;

        for (int level = 0; level <= level_max; level++)
        {
            bool more = false;
            int h = 0;
            int t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);

                int v = 1;
                for (int dim = 0; dim < dim_num; dim++)
                {
                    v *= new_1d[level_1d[dim]];
                }

                point_num += v;

                if (!more)
                {
                    break;
                }
            }
        }

        return point_num;
    }

    public static int f2_se_size(int dim_num, int level_max)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F2_SE_SIZE: Fejer Type 2 Slow Exponential Growth.
        //
        //  Discussion:
        //
        //    The grid is defined as the sum of the product rules whose LEVEL
        //    satisfies:
        //
        //      0 <= LEVEL <= LEVEL_MAX.
        //
        //    This calculation is much faster than a previous method.  It simply
        //    computes the number of new points that are added at each level in the
        //    1D rule, and then counts the new points at a given DIM_NUM dimensional
        //    level vector as the product of the new points added in each dimension.
        //
        //    This approach will work for nested families, and may be extensible
        //    to other families, and to mixed rules.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the maximum value of LEVEL.
        //
        //    Output, int F2_SE_SIZE, the number of points in the grid.
        //
    {
        int point_num;
        switch (level_max)
        {
            //
            //  Special case.
            //
            case < 0:
                point_num = 0;
                return point_num;
            case 0:
                point_num = 1;
                return point_num;
        }

        //
        //  Construct the vector that counts the new points in the 1D rule.
        //
        int[] new_1d = new int[level_max + 1];

        new_1d[0] = 1;

        int o = 1;

        for (int l = 1; l <= level_max; l++)
        {
            int p = 2 * l + 1;
            if (o < p)
            {
                new_1d[l] = o + 1;
                o = 2 * o + 1;
            }
            else
            {
                new_1d[l] = 0;
            }
        }

        //
        //  Count the number of points by counting the number of new points 
        //  associated with each level vector.
        //
        int[] level_1d = new int[dim_num];

        point_num = 0;

        for (int level = 0; level <= level_max; level++)
        {
            bool more = false;
            int h = 0;
            int t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);

                int v = 1;
                for (int dim = 0; dim < dim_num; dim++)
                {
                    v *= new_1d[level_1d[dim]];
                }

                point_num += v;

                if (!more)
                {
                    break;
                }
            }
        }

        return point_num;
    }

    public static int gp_se_size(int dim_num, int level_max)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GP_SE_SIZE: Gauss Patterson Slow Exponential Growth.
        //
        //  Discussion:
        //
        //    The Gauss-Patterson-Slow family assumes that, for the underlying 1D
        //    rules, a precision of 2*L+1 is needed at level L.  Therefore, the
        //    lowest possible order Gauss-Patterson rule is chosen that will achieve
        //    that precision.  This retains a combination of the advantages of
        //    nestedness and high accuracy.
        //
        //    The grid is defined as the sum of the product rules whose LEVEL
        //    satisfies:
        //
        //      0 <= LEVEL <= LEVEL_MAX.
        //
        //    This calculation is much faster than a previous method.  It simply
        //    computes the number of new points that are added at each level in the
        //    1D rule, and then counts the new points at a given DIM_NUM dimensional
        //    level vector as the product of the new points added in each dimension.
        //
        //    This approach will work for nested families, and may be extensible
        //    to other families, and to mixed rules.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the maximum value of LEVEL.
        //
        //    Output, int GP_SE_SIZE, the number of points in the grid.
        //
    {
        int point_num;
        switch (level_max)
        {
            //
            //  Special case.
            //
            case < 0:
                point_num = 0;
                return point_num;
            case 0:
                point_num = 1;
                return point_num;
        }

        //
        //  Count the points in the 1D rule.
        //
        int[] order_1d = new int[level_max + 1];
        order_1d[0] = 1;
        for (int level = 1; level <= level_max; level++)
        {
            int p = 5;
            int o = 3;
            while (p < 2 * level + 1)
            {
                p = 2 * p + 1;
                o = 2 * o + 1;
            }

            order_1d[level] = o;
        }

        //
        //  Count the new points in the 1D rule.
        //
        int[] new_1d = new int[level_max + 1];

        new_1d[0] = 1;
        for (int level = 1; level <= level_max; level++)
        {
            new_1d[level] = order_1d[level] - order_1d[level - 1];
        }

        //
        //  Count the number of points by counting the number of new points 
        //  associated with each level vector.
        //
        int[] level_1d = new int[dim_num];

        point_num = 0;

        for (int level = 0; level <= level_max; level++)
        {
            bool more = false;

            int h = 0;
            int t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);

                int v = 1;
                for (int dim = 0; dim < dim_num; dim++)
                {
                    v *= new_1d[level_1d[dim]];
                }

                point_num += v;

                if (!more)
                {
                    break;
                }
            }
        }

        return point_num;
    }

    public static int ofn_e_size(int dim_num, int level_max)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    OFN_E_SIZE: Open Fully Nested Exponential Growth.
        //
        //  Discussion:
        //
        //    This calculation assumes that an exponential growth rule is being used,
        //    that is, that the 1D rules have orders 1, 3, 7, 15, 31, and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the maximum value of LEVEL.
        //
        //    Output, int OFN_E_SIZE, the number of points in the grid.
        //
    {
        int point_num;
        switch (level_max)
        {
            //
            //  Special case.
            //
            case < 0:
                point_num = 0;
                return point_num;
            case 0:
                point_num = 1;
                return point_num;
        }

        //
        //  Construct the vector that counts the new points in the 1D rule.
        //
        int[] new_1d = new int[level_max + 1];

        new_1d[0] = 1;
        for (int l = 1; l <= level_max; l++)
        {
            new_1d[l] = 2 * new_1d[l - 1];
        }

        //
        //  Count the number of points by counting the number of new points 
        //  associated with each level vector.
        //
        int[] level_1d = new int[dim_num];

        point_num = 0;

        for (int level = 0; level <= level_max; level++)
        {
            bool more = false;
            int h = 0;
            int t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);

                int v = 1;
                for (int dim = 0; dim < dim_num; dim++)
                {
                    v *= new_1d[level_1d[dim]];
                }

                point_num += v;

                if (!more)
                {
                    break;
                }
            }
        }

        return point_num;
    }

    public static int onn_e_size(int dim_num, int level_max)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ONN_E_SIZE Open Non-Nested Exponential Growth.
        //
        //  Discussion:
        //
        //    This calculation assumes that an exponential growth rule is being used,
        //    that is, that the 1D rules have orders 1, 3, 7, 15, 31, and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the maximum value of LEVEL.
        //
        //    Output, int ONN_E_SIZE, the number of points in the grid.
        //
    {
        int point_num;
        switch (level_max)
        {
            //
            //  Special case.
            //
            case < 0:
                point_num = 0;
                return point_num;
            case 0:
                point_num = 1;
                return point_num;
        }

        //
        //  Construct the 1D order vector.
        //
        int[] order_1d = new int[level_max + 1];

        int temp = 2;
        for (int l = 0; l <= level_max; l++)
        {
            order_1d[l] = temp - 1;
            temp *= 2;
        }

        int[] level_1d = new int[dim_num];

        point_num = 0;

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        for (int level = level_min; level <= level_max; level++)
        {
            bool more = false;
            int h = 0;
            int t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);

                int v = 1;
                for (int dim = 0; dim < dim_num; dim++)
                {
                    v *= order_1d[level_1d[dim]];
                }

                point_num += v;

                if (!more)
                {
                    break;
                }
            }
        }

        return point_num;
    }

    public static int onn_l_size(int dim_num, int level_max)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ONN_L_SIZE Open Non-Nested Linear Growth.
        //
        //  Discussion:
        //
        //    This calculation assumes that a linear growth rule is being used,
        //    that is, that the 1D rules have orders 1, 3, 5, 7, 9, and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the maximum value of LEVEL.
        //
        //    Output, int ONN_L_SIZE, the number of points in the grid.
        //
    {
        int point_num;
        switch (level_max)
        {
            //
            //  Special case.
            //
            case < 0:
                point_num = 0;
                return point_num;
            case 0:
                point_num = 1;
                return point_num;
        }

        //
        //  Construct the 1D order vector.
        //
        int[] order_1d = new int[level_max + 1];

        for (int l = 0; l <= level_max; l++)
        {
            order_1d[l] = 2 * l + 1;
        }

        int[] level_1d = new int[dim_num];

        point_num = 0;

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        for (int level = level_min; level <= level_max; level++)
        {
            bool more = false;
            int h = 0;
            int t = 0;

            for (;;)
            {
                Comp.comp_next(level, dim_num, ref level_1d, ref more, ref h, ref t);

                int v = 1;
                for (int dim = 0; dim < dim_num; dim++)
                {
                    v *= order_1d[level_1d[dim]];
                }

                point_num += v;

                if (!more)
                {
                    break;
                }
            }
        }

        return point_num;
    }

    public static int own_e_size(int dim_num, int level_max)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    OWN_E_SIZE: Open Weakly Nested Exponential Growth.
        //
        //  Discussion:
        //
        //    This calculation assumes that an exponential growth rule is being used,
        //    that is, that the 1D rules have orders 1, 3, 7, 15, 31, and so on.
        //
        //    This calculation assumes that the 1D family of quadrature rules 
        //    contains only one repeated point, presumably the value 0.0.
        //    This assumption holds for Gauss-Legendre, Gauss-Hermite and 
        //    Generalized Gauss-Hermite rules.
        //
        //    The routine then counts the number of unique abscissas that will
        //    be generated for a sparse grid of given dimension and level.
        //
        //    The computation is complicated.  It starts by counting just those
        //    abscissas which have no 0.0 in them.  This is relatively easy, since
        //    it is like counting the points in a sparse grid that uses open 
        //    non-nested rules, but for which the order of each rule is reduced by 1.
        //
        //    Then we have to count the abscissas with one 0.0, two 0.0's and so
        //    on to DIM_NUM zeros.  We are assuming this is an isotropic grid,
        //    so for a particular number K of zeroes we only need to count the case
        //    where the first K entries are zero, and multiply by C(DIM_NUM,K).
        //
        //    To count the number of entries with K zeroes, (and assuming 0 < K),
        //    then, we essentially count the number of abscissas in an open 
        //    non-nested rule as before, but modifed so that the minimum level is 0,
        //    rather than LEVEL_MAX - DIM_NUM + 1.
        //
        //    I will mention that this was a rather difficult computation to
        //    figure out!
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the maximum value of LEVEL.
        //
        //    Output, int OWN_E_SIZE, the number of points in the grid.
        //
    {
        int point_num;
        switch (level_max)
        {
            //
            //  Special case.
            //
            case < 0:
                point_num = 0;
                return point_num;
            case 0:
                point_num = 1;
                return point_num;
        }

        //
        //  Construct the vector that counts the new points in the 1D rule.
        //
        int[] new_1d = new int[level_max + 1];

        new_1d[0] = 0;
        int temp = 4;
        for (int l = 1; l <= level_max; l++)
        {
            new_1d[l] = temp - 2;
            temp *= 2;
        }

        //
        //  Count the nonzero points in the full dimensional table with the usual
        //  LEVEL_MIN restriction.
        //
        //  Then count the points with 1, 2, 3, ... DIM_NUM zeroes, by counting
        //  the nonzero points in a DIM_NUM2 table, with LEVEL_MIN set to 0, and
        //  multiplying by the appropriate combinatorial coefficient.
        //
        point_num = 0;

        for (int dim_num2 = dim_num; 0 <= dim_num2; dim_num2--)
        {
            int level_min = dim_num2 == dim_num ? Math.Max(0, level_max - dim_num + 1) : 0;

            int point_num2;
            switch (dim_num2)
            {
                case 0:
                    point_num2 = 1;
                    break;
                default:
                {
                    int[] level_1d = new int[dim_num2];

                    point_num2 = 0;

                    for (int level = level_min; level <= level_max; level++)
                    {
                        bool more = false;
                        int h = 0;
                        int t = 0;

                        for (;;)
                        {
                            Comp.comp_next(level, dim_num2, ref level_1d, ref more, ref h, ref t);

                            int v = 1;
                            for (int dim = 0; dim < dim_num2; dim++)
                            {
                                v *= new_1d[level_1d[dim]];
                            }

                            point_num2 += v;

                            if (!more)
                            {
                                break;
                            }
                        }
                    }

                    break;
                }
            }

            point_num += typeMethods.i4_choose(dim_num, dim_num2) * point_num2;
        }

        return point_num;
    }

    public static int own_l2_size(int dim_num, int level_max)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    OWN_L2_SIZE: Open Weakly Nested Linear 2 Growth.
        //
        //  Discussion:
        //
        //    This calculation assumes that a linear growth rule of size 2 is being used,
        //    that is, that the 1D rules have orders 1, 3, 5, 7, 9, and so on.
        //
        //    This calculation assumes that the 1D family of quadrature rules 
        //    contains only one repeated point, presumably the value 0.0.
        //    This assumption holds for Gauss-Legendre, Gauss-Hermite and 
        //    Generalized Gauss-Hermite rules.
        //
        //    The routine then counts the number of unique abscissas that will
        //    be generated for a sparse grid of given dimension and level.
        //
        //    The computation is complicated.  It starts by counting just those
        //    abscissas which have no 0.0 in them.  This is relatively easy, since
        //    it is like counting the points in a sparse grid that uses open 
        //    non-nested rules, but for which the order of each rule is reduced by 1.
        //
        //    Then we have to count the abscissas with one 0.0, two 0.0's and so
        //    on to DIM_NUM zeros.  We are assuming this is an isotropic grid,
        //    so for a particular number K of zeroes we only need to count the case
        //    where the first K entries are zero, and multiply by C(DIM_NUM,K).
        //
        //    To count the number of entries with K zeroes, (and assuming 0 < K),
        //    then, we essentially count the number of abscissas in an open 
        //    non-nested rule as before, but modifed so that the minimum level is 0,
        //    rather than LEVEL_MAX - DIM_NUM + 1.
        //
        //    I will mention that this was a rather difficult computation to
        //    figure out!
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the maximum value of LEVEL.
        //
        //    Output, int OWN_L_SIZE, the number of points in the grid.
        //
    {
        int point_num;
        switch (level_max)
        {
            //
            //  Special case.
            //
            case < 0:
                point_num = 0;
                return point_num;
            case 0:
                point_num = 1;
                return point_num;
        }

        //
        //  Construct the vector that counts the new points in the 1D rule.
        //
        int[] new_1d = new int[level_max + 1];

        new_1d[0] = 0;
        for (int l = 1; l <= level_max; l++)
        {
            new_1d[l] = 2 * l;
        }

        //
        //  Count the nonzero points in the full dimensional table with the usual
        //  LEVEL_MIN restriction.
        //
        //  Then count the points with 1, 2, 3, ... DIM_NUM zeroes, by counting
        //  the nonzero points in a DIM_NUM2 table, with LEVEL_MIN set to 0, and
        //  multiplying by the appropriate combinatorial coefficient.
        //
        point_num = 0;

        for (int dim_num2 = dim_num; 0 <= dim_num2; dim_num2--)
        {
            int level_min = dim_num2 == dim_num ? Math.Max(0, level_max - dim_num + 1) : 0;

            int point_num2;
            switch (dim_num2)
            {
                case 0:
                    point_num2 = 1;
                    break;
                default:
                {
                    int[] level_1d = new int[dim_num2];

                    point_num2 = 0;

                    for (int level = level_min; level <= level_max; level++)
                    {
                        bool more = false;
                        int h = 0;
                        int t = 0;

                        for (;;)
                        {
                            Comp.comp_next(level, dim_num2, ref level_1d, ref more, ref h, ref t);

                            int v = 1;
                            for (int dim = 0; dim < dim_num2; dim++)
                            {
                                v *= new_1d[level_1d[dim]];
                            }

                            point_num2 += v;

                            if (!more)
                            {
                                break;
                            }
                        }
                    }

                    break;
                }
            }

            point_num += typeMethods.i4_choose(dim_num, dim_num2) * point_num2;
        }

        return point_num;
    }

    public static int own_o_size(int dim_num, int level_max)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    OWN_O_SIZE: Open Weakly Nested Odd Growth.
        //
        //  Discussion:
        //
        //    This calculation assumes that an odd growth rule is being used,
        //    that is, that the 1D rules have orders 1, 3, 3, 5, 5, 7, 7, 9, 9, 
        //    and so on.
        //
        //    This repetition of rules is permissible because the sparse grid rule
        //    only requires that the 1D rule of level L have precision at least 2*L+1.
        //    This is achievable for Gaussian rules.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the maximum value of LEVEL.
        //
        //    Output, int OWN_O_SIZE, the number of points in the grid.
        //
    {
        int point_num;
        switch (level_max)
        {
            //
            //  Special case.
            //
            case < 0:
                point_num = 0;
                return point_num;
            case 0:
                point_num = 1;
                return point_num;
        }

        switch (dim_num)
        {
            case 1:
                point_num = 1 + 2 * ((level_max + 1) / 2);
                return point_num;
        }

        //
        //  Construct the vector that counts the new points in the 1D rule.
        //
        int[] new_1d = new int[level_max + 1];

        for (int l = 0; l <= level_max; l++)
        {
            new_1d[l] = 0;
        }

        for (int l = 1; l <= level_max; l += 2)
        {
            new_1d[l] = l + 1;
        }

        //
        //  Count the nonzero points in the full dimensional table with the usual
        //  LEVEL_MIN restriction.
        //
        //  Then count the points with 1, 2, 3, ... DIM_NUM zeroes, by counting
        //  the nonzero points in a DIM_NUM2 table, with LEVEL_MIN set to 0, and
        //  multiplying by the appropriate combinatorial coefficient.
        //
        point_num = 0;

        for (int dim_num2 = dim_num; 0 <= dim_num2; dim_num2--)
        {
            int level_min = dim_num2 == dim_num ? Math.Max(0, level_max - dim_num + 1) : 0;

            int point_num2;
            switch (dim_num2)
            {
                case 0:
                    point_num2 = 1;
                    break;
                default:
                {
                    int[] level_1d = new int[dim_num2];

                    point_num2 = 0;

                    for (int level = level_min; level <= level_max; level++)
                    {
                        bool more = false;
                        int h = 0;
                        int t = 0;

                        for (;;)
                        {
                            Comp.comp_next(level, dim_num2, ref level_1d, ref more, ref h, ref t);

                            int v = 1;
                            for (int dim = 0; dim < dim_num2; dim++)
                            {
                                v *= new_1d[level_1d[dim]];
                            }

                            point_num2 += v;

                            if (!more)
                            {
                                break;
                            }
                        }
                    }

                    break;
                }
            }

            point_num += typeMethods.i4_choose(dim_num, dim_num2) * point_num2;
        }

        return point_num;
    }
}