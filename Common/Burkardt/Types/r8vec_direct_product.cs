namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public class r8vecDPData
        {
            public int contig = 0;
            public int rep = 0;
            public int skip = 0;
        }
        
        public static void r8vec_direct_product(ref r8vecDPData data, int factor_index, int factor_order,
                double[] factor_value, int factor_num, int point_num, ref double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    To explain what is going on here, suppose we had to construct
            //    a multidimensional quadrature rule as the product of K rules
            //    for 1D quadrature.
            //
            //    The product rule will be represented as a list of points and weights.
            //
            //    The J-th item in the product rule will be associated with
            //      item J1 of 1D rule 1,
            //      item J2 of 1D rule 2,
            //      ...,
            //      item JK of 1D rule K.
            //
            //    In particular,
            //      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
            //    and
            //      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
            //
            //    So we can construct the quadrature rule if we can properly
            //    distribute the information in the 1D quadrature rules.
            //
            //    This routine carries out that task.
            //
            //    Another way to do this would be to compute, one by one, the
            //    set of all possible indices (J1,J2,...,JK), and then index
            //    the appropriate information.  An advantage of the method shown
            //    here is that you can process the K-th set of information and
            //    then discard it.
            //
            //  Example:
            //
            //    Rule 1:
            //      Order = 4
            //      X(1:4) = ( 1, 2, 3, 4 )
            //
            //    Rule 2:
            //      Order = 3
            //      X(1:3) = ( 10, 20, 30 )
            //
            //    Rule 3:
            //      Order = 2
            //      X(1:2) = ( 100, 200 )
            //
            //    Product Rule:
            //      Order = 24
            //      X(1:24) =
            //        ( 1, 10, 100 )
            //        ( 2, 10, 100 )
            //        ( 3, 10, 100 )
            //        ( 4, 10, 100 )
            //        ( 1, 20, 100 )
            //        ( 2, 20, 100 )
            //        ( 3, 20, 100 )
            //        ( 4, 20, 100 )
            //        ( 1, 30, 100 )
            //        ( 2, 30, 100 )
            //        ( 3, 30, 100 )
            //        ( 4, 30, 100 )
            //        ( 1, 10, 200 )
            //        ( 2, 10, 200 )
            //        ( 3, 10, 200 )
            //        ( 4, 10, 200 )
            //        ( 1, 20, 200 )
            //        ( 2, 20, 200 )
            //        ( 3, 20, 200 )
            //        ( 4, 20, 200 )
            //        ( 1, 30, 200 )
            //        ( 2, 30, 200 )
            //        ( 3, 30, 200 )
            //        ( 4, 30, 200 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 April 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int FACTOR_INDEX, the index of the factor being processed.
            //    The first factor processed must be factor 0.
            //
            //    Input, int FACTOR_ORDER, the order of the factor.
            //
            //    Input, double FACTOR_VALUE[FACTOR_ORDER], the factor values
            //    for factor FACTOR_INDEX.
            //
            //    Input, int FACTOR_NUM, the number of factors.
            //
            //    Input, int POINT_NUM, the number of elements in the direct product.
            //
            //    Input/output, double X[FACTOR_NUM*POINT_NUM], the elements of the
            //    direct product, which are built up gradually.
            //
            //  Local Parameters:
            //
            //    Local, int START, the first location of a block of values to set.
            //
            //    Local, int CONTIG, the number of consecutive values to set.
            //
            //    Local, int SKIP, the distance from the current value of START
            //    to the next location of a block of values to set.
            //
            //    Local, int REP, the number of blocks of values to set.
            //
        {

            int i;
            int j;
            int k;
            int start;

            if (factor_index == 0)
            {
                data.contig = 1;
                data.skip = 1;
                data.rep = point_num;
                for (j = 0; j < point_num; j++)
                {
                    for (i = 0; i < factor_num; i++)
                    {
                        x[i + j * factor_num] = 0.0;
                    }
                }
            }

            data.rep = data.rep / factor_order;data.
            skip = data.skip * factor_order;

            for (i = 0; i < factor_order; i++)
            {
                start = 0 + i * data.contig;

                for (k = 1; k <= data.rep; k++)
                {
                    for (j = start; j < start + data.contig; j++)
                    {
                        x[factor_index + j * factor_num] = factor_value[i];
                    }

                    start = start + data.skip;
                }
            }

            data.contig = data.contig * factor_order;
        }

        public static void r8vec_direct_product2(ref r8vecDPData data, int factor_index, int factor_order,
                double[] factor_value, int factor_num, int point_num, ref double[] w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    To explain what is going on here, suppose we had to construct
            //    a multidimensional quadrature rule as the product of K rules
            //    for 1D quadrature.
            //
            //    The product rule will be represented as a list of points and weights.
            //
            //    The J-th item in the product rule will be associated with
            //      item J1 of 1D rule 1,
            //      item J2 of 1D rule 2,
            //      ...,
            //      item JK of 1D rule K.
            //
            //    In particular,
            //      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
            //    and
            //      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
            //
            //    So we can construct the quadrature rule if we can properly
            //    distribute the information in the 1D quadrature rules.
            //
            //    This routine carries out that task for the weights W.
            //
            //    Another way to do this would be to compute, one by one, the
            //    set of all possible indices (J1,J2,...,JK), and then index
            //    the appropriate information.  An advantage of the method shown
            //    here is that you can process the K-th set of information and
            //    then discard it.
            //
            //  Example:
            //
            //    Rule 1:
            //      Order = 4
            //      W(1:4) = ( 2, 3, 5, 7 )
            //
            //    Rule 2:
            //      Order = 3
            //      W(1:3) = ( 11, 13, 17 )
            //
            //    Rule 3:
            //      Order = 2
            //      W(1:2) = ( 19, 23 )
            //
            //    Product Rule:
            //      Order = 24
            //      W(1:24) =
            //        ( 2 * 11 * 19 )
            //        ( 3 * 11 * 19 )
            //        ( 4 * 11 * 19 )
            //        ( 7 * 11 * 19 )
            //        ( 2 * 13 * 19 )
            //        ( 3 * 13 * 19 )
            //        ( 5 * 13 * 19 )
            //        ( 7 * 13 * 19 )
            //        ( 2 * 17 * 19 )
            //        ( 3 * 17 * 19 )
            //        ( 5 * 17 * 19 )
            //        ( 7 * 17 * 19 )
            //        ( 2 * 11 * 23 )
            //        ( 3 * 11 * 23 )
            //        ( 5 * 11 * 23 )
            //        ( 7 * 11 * 23 )
            //        ( 2 * 13 * 23 )
            //        ( 3 * 13 * 23 )
            //        ( 5 * 13 * 23 )
            //        ( 7 * 13 * 23 )
            //        ( 2 * 17 * 23 )
            //        ( 3 * 17 * 23 )
            //        ( 5 * 17 * 23 )
            //        ( 7 * 17 * 23 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 April 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int FACTOR_INDEX, the index of the factor being processed.
            //    The first factor processed must be factor 0.
            //
            //    Input, int FACTOR_ORDER, the order of the factor.
            //
            //    Input, double FACTOR_VALUE[FACTOR_ORDER], the factor values for
            //    factor FACTOR_INDEX.
            //
            //    Input, int FACTOR_NUM, the number of factors.
            //
            //    Input, int POINT_NUM, the number of elements in the direct product.
            //
            //    Input/output, double W[POINT_NUM], the elements of the
            //    direct product, which are built up gradually.
            //
            //  Local Parameters:
            //
            //    Local, integer START, the first location of a block of values to set.
            //
            //    Local, integer CONTIG, the number of consecutive values to set.
            //
            //    Local, integer SKIP, the distance from the current value of START
            //    to the next location of a block of values to set.
            //
            //    Local, integer REP, the number of blocks of values to set.
            //
        {
            int i;
            int j;
            int k;
            int start;

            if (factor_index == 0)
            {
                data.contig = 1;
                data.skip = 1;
                data.rep = point_num;
                for (i = 0; i < point_num; i++)
                {
                    w[i] = 1.0;
                }
            }

            data.rep = data.rep / factor_order;
            data.skip = data.skip * factor_order;

            for (j = 0; j < factor_order; j++)
            {
                start = 0 + j * data.contig;

                for (k = 1; k <= data.rep; k++)
                {
                    for (i = start; i < start + data.contig; i++)
                    {
                        w[i] = w[i] * factor_value[j];
                    }

                    start = start + data.skip;
                }
            }

            data.contig = data.contig * factor_order;
        }

    }
}