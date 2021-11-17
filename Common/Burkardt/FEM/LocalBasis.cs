namespace Burkardt.FEM;

public static class LocalBasis
{
    public static double[] local_basis_1d(int order, double[] node_x, double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOCAL_BASIS_1D evaluates the basis functions in an element.
        //
        //  Discussion:
        //
        //    PHI(I)(X) = product ( J ~= I ) ( X         - NODE_X(I) ) 
        //                                 / ( NODE_X(J) - NODE_X(I) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the element.
        //    0 <= ORDER.  ORDER = 1 means piecewise linear.
        //
        //    Input, double NODE_X[ORDER], the element nodes.  
        //    These must be distinct.  Basis function I is 1 when X = NODE_X(I) 
        //    and 0 when X is equal to any other node.
        //
        //    Input, double X, the point at which the basis functions are to 
        //    be evaluated.
        //
        //    Output, double LOCAL_BASIS_1D[ORDER], the basis functions.
        //
    {
        double[] phi = new double[order];

        for (int j = 0; j < order; j++)
        {
            phi[j] = 1.0;
            for (int i = 0; i < order; i++)
            {
                if (j != i)
                {
                    phi[j] = phi[j] * (x - node_x[i]) / (node_x[j] - node_x[i]);
                }
            }
        }

        return phi;
    }

    public static double[] local_basis_prime_1d(int order, double[] node_x, double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOCAL_BASIS_PRIME_1D evaluates the basis function derivatives in an element.
        //
        //  Discussion:
        //
        //    PHI(I)(X) = product ( J ~= I ) ( X - NODE_X(I) ) 
        //                                 / ( NODE_X(J) - NODE_X(I) )
        //
        //    dPHIdx(I)(X) = sum ( J ~= I ) ( 1 / ( NODE_X(J) - NODE_X(I) ) *
        //      product ( K ~= ( J, I ) ) ( X - NODE_X(I) ) / ( NODE_X(J) - NODE_X(I) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 June 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the element.
        //    0 <= ORDER.  ORDER = 1 means piecewise linear.
        //
        //    Input, double NODE_X[ORDER], the element nodes.  
        //    These must be distinct.  Basis function I is 1 when X = NODE_X(I) 
        //    and 0 when X is equal to any other node.
        //
        //    Input, double X, the point at which the basis functions are to 
        //    be evaluated.
        //
        //    Output, double LOCAL_BASIS_PRIME_1D[ORDER], the basis functions.
        //
    {
        double[] dphidx = new double[order];

        for (int i = 0; i < order; i++)
        {
            dphidx[i] = 0.0;
            for (int j = 0; j < order; j++)
            {
                if (j != i)
                {
                    double term = 1.0 / (node_x[j] - node_x[i]);
                    for (int k = 0; k < order; k++)
                    {
                        if (k != i && k != j)
                        {
                            term = term * (x - node_x[i]) / (node_x[k] - node_x[i]);
                        }
                    }

                    dphidx[i] += term;
                }
            }
        }

        return dphidx;
    }

}