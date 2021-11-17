using System;
using Burkardt.Types;

namespace Burkardt.IntegralNS;

public static class Lattice
{
    public static double lattice(int dim_num, int m, int[] z,
            Func<int, double[], double> f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LATTICE applies a lattice integration rule.
        //
        //  Discussion:
        //
        //    Because this is a standard lattice rule, it is really only suited
        //    for functions which are periodic, of period 1, in both X and Y.
        //
        //    For a suitable F, and a given value of M (the number of lattice points),
        //    the performance of the routine is affected by the choice of the
        //    generator vector Z.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 November 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Ian Sloan, Stephen Joe,
        //    Lattice Methods for Multiple Integration,
        //    Oxford, 1994,
        //    ISBN: 0198534728,
        //    LC: QA311.S56
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int M, the order (number of points) of the rule.
        //
        //    Input, int Z[DIM_NUM], the generator vector.  Typically, the elements
        //    of Z satisfy 1 <= Z[*] < M, and are relatively prime to M.
        //    This is easy to guarantee if M is itself a prime number.
        //
        //    Input, Func < int, double[], double > f, the name of the 
        //    user-supplied routine which evaluates the function.
        //
        //    Output, double LATTICE, the estimated integral.
        //
    {
        int i;
        int j;
        double quad;
        double[] x;

        x = new double[dim_num];

        quad = 0.0;

        for (j = 0; j <= m - 1; j++)
        {
            for (i = 0; i < dim_num; i++)
            {
                x[i] = j * z[i] / (double) m % 1.0;
            }

            quad += f(dim_num, x);
        }

        quad /= m;

        return quad;
    }

    public static double lattice_np0(int dim_num, int m, int[] z,
            Func<int, double[], double> f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LATTICE_NP0 applies a lattice integration rule to a nonperiodic function.
        //
        //  Discussion:
        //
        //    This routine attempts to modify a lattice rule, suitable for use
        //    with a periodic function, for use with a nonperiodic function F(X),
        //    essentially by applying the lattice rule to the function
        //
        //      G(X) = ( F(X) + F(1-X) ) / 2
        //
        //    This is the rule in 1 dimension.  In two dimensions, we have
        //
        //      G(X,Y) = ( F(X,Y) + F(X,1-Y) + F(1-X,Y) + F(1-X,1-Y) ) / 4
        //
        //    with the obvious generalizations to higher dimensions.  
        //
        //    Drawbacks of this approach include:
        //
        //    * in dimension DIM_NUM, we must evaluate the function F at 
        //      2**DIM_NUM points for every single evaluation of G;
        //
        //    * the function G, regarded as a periodic function, is continuous,
        //      but not generally differentiable, at 0 and 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 November 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Seymour Haber,
        //    Parameters for Integrating Periodic Functions of Several Variables,
        //    Mathematics of Computation,
        //    Volume 41, Number 163, July 1983, pages 115-129.
        //
        //    Ian Sloan, Stephen Joe,
        //    Lattice Methods for Multiple Integration,
        //    Oxford, 1994,
        //    ISBN: 0198534728,
        //    LC: QA311.S56
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int M, the order (number of points) of the rule.
        //
        //    Input, int Z[DIM_NUM], the generator vector.  Typically, the elements
        //    of Z satisfy 1 <= Z[*] < M, and are relatively prime to M.
        //    This is easy to guarantee if M is itself a prime number.
        //
        //    Input, Func < int, double[], double > f, the name of the 
        //    user-supplied routine which evaluates the function.
        //
        //    Output, double LATTICE_NP0, the estimated integral.
        //
    {
        int dim;
        int gray_bit;
        int j;
        double quad;
        double[] x;
        double[] y;
        typeMethods.GrayData data = new();

        x = new double[dim_num];
        y = new double[dim_num];

        quad = 0.0;

        for (j = 0; j <= m - 1; j++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                x[dim] = j * z[dim] / (double) m % 1.0;
            }

            //
            //  Generate all DIM_NUM-tuples for which the I-th element is X(I) or 1-X(I).
            //
            gray_bit = -dim_num;

            for (;;)
            {
                typeMethods.gray_next(ref data, dim_num, ref gray_bit);

                if (gray_bit == -dim_num)
                {
                    break;
                }

                switch (gray_bit)
                {
                    case 0:
                    {
                        for (dim = 0; dim < dim_num; dim++)
                        {
                            y[dim] = x[dim];
                        }

                        break;
                    }
                    default:
                        dim = Math.Abs(gray_bit) - 1;
                        y[dim] = 1.0 - y[dim];
                        break;
                }

                quad += f(dim_num, y);
            }
        }

        quad /= (int) Math.Pow(2, dim_num) * m;

        return quad;
    }

    public static double lattice_np1(int dim_num, int m, int[] z,
            Func<int, double[], double> f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LATTICE_NP1 applies a lattice integration rule to a nonperiodic function.
        //
        //  Discussion:
        //
        //    This routine applies the transformation function
        //
        //      PHI(T) = 3*T^2 - 2*T^3
        //
        //    to a nonperiodic integrand to make it suitable for treatment
        //    by a lattice rule.
        //
        //    For a suitable F, and a given value of M (the number of lattice points),
        //    the performance of the routine is affected by the choice of the
        //    generator vector Z.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 November 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Ian Sloan, Stephen Joe,
        //    Lattice Methods for Multiple Integration,
        //    Oxford, 1994,
        //    ISBN: 0198534728,
        //    LC: QA311.S56
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int M, the order (number of points) of the rule.
        //
        //    Input, int Z[DIM_NUM], the generator vector.  Typically, the elements
        //    of Z satisfy 1 <= Z(1:DIM_NUM) < M, and are relatively prime to M.
        //    This is easy to guarantee if M is itself a prime number.
        //
        //    Input, Func < int, double[], double > f, the name of the 
        //    user-supplied routine which evaluates the function.
        //
        //    Output, double LATTICE_NP1, the estimated integral.
        //
    {
        int dim;
        double dphi;
        int j;
        double quad;
        double[] x;
        double[] y;

        x = new double[dim_num];
        y = new double[dim_num];

        quad = 0.0;

        for (j = 0; j <= m - 1; j++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                x[dim] = j * z[dim] / (double) m % 1.0;
            }

            dphi = 1.0;
            for (dim = 0; dim < dim_num; dim++)
            {
                y[dim] = (3.0 - 2.0 * x[dim]) * x[dim] * x[dim];
                dphi = dphi * 6.0 * (1.0 - x[dim]) * x[dim];
            }

            quad += f(dim_num, y) * dphi;
        }

        quad /= m;

        return quad;
    }

    public static void lattice_print(int dim_num, int m, int[] z, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LATTICE_PRINT prints the points in a lattice rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 November 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Ian Sloan, Stephen Joe,
        //    Lattice Methods for Multiple Integration,
        //    Oxford, 1994,
        //    ISBN: 0198534728,
        //    LC: QA311.S56
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int M, the number of points to use.
        //
        //    Input, int Z[DIM_NUM], the generator vector.
        //
        //    Input, string TITLE, a title.
        //
    {
        int dim;
        int i;
        int[] y;

        y = new int[dim_num];

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");

        for (i = 0; i <= m - 1; i++)
        {
            string cout = (i + 1).ToString().PadLeft(4) + "    ";
            for (dim = 0; dim < dim_num; dim++)
            {
                y[dim] = i * z[dim] % m;
                cout += y[dim].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }
}