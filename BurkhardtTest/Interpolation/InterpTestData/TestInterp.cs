﻿using System;

namespace InterpTest;

public static class TestInterp
{
    public static double[] data_copy_new(int m, int n, double[] a1)

//****************************************************************************80
//
//  Purpose:
//
//    DATA_COPY_NEW copies one R8MAT to a "new" R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's, which
//    may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A1[M*N], the matrix to be copied.
//
//    Output, double DATA_COPY_NEW[M*N], the copy of A1.
//
    {
        int j;

        double[] a2 = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a2[i + j * m] = a1[i + j * m];
            }
        }

        return a2;
    }

    public static double[] p00_data(int prob, int dim_num, int data_num)

//****************************************************************************80
//
//  Purpose:
//
//    P00_DATA returns the data for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int DIM_NUM, the spatial dimension of the dependent
//    variables.
//
//    Input, int DATA_NUM, the number of data points.
//
//    Output, double P00_DATA[DIM_NUM*DATA_NUM], the data.
//
    {
        double[] p_data;

        switch (prob)
        {
            case 1:
                p_data = p01_data(dim_num, data_num);
                break;
            case 2:
                p_data = p02_data(dim_num, data_num);
                break;
            case 3:
                p_data = p03_data(dim_num, data_num);
                break;
            case 4:
                p_data = p04_data(dim_num, data_num);
                break;
            case 5:
                p_data = p05_data(dim_num, data_num);
                break;
            case 6:
                p_data = p06_data(dim_num, data_num);
                break;
            case 7:
                p_data = p07_data(dim_num, data_num);
                break;
            case 8:
                p_data = p08_data(dim_num, data_num);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_DATA - Fatal error!");
                Console.WriteLine("  Unexpected input value of PROB.");
                return new double[1];
        }

        return p_data;
    }

    public static int p00_data_num(int prob)

//****************************************************************************80
//
//  Purpose:
//
//    P00_DATA_NUM returns the number of data points for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Output, int P00_DATA_NUM, the number of data points.
//
    {
        int data_num;

        switch (prob)
        {
            case 1:
                data_num = p01_data_num();
                break;
            case 2:
                data_num = p02_data_num();
                break;
            case 3:
                data_num = p03_data_num();
                break;
            case 4:
                data_num = p04_data_num();
                break;
            case 5:
                data_num = p05_data_num();
                break;
            case 6:
                data_num = p06_data_num();
                break;
            case 7:
                data_num = p07_data_num();
                break;
            case 8:
                data_num = p08_data_num();
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_DATA_NUM - Fatal error!");
                Console.WriteLine("  Unexpected input value of PROB.");
                return 1;
        }

        return data_num;
    }

    public static int p00_dim_num(int prob)

//****************************************************************************80
//
//  Purpose:
//
//    P00_DIM_NUM returns the spatial dimension for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Output, int P00_DIM_NUM, the spatial dimension of the
//    dependent variables.
//
    {
        int dim_num;

        switch (prob)
        {
            case 1:
                dim_num = p01_dim_num();
                break;
            case 2:
                dim_num = p02_dim_num();
                break;
            case 3:
                dim_num = p03_dim_num();
                break;
            case 4:
                dim_num = p04_dim_num();
                break;
            case 5:
                dim_num = p05_dim_num();
                break;
            case 6:
                dim_num = p06_dim_num();
                break;
            case 7:
                dim_num = p07_dim_num();
                break;
            case 8:
                dim_num = p08_dim_num();
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_DIM_NUM - Fatal error!");
                Console.WriteLine("  Unexpected input value of PROB.");
                return 1;
        }

        return dim_num;
    }

    public static int p00_prob_num()

//****************************************************************************80
//
//  Purpose:
//
//    P00_PROB_NUM returns the number of test problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P00_PROB_NUM, the number of test problems.
//
    {
        const int prob_num = 8;

        return prob_num;
    }
//****************************************************************************80

    public static void p00_story(int prob)

//****************************************************************************80
//
//  Purpose:
//
//    P00_STORY prints the "story" for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
    {
        switch (prob)
        {
            case 1:
                p01_story();
                break;
            case 2:
                p02_story();
                break;
            case 3:
                p03_story();
                break;
            case 4:
                p04_story();
                break;
            case 5:
                p05_story();
                break;
            case 6:
                p06_story();
                break;
            case 7:
                p07_story();
                break;
            case 8:
                p08_story();
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("P00_STORY - Fatal error!");
                Console.WriteLine("  Unexpected input value of PROB.");
                break;
        }
    }
//****************************************************************************80

    public static double[] p01_data(int dim_num, int data_num)

//****************************************************************************80
//
//  Purpose:
//
//    P01_DATA returns the data for problem p01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension of the dependent
//    variables.
//
//    Input, int DATA_NUM, the number of data points.
//
//    Output, double P01_DATA[DIM_NUM*DATA_NUM], the data.
//
    {
        double[] p_data_save = {
                0.0, 4.0,
                1.0, 5.0,
                2.0, 6.0,
                4.0, 6.0,
                5.0, 5.0,
                6.0, 3.0,
                7.0, 1.0,
                8.0, 1.0,
                9.0, 1.0,
                10.0, 3.0,
                11.0, 4.0,
                12.0, 4.0,
                13.0, 3.0,
                14.0, 3.0,
                15.0, 4.0,
                16.0, 4.0,
                17.0, 3.0,
                18.0, 0.0
            }
            ;

        double[] p_data = data_copy_new(dim_num, data_num, p_data_save);

        return p_data;
    }

    public static int p01_data_num()

//****************************************************************************80
//
//  Purpose:
//
//    P01_DATA_NUM returns the number of data points for problem p01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P01_DATA_NUM, the number of data points.
//
    {
        const int data_num = 18;

        return data_num;
    }

    public static int p01_dim_num()

//****************************************************************************80
//
//  Purpose:
//
//    P01_DIM_NUM returns the spatial dimension for problem p01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P01_DIM_NUM, the spatial dimension of the 
//    dependent variables.
//
    {
        const int dim_num = 2;

        return dim_num;
    }

    public static void p01_story()

//****************************************************************************80
//
//  Purpose:
//
//    P01_STORY prints the "story" for problem p01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Hans-Joerg Wenz,
//    Interpolation of Curve Data by Blended Generalized Circles,
//    Computer Aided Geometric Design,
//    Volume 13, Number 8, November 1996, pages 673-680.
//
//  Parameters:
//
//    None
//
    {
        Console.WriteLine("");
        Console.WriteLine("  This example is due to Hans-Joerg Wenz.");
        Console.WriteLine("  It is an example of good data, which is dense enough in areas");
        Console.WriteLine("  where the expected curvature of the interpolant is large.");
        Console.WriteLine("  Good results can be expected with almost any reasonable");
        Console.WriteLine("  interpolation method.");

    }

    public static double[] p02_data(int dim_num, int data_num)

//****************************************************************************80
//
//  Purpose:
//
//    P02_DATA returns the data for problem p02.
//
//  Discussion:
//
//    Two pairs of identical X values have now been slightly separated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension of the dependent
//    variables.
//
//    Input, int DATA_NUM, the number of data points.
//
//    Output, double P02_DATA[DIM_NUM*DATA_NUM], the data.
//
    {
        double[] p_data_save = {
                0.00, 0.00,
                1.34, 5.00,
                5.00, 8.66,
                10.00, 10.00,
                10.60, 10.40,
                10.70, 12.00,
                10.705, 28.60,
                10.80, 30.20,
                11.40, 30.60,
                19.60, 30.60,
                20.20, 30.20,
                20.295, 28.60,
                20.30, 12.00,
                20.40, 10.40,
                21.00, 10.00,
                26.00, 8.66,
                29.66, 5.00,
                31.00, 0.00
            }
            ;

        double[] p_data = data_copy_new(dim_num, data_num, p_data_save);

        return p_data;
    }

    public static int p02_data_num()

//****************************************************************************80
//
//  Purpose:
//
//    P02_DATA_NUM returns the number of data points for problem p02.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P02_DATA_NUM, the number of data points.
//
    {
        const int data_num = 18;

        return data_num;
    }

    public static int p02_dim_num()

//****************************************************************************80
//
//  Purpose:
//
//    P02_DIM_NUM returns the spatial dimension for problem p02.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P02_DIM_NUM, the spatial dimension of the 
//    dependent variables.
//
    {
        const int dim_num = 2;

        return dim_num;
    }

    public static void p02_story()

//****************************************************************************80
//
//  Purpose:
//
//    P02_STORY prints the "story" for problem p02.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ETY Lee,
//    Choosing Nodes in Parametric Curve Interpolation,
//    Computer-Aided Design,
//    Volume 21, Number 6, July/August 1989, pages 363-370.
//
//  Parameters:
//
//    None
//
    {
        Console.WriteLine("");
        Console.WriteLine("  This example is due to ETY Lee of Boeing.");
        Console.WriteLine("  Data near the corners is more dense than in regions of small curvature.");
        Console.WriteLine("  A local interpolation method will produce a more plausible");
        Console.WriteLine("  interpolant than a nonlocal interpolation method, such as");
        Console.WriteLine("  cubic splines.");

    }

    public static double[] p03_data(int dim_num, int data_num)

//****************************************************************************80
//
//  Purpose:
//
//    P03_DATA returns the data for problem p03.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension of the dependent
//    variables.
//
//    Input, int DATA_NUM, the number of data points.
//
//    Output, double P03_DATA[DIM_NUM*DATA_NUM], the data.
//
    {
        double[] p_data_save = {
                0.0, 0.0,
                2.0, 10.0,
                3.0, 10.0,
                5.0, 10.0,
                6.0, 10.0,
                8.0, 10.0,
                9.0, 10.5,
                11.0, 15.0,
                12.0, 50.0,
                14.0, 60.0,
                15.0, 85.0
            }
            ;

        double[] p_data = data_copy_new(dim_num, data_num, p_data_save);

        return p_data;
    }

    public static int p03_data_num()

//****************************************************************************80
//
//  Purpose:
//
//    P03_DATA_NUM returns the number of data points for problem p03.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P03_DATA_NUM, the number of data points.
//
    {
        const int data_num = 11;

        return data_num;
    }

    public static int p03_dim_num()

//****************************************************************************80
//
//  Purpose:
//
//    P03_DIM_NUM returns the spatial dimension for problem p03.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P03_DIM_NUM, the spatial dimension of the 
//    dependent variables.
//
    {
        const int dim_num = 2;

        return dim_num;
    }

    public static void p03_story()

//****************************************************************************80
//
//  Purpose:
//
//    P03_STORY prints the "story" for problem p03.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    None
//
    {
        Console.WriteLine("");
        Console.WriteLine("  This example is due to Fred Fritsch and Ralph Carlson.");
        Console.WriteLine("  This data can cause problems for interpolation methods.");
        Console.WriteLine("  There are sudden changes in direction, and at the same time,");
        Console.WriteLine("  sparsely-placed data.  This can cause an interpolant to overshoot");
        Console.WriteLine("  the data in a way that seems implausible.");
    }

    public static double[] p04_data(int dim_num, int data_num)

//****************************************************************************80
//
//  Purpose:
//
//    P04_DATA returns the data for problem p04.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension of the dependent
//    variables.
//
//    Input, int DATA_NUM, the number of data points.
//
//    Output, double P04_DATA[DIM_NUM*DATA_NUM], the data.
//
    {
        double[] p_data_save = {
                0.00, 0.00,
                0.05, 0.70,
                0.10, 1.00,
                0.20, 1.00,
                0.80, 0.30,
                0.85, 0.05,
                0.90, 0.10,
                1.00, 1.00
            }
            ;

        double[] p_data = data_copy_new(dim_num, data_num, p_data_save);

        return p_data;
    }

    public static int p04_data_num()

//****************************************************************************80
//
//  Purpose:
//
//    P04_DATA_NUM returns the number of data points for problem p04.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P04_DATA_NUM, the number of data points.
//
    {
        const int data_num = 8;

        return data_num;
    }

    public static int p04_dim_num()

//****************************************************************************80
//
//  Purpose:
//
//    P04_DIM_NUM returns the spatial dimension for problem p04.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P04_DIM_NUM, the spatial dimension of the 
//    dependent variables.
//
    {
        const int dim_num = 2;

        return dim_num;
    }

    public static void p04_story()

//****************************************************************************80
//
//  Purpose:
//
//    P04_STORY prints the "story" for problem p04.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Larry Irvine, Samuel Marin, Philip Smith,
//    Constrained Interpolation and Smoothing,
//    Constructive Approximation,
//    Volume 2, Number 1, December 1986, pages 129-151.
//
//  Parameters:
//
//    None
//
    {
        Console.WriteLine("");
        Console.WriteLine("  This example is due to Larry Irvine, Samuel Marin and Philip Smith.");
        Console.WriteLine("  This data can cause problems for interpolation methods.");
        Console.WriteLine("  There are sudden changes in direction, and at the same time,");
        Console.WriteLine("  sparsely-placed data.  This can cause an interpolant to overshoot");
        Console.WriteLine("  the data in a way that seems implausible.");
    }

    public static double[] p05_data(int dim_num, int data_num)

//****************************************************************************80
//
//  Purpose:
//
//    P05_DATA returns the data for problem p05.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension of the dependent
//    variables.
//
//    Input, int DATA_NUM, the number of data points.
//
//    Output, double P05_DATA[DIM_NUM*DATA_NUM], the data.
//
    {
        double[] p_data_save = {
                0.00, 0.00,
                0.10, 0.90,
                0.20, 0.95,
                0.30, 0.90,
                0.40, 0.10,
                0.50, 0.05,
                0.60, 0.05,
                0.80, 0.20,
                1.00, 1.00
            }
            ;

        double[] p_data = data_copy_new(dim_num, data_num, p_data_save);

        return p_data;
    }

    public static int p05_data_num()

//****************************************************************************80
//
//  Purpose:
//
//    P05_DATA_NUM returns the number of data points for problem p05.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P05_DATA_NUM, the number of data points.
//
    {
        const int data_num = 9;

        return data_num;
    }

    public static int p05_dim_num()

//****************************************************************************80
//
//  Purpose:
//
//    P05_DIM_NUM returns the spatial dimension for problem p05.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P05_DIM_NUM, the spatial dimension of the 
//    dependent variables.
//
    {
        const int dim_num = 2;

        return dim_num;
    }

    public static void p05_story()

//****************************************************************************80
//
//  Purpose:
//
//    P05_STORY prints the "story" for problem p05.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Larry Irvine, Samuel Marin, Philip Smith,
//    Constrained Interpolation and Smoothing,
//    Constructive Approximation,
//    Volume 2, Number 1, December 1986, pages 129-151.
//
//  Parameters:
//
//    None
//
    {
        Console.WriteLine("");
        Console.WriteLine("  This example is due to Larry Irvine, Samuel Marin and Philip Smith.");
        Console.WriteLine("  This data can cause problems for interpolation methods.");
        Console.WriteLine("  There are sudden changes in direction, and at the same time,");
        Console.WriteLine("  sparsely-placed data.  This can cause an interpolant to overshoot");
        Console.WriteLine("  the data in a way that seems implausible.");

    }

    public static double[] p06_data(int dim_num, int data_num)

//****************************************************************************80
//
//  Purpose:
//
//    P06_DATA returns the data for problem p06.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension of the dependent
//    variables.
//
//    Input, int DATA_NUM, the number of data points.
//
//    Output, double P06_DATA[DIM_NUM*DATA_NUM], the data.
//
    {
        double[] p_data_save = {
                595.0, 0.644,
                605.0, 0.622,
                615.0, 0.638,
                625.0, 0.649,
                635.0, 0.652,
                645.0, 0.639,
                655.0, 0.646,
                665.0, 0.657,
                675.0, 0.652,
                685.0, 0.655,
                695.0, 0.644,
                705.0, 0.663,
                715.0, 0.663,
                725.0, 0.668,
                735.0, 0.676,
                745.0, 0.676,
                755.0, 0.686,
                765.0, 0.679,
                775.0, 0.678,
                785.0, 0.683,
                795.0, 0.694,
                805.0, 0.699,
                815.0, 0.710,
                825.0, 0.730,
                835.0, 0.763,
                845.0, 0.812,
                855.0, 0.907,
                865.0, 1.044,
                875.0, 1.336,
                885.0, 1.881,
                895.0, 2.169,
                905.0, 2.075,
                915.0, 1.598,
                925.0, 1.211,
                935.0, 0.916,
                945.0, 0.746,
                955.0, 0.672,
                965.0, 0.627,
                975.0, 0.615,
                985.0, 0.607,
                995.0, 0.606,
                1005.0, 0.609,
                1015.0, 0.603,
                1025.0, 0.601,
                1035.0, 0.603,
                1045.0, 0.601,
                1055.0, 0.611,
                1065.0, 0.601,
                1075.0, 0.608
            }
            ;

        double[] p_data = data_copy_new(dim_num, data_num, p_data_save);

        return p_data;
    }

    public static int p06_data_num()

//****************************************************************************80
//
//  Purpose:
//
//    P06_DATA_NUM returns the number of data points for problem p06.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P06_DATA_NUM, the number of data points.
//
    {
        const int data_num = 49;

        return data_num;
    }

    public static int p06_dim_num()

//****************************************************************************80
//
//  Purpose:
//
//    P06_DIM_NUM returns the spatial dimension for problem p06.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P06_DIM_NUM, the spatial dimension of the 
//    dependent variables.
//
    {
        const int dim_num = 2;

        return dim_num;
    }

    public static void p06_story()

//****************************************************************************80
//
//  Purpose:
//
//    P06_STORY prints the "story" for problem p06.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl DeBoor, John Rice,
//    Least-squares cubic spline approximation II - variable knots. 
//    Technical Report CSD TR 21, 
//    Purdue University, Lafayette, Indiana, 1968.
//
//    Carl DeBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    None
//
    {
        Console.WriteLine("");
        Console.WriteLine("  The data is due to deBoor and Rice.");
        Console.WriteLine("  The data represents a temperature dependent property of titanium.");
        Console.WriteLine("  The data has been used extensively as an example in spline");
        Console.WriteLine("  approximation with variably-spaced knots.");
        Console.WriteLine("  DeBoor considers two sets of knots:");
        Console.WriteLine("  (595,675,755,835,915,995,1075)");
        Console.WriteLine("  and");
        Console.WriteLine("  (595,725,850,910,975,1040,1075).");
    }

    public static double[] p07_data(int dim_num, int data_num)

//****************************************************************************80
//
//  Purpose:
//
//    P07_DATA returns the data for problem p07.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension of the dependent
//    variables.
//
//    Input, int DATA_NUM, the number of data points.
//
//    Output, double P07_DATA[DIM_NUM*DATA_NUM], the data.
//
    {
        double[] p_data_save = {
                0.0, 1.0,
                1.0, 2.0,
                4.0, 2.0,
                5.0, 1.0
            }
            ;

        double[] p_data = data_copy_new(dim_num, data_num, p_data_save);

        return p_data;
    }

    public static int p07_data_num()

//****************************************************************************80
//
//  Purpose:
//
//    P07_DATA_NUM returns the number of data points for problem p07.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P07_DATA_NUM, the number of data points.
//
    {
        const int data_num = 4;

        return data_num;
    }

    public static int p07_dim_num()

//****************************************************************************80
//
//  Purpose:
//
//    P07_DIM_NUM returns the spatial dimension for problem p07.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P07_DIM_NUM, the spatial dimension of the 
//    dependent variables.
//
    {
        const int dim_num = 2;

        return dim_num;
    }

    public static void p07_story()

//****************************************************************************80
//
//  Purpose:
//
//    P07_STORY prints the "story" for problem p07.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
    {
        Console.WriteLine("");
        Console.WriteLine("  This data is a simple symmetric set of 4 points,");
        Console.WriteLine("  for which it is interesting to develop the Shepard");
        Console.WriteLine("  interpolants for varying values of the exponent p.");

    }

    public static double[] p08_data(int dim_num, int data_num)

//****************************************************************************80
//
//  Purpose:
//
//    P08_DATA returns the data for problem p08.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension of the dependent
//    variables.
//
//    Input, int DATA_NUM, the number of data points.
//
//    Output, double P08_DATA[DIM_NUM*DATA_NUM], the data.
//
    {
        double[] p_data_save = {
                -1.0, 1.00,
                -0.8, 0.64,
                -0.6, 0.36,
                -0.4, 0.16,
                -0.2, 0.04,
                0.0, 0.00,
                0.2, 0.04,
                0.20001, 0.05,
                0.4, 0.16,
                0.6, 0.36,
                0.8, 0.64,
                1.0, 1.00
            }
            ;

        double[] p_data = data_copy_new(dim_num, data_num, p_data_save);

        return p_data;
    }

    public static int p08_data_num()

//****************************************************************************80
//
//  Purpose:
//
//    P08_DATA_NUM returns the number of data points for problem p08.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P08_DATA_NUM, the number of data points.
//
    {
        const int data_num = 12;

        return data_num;
    }

    public static int p08_dim_num()

//****************************************************************************80
//
//  Purpose:
//
//    P08_DIM_NUM returns the spatial dimension for problem p08.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P08_DIM_NUM, the spatial dimension of the 
//    dependent variables.
//
    {
        const int dim_num = 2;

        return dim_num;
    }

    public static void p08_story()

//****************************************************************************80
//
//  Purpose:
//
//    P08_STORY prints the "story" for problem p08.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
    {
        Console.WriteLine("");
        Console.WriteLine("  This is equally spaced data for y = x^2,");
        Console.WriteLine("  except for one extra point whose x value is");
        Console.WriteLine("  close to another, but whose y value is not so close.");
        Console.WriteLine("  A small disagreement in nearby data can be disaster.");
    }
}