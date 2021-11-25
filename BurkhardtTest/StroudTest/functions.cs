using System;

namespace StroudTest;

public static class functions
{
    public static double function_1d(int function_1d_index, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FUNCTION_1D evaluates the current 1D function.
        //
        //  Discussion:
        //
        //    This routine assumes that the global variable FUNCTION_1D_INDEX has been
        //    set, and is accessible.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the value of the variable.
        //
        //    Output, double FUNCTION_1D, the value of the function.
        //
    {
        double value = 0;

        switch (function_1d_index)
        {
            case 1:
                value = 1.0;
                break;
            case 2:
                value = x;
                break;
            case 3:
                value = x * x;
                break;
            case 4:
            case 5:
                value = x * x * x * x;
                break;
            case 6:
                value = x * x * x * x * x;
                break;
            case 7:
                value = x * x * x * x * x * x;
                break;
            case 8:
                value = Math.Abs(x);
                break;
            case 9:
                value = Math.Sin(x);
                break;
            case 10:
                value = Math.Exp(x);
                break;
            case 11:
                value = 1.0 / (1.0 + Math.Abs(x));
                break;
            case 12:
                value = Math.Sqrt(Math.Abs(x));
                break;
            default:
                value = 0.0;
                break;
        }

        return value;
    }

    public static void function_1d_name(int function_1d_index, ref string name)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FUNCTION_1D_NAME returns the name of the current 1D function.
        //
        //  Discussion:
        //
        //    This routine assumes that the global variable FUNCTION_1D_INDEX has been
        //    set, and is accessible.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, char NAME[8], the name of the current 1D function.
        //
    {
        name = function_1d_index switch
        {
            1 => "      1",
            2 => "      X",
            3 => "    X^2",
            4 => "    X^3",
            5 => "    X^4",
            6 => "    X^5",
            7 => "    X^6",
            8 => "      R",
            9 => " SIN(X)",
            10 => " EXP(X)",
            11 => "1/(1+R)",
            12 => "SQRT(R)",
            _ => "???????"
        };
    }

    public static int function_1d_num()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FUNCTION_1D_NUM returns the number of 1D functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int FUNCTION_1D_NUM, the number of 1D functions.
        //
    {
        const int value = 12;

        return value;
    }

    public static double function_2d(int function_2d_index, double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FUNCTION_2D evaluates the current 2D function.
        //
        //  Discussion:
        //
        //    This routine assumes that the global variable FUNCTION_2D_INDEX has been
        //    set, and is accessible.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, double Y, the value of the variables.
        //
        //    Output, double FUNCTION_2D, the value of the function.
        //
    {
        double value = function_2d_index switch
        {
            1 => 1.0,
            2 => x,
            3 => x * x,
            4 => x * x * x,
            5 => x * x * x * x,
            6 => x * x * x * x * x,
            7 => x * x * x * x * x * x,
            8 => Math.Sqrt(x * x + y * y),
            9 => Math.Sin(x),
            10 => Math.Exp(x),
            11 => 1.0 / (1.0 + Math.Sqrt(x * x + y * y)),
            12 => Math.Sqrt(Math.Sqrt(x * x + y * y)),
            _ => 0.0
        };

        return value;
    }

    public static void function_2d_name(int function_2d_index, ref string name)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FUNCTION_2D_NAME returns the name of the current 2D function.
        //
        //  Discussion:
        //
        //    This routine assumes that the global variable FUNCTION_2D_INDEX has been
        //    set, and is accessible.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, char NAME[8], the name of the current 2D function.
        //
    {
        name = function_2d_index switch
        {
            1 => "      1",
            2 => "      X",
            3 => "    X^2",
            4 => "    X^3",
            5 => "    X^4",
            6 => "    X^5",
            7 => "    X^6",
            8 => "      R",
            9 => " SIN(X)",
            10 => " EXP(X)",
            11 => "1/(1+R)",
            12 => "SQRT(R)",
            _ => "???????"
        };
    }

    public static int function_2d_num()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FUNCTION_2D_NUM returns the number of 2D functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int FUNCTION_2D_NUM, the number of 2D functions.
        //
    {
        const int value = 12;

        return value;
    }

    public static double function_3d(int function_3d_index, double x, double y, double z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FUNCTION_3D evaluates a function F(X,Y,Z) of 3 variables.
        //
        //  Discussion:
        //
        //    This routine assumes that the global variable FUNCTION_3D_INDEX has been
        //    set, and is accessible.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, Z, the value of the variables.
        //
        //    Output, double FUNC_3D, the value of the function.
        //
    {
        double value = function_3d_index switch
        {
            1 => 1.0,
            2 => x,
            3 => y,
            4 => z,
            5 => x * x,
            6 => x * y,
            7 => x * z,
            8 => y * y,
            9 => y * z,
            10 => z * z,
            11 => x * x * x,
            12 => x * y * z,
            13 => z * z * z,
            14 => x * x * x * x,
            15 => x * x * z * z,
            16 => z * z * z * z,
            17 => x * x * x * x * x,
            18 => Math.Pow(x, 6),
            19 => Math.Sqrt(x * x + y * y + z * z),
            20 => Math.Sin(x),
            21 => Math.Exp(x),
            22 => 1.0 / Math.Sqrt(1.0 + x * x + y * y + z * z),
            23 => Math.Sqrt(Math.Sqrt(x * x + y * y + z * z)),
            _ => 0.0
        };

        return value;
    }

    public static void function_3d_name(int function_3d_index, ref string name)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FUNCTION_3D_NAME returns the name of the current 3D function.
        //
        //  Discussion:
        //
        //    This routine assumes that the global variable FUNCTION_3D_INDEX has been
        //    set, and is accessible.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, char NAME[8], the name of the current 3D function.
        //
    {
        name = function_3d_index switch
        {
            1 => "      1",
            2 => "      X",
            3 => "      Y",
            4 => "      Z",
            5 => "    X*X",
            6 => "    X*Y",
            7 => "    X*Z",
            8 => "    Y*Y",
            9 => "    Y*Z",
            10 => "    Z*Z",
            11 => "    X^3",
            12 => "  X*Y*Z",
            13 => "  Z*Z*Z",
            14 => "    X^4",
            15 => "X^2 Z^2",
            16 => "    Z^4",
            17 => "    X^5",
            18 => "    X^6",
            19 => "      R",
            20 => " SIN(X)",
            21 => " EXP(X)",
            22 => "1/(1+R)",
            23 => "SQRT(R)",
            _ => "???????"
        };
    }

    public static int function_3d_num()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FUNCTION_3D_NUM returns the number of 3D functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int FUNCTION_3D_NUM, the number of 3D functions.
        //
    {
        const int value = 23;

        return value;
    }

    public static double function_nd(int function_nd_index, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FUNCTION_ND evaluates the current ND function.
        //
        //  Discussion:
        //
        //    This routine assumes that the global variable FUNCTION_ND_INDEX has been
        //    set, and is accessible.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of variables.
        //
        //    Input, double X[N], the value of the variables.
        //
        //    Output, double FUNCTION_ND, the value of the function.
        //
    {
        int i;
        double temp;
        double value;

        switch (function_nd_index)
        {
            case 1:
                value = 1.0;
                break;
            case 2:
                value = x[0];
                break;
            case 3:
                value = Math.Pow(x[0], 2);
                break;
            case 4:
                value = Math.Pow(x[0], 3);
                break;
            case 5:
                value = Math.Pow(x[0], 4);
                break;
            case 6:
                value = Math.Pow(x[0], 5);
                break;
            case 7:
                value = Math.Pow(x[0], 6);
                break;
            case 8:
            {
                temp = 0.0;
                for (i = 0; i < n; i++)
                {
                    temp += x[i] * x[i];
                }

                value = Math.Sqrt(temp);
                break;
            }
            case 9:
                value = Math.Sin(x[0]);
                break;
            case 10:
                value = Math.Exp(x[0]);
                break;
            case 11:
            {
                temp = 0.0;
                for (i = 0; i < n; i++)
                {
                    temp += x[i] * x[i];
                }

                value = 1.0 / (1.0 + Math.Sqrt(temp));
                break;
            }
            case 12:
            {
                temp = 0.0;
                for (i = 0; i < n; i++)
                {
                    temp += x[i] * x[i];
                }

                value = Math.Sqrt(Math.Sqrt(temp));
                break;
            }
            default:
                value = 0.0;
                break;
        }

        return value;
    }

    public static void function_nd_name(int function_nd_index, ref string name)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FUNCTION_ND_NAME returns the name of the current ND function.
        //
        //  Discussion:
        //
        //    This routine assumes that the global variable FUNCTION_ND_INDEX has been
        //    set, and is accessible.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, char NAME[8], the name of the current ND function.
        //
    {
        name = function_nd_index switch
        {
            1 => "      1",
            2 => "      X",
            3 => "    X^2",
            4 => "    X^3",
            5 => "    X^4",
            6 => "    X^5",
            7 => "    X^6",
            8 => "      R",
            9 => " SIN(X)",
            10 => " EXP(X)",
            11 => "1/(1+R)",
            12 => "SQRT(R)",
            _ => "???????"
        };
    }

    public static int function_nd_num()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FUNCTION_ND_NUM returns the number of ND functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int FUNCTION_ND_NUM, the number of ND functions.
        //
    {
        const int value = 12;

        return value;
    }

    public static double f_1_2d(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_1_2D evaluates the function 1 in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the arguments.
        //
        //    Output, double F_1_2D, the value of the function.
        //
    {
        const double value = 1.0;

        return value;
    }

    public static double f_x_2d(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_X_2D evaluates the function X in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the arguments.
        //
        //    Output, double F_X_2D, the value of the function.
        //
    {
        return x;
    }

    public static double f_r_2d(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_R_2D evaluates the function Math.Sqrt ( X**2 + Y**2 ) in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the arguments.
        //
        //    Output, double F_R_2D, the value of the function.
        //
    {
        double value = Math.Sqrt(x * x + y * y);

        return value;
    }


}