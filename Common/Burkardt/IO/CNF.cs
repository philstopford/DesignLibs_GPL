using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using Burkardt.Types;
using static System.String;

namespace Burkardt.IO;

public static class CNF
{
    public static bool cnf_data_read(string cnf_file_name, int v_num, int c_num,
            int l_num, ref int[] l_c_num, ref int[] l_val )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CNF_DATA_READ reads the data of a CNF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string CNF_FILE_NAME, the name of the CNF file.
        //
        //    Input, int V_NUM, the number of variables.
        //
        //    Input, int C_NUM, the number of clauses.
        //
        //    Input, int L_NUM, the number of signed literals.
        //
        //    Output, int L_C_NUM[C_NUM], the number of signed
        //    literals occuring in each clause.
        //
        //    Output, int L_VAL[L_NUM], a list of all the signed 
        //    literals in all the clauses, ordered by clause.
        //
        //    Output, bool CNF_DATA_READ, is TRUE if there was an error during 
        //    the read.
        //
    {
        string[] input;

        try
        {
            input = File.ReadAllLines(cnf_file_name);

        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("CNF_DATA_READ - Fatal error!");
            Console.WriteLine("  Could not open file.");
            return true;
        }
            
        //
        //  Read lines until you find one that is not blank and does not begin
        //  with a "c".  This should be the header line.
        //
        int index = 0;
        string line = "";
        foreach (string tmp in input)
        {
            index++;
            switch (tmp[0])
            {
                case 'c':
                case 'C':
                    continue;
            }

            if (0 >= typeMethods.s_len_trim(tmp))
            {
                continue;
            }

            line = tmp;
            break;
        }

        switch (line)
        {
            case "":
                Console.WriteLine("");
                Console.WriteLine("CNF_DATA_READ - Fatal error!");
                Console.WriteLine("  Error3 while reading the file.");
                return true;
        }

        switch (line.StartsWith("p cnf "))
        {
            //
            //  We expect to be reading the line "p cnf V_NUM C_NUM"
            //
            case false:
                Console.WriteLine("");
                Console.WriteLine("CNF_DATA_READ - Fatal error!");
                Console.WriteLine("  First non-comment non-blank line does not start");
                Console.WriteLine("  with 'p cnf' marker.");
                return true;
        }

        //
        //  Remove the first four characters and shift left.
        //
        line = line.Replace("p cnf ", "");

        //
        //  Extract the next word, which is the number of variables.
        //  You can compare this to V_NUM for an extra check.
        //
        try
        {
            typeMethods.s_to_i4vec(line, 2);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("CNF_DATA_READ - Fatal error!");
            Console.WriteLine("  Unexpected End of input.");
            return true;
        }

        //
        //  Read remaining lines, counting the literals while ignoring occurrences of '0'.
        //
        int l_num2 = 0;
        int c_num2 = 0;
        int l_c_num2 = 0;

        foreach (string tmp in input.Skip(index))
        {
            line = tmp;

            switch (line[0])
            {
                case 'c':
                case 'C':
                    continue;
            }

            if (typeMethods.s_len_trim(line) == 0)
            {
                continue;
            }

            while (true)
            {
                string tmp2 = tmp.Replace("       ", " ");
                string[] tokens = tmp2.Split(' ');
                string word = tokens[0];
                Join( " ", tokens.Skip(1) );

                if (typeMethods.s_len_trim(word) <= 0)
                {
                    break;
                }

                int l_val2 = typeMethods.s_to_i4(word).val;
                l_val = new int[tokens.Length];
                l_c_num = new int[tokens.Length];

                if (l_val2 != 0)
                {
                    l_val[l_num2] = l_val2;
                    l_num2 += 1;
                    l_c_num2 += 1;
                }
                else
                {
                    l_c_num[c_num2] = l_c_num2;
                    c_num2 += 1;
                    l_c_num2 = 0;
                }
            }
        }

        //
        //  At the end:
        //
        //    C_NUM2 should equal C_NUM.
        //    L_NUM2 should equal L_NUM.
        //
        //  Close file and return.
        //

        return false;
    }

    public static bool cnf_data_write(int c_num, int l_num, int[] l_c_num, int[] l_val,
            ref List<string> output_unit )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CNF_DATA_WRITE writes data to a CNF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int C_NUM, the number of clauses.
        //
        //    Input, int L_NUM, the total number of signed literals.
        //
        //    Input, int L_C_NUM[C_NUM], the number of signed
        //    literals occuring in each clause.
        //
        //    Input, int L_VAL[L_NUM], a list of all the signed 
        //    literals in all the clauses, ordered by clause.
        //
        //    Input, ofstream &OUTPUT_UNIT, the output unit.
        //
    {
        int c;
        
        int l = 0;

        string line = "";

        for (c = 0; c < c_num; c++)
        {
            int l_c;
            for (l_c = 0; l_c < l_c_num[c]; l_c++)
            {
                line += " " + l_val[l].ToString(CultureInfo.InvariantCulture).PadLeft(7);
                l += 1;

                switch ((l_c + 1) % 10)
                {
                    case 0:
                        output_unit.Add(line);
                        break;
                }
            }

            line += " " + 0.ToString(CultureInfo.InvariantCulture).PadLeft(7);
            output_unit.Add(line);
        }

        return false;
    }

    public static bool cnf_evaluate(int v_num, int c_num, int l_num, int[] l_c_num, int[] l_val,
            bool[] v_val )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CNF_EVALUATE evaluates a formula in CNF form.
        //
        //  Discussion:
        //
        //    The formula is in conjunctive normal form.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int V_NUM, the number of variables.
        //
        //    Input, int C_NUM, the number of clauses.
        //
        //    Input, int L_NUM, the total number of signed literals.
        //
        //    Input, int L_C_NUM[C_NUM], the number of signed
        //    literals occuring in each clause.
        //
        //    Input, int L_VAL[L_NUM], a list of all the signed 
        //    literals in all the clauses, ordered by clause.
        //
        //    Input, bool V_VAL[V_NUM], the values assigned to the variables.
        //
        //    Output, bool CNF_EVALUATE, the value of the CNF formula for the
        //    given variable values.
        //
    {
        int c;

        bool f_val = true;

        int l = 0;

        for (c = 0; c < c_num; c++)
        {
            //
            //  The clause is false unless some signed literal is true.
            //
            bool c_val = false;
            int l_c;
            for (l_c = 0; l_c < l_c_num[c]; l_c++)
            {
                bool s_val = 0 < l_val[l];
                int v_index = Math.Abs(l_val[l]);
                l += 1;
                //
                //  The signed literal is true if the sign "equals" the value.
                //  Note that we CAN'T exit the loop because we need to run out the 
                //  L index!
                //
                if (v_val[v_index - 1] == s_val)
                {
                    c_val = true;
                }
            }

            //
            //  The formula is false if any clause is false.
            //
            if (c_val)
            {
                continue;
            }

            f_val = false;
            break;
        }

        return f_val;
    }

    public static bool cnf_header_read(string cnf_file_name, ref int v_num, ref int c_num,
            ref int l_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CNF_HEADER_READ reads the header of a CNF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string CNF_FILE_NAME, the name of the CNF file.
        //
        //    Output, int *V_NUM, the number of variables.
        //
        //    Output, int *C_NUM, the number of clauses.
        //
        //    Output, int *L_NUM, the number of signed literals.
        //
        //    Output, bool CNF_HEADER_READ, is TRUE if there was an error during 
        //    the read.
        //
    {
        string[] input;
        string line = "";
        
        try
        {
            input = File.ReadAllLines(cnf_file_name);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("CNF_HEADER_READ - Fatal error!");
            Console.WriteLine("  Could not open file.");
            return true;
        }

        //
        //  Read lines until you find one that is not blank and does not begin
        //  with a "c".  This should be the header line.
        //
        int index = 0;
        foreach (string tmp in input)
        {
            index++;
                
            switch (tmp[0])
            {
                case 'c':
                case 'C':
                    continue;
            }

            if (0 >= typeMethods.s_len_trim(tmp))
            {
                continue;
            }

            line = tmp;
            break;
        }

        switch (line)
        {
            case "":
                Console.WriteLine("");
                Console.WriteLine("CNF_HEADER_READ - Fatal error!");
                Console.WriteLine("  Error3 while reading the file.");
                return true;
        }

        switch (line.StartsWith("p cnf "))
        {
            //
            //  We expect to be reading the line "p cnf V_NUM C_NUM"
            //
            case false:
                Console.WriteLine("");
                Console.WriteLine("CNF_DATA_READ - Fatal error!");
                Console.WriteLine("  First non-comment non-blank line does not start");
                Console.WriteLine("  with 'p cnf' marker.");
                return true;
        }

        //
        //  Remove the first four characters and shift left.
        //
        line = line.Replace("p cnf ", "");

        //
        //  Extract the next word, which is the number of variables.
        //
        try
        {
            i4vec ti = typeMethods.s_to_i4vec(line, 2);
            v_num = ti.ivec[0];
            c_num = ti.ivec[1];
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("CNF_HEADER_READ - Fatal error!");
            Console.WriteLine("  Unexpected End of input.");
            return true;
        }

        //
        //  Read remaining lines, counting the literals while ignoring occurrences of '0'.
        //
        l_num = 0;

        foreach (string tmp in input.Skip(index))
        {
            line = tmp;

            switch (line[0])
            {
                case 'c':
                case 'C':
                    continue;
            }

            if (typeMethods.s_len_trim(line) < 0)
            {
                break;
            }

            if (typeMethods.s_len_trim(line) == 0)
            {
                continue;
            }

            while (true)
            {
                string tmp2 = tmp.Replace("       ", " ");
                string[] tokens = tmp2.Split(' ');
                string word = tokens[0];
                Join(" ", tokens.Skip(1));

                if (typeMethods.s_len_trim(word) <= 0)
                {
                    break;
                }

                int l_val = typeMethods.s_to_i4(word).val;
                
                if (l_val != 0)
                {
                    l_num += 1;
                }
            }
        }

        return false;
    }

    public static bool cnf_header_write(int v_num, int c_num, string output_name,
            ref List<string> output_unit)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CNF_HEADER_WRITE writes the header for a CNF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int V_NUM, the number of variables.
        //
        //    Input, int C_NUM, the number of clauses.
        //
        //    Input, string OUTPUT_NAME, the name of the output file.
        //
        //    Input, ofstream &OUTPUT_UNIT, the output unit.
        //
    {
        const bool error = false;

        output_unit.Add("c " + output_name + "");
        output_unit.Add("c");
        output_unit.Add("p cnf " + v_num + " " + c_num + "");

        return error;
    }

    public static void cnf_print(int v_num, int c_num, int l_num, int[] l_c_num, int[] l_val )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CNF_PRINT prints CNF information.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int V_NUM, the number of variables.
        //
        //    Input, int C_NUM, the number of clauses.
        //
        //    Input, int L_NUM, the total number of signed literals.
        //
        //    Input, int L_C_NUM[C_NUM], the number of signed
        //    literals occuring in each clause.
        //
        //    Input, int L_VAL[L_NUM], a list of all the signed 
        //    literals in all the clauses, ordered by clause.
        //
    {
        int c;

        Console.WriteLine("");
        Console.WriteLine("CNF data printout:");
        Console.WriteLine("");
        Console.WriteLine("  The number of variables       V_NUM  = " + v_num + "");
        Console.WriteLine("  The number of clauses         C_NUM  = " + c_num + "");
        Console.WriteLine("  The number of signed literals L_NUM  = " + l_num + "");
        int l = 0;
        for (c = 0; c < c_num; c++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Clause " + c
                                          + " includes " + l_c_num[c] + " signed literals:");
            int l_c;
            for (l_c = 0; l_c < l_c_num[c]; l_c++)
            {
                Console.WriteLine(l_val[l].ToString(CultureInfo.InvariantCulture).PadLeft(4) + "");
                l += 1;
            }
        }
    }

    public static bool cnf_write(int v_num, int c_num, int l_num, int[] l_c_num, int[] l_val,
            string output_name )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CNF_WRITE writes the header and data of a CNF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int V_NUM, the number of variables.
        //
        //    Input, int C_NUM, the number of clauses.
        //
        //    Input, int L_NUM, the total number of signed literals.
        //
        //    Input, int L_C_NUM[C_NUM], the number of signed
        //    literals occuring in each clause.
        //
        //    Input, int L_VAL[L_NUM], a list of all the signed 
        //    literals in all the clauses, ordered by clause.
        //
        //    Input, string OUTPUT_NAME, the name of the output file.
        //
    {
        List<string> output_unit = new();

        bool error = false;
        //
        //  Write the header.
        //
        error = cnf_header_write(v_num, c_num, output_name, ref output_unit);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("CNF_WRITE - Fatal error!");
                Console.WriteLine("  Cannot write the header for the output file \"" + output_name + "\".");
                return true;
        }

        //
        //  Write the data.
        //
        error = cnf_data_write(c_num, l_num, l_c_num, l_val, ref output_unit);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("CNF_WRITE - Fatal error!");
                Console.WriteLine("  Cannot write the data for the output file \"" + output_name + "\".");
                return true;
            default:
                try
                {
                    File.WriteAllLines(output_name, output_unit);
                }
                catch (Exception)
                {
                    Console.WriteLine("");
                    Console.WriteLine("CNF_WRITE - Fatal error!");
                    Console.WriteLine("  Cannot open the output file \"" + output_name + "\".");
                    error = true;
                }
            
                return error;
        }
    }
}