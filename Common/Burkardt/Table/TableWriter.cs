using System;
using System.IO;

namespace Burkardt.Table
{
    public static partial class TableWriter
    {
        public static string file_name_ext_swap(string filename, string extension)
        {
            string[] tokens = filename.Split('.');
            if (tokens.Length == 1)
            {
                return String.Join('.', filename, extension);
            }

            string ret = "";
            for (int i = 0; i < tokens.Length - 1; i++)
            {
                ret += tokens[i] + ".";
            }

            ret += extension;

            return ret;
        }

    }
}