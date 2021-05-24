namespace Burkardt.Table
{
    public static partial class TableReader
    {
        static TableHeader readHeader(string input_filename)
        {
            TableHeader ret = new TableHeader {m = TableMisc.file_column_count(input_filename), n = TableMisc.file_row_count ( input_filename )};

            return ret;
        }

    }
}