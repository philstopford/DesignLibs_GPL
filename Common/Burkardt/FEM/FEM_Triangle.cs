using System.Collections.Generic;
using System.IO;

namespace Burkardt.FEM;

public static class Triangle
{
    public static void triangle_element_write(string element_file, int element_num, int element_order,
            int element_att_num, int[] element_node, double[] element_att )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_ELEMENT_WRITE writes a TRIANGLE ".ele" file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 December 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string ELEMENT_FILE, the name of the file to be written.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_ORDER, the order of the elements.
        //
        //    Input, int ELEMENT_ATT_NUM, the number of element attributes.
        //
        //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the indices of the
        //    nodes that make up each element.
        //
        //    Input, double ELEMENT_ATT[ELEMENT_ATT_NUM*ELEMENT_NUM], the attributes
        //    of each element.
        //
    {
        List<string> lines = new()
        {
            "  " + element_num
                 + "  " + element_order
                 + "  " + element_att_num
        };

        for (int element = 0; element < element_num; element++)
        {
            string line =  "  " + (element + 1);
            for (int order = 0; order < element_order; order++)
            {
                line += "  " + element_node[order + element * element_order];
            }

            for (int att = 0; att < element_att_num; att++)
            {
                line += element_att[att + element * element_att_num];
            }

            lines.Add(line);
        }

        File.WriteAllLines(element_file, lines);
    }

    public static void triangle_node_write(string node_file, int node_num, int node_dim,
            int node_att_num, int node_marker_num, double[] node_coord,
            double[] node_att, int[] node_marker )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NODE_WRITE writes a TRIANGLE ".node" file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 December 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string NODE_FILE, the name of the node file to be written.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int NODE_DIM, the spatial dimension.
        //
        //    Input, int NODE_ATT_NUM, number of node attributes listed on each 
        //    node record.
        //
        //    Input, int NODE_MARKER_NUM, 1 if every node record includes a final
        //    boundary marker value.
        //
        //    Output, double NODE_COORD[NODE_DIM*NODE_NUM], the nodal coordinates.
        //
        //    Output, double NODE_ATT[NODE_ATT_NUM*NODE_NUM], the nodal attributes.
        //
        //    Output, int NODE_MARKER[NODE_MARKER_NUM*NODE_NUM], the node markers.
        //
    {
        List<string> lines = new()
        {
            "  " + node_num
                 + "  " + node_dim
                 + "  " + node_att_num
                 + "  " + node_marker_num
        };

        for (int node = 0; node < node_num; node++)
        {
            string line = "  " + (node + 1);
            for (int dim = 0; dim < node_dim; dim++)
            {
                line += "  " + node_coord[dim + node * node_dim];
            }

            for (int att = 0; att < node_att_num; att++)
            {
                line += "  " + node_att[att + node * node_att_num];
            }

            switch (node_marker_num)
            {
                case 1:
                    line += "  " + node_marker[node];
                    break;
            }

            lines.Add(line);
        }

        File.WriteAllLines(node_file, lines);
    }
}