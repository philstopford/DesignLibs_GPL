using System.Collections.Generic;
using System.IO;

namespace Burkardt
{
    public static class XML
    {
        public static void xml_mesh1d_write(string xml_filename, int m, int node_num,
            double[] node_x, int element_order, int element_num, int[] element_node )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XML_MESH1D_WRITE writes a 1D mesh as a DOLFIN XML file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Anders Logg, Kent-Andre Mardal, Garth Wells,
        //    Automated Solution of Differential Equations by the Finite Element
        //    Method: The FEniCS Book,
        //    Lecture Notes in Computational Science and Engineering,
        //    Springer, 2011,
        //    ISBN13: 978-364223098
        //
        //  Parameters:
        //
        //    Input, string XML_FILENAME, the name of the XML file 
        //    to create.
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, inte NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_X[M*NODE_NUM], the node coordinates.
        //
        //    Input, int ELEMENT_ORDER, the order of the elements.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
        //    the nodes that make up each element.
        //
        {
            int element;
            int node;
            //
            //  Force 0-based indexing.
            //
            Mesh.mesh_base_zero(node_num, element_order, element_num, ref element_node);

            List<string> lines = new List<string>();
            //
            //  Write the data.
            //
            lines.Add("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            lines.Add("");
            lines.Add("<dolfin xmlns:dolfin=\"http://www.fenics.org/dolfin/\">");
            lines.Add("  <mesh celltype=\"interval\" dim=\"" + m + "\">");

            lines.Add("    <vertices size=\"" + node_num + "\">");
            for (node = 0; node < node_num; node++)
            {
                lines.Add("      <vertex index =\"" + node
                    + "\" x =\"" + node_x[0 + node * m] + "\"/>");
            }

            lines.Add("    </vertices>");

            lines.Add("    <cells size=\"" + element_num + "\">");
            for (element = 0; element < element_num; element++)
            {
                lines.Add("      <interval index =\"" + element
                    + "\" v0 =\"" + element_node[0 + element * element_order]
                    + "\" v1 =\"" + element_node[1 + element * element_order] + "\"/>");
            }

            lines.Add("    </cells>");
            lines.Add("  </mesh>");
            lines.Add("</dolfin>");

            File.WriteAllLines(xml_filename, lines);
        }

        public static void xml_mesh2d_write(string xml_filename, int m, int node_num,
            double[] node_x, int element_order, int element_num, int[] element_node )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XML_MESH2D_WRITE writes a 2D mesh as a DOLFIN XML file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Anders Logg, Kent-Andre Mardal, Garth Wells,
        //    Automated Solution of Differential Equations by the Finite Element
        //    Method: The FEniCS Book,
        //    Lecture Notes in Computational Science and Engineering,
        //    Springer, 2011,
        //    ISBN13: 978-364223098
        //
        //  Parameters:
        //
        //    Input, string XML_FILENAME, the name of the XML file 
        //    to create.
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, inte NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_X[M*NODE_NUM], the node coordinates.
        //
        //    Input, int ELEMENT_ORDER, the order of the elements.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
        //    the nodes that make up each element.
        //
        {
            int element;
            int node;
            //
            //  Force 0-based indexing.
            //
            Mesh.mesh_base_zero(node_num, element_order, element_num, ref element_node);

            List<string> lines = new List<string>();
            
            //
            //  Write the data.
            //
            lines.Add("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            lines.Add("");
            lines.Add("<dolfin xmlns:dolfin=\"http://www.fenics.org/dolfin/\">");
            lines.Add("  <mesh celltype=\"triangle\" dim=\"" + m + "\">");

            lines.Add("    <vertices size=\"" + node_num + "\">");
            for (node = 0; node < node_num; node++)
            {
                lines.Add("      <vertex index =\"" + node
                    + "\" x =\"" + node_x[0 + node * m]
                    + "\" y =\"" + node_x[1 + node * m] + "/>");
            }

            lines.Add("    </vertices>");

            lines.Add("    <cells size=\"" + element_num + "\">");
            for (element = 0; element < element_num; element++)
            {
                lines.Add("      <triangle index =\"" + element
                    + "\" v0 =\"" + element_node[0 + element * element_order]
                    + "\" v1 =\"" + element_node[1 + element * element_order]
                    + "\" v2 =\"" + element_node[2 + element * element_order] + "\"/>");
            }

            lines.Add("    </cells>");
            lines.Add("  </mesh>");
            lines.Add("</dolfin>");

            File.WriteAllLines(xml_filename, lines);
        }

        public static void xml_mesh3d_write(string xml_filename, int m, int node_num,
            double[] node_x, int element_order, int element_num, int[] element_node )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XML_MESH3D_WRITE writes a 3D mesh as a DOLFIN XML file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Anders Logg, Kent-Andre Mardal, Garth Wells,
        //    Automated Solution of Differential Equations by the Finite Element
        //    Method: The FEniCS Book,
        //    Lecture Notes in Computational Science and Engineering,
        //    Springer, 2011,
        //    ISBN13: 978-364223098
        //
        //  Parameters:
        //
        //    Input, string XML_FILENAME, the name of the XML file 
        //    to create.
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, inte NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_X[M*NODE_NUM], the node coordinates.
        //
        //    Input, int ELEMENT_ORDER, the order of the elements.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
        //    the nodes that make up each element.
        //
        {
            int element;
            int node;
            //
            //  Force 0-based indexing.
            //
            Mesh.mesh_base_zero(node_num, element_order, element_num, ref element_node);

            List<string> lines = new List<string>();
            
            //
            //  Write the data.
            //
            lines.Add("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            lines.Add("");
            lines.Add("<dolfin xmlns:dolfin=\"http://www.fenics.org/dolfin/\">");
            lines.Add("  <mesh celltype=\"tetrahedron\" dim=\"" + m + "\">");

            lines.Add("    <vertices size=\"" + node_num + "\">");
            for (node = 0; node < node_num; node++)
            {
                lines.Add("      <vertex index =\"" + node
                    + "\" x =\"" + node_x[0 + node * m]
                    + "\" y =\"" + node_x[1 + node * m]
                    + "\" z =\"" + node_x[2 + node * m] + "\"/>");
            }

            lines.Add("    </vertices>");

            lines.Add("    <cells size=\"" + element_num + "\">");
            for (element = 0; element < element_num; element++)
            {
                lines.Add("      <tetrahedron index =\"" + element
                    + "\" v0 =\"" + element_node[0 + element * element_order]
                    + "\" v1 =\"" + element_node[1 + element * element_order]
                    + "\" v2 =\"" + element_node[1 + element * element_order]
                    + "\" v3 =\"" + element_node[2 + element * element_order] + "\"/>");
            }

            lines.Add("    </cells>");
            lines.Add("  </mesh>");
            lines.Add("</dolfin>");

            File.WriteAllLines(xml_filename, lines);
        }

        public static void xml_write ( string xml_filename, int dim_num, int node_num, 
        double[] node_xyz, int element_order, int element_num, int[] element_node )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XML_WRITE writes the triangulation data as a DOLFIN XML mesh file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Anders Logg, Kent-Andre Mardal, Garth Wells,
        //    Automated Solution of Differential Equations by the Finite Element
        //    Method: The FEniCS Book,
        //    Lecture Notes in Computational Science and Engineering,
        //    Springer, 2011,
        //    ISBN13: 978-364223098
        //
        //  Parameters:
        //
        //    Input, string XML_FILENAME, the name of the XML file 
        //    to create.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, inte NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XYZ[3*NODE_NUM], the node coordinates.
        //
        //    Input, int ELEMENT_ORDER, the order of the elements.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
        //    the nodes that make up each element.
        //
        {
            int element;
            List<string> xml = new List<string>();
            int node;
            //  Write the data.
            //
            xml.Add("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            xml.Add("");
            xml.Add("<dolfin xmlns:dolfin=\"http://www.fenics.org/dolfin/\">");
            xml.Add("  <mesh celltype=\"tetrahedron\" dim=\"3\">");

            xml.Add("    <vertices size=\"" + node_num + "\">");
            for (node = 0; node < node_num; node++)
            {
                xml.Add("      <vertex index =\"" + node
                    + "\" x =\"" + node_xyz[0 + node * 3]
                    + "\" y =\"" + node_xyz[1 + node * 3]
                    + "\" z =\"" + node_xyz[2 + node * 3] + "\"/>");
            }

            xml.Add("    </vertices>");

            xml.Add("    <cells size=\"" + element_num + "\">");
            for (element = 0; element < element_num; element++)
            {
                xml.Add("      <tetrahedron index =\"" + element
                    + "\" v0 =\"" + element_node[0 + element * 4]
                    + "\" v1 =\"" + element_node[1 + element * 4]
                    + "\" v2 =\"" + element_node[2 + element * 4]
                    + "\" v3 =\"" + element_node[3 + element * 4] + "\"/>");
            }

            xml.Add("    </cells>");
            xml.Add("  </mesh>");
            xml.Add("</dolfin>");

            File.WriteAllLines(xml_filename, xml);
        }
    }
}