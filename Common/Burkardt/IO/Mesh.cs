using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.IO;

public static class Mesh
{
    public static void cyl248_data(int dim, int vertices, int edges, int triangles,
            int quadrilaterals, int tetrahedrons, int hexahedrons,
            ref double[] vertex_coordinate, ref int[] vertex_label, ref int[] edge_vertex,
            ref int[] edge_label, ref int[] triangle_vertex, ref int[] triangle_label,
            ref int[] quadrilateral_vertex, ref int[] quadrilateral_label,
            ref int[] tetrahedron_vertex, ref int[] tetrahedron_label, ref int[] hexahedron_vertex,
            ref int[] hexahedron_label)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CYL248_DATA defines the data for a 3D tetrahedral mesh.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Pascal Frey,
        //    MEDIT: An interactive mesh visualization software,
        //    Technical Report RT-0253,
        //    Institut National de Recherche en Informatique et en Automatique,
        //    03 December 2001.
        //
        //  Parameters:
        //
        //    Input, int DIM, the spatial dimension, which should be 2 or 3.
        //
        //    Input, int VERTICES, the number of vertices.
        //
        //    Input, int EDGES, the number of edges (may be 0).
        //
        //    Input, int TRIANGLES, the number of triangles (may be 0).
        //
        //    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).
        //
        //    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).
        //
        //    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
        //
        //    Output, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
        //    of each vertex.
        //
        //    Output, int VERTEX_LABEL[VERTICES], a label for each vertex.
        //
        //    Output, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.
        //
        //    Output, int EDGE_LABEL[EDGES], a label for each edge.
        //
        //    Output, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
        //    each triangle.
        //
        //    Output, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.
        //
        //    Output, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
        //    form each quadrilateral.
        //
        //    Output, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
        //    each quadrilateral.
        //
        //    Output, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
        //    form each tetrahedron.
        //
        //    Output, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
        //    each tetrahedron.
        //
        //    Output, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
        //    each hexahedron.
        //
        //    Output, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
        //
    {
        int i;

        int[] tetrahedron_vertex_save =
        {
            23, 1, 9, 8,
            27, 9, 23, 1,
            26, 8, 23, 9,
            26, 9, 7, 8,
            2, 9, 27, 1,
            26, 9, 10, 7,
            26, 28, 7, 10,
            11, 29, 3, 2,
            7, 6, 10, 28,
            10, 6, 31, 28,
            11, 29, 30, 3,
            11, 30, 4, 3,
            11, 30, 32, 4,
            10, 6, 5, 31,
            11, 5, 4, 32,
            19, 33, 34, 20,
            39, 22, 40, 16,
            39, 17, 36, 22,
            39, 22, 16, 17,
            40, 22, 15, 16,
            12, 19, 20, 33,
            19, 20, 34, 18,
            12, 33, 20, 35,
            38, 37, 14, 21,
            36, 22, 17, 18,
            38, 14, 15, 21,
            13, 14, 37, 21,
            12, 20, 13, 35,
            80, 32, 11, 30,
            80, 28, 10, 31,
            80, 31, 59, 28,
            80, 58, 57, 26,
            80, 28, 58, 26,
            80, 59, 58, 28,
            80, 28, 26, 10,
            80, 10, 26, 9,
            80, 9, 11, 10,
            80, 9, 26, 23,
            80, 23, 26, 57,
            80, 23, 27, 9,
            80, 23, 56, 27,
            80, 30, 11, 29,
            80, 5, 10, 11,
            80, 5, 11, 32,
            80, 5, 32, 31,
            80, 31, 10, 5,
            80, 2, 11, 9,
            80, 29, 11, 2,
            80, 2, 9, 27,
            80, 27, 29, 2,
            81, 40, 39, 22,
            81, 22, 39, 36,
            81, 18, 36, 34,
            81, 34, 20, 18,
            81, 22, 36, 18,
            81, 20, 22, 18,
            81, 37, 38, 21,
            81, 20, 33, 35,
            81, 13, 21, 20,
            81, 13, 20, 35,
            81, 13, 37, 21,
            81, 35, 37, 13,
            81, 20, 21, 22,
            81, 34, 33, 20,
            81, 21, 38, 15,
            81, 38, 40, 15,
            81, 22, 21, 15,
            81, 15, 40, 22,
            82, 60, 74, 59,
            82, 74, 25, 59,
            82, 73, 72, 58,
            82, 25, 73, 58,
            82, 59, 25, 58,
            82, 58, 72, 57,
            82, 57, 80, 58,
            82, 58, 80, 59,
            83, 71, 79, 70,
            83, 70, 76, 78,
            83, 79, 76, 70,
            83, 79, 60, 76,
            83, 82, 60, 74,
            84, 54, 64, 55,
            84, 64, 65, 55,
            84, 65, 63, 55,
            84, 65, 71, 63,
            85, 29, 62, 30,
            85, 80, 29, 30,
            85, 29, 61, 62,
            85, 78, 83, 76,
            85, 78, 76, 30,
            85, 62, 78, 30,
            85, 76, 83, 60,
            85, 76, 32, 30,
            85, 32, 80, 30,
            85, 32, 76, 60,
            85, 27, 61, 29,
            85, 80, 27, 29,
            85, 83, 82, 60,
            85, 77, 78, 62,
            85, 60, 82, 59,
            85, 59, 82, 80,
            85, 32, 60, 31,
            85, 80, 32, 31,
            85, 60, 59, 31,
            85, 59, 80, 31,
            86, 51, 68, 52,
            86, 69, 68, 51,
            86, 68, 67, 52,
            86, 52, 67, 53,
            86, 67, 66, 53,
            86, 53, 66, 54,
            87, 50, 70, 49,
            87, 71, 70, 50,
            87, 63, 71, 50,
            87, 63, 84, 71,
            87, 70, 69, 49,
            87, 71, 83, 70,
            87, 49, 69, 51,
            87, 69, 86, 51,
            88, 64, 66, 73,
            88, 72, 73, 66,
            88, 72, 82, 73,
            88, 24, 72, 66,
            88, 64, 73, 25,
            88, 73, 82, 25,
            88, 66, 64, 54,
            88, 84, 54, 64,
            88, 87, 86, 84,
            88, 67, 24, 66,
            88, 66, 86, 67,
            88, 64, 25, 65,
            88, 65, 84, 64,
            88, 25, 74, 65,
            88, 25, 82, 74,
            88, 83, 87, 71,
            88, 71, 87, 84,
            88, 82, 83, 74,
            88, 74, 83, 71,
            88, 65, 74, 71,
            88, 71, 84, 65,
            89, 86, 87, 84,
            89, 39, 48, 44,
            89, 44, 49, 43,
            89, 44, 43, 36,
            89, 44, 48, 50,
            89, 48, 63, 50,
            89, 86, 84, 54,
            89, 51, 87, 86,
            89, 44, 50, 49,
            89, 50, 87, 49,
            89, 43, 49, 51,
            89, 49, 87, 51,
            89, 39, 44, 36,
            89, 36, 81, 39,
            89, 63, 48, 47,
            89, 47, 48, 40,
            89, 46, 55, 47,
            89, 38, 46, 47,
            89, 55, 63, 47,
            89, 55, 84, 63,
            89, 43, 42, 34,
            89, 43, 51, 42,
            89, 45, 53, 54,
            89, 53, 86, 54,
            89, 45, 54, 46,
            89, 42, 52, 41,
            89, 41, 52, 53,
            89, 52, 86, 53,
            89, 42, 51, 52,
            89, 51, 86, 52,
            89, 46, 54, 55,
            89, 54, 84, 55,
            90, 56, 75, 61,
            90, 24, 75, 56,
            90, 27, 56, 61,
            90, 61, 85, 27,
            90, 75, 77, 61,
            90, 80, 82, 57,
            90, 85, 82, 80,
            90, 57, 24, 56,
            90, 72, 24, 57,
            90, 57, 82, 72,
            90, 80, 56, 27,
            90, 85, 80, 27,
            91, 85, 90, 77,
            91, 86, 87, 69,
            91, 78, 77, 69,
            91, 83, 88, 82,
            91, 90, 82, 88,
            91, 67, 88, 86,
            91, 88, 87, 86,
            91, 87, 88, 83,
            91, 83, 85, 78,
            91, 78, 85, 77,
            91, 77, 75, 68,
            91, 77, 90, 75,
            91, 69, 77, 68,
            91, 68, 86, 69,
            91, 68, 75, 67,
            91, 67, 86, 68,
            91, 24, 88, 67,
            91, 90, 88, 24,
            91, 69, 87, 70,
            91, 87, 83, 70,
            91, 75, 24, 67,
            91, 75, 90, 24,
            92, 89, 46, 45,
            92, 41, 53, 45,
            92, 89, 45, 53,
            92, 89, 53, 41,
            92, 89, 41, 42,
            92, 35, 41, 45,
            92, 33, 41, 35,
            92, 35, 81, 33,
            92, 35, 45, 37,
            92, 81, 35, 37,
            92, 34, 89, 42,
            92, 81, 89, 34,
            92, 33, 42, 41,
            92, 37, 45, 46,
            92, 37, 46, 38,
            92, 81, 37, 38,
            92, 33, 34, 42,
            92, 33, 81, 34,
            83, 74, 60, 71,
            83, 60, 79, 71,
            89, 39, 40, 48,
            89, 39, 81, 40,
            89, 36, 43, 34,
            89, 34, 81, 36,
            89, 63, 87, 50,
            89, 84, 87, 63,
            54, 88, 66, 86,
            54, 88, 86, 84,
            90, 72, 88, 24,
            90, 82, 88, 72,
            38, 47, 89, 40,
            38, 89, 81, 40,
            92, 46, 89, 38,
            92, 89, 81, 38,
            80, 23, 57, 56,
            80, 57, 90, 56,
            61, 85, 62, 77,
            61, 90, 85, 77,
            82, 85, 91, 83,
            82, 90, 91, 85,
            70, 91, 78, 83,
            70, 78, 91, 69
        };
        int[] triangle_label_save =
        {
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3
        };
        int[] triangle_vertex_save =
        {
            12, 20, 19,
            12, 13, 20,
            19, 20, 18,
            20, 22, 18,
            22, 17, 18,
            22, 16, 17,
            13, 21, 20,
            13, 14, 21,
            14, 15, 21,
            22, 15, 16,
            22, 21, 15,
            20, 21, 22,
            1, 9, 8,
            2, 9, 1,
            9, 7, 8,
            2, 11, 9,
            11, 2, 3,
            11, 3, 4,
            9, 10, 7,
            7, 10, 6,
            10, 5, 6,
            11, 4, 5,
            5, 10, 11,
            9, 11, 10,
            23, 1, 8,
            26, 23, 8,
            26, 8, 7,
            27, 1, 23,
            2, 1, 27,
            26, 7, 28,
            7, 6, 28,
            27, 29, 2,
            29, 3, 2,
            29, 30, 3,
            30, 4, 3,
            6, 31, 28,
            6, 5, 31,
            5, 32, 31,
            5, 4, 32,
            12, 19, 33,
            19, 34, 33,
            19, 18, 34,
            12, 33, 35,
            12, 35, 13,
            18, 36, 34,
            36, 18, 17,
            35, 37, 13,
            13, 37, 14,
            38, 14, 37,
            38, 15, 14,
            39, 36, 17,
            39, 17, 16,
            38, 40, 15,
            40, 16, 15,
            39, 16, 40,
            33, 41, 35,
            33, 42, 41,
            33, 34, 42,
            36, 43, 34,
            43, 42, 34,
            39, 44, 36,
            44, 43, 36,
            35, 45, 37,
            35, 41, 45,
            37, 46, 38,
            37, 45, 46,
            38, 47, 40,
            38, 46, 47,
            39, 48, 44,
            39, 40, 48,
            47, 48, 40,
            44, 49, 43,
            44, 50, 49,
            44, 48, 50,
            43, 51, 42,
            43, 49, 51,
            42, 52, 41,
            42, 51, 52,
            41, 53, 45,
            41, 52, 53,
            45, 54, 46,
            45, 53, 54,
            46, 55, 47,
            46, 54, 55,
            30, 32, 4,
            23, 56, 27,
            23, 57, 56,
            23, 26, 57,
            28, 58, 26,
            58, 57, 26,
            31, 59, 28,
            59, 58, 28,
            32, 60, 31,
            60, 59, 31,
            27, 61, 29,
            27, 56, 61,
            29, 62, 30,
            29, 61, 62,
            55, 63, 47,
            63, 48, 47,
            48, 63, 50,
            54, 64, 55,
            64, 65, 55,
            65, 63, 55,
            53, 66, 54,
            66, 64, 54,
            52, 67, 53,
            67, 66, 53,
            51, 68, 52,
            68, 67, 52,
            49, 69, 51,
            69, 68, 51,
            50, 70, 49,
            70, 69, 49,
            63, 71, 50,
            71, 70, 50,
            65, 71, 63,
            64, 25, 65,
            64, 73, 25,
            64, 66, 73,
            67, 24, 66,
            24, 72, 66,
            72, 73, 66,
            68, 75, 67,
            75, 24, 67,
            69, 77, 68,
            77, 75, 68,
            70, 78, 69,
            78, 77, 69,
            62, 78, 30,
            78, 76, 30,
            76, 32, 30,
            32, 76, 60,
            61, 77, 62,
            77, 78, 62,
            56, 75, 61,
            75, 77, 61,
            57, 24, 56,
            24, 75, 56,
            58, 72, 57,
            72, 24, 57,
            59, 25, 58,
            25, 73, 58,
            73, 72, 58,
            60, 74, 59,
            74, 25, 59,
            25, 74, 65,
            65, 74, 71,
            70, 76, 78,
            71, 79, 70,
            79, 76, 70,
            79, 60, 76,
            74, 60, 71,
            60, 79, 71
        };
        int[] vertex_label_save =
        {
            3, 3, 3, 3, 3, 3, 3, 3, 2, 2,
            2, 3, 3, 3, 3, 3, 3, 3, 3, 4,
            4, 4, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0
        };
        double[] vertex_coordinate_save =
        {
            1.0, 0.2, 0.0,
            1.0, 0.141421, 0.141421,
            1.0, 0.0, 0.2,
            1.0, -0.141421, 0.141421,
            1.0, -0.2, 0.0,
            1.0, -0.141421, -0.141421,
            1.0, 0.0, -0.2,
            1.0, 0.141421, -0.141421,
            1.0, 0.066163, -0.0302872,
            1.0, -0.0615154, -0.0610739,
            1.0, -0.0306985, 0.0668017,
            0.0, 0.2, 0.0,
            0.0, 0.141421, -0.141421,
            0.0, 0.0, -0.2,
            0.0, -0.141421, -0.141421,
            0.0, -0.2, 0.0,
            0.0, -0.141421, 0.141421,
            0.0, 0.0, 0.2,
            0.0, 0.141421, 0.141421,
            0.0, 0.0686748, 0.0255359,
            0.0, 0.0, -0.0865993,
            0.0, -0.0686749, 0.0255359,
            0.8816, 0.185522, -0.0747102,
            0.642415, 0.187806, -0.0687668,
            0.627606, -0.0696445, -0.187482,
            0.876431, 0.0811908, -0.182779,
            0.881613, 0.186118, 0.0732131,
            0.872048, -0.0699008, -0.187387,
            0.878318, 0.0844232, 0.181308,
            0.845861, -0.0716063, 0.186742,
            0.866503, -0.182493, -0.0818307,
            0.859402, -0.186751, 0.0715813,
            0.131355, 0.18477, 0.0765501,
            0.13317, 0.077694, 0.184292,
            0.130862, 0.185301, -0.0752567,
            0.135181, -0.0749468, 0.185426,
            0.130839, 0.0781729, -0.18409,
            0.131856, -0.0754694, -0.185214,
            0.135683, -0.184121, 0.0780993,
            0.134207, -0.184959, -0.0760928,
            0.261923, 0.199982, 0.00264585,
            0.263928, 0.144161, 0.138627,
            0.268645, 0.00535339, 0.199928,
            0.272346, -0.137646, 0.145098,
            0.26108, 0.144683, -0.138082,
            0.260772, 0.00498797, -0.199938,
            0.264253, -0.139152, -0.143655,
            0.270288, -0.199962, 0.00389323,
            0.408181, -0.0730357, 0.186187,
            0.411818, -0.184374, 0.0774991,
            0.397539, 0.080738, 0.182979,
            0.39192, 0.185619, 0.0744699,
            0.392192, 0.184438, -0.0773479,
            0.389194, 0.0770141, -0.184577,
            0.38786, -0.0747817, -0.185493,
            0.762413, 0.199986, -0.0023425,
            0.762987, 0.151152, -0.13097,
            0.741526, 0.0187858, -0.199116,
            0.746899, -0.128364, -0.153371,
            0.720076, -0.19917, -0.0182053,
            0.7628, 0.152219, 0.129728,
            0.763882, 0.0434475, 0.195224,
            0.399903, -0.1841, -0.0781489,
            0.506331, -0.00579066, -0.199916,
            0.514514, -0.133894, -0.148568,
            0.526121, 0.135152, -0.147424,
            0.517967, 0.199953, -0.0043215,
            0.520585, 0.147847, 0.13469,
            0.533956, 0.0124181, 0.199614,
            0.558316, -0.136902, 0.145801,
            0.549126, -0.199624, -0.0122659,
            0.657307, 0.117735, -0.161674,
            0.611189, 0.041829, -0.195577,
            0.631917, -0.164669, -0.113508,
            0.641444, 0.187001, 0.0709267,
            0.720251, -0.155557, 0.125706,
            0.647345, 0.0932963, 0.176906,
            0.677484, -0.0430068, 0.195321,
            0.635293, -0.188734, 0.0661777,
            0.888023, -0.00868364, -0.00818647,
            0.112146, 0.0, -0.0118425,
            0.676228, 0.0124197, -0.0856487,
            0.638436, -0.0639898, 0.0525795,
            0.452586, -0.0410297, -0.0704842,
            0.762004, -0.0188614, 0.0693717,
            0.463368, 0.0649048, 0.0262133,
            0.473921, -0.0356443, 0.0388516,
            0.557002, 0.0123705, -0.0932599,
            0.290986, -0.0200898, 0.00857934,
            0.7038, 0.0856777, 0.0182744,
            0.576134, 0.0436218, 0.0828782,
            0.215187, 0.080855, -0.0314946
        };

        typeMethods.r8vec_copy(dim * vertices, vertex_coordinate_save, ref vertex_coordinate);

        typeMethods.i4vec_copy(vertices, vertex_label_save, ref vertex_label);

        typeMethods.i4vec_copy(3 * triangles, triangle_vertex_save, ref triangle_vertex);

        typeMethods.i4vec_copy(triangles, triangle_label_save, ref triangle_label);

        typeMethods.i4vec_copy(4 * tetrahedrons, tetrahedron_vertex_save, ref tetrahedron_vertex);

        for (i = 0; i < tetrahedrons; i++)
        {
            tetrahedron_label[i] = 1;
        }
    }

    public static void cyl248_size(ref int dim, ref int vertices, ref int edges, ref int triangles,
            ref int quadrilaterals, ref int tetrahedrons, ref int hexahedrons)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CYL248_SIZE defines the sizes for a 3D tetrahedral mesh.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Pascal Frey,
        //    MEDIT: An interactive mesh visualization software,
        //    Technical Report RT-0253,
        //    Institut National de Recherche en Informatique et en Automatique,
        //    03 December 2001.
        //
        //  Parameters:
        //
        //    Output, int *DIM, the spatial dimension, which should be 2 or 3.
        //
        //    Output, int *VERTICES, the number of vertices.
        //
        //    Output, int *EDGES, the number of edges (may be 0).
        //
        //    Output, int *TRIANGLES, the number of triangles (may be 0).
        //
        //    Output, int *QUADRILATERALS, the number of quadrilaterals (may be 0).
        //
        //    Output, int *TETRAHEDRONS, the number of tetrahedrons (may be 0).
        //
        //    Output, int *HEXAHEDRONS, the number of hexahedrons (may be 0).
        //
    {
        dim = 3;
        vertices = 92;
        edges = 0;
        triangles = 154;
        quadrilaterals = 0;
        tetrahedrons = 248;
        hexahedrons = 0;

    }

    public static void hexahexa_2x2x2_data(int dim, int vertices, int edges, int triangles,
            int quadrilaterals, int tetrahedrons, int hexahedrons,
            ref double[] vertex_coordinate, ref int[] vertex_label, ref int[] edge_vertex,
            ref int[] edge_label, ref int[] triangle_vertex, ref int[] triangle_label,
            ref int[] quadrilateral_vertex, ref int[] quadrilateral_label,
            ref int[] tetrahedron_vertex, ref int[] tetrahedron_label, ref int[] hexahedron_vertex,
            ref int[] hexahedron_label)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HEXAHEXA_2X2X2_DATA defines the data for a 3D hexahedral mesh.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Pascal Frey,
        //    MEDIT: An interactive mesh visualization software,
        //    Technical Report RT-0253,
        //    Institut National de Recherche en Informatique et en Automatique,
        //    03 December 2001.
        //
        //  Parameters:
        //
        //    Input, int DIM, the spatial dimension, which should be 2 or 3.
        //
        //    Input, int VERTICES, the number of vertices.
        //
        //    Input, int EDGES, the number of edges (may be 0).
        //
        //    Input, int TRIANGLES, the number of triangles (may be 0).
        //
        //    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).
        //
        //    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).
        //
        //    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
        //
        //    Output, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
        //    of each vertex.
        //
        //    Output, int VERTEX_LABEL[VERTICES], a label for each vertex.
        //
        //    Output, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.
        //
        //    Output, int EDGE_LABEL[EDGES], a label for each edge.
        //
        //    Output, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
        //    each triangle.
        //
        //    Output, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.
        //
        //    Output, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
        //    form each quadrilateral.
        //
        //    Output, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
        //    each quadrilateral.
        //
        //    Output, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
        //    form each tetrahedron.
        //
        //    Output, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
        //    each tetrahedron.
        //
        //    Output, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
        //    each hexahedron.
        //
        //    Output, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
        //
    {
        int[] hexahedron_label_save =
        {
            1, 1, 1, 1, 1, 1, 1, 1
        };
        int[] hexahedron_vertex_save =
        {
            1, 2, 5, 4, 10, 11, 14, 13,
            2, 3, 6, 5, 11, 12, 15, 14,
            4, 5, 8, 7, 13, 14, 17, 16,
            5, 6, 9, 8, 14, 15, 18, 17,
            10, 11, 14, 13, 19, 20, 23, 22,
            11, 12, 15, 14, 20, 21, 24, 23,
            13, 14, 17, 16, 22, 23, 26, 25,
            14, 15, 18, 17, 23, 24, 27, 26
        };
        int[] quadrilateral_label_save =
        {
            1, 1, 1, 1, 2, 2, 2, 2, 3, 3,
            3, 3, 4, 4, 4, 4, 5, 5, 5, 5,
            6, 6, 6, 6
        };
        int[] quadrilateral_vertex_save =
        {
            1, 4, 5, 2,
            2, 5, 6, 3,
            4, 7, 8, 5,
            5, 8, 9, 6,
            1, 2, 11, 10,
            2, 3, 12, 11,
            10, 11, 20, 19,
            11, 12, 21, 20,
            3, 6, 15, 12,
            6, 9, 18, 15,
            12, 15, 24, 21,
            15, 18, 27, 24,
            7, 16, 17, 8,
            8, 17, 18, 9,
            16, 25, 26, 17,
            17, 26, 27, 18,
            1, 10, 13, 4,
            4, 13, 16, 7,
            10, 19, 22, 13,
            13, 22, 25, 16,
            19, 20, 23, 22,
            20, 21, 24, 23,
            22, 23, 26, 25,
            23, 24, 27, 26
        };
        int[] vertex_label_save =
        {
            5, 2, 3, 5, 1, 3, 5, 4, 4, 5,
            2, 3, 5, 0, 3, 5, 4, 4, 6, 6,
            6, 6, 6, 6, 6, 6, 6
        };
        double[] vertex_coordinate_save =
        {
            0.0, 0.0, 0.0,
            0.5, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 0.5, 0.0,
            0.5, 0.5, 0.0,
            1.0, 0.5, 0.0,
            0.0, 1.0, 0.0,
            0.5, 1.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 0.0, 0.5,
            0.5, 0.0, 0.5,
            1.0, 0.0, 0.5,
            0.0, 0.5, 0.5,
            0.5, 0.5, 0.5,
            1.0, 0.5, 0.5,
            0.0, 1.0, 0.5,
            0.5, 1.0, 0.5,
            1.0, 1.0, 0.5,
            0.0, 0.0, 1.0,
            0.5, 0.0, 1.0,
            1.0, 0.0, 1.0,
            0.0, 0.5, 1.0,
            0.5, 0.5, 1.0,
            1.0, 0.5, 1.0,
            0.0, 1.0, 1.0,
            0.5, 1.0, 1.0,
            1.0, 1.0, 1.0
        };

        typeMethods.r8vec_copy(dim * vertices, vertex_coordinate_save, ref vertex_coordinate);
        typeMethods.i4vec_copy(vertices, vertex_label_save, ref vertex_label);
        typeMethods.i4vec_copy(4 * quadrilaterals, quadrilateral_vertex_save, ref quadrilateral_vertex);
        typeMethods.i4vec_copy(quadrilaterals, quadrilateral_label_save, ref quadrilateral_label);
        typeMethods.i4vec_copy(8 * hexahedrons, hexahedron_vertex_save, ref hexahedron_vertex);
        typeMethods.i4vec_copy(hexahedrons, hexahedron_label_save, ref hexahedron_label);

    }

    public static void hexahexa_2x2x2_size(ref int dim, ref int vertices, ref int edges, ref int triangles,
            ref int quadrilaterals, ref int tetrahedrons, ref int hexahedrons)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HEXAHEXA_2X2X2_SIZE defines the sizes for a 3D hexahedral mesh.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Pascal Frey,
        //    MEDIT: An interactive mesh visualization software,
        //    Technical Report RT-0253,
        //    Institut National de Recherche en Informatique et en Automatique,
        //    03 December 2001.
        //
        //  Parameters:
        //
        //    Output, int *DIM, the spatial dimension, which should be 2 or 3.
        //
        //    Output, int *VERTICES, the number of vertices.
        //
        //    Output, int *EDGES, the number of edges (may be 0).
        //
        //    Output, int *TRIANGLES, the number of triangles (may be 0).
        //
        //    Output, int *QUADRILATERALS, the number of quadrilaterals (may be 0).
        //
        //    Output, int *TETRAHEDRONS, the number of tetrahedrons (may be 0).
        //
        //    Output, int *HEXAHEDRONS, the number of hexahedrons (may be 0).
        //
    {
        dim = 3;
        vertices = 27;
        edges = 0;
        triangles = 0;
        quadrilaterals = 24;
        tetrahedrons = 0;
        hexahedrons = 8;
    }

    public static void mesh_data_print(int dim, int vertices, int edges, int triangles,
            int quadrilaterals, int tetrahedrons, int hexahedrons,
            double[] vertex_coordinate, int[] vertex_label, int[] edge_vertex,
            int[] edge_label, int[] triangle_vertex, int[] triangle_label,
            int[] quadrilateral_vertex, int[] quadrilateral_label,
            int[] tetrahedron_vertex, int[] tetrahedron_label,
            int[] hexahedron_vertex, int[] hexahedron_label)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MESH_DATA_PRINT prints the data of a MESH dataset.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Pascal Frey,
        //    MEDIT: An interactive mesh visualization software,
        //    Technical Report RT-0253,
        //    Institut National de Recherche en Informatique et en Automatique,
        //    03 December 2001.
        //
        //  Parameters:
        //
        //    Input, int DIM, the spatial dimension, which should be 2 or 3.
        //
        //    Input, int VERTICES, the number of vertices.
        //
        //    Input, int EDGES, the number of edges (may be 0).
        //
        //    Input, int TRIANGLES, the number of triangles (may be 0).
        //
        //    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).
        //
        //    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).
        //
        //    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
        //
        //    Input, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
        //    of each vertex.
        //
        //    Input, int VERTEX_LABEL[VERTICES], a label for each vertex.
        //
        //    Input, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.
        //
        //    Input, int EDGE_LABEL[EDGES], a label for each edge.
        //
        //    Input, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
        //    each triangle.
        //
        //    Input, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.
        //
        //    Input, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
        //    form each quadrilateral.
        //
        //    Input, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
        //    each quadrilateral.
        //
        //    Input, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
        //    form each tetrahedron.
        //
        //    Input, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
        //    each tetrahedron.
        //
        //    Input, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
        //    each hexahedron.
        //
        //    Input, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
        //
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("  Vertices:");
        Console.WriteLine("");
        for (j = 0; j < vertices; j++)
        {
            string cout = "";
            for (i = 0; i < dim; i++)
            {
                cout += "  " + vertex_coordinate[i + j * dim];
            }

            Console.WriteLine(cout + "  (" + vertex_label[j] + ")");
        }

        switch (edges)
        {
            case > 0:
            {
                Console.WriteLine("");
                Console.WriteLine("  Edges:");
                Console.WriteLine("");
                for (j = 0; j < edges; j++)
                {
                    string cout = "";
                    for (i = 0; i < 2; i++)
                    {
                        cout += "  " + edge_vertex[i + j * 2];
                    }

                    Console.WriteLine(cout + "  (" + edge_label[j] + ")");
                }

                break;
            }
        }

        switch (triangles)
        {
            case > 0:
            {
                Console.WriteLine("");
                Console.WriteLine("  Triangles:");
                Console.WriteLine("");
                for (j = 0; j < triangles; j++)
                {
                    string cout = "";
                    for (i = 0; i < 3; i++)
                    {
                        cout += "  " + triangle_vertex[i + j * 3];
                    }

                    Console.WriteLine(cout + "  (" + triangle_label[j] + ")");
                }

                break;
            }
        }

        switch (quadrilaterals)
        {
            case > 0:
            {
                Console.WriteLine("");
                Console.WriteLine("  Quadrilaterals:");
                Console.WriteLine("");
                for (j = 0; j < quadrilaterals; j++)
                {
                    string cout = "";
                    for (i = 0; i < 4; i++)
                    {
                        cout += "  " + quadrilateral_vertex[i + j * 4];
                    }

                    Console.WriteLine(cout + "  (" + quadrilateral_label[j] + ")");
                }

                break;
            }
        }

        switch (tetrahedrons)
        {
            case > 0:
            {
                Console.WriteLine("");
                Console.WriteLine("  Tetrahedrons:");
                Console.WriteLine("");
                for (j = 0; j < tetrahedrons; j++)
                {
                    string cout = "";
                    for (i = 0; i < 4; i++)
                    {
                        cout += "  " + tetrahedron_vertex[i + j * 4];
                    }

                    Console.WriteLine(cout + "  (" + tetrahedron_label[j] + ")");
                }

                break;
            }
        }

        switch (hexahedrons)
        {
            case > 0:
            {
                Console.WriteLine("");
                Console.WriteLine("  Hexahedrons:");
                Console.WriteLine("");
                for (j = 0; j < hexahedrons; j++)
                {
                    string cout = "";
                    for (i = 0; i < 8; i++)
                    {
                        cout += "  " + hexahedron_vertex[i + j * 8];
                    }

                    Console.WriteLine(cout + "  (" + hexahedron_label[j] + ")");
                }

                break;
            }
        }
    }

    public static void mesh_data_read(string filename, int dim, int vertices, int edges,
            int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
            ref double[] vertex_coordinate, ref int[] vertex_label, ref int[] edge_vertex,
            ref int[] edge_label, ref int[] triangle_vertex, ref int[] triangle_label,
            ref int[] quadrilateral_vertex, ref int[] quadrilateral_label,
            ref int[] tetrahedron_vertex, ref int[] tetrahedron_label,
            ref int[] hexahedron_vertex, ref int[] hexahedron_label)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MESH_DATA_READ reads data from a MESH file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Pascal Frey,
        //    MEDIT: An interactive mesh visualization software,
        //    Technical Report RT-0253,
        //    Institut National de Recherche en Informatique et en Automatique,
        //    03 December 2001.
        //
        //  Parameters:
        //
        //    Input, string FILENAME, the name of the MESH file.
        //
        //    Input, int DIM, the spatial dimension, which should be 2 or 3.
        //
        //    Input, int VERTICES, the number of vertices.
        //
        //    Input, int EDGES, the number of edges (may be 0).
        //
        //    Input, int TRIANGLES, the number of triangles (may be 0).
        //
        //    Input, int QUADRILATERALS, the number of quadrilaterals
        //    (may be 0).
        //
        //    Input, int TETRAHEDRONS, the number of tetrahedrons
        //    (may be 0).
        //
        //    Input, int HEXAHEDRONS, the number of hexahedrons
        //    (may be 0).
        //
        //    Output, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
        //    of each vertex.
        //
        //    Output, int VERTEX_LABEL[VERTICES], a label for each vertex.
        //
        //    Output, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.
        //
        //    Output, int EDGE_LABEL[EDGES], a label for each edge.
        //
        //    Output, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
        //    each triangle.
        //
        //    Output, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.
        //
        //    Output, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
        //    form each quadrilateral.
        //
        //    Output, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
        //    each quadrilateral.
        //
        //    Output, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
        //    form each tetrahedron.
        //
        //    Output, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
        //    each tetrahedron.
        //
        //    Output, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
        //    each hexahedron.
        //
        //    Output, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
        //
    {
        int edge = 0;
        int hexahedron = 0;
        string[] input;
        int quadrilateral = 0;
        int tetrahedron = 0;
        int triangle = 0;
        int vertex = 0;
        //
        //  Initialize everything to nothing.
        //
        typeMethods.i4vec_zero(edges, ref edge_label);
        typeMethods.i4vec_zero(2 * edges, ref edge_vertex);
        typeMethods.i4vec_zero(hexahedrons, ref hexahedron_label);
        typeMethods.i4vec_zero(8 * hexahedrons, ref hexahedron_vertex);
        typeMethods.i4vec_zero(quadrilaterals, ref quadrilateral_label);
        typeMethods.i4vec_zero(4 * quadrilaterals, ref quadrilateral_vertex);
        typeMethods.i4vec_zero(tetrahedrons, ref tetrahedron_label);
        typeMethods.i4vec_zero(4 * tetrahedrons, ref tetrahedron_vertex);
        typeMethods.i4vec_zero(triangles, ref triangle_label);
        typeMethods.i4vec_zero(3 * triangles, ref triangle_vertex);
        typeMethods.r8vec_zero(dim * vertices, ref vertex_coordinate);
        typeMethods.i4vec_zero(vertices, ref vertex_label);

        try
        {
            input = File.ReadAllLines(filename);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("MESH_DATA_READ - Fatal error!");
            Console.WriteLine("  Could not open file.");
            return;
        }

        //
        //  Read lines til you get alphanumerics and determine a "mode"
        //
        int line_num = 0;
        string keyword = "NONE";

        for (;;)
        {
            string text;
            try
            {
                text = input[line_num];
            }
            catch (Exception)
            {
                break;
            }

            line_num += 1;

            if (typeMethods.s_len_trim(text) == 0)
            {
                keyword = "NONE";
                continue;
            }

            switch (text[0])
            {
                case '#':
                    continue;
            }
            //
            //  Remove initial blanks.
            //

            //
            //  Expecting a keyword.
            //
            if (typeMethods.s_eqi(text, "CORNERS"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "DIMENSION"))
            {
                keyword = "DIMENSION";
            }
            else if (typeMethods.s_eqi(text, "EDGES"))
            {
                keyword = "EDGES";
            }
            else if (typeMethods.s_eqi(text, "END"))
            {
                Console.WriteLine("");
                Console.WriteLine("  END statement encountered.");
                break;
            }
            else if (typeMethods.s_eqi(text, "HEXAHEDRA") ||
                     typeMethods.s_eqi(text, "HEXAHEDRONS"))
            {
                keyword = "HEXAHEDRONS";
            }
            else if (text.StartsWith("MESHVERSIONFORMATTED"))
            {
            }
            else if (typeMethods.s_eqi(text, "NORMALATQUADRILATERALVERTICES"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "NORMALATTRIANGLEVERTICES"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "NORMALATVERTICES"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "NORMALS"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "QUADRILATERALS"))
            {
                keyword = "QUADRILATERALS";
            }
            else if (typeMethods.s_eqi(text, "REQUIREDEDGES"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "REQUIREDVERTICES"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "RIDGES"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "TANGENTATEDGES"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "TANGENTS"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "TETRAHEDRA") ||
                     typeMethods.s_eqi(text, "TETRAHEDRONS"))
            {
                keyword = "TETRAHEDRONS";
            }
            else if (typeMethods.s_eqi(text, "TRIANGLES"))
            {
                keyword = "TRIANGLES";
            }
            else if (typeMethods.s_eqi(text, "VERTICES"))
            {
                keyword = "VERTICES";
            }
            //
            //  Presumably, numeric data to be processed by keyword.
            //
            else if (typeMethods.s_eqi(keyword, "DIMENSION"))
            {
                keyword = "NONE";
            }
            else if (typeMethods.s_eqi(keyword, "EDGES"))
            {
                keyword = "EDGE_VERTEX";
                edge = 0;
            }
            else
            {
                int i;
                int[] i4vec_;
                if (typeMethods.s_eqi(keyword, "EDGE_VERTEX"))
                {
                    i4vec_ = typeMethods.s_to_i4vec(text, 3).ivec;
                    for (i = 0; i < 2; i++)
                    {
                        edge_vertex[i + edge * 2] = i4vec_[i];
                    }

                    edge_label[edge] = i4vec_[2];
                    edge += 1;
                }
                else if (typeMethods.s_eqi(keyword, "HEXAHEDRONS"))
                {
                    keyword = "HEXAHEDRON_VERTEX";
                    hexahedron = 0;
                }
                else if (typeMethods.s_eqi(keyword, "HEXAHEDRON_VERTEX"))
                {
                    i4vec_ = typeMethods.s_to_i4vec(text, 9).ivec;
                    for (i = 0; i < 8; i++)
                    {
                        hexahedron_vertex[i + hexahedron * 8] = i4vec_[i];
                    }

                    hexahedron_label[hexahedron] = i4vec_[8];
                    hexahedron += 1;
                }
                else if (typeMethods.s_eqi(keyword, "QUADRILATERALS"))
                {
                    keyword = "QUADRILATERAL_VERTEX";
                    quadrilateral = 0;
                }
                else if (typeMethods.s_eqi(keyword, "QUADRILATERAL_VERTEX"))
                {
                    i4vec_ = typeMethods.s_to_i4vec(text, 5).ivec;
                    for (i = 0; i < 4; i++)
                    {
                        quadrilateral_vertex[i + quadrilateral * 4] = i4vec_[i];
                    }

                    quadrilateral_label[quadrilateral] = i4vec_[4];
                    quadrilateral += 1;
                }
                else if (typeMethods.s_eqi(keyword, "TETRAHEDRONS"))
                {
                    keyword = "TETRAHEDRON_VERTEX";
                    tetrahedron = 0;
                }
                else if (typeMethods.s_eqi(keyword, "TETRAHEDRON_VERTEX"))
                {
                    i4vec_ = typeMethods.s_to_i4vec(text, 5).ivec;
                    for (i = 0; i < 4; i++)
                    {
                        tetrahedron_vertex[i + tetrahedron * 4] = i4vec_[i];
                    }

                    tetrahedron_label[tetrahedron] = i4vec_[4];
                    tetrahedron += 1;
                }
                else if (typeMethods.s_eqi(keyword, "TRIANGLES"))
                {
                    keyword = "TRIANGLE_VERTEX";
                    triangle = 0;
                }
                else if (typeMethods.s_eqi(keyword, "TRIANGLE_VERTEX"))
                {
                    i4vec_ = typeMethods.s_to_i4vec(text, 4).ivec;
                    for (i = 0; i < 3; i++)
                    {
                        triangle_vertex[i + triangle * 3] = i4vec_[i];
                    }

                    triangle_label[triangle] = i4vec_[3];
                    triangle += 1;
                }
                else if (typeMethods.s_eqi(keyword, "VERTICES"))
                {
                    keyword = "VERTEX_COORDINATE";
                    vertex = 0;
                }
                else if (typeMethods.s_eqi(keyword, "VERTEX_COORDINATE"))
                {
                    double[] r8vec_ = typeMethods.s_to_r8vec(text, dim + 1).rvec;
                    for (i = 0; i < dim; i++)
                    {
                        vertex_coordinate[i + vertex * dim] = r8vec_[i];
                    }

                    vertex_label[vertex] = (int) r8vec_[dim];
                    vertex += 1;
                }
                else if (typeMethods.s_eqi(keyword, "SKIP"))
                {
                }
                else
                {
                    Console.WriteLine("");
                    Console.WriteLine("MESH_DATA_READ - Fatal error!");
                    Console.WriteLine("  Could not find keyword while reading line "
                                      + line_num + "");
                    Console.WriteLine("\"" + text + "\".");
                    return;
                }
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Read " + line_num + " lines from \"" + filename + "\".");

    }

    public static void mesh_size_print(int dim, int vertices, int edges, int triangles,
            int quadrilaterals, int tetrahedrons, int hexahedrons)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MESH_SIZE_PRINT prints the sizes of an ICE dataset.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Pascal Frey,
        //    MEDIT: An interactive mesh visualization software,
        //    Technical Report RT-0253,
        //    Institut National de Recherche en Informatique et en Automatique,
        //    03 December 2001.
        //
        //  Parameters:
        //
        //    Input, int DIM, the spatial dimension, which should be 2 or 3.
        //
        //    Input, int VERTICES, the number of vertices.
        //
        //    Input, int EDGES, the number of edges (may be 0).
        //
        //    Input, int TRIANGLES, the number of triangles (may be 0).
        //
        //    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).
        //
        //    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).
        //
        //    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
        //
    {
        Console.WriteLine("");
        Console.WriteLine("  Number of dimensions = " + dim + "");
        Console.WriteLine("  Number of vertices = " + vertices + "");
        Console.WriteLine("  Number of edges = " + edges + "");
        Console.WriteLine("  Number of triangles = " + triangles + "");
        Console.WriteLine("  Number of quadrilaterals = " + quadrilaterals + "");
        Console.WriteLine("  Number of tetrahedrons = " + tetrahedrons + "");
        Console.WriteLine("  Number of hexahedrons = " + hexahedrons + "");
    }

    public static void mesh_size_read(string filename, ref int dim, ref int vertices, ref int edges,
            ref int triangles, ref int quadrilaterals, ref int tetrahedrons, ref int hexahedrons)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MESH_SIZE_READ reads sizes from a MESH file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Pascal Frey,
        //    MEDIT: An interactive mesh visualization software,
        //    Technical Report RT-0253,
        //    Institut National de Recherche en Informatique et en Automatique,
        //    03 December 2001.
        //
        //  Parameters:
        //
        //    Input, string FILENAME, the name of the MESH file.
        //
        //    Output, int *DIM, the spatial dimension, which should be 2 or 3.
        //
        //    Output, int *VERTICES, the number of vertices.
        //
        //    Output, int *EDGES, the number of edges (may be 0).
        //
        //    Output, int *TRIANGLES, the number of triangles (may be 0).
        //
        //    Output, int *QUADRILATERALS, the number of quadrilaterals
        //    (may be 0).
        //
        //    Output, int *TETRAHEDRONS, the number of tetrahedrons
        //    (may be 0).
        //
        //    Output, int *HEXAHEDRONS, the number of hexahedrons
        //    (may be 0).
        //
    {
        string[] input;
        //
        //  Initialize everything to nothing.
        //
        dim = 0;
        vertices = 0;
        edges = 0;
        triangles = 0;
        quadrilaterals = 0;
        tetrahedrons = 0;
        hexahedrons = 0;

        try
        {
            input = File.ReadAllLines(filename);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("MESH_SIZE_READ - Fatal error!");
            Console.WriteLine("  Could not open file.");
            return;
        }

        //
        //  Read lines til you get alphanumerics and determine a "mode"
        //
        int line_num = 0;
        string keyword = "NONE";

        for (;;)
        {
            string text;
            try
            {
                text = input[line_num];
            }
            catch (Exception)
            {
                break;
            }

            line_num += 1;

            if (typeMethods.s_len_trim(text) == 0)
            {
                keyword = "NONE";
                continue;
            }

            switch (text[0])
            {
                case '#':
                    continue;
            }
            //
            //  Remove initial blanks.
            //

            //
            //  Expecting a keyword.
            //
            if (typeMethods.s_eqi(text, "CORNERS"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "DIMENSION"))
            {
                keyword = "DIMENSION";
            }
            else if (typeMethods.s_eqi(text, "EDGES"))
            {
                keyword = "EDGES";
            }
            else if (typeMethods.s_eqi(text, "END"))
            {
                Console.WriteLine("");
                Console.WriteLine("  END statement encountered.");
                break;
            }
            else if (typeMethods.s_eqi(text, "HEXAHEDRA") ||
                     typeMethods.s_eqi(text, "HEXAHEDRONS"))
            {
                keyword = "HEXAHEDRONS";
            }
            else if (text.StartsWith("MESHVERSIONFORMATTED"))
            {
            }
            else if (typeMethods.s_eqi(text, "NORMALATQUADRILATERALVERTICES"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "NORMALATTRIANGLEVERTICES"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "NORMALATVERTICES"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "NORMALS"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "QUADRILATERALS"))
            {
                keyword = "QUADRILATERALS";
            }
            else if (typeMethods.s_eqi(text, "REQUIREDEDGES"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "REQUIREDVERTICES"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "RIDGES"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "TANGENTATEDGES"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "TANGENTS"))
            {
                keyword = "SKIP";
            }
            else if (typeMethods.s_eqi(text, "TETRAHEDRA") ||
                     typeMethods.s_eqi(text, "TETRAHEDRONS"))
            {
                keyword = "TETRAHEDRONS";
            }
            else if (typeMethods.s_eqi(text, "TRIANGLES"))
            {
                keyword = "TRIANGLES";
            }
            else if (typeMethods.s_eqi(text, "VERTICES"))
            {
                keyword = "VERTICES";
            }
            //
            //  Presumably, numeric data to be processed by keyword.
            //
            else if (typeMethods.s_eqi(keyword, "DIMENSION"))
            {
                dim = Convert.ToInt32(text);
                keyword = "NONE";
            }
            else if (typeMethods.s_eqi(keyword, "EDGES"))
            {
                edges = Convert.ToInt32(text);
                keyword = "EDGE_VERTEX";
            }
            else if (typeMethods.s_eqi(keyword, "EDGE_VERTEX"))
            {
            }
            else if (typeMethods.s_eqi(keyword, "HEXAHEDRONS"))
            {
                hexahedrons = Convert.ToInt32(text);
                keyword = "HEXAHEDRON_VERTEX";
            }
            else if (typeMethods.s_eqi(keyword, "HEXAHEDRON_VERTEX"))
            {
            }
            else if (typeMethods.s_eqi(keyword, "QUADRILATERALS"))
            {
                quadrilaterals = Convert.ToInt32(text);
                keyword = "QUADRILATERAL_VERTEX";
            }
            else if (typeMethods.s_eqi(keyword, "QUADRILATERAL_VERTEX"))
            {
            }
            else if (typeMethods.s_eqi(keyword, "TETRAHEDRONS"))
            {
                tetrahedrons = Convert.ToInt32(text);
                keyword = "TETRAHEDRON_VERTEX";
            }
            else if (typeMethods.s_eqi(keyword, "TETRAHEDRON_VERTEX"))
            {
            }
            else if (typeMethods.s_eqi(keyword, "TRIANGLES"))
            {
                triangles = Convert.ToInt32(text);
                keyword = "TRIANGLE_VERTEX";
            }
            else if (typeMethods.s_eqi(keyword, "TRIANGLE_VERTEX"))
            {
            }
            else if (typeMethods.s_eqi(keyword, "VERTICES"))
            {
                vertices = Convert.ToInt32(text);
                keyword = "VERTEX_COORDINATE";
            }
            else if (typeMethods.s_eqi(keyword, "VERTEX_COORDINATE"))
            {
            }
            else if (typeMethods.s_eqi(keyword, "SKIP"))
            {
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("MESH_SIZE_READ - Fatal error!");
                Console.WriteLine("  Could not find keyword while reading line "
                                  + line_num + "");
                Console.WriteLine("\"" + text + "\".");
                return;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Read " + line_num + " lines from \"" + filename + "\".");

    }

    public static void mesh_write(string filename, int dim, int vertices, int edges,
            int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
            double[] vertex_coordinate, int[] vertex_label, int[] edge_vertex,
            int[] edge_label, int[] triangle_vertex, int[] triangle_label,
            int[] quadrilateral_vertex, int[] quadrilateral_label,
            int[] tetrahedron_vertex, int[] tetrahedron_label,
            int[] hexahedron_vertex, int[] hexahedron_label)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MESH_WRITE writes mesh data to a MESH file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Pascal Frey,
        //    MEDIT: An interactive mesh visualization software,
        //    Technical Report RT-0253,
        //    Institut National de Recherche en Informatique et en Automatique,
        //    03 December 2001.
        //
        //  Parameters:
        //
        //    Input, string FILENAME, the name of the file to be created.
        //    Ordinarily, the name should include the extension ".mesh".
        //
        //    Input, int DIM, the spatial dimension, which should be 2 or 3.
        //
        //    Input, int VERTICES, the number of vertices.
        //
        //    Input, int EDGES, the number of edges (may be 0).
        //
        //    Input, int TRIANGLES, the number of triangles (may be 0).
        //
        //    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).
        //
        //    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).
        //
        //    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
        //
        //    Input, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
        //    of each vertex.
        //
        //    Input, int VERTEX_LABEL[VERTICES], a label for each vertex.
        //
        //    Input, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.
        //
        //    Input, int EDGE_LABEL[EDGES], a label for each edge.
        //
        //    Input, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
        //    each triangle.
        //
        //    Input, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.
        //
        //    Input, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
        //    form each quadrilateral.
        //
        //    Input, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
        //    each quadrilateral.
        //
        //    Input, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
        //    form each tetrahedron.
        //
        //    Input, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
        //    each tetrahedron.
        //
        //    Input, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
        //    each hexahedron.
        //
        //    Input, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
        //
    {
        int i;
        int j;
        List<string> output = new()
        {
            "MeshVersionFormatted 1",
            "#  Created by mesh_write.C",
            //
            //  Vertices.
            //
            "",
            "Vertices",
            vertices + ""
        };

        for (j = 0; j < vertices; j++)
        {
            string cout = "";

            for (i = 0; i < dim; i++)
            {
                cout += "  " + vertex_coordinate[i + j * dim];
            }

            output.Add(cout + "  " + vertex_label[j] + "");
        }

        switch (edges)
        {
            //
            //  Edges.
            //
            case > 0:
            {
                output.Add("");
                output.Add("Edges");
                output.Add(edges + "");
                for (j = 0; j < edges; j++)
                {
                    string cout = "";

                    for (i = 0; i < 2; i++)
                    {
                        cout += "  " + edge_vertex[i + j * 2];
                    }

                    output.Add(cout + "  " + edge_label[j] + "");
                }

                break;
            }
        }

        switch (triangles)
        {
            //
            //  Triangles.
            //
            case > 0:
            {
                output.Add("");
                output.Add("Triangles");
                output.Add(triangles + "");
                for (j = 0; j < triangles; j++)
                {
                    string cout = "";

                    for (i = 0; i < 3; i++)
                    {
                        cout += "  " + triangle_vertex[i + j * 3];
                    }

                    output.Add(cout + "  " + triangle_label[j] + "");
                }

                break;
            }
        }

        switch (quadrilaterals)
        {
            //
            //  Quadrilaterals.
            //
            case > 0:
            {
                output.Add("");
                output.Add("Quadrilaterals");
                output.Add(quadrilaterals + "");
                for (j = 0; j < quadrilaterals; j++)
                {
                    string cout = "";
                    for (i = 0; i < 4; i++)
                    {
                        cout += "  " + quadrilateral_vertex[i + j * 4];
                    }

                    output.Add(cout + "  " + quadrilateral_label[j] + "");
                }

                break;
            }
        }

        switch (tetrahedrons)
        {
            //
            //  Tetrahedra.
            //
            case > 0:
            {
                output.Add("");
                output.Add("Tetrahedra");
                output.Add(tetrahedrons + "");
                for (j = 0; j < tetrahedrons; j++)
                {
                    string cout = "";
                    for (i = 0; i < 4; i++)
                    {
                        cout += "  " + tetrahedron_vertex[i + j * 4];
                    }

                    output.Add(cout + "  " + tetrahedron_label[j] + "");
                }

                break;
            }
        }

        switch (hexahedrons)
        {
            //
            //  Hexahedra.
            //
            case > 0:
            {
                output.Add("");
                output.Add("Hexahedra");
                output.Add(hexahedrons + "");
                for (j = 0; j < hexahedrons; j++)
                {
                    string cout = "";
                    for (i = 0; i < 8; i++)
                    {
                        cout += "  " + hexahedron_vertex[i + j * 8];
                    }

                    output.Add(cout + "  " + hexahedron_label[j] + "");
                }

                break;
            }
        }

        //
        //  End
        //
        output.Add("");
        output.Add("End");

        try
        {
            File.WriteAllLines(filename, output);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("MESH_WRITE - Fatal error!");
            Console.WriteLine("  Unable to open output file.");
        }

    }
}