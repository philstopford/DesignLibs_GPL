using geoCoreLib;
using geoLib;
using MiscUtil.Conversion;
using MiscUtil.IO;
using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Text;

namespace oasis
{
    static class Utils
    {
        public static string Left(this string str, int length)
        {
            return str.Substring(0, Math.Min(length, str.Length));
        }
    }

    partial class oasReader
    {
        public delegate void StatusUpdateUI(string text);
        public StatusUpdateUI statusUpdateUI { get; set; }
        public delegate void ProgressUpdateUI(double progress);
        public ProgressUpdateUI progressUpdateUI { get; set; }

        public List<string> error_msgs;
        public bool valid { get; set; }

        GCDrawingfield drawing_;

        enum elementType { boxElement, polygonElement, pathElement, cellrefElement, textElement, circleElement, trapezoidElement, ctrapezoidElement };
        public struct modals
        {
            public bool absoluteMode { get; set; }
            public Int32 placement_x { get; set; }
            public Int32 placement_y { get; set; }
            public Int32 layer { get; set; }
            public Int32 datatype { get; set; }
            public string placement_cell { get; set; }
            public double mag { get; set; }
            public double angle { get; set; }
            public bool mirror_x { get; set; }
            public Int32 textlayer { get; set; }
            public Int32 texttype{ get; set; }
            public Int32 text_x { get; set; }
            public Int32 text_y { get; set; }
            public string text_string { get; set; }
            public Int32 geometry_x { get; set; }
            public Int32 geometry_y { get; set; }
            public Int32 xy_model { get; set; }
            public Int32 geometry_w { get; set; }
            public Int32 geometry_h { get; set; }
            public List<GeoLibPoint> polygon_point_list { get; set; }
            public Int32 path_halfwidth { get; set; }
            public Int32 path_point_list { get; set; }
            public Int32 path_start_extension { get; set; }
            public Int32 path_end_extension { get; set; }
            public Int32 path_start_extension_value { get; set; }
            public Int32 path_end_extension_value { get; set; }
            public Int32 ctrapezoid_type { get; set; }
            public double circle_radius { get; set; }
            public Int32 last_property_name { get; set; }
            public Int32 last_value_list { get; set; }
            //  repetition;
            public Int32 repetition { get; set; }
            public Int32 x_dimension { get; set; }
            public Int32 y_dimension { get; set; }
            public Int32 x_space { get; set; }
            public Int32 y_space { get; set; }
            public List<GeoLibPoint> repArray { get; set; }
            public Int32 trapezoid_delta_a { get; set; }
            public Int32 trapezoid_delta_b { get; set; }
            public bool trapezoid_orientation { get; set; }
        }

        GCCell cell_;
        string[] cellNames = new string[1024000]; // maxLayers
        string[] textNames = new string[1024000]; // maxLayers
        public Dictionary<string, string> layerNames { get; set; }
        modals modal;

        bool zLibUsed; // ZLib data in use.

        byte[] zLibOut; // output buffer of decompressed ZLib
        uint zlibOutPos; // position in output buffer

        string filename;
        EndianBinaryReader br;

        public oasReader(String filename_)
        {
            pOASReader(filename_);
        }

        void pOASReader(String filename_)
        {
            drawing_ = new GCDrawingfield(filename_);
            error_msgs = new List<string>();
            modal.repArray = new List<GeoLibPoint>();
            modal.polygon_point_list = new List<GeoLibPoint>();
            filename = filename_;
        }

        public void reset()
        {
            pReset();
        }

        void pReset()
        {
            drawing_.reset();
            cell_ = null;
            resetModal();
            valid = false;
            error_msgs.Clear();
        }

        void resetModal()
        {
            modal.placement_x = 0;
            modal.placement_y = 0;
            modal.geometry_x = 0;
            modal.geometry_y = 0;
            modal.text_x = 0;
            modal.text_y = 0;
            modal.absoluteMode = true;
            // undefined
            modal.layer = -1;
            modal.datatype = -1;
            modal.text_string = "";
            modal.textlayer = -1;
            modal.texttype = -1;
            modal.circle_radius = -1;
            modal.repetition = -1;
            modal.polygon_point_list = new List<GeoLibPoint>();
            modal.repArray = new List<GeoLibPoint>();
        }

        public bool load(ref GCDrawingfield drawing)
        {
            return pLoad(ref drawing);
        }

        bool pLoad(ref GCDrawingfield drawing)
        {
            drawing_ = drawing;
            valid = true;
            error_msgs.Clear();

            layerNames = new Dictionary<string, string>();
            try
            {
                statusUpdateUI?.Invoke("Loading");
                progressUpdateUI?.Invoke(0);

                Stream stream = File.OpenRead(filename);
                if (filename.ToLower().EndsWith("gz"))
                {
                    using (GZipStream gzs = new GZipStream(stream, CompressionMode.Decompress))
                    {
                        MemoryStream ms = new MemoryStream();
                        gzs.CopyTo(ms);
                        ms.Seek(0, SeekOrigin.Begin);
                        br = new EndianBinaryReader(EndianBitConverter.Little, ms);
                    }
                }
                else
                {
                    br = new EndianBinaryReader(EndianBitConverter.Little, stream);
                }

                string s;
                Int32 i;
                string s1 = "";
                zLibUsed = false;
                for (i = 0; i < 13; i++)
                {
                    byte help = readRaw();
                    if (help != 0)
                    {
                        s = Encoding.UTF8.GetString(new [] { help });
                        s1 += s;
                    }
                }
                if (s1 != "%SEMI-OASIS\r\n")
                {
                    string err = "Invalid Format.";
                    error_msgs.Add(err);
                    throw new Exception(err);
                }
                Int32 record;
                bool tableAtEnd = false;
                cell_ = null;
                int cellNameCount = 0;
                int textNameCount = 0;
                do
                {
                    if (cellNameCount > cellNames.Length)
                    {
                        string err = "More cells (" + cellNameCount + ") than are able to be supported (" + cellNames.Length + ")!";
                        error_msgs.Add(err);
                        throw new Exception(err);
                    }
                    if (textNameCount > textNames.Length)
                    {
                        string err = "More text names (" + textNameCount + ") than are able to be supported (" + textNames.Length + ")!";
                        error_msgs.Add(err);
                        throw new Exception(err);
                    }
                    record = readUnsignedInteger();
                    byte info_byte;
                    switch (record)
                    {
                        case 0: //pad
                            break;
                        case 1: //start
                            s = readString();
                            if (s != "1.0")
                            {
                                string err2 = "Unknown/unsupported version of OASIS: " + s;
                                error_msgs.Add(err2);
                                throw new Exception(err2);
                            }
                            drawing_.databaseunits = readReal();
                            i = readUnsignedInteger();
                            if (i == 0)
                            {
                                tableAtEnd = false;
                                for (i = 0; i < 12; i++)
                                {
                                    readUnsignedInteger();
                                }
                            }
                            else tableAtEnd = true;
                            break;
                        case 2: //end
                            if (tableAtEnd)
                            {
                                for (i = 0; i < 12; i++)
                                {
                                    readUnsignedInteger();
                                }
                            }
                            break;
                        case 3: //cellname
                            cellNames[cellNameCount] = readString();
                            cellNameCount++;
                            break;
                        case 4: //cellname
                            s = readString();
                            i = readUnsignedInteger();
                            cellNames[i] = s;
                            break;
                        case 5: //textname
                            textNames[textNameCount] = readString();
                            textNameCount++;
                            break;
                        case 6: //textname
                            s = readString();
                            i = readUnsignedInteger();
                            textNames[i] = s;
                            break;
                        case 7: //property
                            s = readString();
                            break;
                        case 8: //property
                            s = readString();
                            i = readUnsignedInteger();
                            break;
                        case 9: //property string
                            s = readString();
                            break;
                        case 10: //property string
                            s = readString();
                            i = readUnsignedInteger();
                            break;
                        case 11: //layername textlayername 
                        case 12:
                            s = readString();
                            i = readUnsignedInteger();
                            switch (i)
                            {
                                case 0:
                                    break;
                                case 1:
                                    i = readUnsignedInteger();
                                    break;
                                case 3:
                                    i = readUnsignedInteger();
                                    break;
                                case 2:
                                    i = readUnsignedInteger();
                                    break;
                                case 4:
                                    i = readUnsignedInteger();
                                    readUnsignedInteger();
                                    break;
                                default:
                                    //"Error in layername/textlayername"
                                    break;
                            }
                            int l = i;
                            //datatype
                            i = readUnsignedInteger();
                            switch (i)
                            {
                                case 0:
                                    break;
                                case 1:
                                    i = readUnsignedInteger();
                                    break;
                                case 2:
                                    i = readUnsignedInteger();
                                    break;
                                case 3:
                                    i = readUnsignedInteger();
                                    break;
                                case 4:
                                    i = readUnsignedInteger();
                                    readUnsignedInteger();
                                    break;
                                default:
                                    // "Error in layername/textlayername"
                                    break;
                            }
                            try
                            {
                                layerNames.Add("L" + l + "D" + i, s);
                            }
                            catch (Exception)
                            {
                            }
                            break;
                        case 13: // cellrecord
                            cell_ = drawing_.addCell();
                            i = readUnsignedInteger();
                            if (cellNames[i] == "")
                            {
                                cellNames[i] = "layout#cell~" + i;
                            }
                            cell_.cellName = cellNames[i];
                            resetModal();
                            break;
                        case 14: // cellrecord
                            cell_ = drawing_.addCell();
                            cell_.cellName = readString();
                            resetModal();
                            break;
                        case 15: //xyabsolute
                            modal.absoluteMode = true;
                            break;
                        case 16: //xyrelative
                            modal.absoluteMode = false;
                            break;
                        case 17: //cellref
                            modal.mag = 1;
                            modal.angle = 0;
                            modal.mirror_x = false;
                            info_byte = readRaw();
                            if ((info_byte & 128) != 0)
                            {
                                if ((info_byte & 64) != 0)
                                {
                                    i = readUnsignedInteger();
                                    if (cellNames[i] == "")
                                    {
                                        cellNames[i] = "layout#cell~" + i;
                                    }
                                    modal.placement_cell = cellNames[i];
                                }
                                else modal.placement_cell = readString();
                            }
                            if ((info_byte & 1) != 0)
                            {
                                modal.mirror_x = true;
                            }
                            if ((info_byte & 6) == 2)
                            {
                                modal.angle = 90;
                            }
                            if ((info_byte & 6) == 4)
                            {
                                modal.angle = 180;
                            }
                            if ((info_byte & 6) == 6)
                            {
                                modal.angle = 270;
                            }
                            if (modal.absoluteMode)
                            {
                                if ((info_byte & 32) != 0)
                                {
                                    modal.placement_x = readSignedInteger();
                                }
                                if ((info_byte & 16) != 0)
                                {
                                    modal.placement_y = readSignedInteger();
                                }
                            }
                            else
                            {
                                if ((info_byte & 32) != 0)
                                {
                                    modal.placement_x += readSignedInteger();
                                }
                                if ((info_byte & 16) != 0)
                                {
                                    modal.placement_y += readSignedInteger();
                                }
                            }
                            if ((info_byte & 8) != 0)
                            {
                                readRepetition();
                                processRepetition(elementType.cellrefElement);
                            }
                            else addCellref();
                            break;
                        case 18: //cellref
                            modal.mag = 1;
                            modal.angle = 0;
                            modal.mirror_x = false;
                            info_byte = readRaw();
                            if ((info_byte & 128) != 0)
                            {
                                if ((info_byte & 64) != 0)
                                {
                                    i = readUnsignedInteger();
                                    if (cellNames[i] == "")
                                    {
                                        cellNames[i] = "layout#cell~" + i;
                                    }
                                    modal.placement_cell = cellNames[i];
                                }
                                else modal.placement_cell = readString();
                            }
                            if ((info_byte & 1) != 0)
                            {
                                modal.mirror_x = true;
                            }
                            if ((info_byte & 4) != 0)
                            {
                                modal.mag = readReal();
                            }
                            if ((info_byte & 2) != 0)
                            {
                                modal.angle = readReal();
                            }
                            if (modal.absoluteMode)
                            {
                                if ((info_byte & 32) != 0)
                                {
                                    modal.placement_x = readSignedInteger();
                                }
                                if ((info_byte & 16) != 0)
                                {
                                    modal.placement_y = readSignedInteger();
                                }
                            }
                            else
                            {
                                if ((info_byte & 32) != 0)
                                {
                                    modal.placement_x += readSignedInteger();
                                }
                                if ((info_byte & 16) != 0)
                                {
                                    modal.placement_y += readSignedInteger();
                                }
                            }
                            if ((info_byte & 8) != 0)
                            {
                                readRepetition();
                                processRepetition(elementType.cellrefElement);
                            }
                            else addCellref();
                            break;
                        case 19: //text
                            info_byte = readRaw();
                            if ((info_byte & 64) != 0)
                            {
                                if ((info_byte & 32) != 0)
                                {
                                    i = readUnsignedInteger();
                                    if (textNames[i] == "")
                                    {
                                        textNames[i] = "layout#text~" + i;
                                    }
                                    modal.text_string = textNames[i];
                                }
                                else modal.text_string = readString();
                            }
                            if ((info_byte & 1) != 0)
                            {
                                modal.textlayer = readUnsignedInteger();
                            }
                            if ((info_byte & 2) != 0)
                            {
                                modal.texttype = readUnsignedInteger();
                                modal.datatype = modal.texttype; // we don't treat text differently.
                            }
                            if (modal.absoluteMode)
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.text_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.text_y = readSignedInteger();
                                }
                            }
                            else
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.text_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.text_y += readSignedInteger();
                                }
                            }
                            if ((info_byte & 4) != 0)
                            {
                                readRepetition();
                                processRepetition(elementType.textElement);
                            }
                            else addText();
                            break;
                        case 20: //rectangle/box
                            info_byte = readRaw();
                            if ((info_byte & 1) != 0)
                            {
                                modal.layer = readUnsignedInteger();
                            }
                            if ((info_byte & 2) != 0)
                            {
                                modal.datatype = readUnsignedInteger();
                            }
                            if ((info_byte & 64) != 0)
                            {
                                modal.geometry_w = readUnsignedInteger();
                            }
                            if ((info_byte & 32) != 0)
                            {
                                modal.geometry_h = readUnsignedInteger();
                            }
                            if ((info_byte & 128) != 0)
                            {
                                modal.geometry_h = modal.geometry_w;
                            }
                            if (modal.absoluteMode)
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }
                            }
                            else
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }
                            }
                            if ((info_byte & 4) != 0)
                            {
                                readRepetition();
                                processRepetition(elementType.boxElement);
                            }
                            else addBox();
                            break;
                        case 21: //polygon
                            info_byte = readRaw();
                            if ((info_byte & 1) != 0)
                            {
                                modal.layer = readUnsignedInteger();
                            }
                            if ((info_byte & 2) != 0)
                            {
                                modal.datatype = readUnsignedInteger();
                            }
                            if ((info_byte & 32) != 0)
                            {
                                readPointList(true);
                            }
                            if (modal.absoluteMode)
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }
                            }
                            else
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }
                            }
                            if ((info_byte & 4) != 0)
                            {
                                readRepetition();
                                processRepetition(elementType.polygonElement);
                            }
                            else
                            {
                                addPolygon();
                            }
                            break;
                        case 22: //path
                            info_byte = readRaw();
                            if ((info_byte & 1) != 0)
                            {
                                modal.layer = readUnsignedInteger();
                            }
                            if ((info_byte & 2) != 0)
                            {
                                modal.datatype = readUnsignedInteger();
                            }
                            if ((info_byte & 64) != 0)
                            {
                                modal.geometry_w = readUnsignedInteger();
                                // Avoid zero width path.
                                if (modal.geometry_w == 0)
                                {
                                    modal.geometry_w = 10;
                                }
                            }
                            if ((info_byte & 128) != 0)
                            {
                                readExtension();
                            }
                            if ((info_byte & 32) != 0)
                            {
                                readPointList(false);
                            }
                            if (modal.absoluteMode)
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }
                            }
                            else
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }
                            }
                            if ((info_byte & 4) != 0)
                            {
                                readRepetition();
                                processRepetition(elementType.pathElement);
                            }
                            else addPath();
                            break;
                        case 23: //trapezoid
                            modal.trapezoid_delta_a = 0;
                            modal.trapezoid_delta_b = 0;
                            modal.trapezoid_orientation = false;
                            info_byte = readRaw();
                            if ((info_byte & 1) != 0)
                            {
                                modal.layer = readUnsignedInteger();
                            }
                            if ((info_byte & 2) != 0)
                            {
                                modal.datatype = readUnsignedInteger();
                            }
                            if ((info_byte & 64) != 0)
                            {
                                modal.geometry_w = readUnsignedInteger();
                            }
                            if ((info_byte & 32) != 0)
                            {
                                modal.geometry_h = readUnsignedInteger();
                            }
                            if ((info_byte & 128) != 0)
                            {
                                modal.trapezoid_orientation = true;
                            }
                            modal.trapezoid_delta_a = read1Delta(false).X;
                            modal.trapezoid_delta_b = read1Delta(false).X;
                            if (modal.absoluteMode)
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }
                            }
                            else
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }
                            }
                            if ((info_byte & 4) != 0)
                            {
                                readRepetition();
                                processRepetition(elementType.trapezoidElement);
                            }
                            else
                            {
                                addTrapezoid();
                            }
                            break;
                        case 24:
                        case 25:
                            // trapezoid
                            info_byte = readRaw();
                            modal.trapezoid_delta_a = 0;
                            modal.trapezoid_delta_b = 0;
                            modal.trapezoid_orientation = false;
                            if ((info_byte & 1) != 0)
                            {
                                modal.layer = readUnsignedInteger();
                            }
                            if ((info_byte & 2) != 0)
                            {
                                modal.datatype = readUnsignedInteger();
                            }
                            if ((info_byte & 64) != 0)
                            {
                                modal.geometry_w = readUnsignedInteger();
                            }
                            if ((info_byte & 32) != 0)
                            {
                                modal.geometry_h = readUnsignedInteger();
                            }
                            if ((info_byte & 128) != 0)
                            {
                                modal.trapezoid_orientation = true;
                            }
                            if (record == 24)
                            {
                                modal.trapezoid_delta_a = read1Delta(false).X;
                            }
                            else
                            {
                                modal.trapezoid_delta_b = read1Delta(false).X;
                            }
                            if (modal.absoluteMode)
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }
                            }
                            else
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }
                            }
                            if ((info_byte & 4) != 0)
                            {
                                readRepetition();
                                processRepetition(elementType.trapezoidElement);
                            }
                            else
                            {
                                addTrapezoid();
                            }
                            break;
                        case 26:
                            // ctrapezoid
                            info_byte = readRaw();
                            if ((info_byte & 1) != 0)
                            {
                                modal.layer = readUnsignedInteger();
                            }
                            if ((info_byte & 2) != 0)
                            {
                                modal.datatype = readUnsignedInteger();
                            }
                            if ((info_byte & 128) != 0)
                            {
                                modal.ctrapezoid_type = readUnsignedInteger();
                            }
                            if ((info_byte & 64) != 0)
                            {
                                modal.geometry_w = readUnsignedInteger();
                            }
                            if ((info_byte & 32) != 0)
                            {
                                modal.geometry_h = readUnsignedInteger();
                            }
                            if (modal.absoluteMode)
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }
                            }
                            else
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }
                            }
                            if ((info_byte & 4) != 0)
                            {
                                readRepetition();
                                processRepetition(elementType.ctrapezoidElement);
                            }
                            else
                            {
                                addCtrapezoid();
                            }
                            break;
                        case 27:
                            // circle
                            info_byte = readRaw();
                            if ((info_byte & 1) != 0)
                            {
                                modal.layer = readUnsignedInteger();
                            }
                            if ((info_byte & 2) != 0)
                            {
                                modal.datatype = readUnsignedInteger();
                            }
                            if ((info_byte & 32) != 0)
                            {
                                modal.circle_radius = readUnsignedInteger(); // unsure : need to check if scaling needed here.
                            }
                            if (modal.absoluteMode)
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }
                            }
                            else
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }
                            }
                            if ((info_byte & 4) != 0)
                            {
                                readRepetition();
                                processRepetition(elementType.circleElement);
                            }
                            else
                            {
                                addCircle();
                            }
                            break;
                        case 28: //property
                        case 29:
                            info_byte = readRaw();
                            if (((info_byte & 4) != 0) && ((info_byte & 2) != 0))
                            {
                                readUnsignedInteger();
                            }
                            if (((info_byte & 4) != 0) && ((info_byte & 2) == 0))
                            {
                                readString();
                            }
                            if ((info_byte & 8) == 0)
                            {
                                int count = info_byte >> 4;
                                if (count == 15)
                                {
                                    count = readUnsignedInteger();
                                }
                                for (int j = 0; j < count; j++)
                                {
                                    readProperty();
                                }
                            }
                            break;
                        case 30: //Xname record
                            readUnsignedInteger();
                            readString();
                            readUnsignedInteger();
                            break;
                        case 31: //Xname record
                            readUnsignedInteger();
                            readString();
                            break;
                        case 32: //Xelement record
                            readUnsignedInteger();
                            readString();
                            break;
                        case 33: //xgeometry record
                            info_byte = readRaw();
                            readUnsignedInteger();
                            if ((info_byte & 1) != 0)
                            {
                                modal.layer = readUnsignedInteger();
                            }
                            if ((info_byte & 2) != 0)
                            {
                                modal.datatype = readUnsignedInteger();
                            }
                            readString();
                            if ((info_byte & 32) != 0)
                            {
                                modal.circle_radius = readUnsignedInteger(); // unsure : need to check if scaling needed here.
                            }
                            if (modal.absoluteMode)
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }
                            }
                            else
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }
                            }
                            if ((info_byte & 4) != 0)
                            {
                                readRepetition();
                            }
                            break;
                        case 34: //compression
                                 // throw new Exception("Compression not supported at this time");
                            readUnsignedInteger();
                            uint before = (uint)readUnsignedInteger();
                            uint after = (uint)readUnsignedInteger();
                            // Arguments are switched around here. Not sure whether this is correct.
                            zLibInit(after, before);
                            break;
                        default:
                            string err = "Unknown/unsupported Record." + record;
                            error_msgs.Add(err);
                            throw new Exception(err);
                    }
                }

                while (record != 2);

                // update cellref/text, if table at end
                foreach (var t in drawing_.cellList)
                {
                    if (t != null)
                    {
                        s1 = t.cellName;
                        if (s1.Left(12) == "layout#cell~")
                        {
                            s1 = s1.Substring(12, s1.Length - 12);
                            t.cellName = cellNames[Convert.ToInt32(s1)];
                        }
                    }
                }
                foreach (var t in drawing_.cellList)
                {
                    if (t != null)
                    {
                        foreach (var t1 in t.elementList)
                        {
                            if (t1.isCellref() || t1.isCellrefArray())
                            {
                                if (t1.depend() == null)
                                {
                                    s1 = t1.getName();
                                    if ((s1 != null) && (s1.Left(12) == "layout#cell~"))
                                    {
                                        s1 = s1.Substring(12, s1.Length - 12);
                                        t1.setName(cellNames[Convert.ToInt32(s1)]);
                                        t1.setCellRef(drawing_.findCell(cellNames[Convert.ToInt32(s1)]));
                                    }
                                }
                            }
                            if (t1.isText())
                            {
                                s1 = t1.getName();
                                if (s1.Left(12) == "layout#text~")
                                {
                                    s1 = s1.Substring(12, s1.Length - 12);
                                    t1.setName(textNames[Convert.ToInt32(s1)]);
                                }
                            }
                        }
                    }
                }

                try
                {
                    drawing_.active_cell = drawing.findCellIndex(cell_.cellName);
                    drawing_.resize(1000.0 / drawing_.databaseunits);
                }
                catch (Exception)
                {
                    string err = "Unable to find any cells. This library only supports Oasis saved in strict mode.";
                    error_msgs.Add(err);
                    throw new Exception(err);
                }

                statusUpdateUI?.Invoke("Done");
                progressUpdateUI?.Invoke(1.0f);
            }
            catch (Exception e)
            {
                valid = false;
                if (error_msgs.IndexOf(e.Message) == -1)
                {
                    error_msgs.Add(e.Message);
                }
            }
            return valid;
        }
    }
}
