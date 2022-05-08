using geoCoreLib;
using geoLib;
using MiscUtil.Conversion;
using MiscUtil.IO;
using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text;

namespace oasis;

internal static class Utils
{
    public static string Left(this string str, int length)
    {
        return str.Substring(0, Math.Min(length, str.Length));
    }
}

internal partial class oasReader
{
    public delegate void StatusUpdateUI(string text);
    public StatusUpdateUI statusUpdateUI { get; set; }
    public delegate void ProgressUpdateUI(double progress);
    public ProgressUpdateUI progressUpdateUI { get; set; }

    public List<string> error_msgs;
    public bool valid { get; set; }

    private GCDrawingfield drawing_;

    private enum elementType { boxElement, polygonElement, pathElement, cellrefElement, textElement, circleElement, trapezoidElement, ctrapezoidElement }
    public struct modals
    {
        public bool absoluteMode { get; set; }
        public int placement_x { get; set; }
        public int placement_y { get; set; }
        public int layer { get; set; }
        public int datatype { get; set; }
        public string placement_cell { get; set; }
        public double mag { get; set; }
        public double angle { get; set; }
        public bool mirror_x { get; set; }
        public int textlayer { get; set; }
        public int texttype{ get; set; }
        public int text_x { get; set; }
        public int text_y { get; set; }
        public string text_string { get; set; }
        public int geometry_x { get; set; }
        public int geometry_y { get; set; }
        public int xy_model { get; set; }
        public int geometry_w { get; set; }
        public int geometry_h { get; set; }
        public List<GeoLibPoint> polygon_point_list { get; set; }
        public int path_halfwidth { get; set; }
        public int path_point_list { get; set; }
        public int path_start_extension { get; set; }
        public int path_end_extension { get; set; }
        public int path_start_extension_value { get; set; }
        public int path_end_extension_value { get; set; }
        public int ctrapezoid_type { get; set; }
        public double circle_radius { get; set; }
        public int last_property_name { get; set; }
        public int last_value_list { get; set; }
        //  repetition;
        public int repetition { get; set; }
        public int x_dimension { get; set; }
        public int y_dimension { get; set; }
        public int x_space { get; set; }
        public int y_space { get; set; }
        public List<GeoLibPoint> repArray { get; set; }
        public int trapezoid_delta_a { get; set; }
        public int trapezoid_delta_b { get; set; }
        public bool trapezoid_orientation { get; set; }
        
        public string s { get; set; }
    }

    private GCCell cell_;
    private string[] cellNames = new string[1024000]; // maxLayers
    private string[] textNames = new string[1024000]; // maxLayers
    public Dictionary<string, string> layerNames { get; set; }
    private modals modal;

    private bool zLibUsed; // ZLib data in use.

    private byte[] zLibOut; // output buffer of decompressed ZLib
    private uint zlibOutPos; // position in output buffer

    private string filename;
    private EndianBinaryReader br;

    public oasReader(string filename_)
    {
        pOASReader(filename_);
    }

    private void pOASReader(string filename_)
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

    private void pReset()
    {
        drawing_.reset();
        cell_ = null;
        resetModal();
        valid = false;
        error_msgs.Clear();
    }

    private void resetModal()
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
        modal.s = "";
    }

    public bool load(ref GCDrawingfield drawing)
    {
        return pLoad(ref drawing);
    }

    private bool pLoad(ref GCDrawingfield drawing)
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
                using GZipStream gzs = new(stream, CompressionMode.Decompress);
                MemoryStream ms = new();
                gzs.CopyTo(ms);
                ms.Seek(0, SeekOrigin.Begin);
                br = new EndianBinaryReader(EndianBitConverter.Little, ms);
            }
            else
            {
                br = new EndianBinaryReader(EndianBitConverter.Little, stream);
            }

            int i;
            string s1 = "";
            zLibUsed = false;
            for (i = 0; i < 13; i++)
            {
                byte help = readRaw();
                if (help == 0)
                {
                    continue;
                }

                modal.s = Encoding.UTF8.GetString(new [] { help });
                s1 += modal.s;
            }
            if (s1 != "%SEMI-OASIS\r\n")
            {
                const string err = "Invalid Format.";
                error_msgs.Add(err);
                throw new Exception(err);
            }
            int record;
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
                        modal.s = readString();
                        if (modal.s != "1.0")
                        {
                            string err2 = "Unknown/unsupported version of OASIS: " + modal.s;
                            error_msgs.Add(err2);
                            throw new Exception(err2);
                        }
                        drawing_.databaseunits = readReal();
                        i = readUnsignedInteger();
                        switch (i)
                        {
                            case 0:
                            {
                                tableAtEnd = false;
                                for (i = 0; i < 12; i++)
                                {
                                    readUnsignedInteger();
                                }

                                break;
                            }
                            default:
                                tableAtEnd = true;
                                break;
                        }

                        break;
                    case 2: //end
                        switch (tableAtEnd)
                        {
                            case true:
                            {
                                for (i = 0; i < 12; i++)
                                {
                                    readUnsignedInteger();
                                }

                                break;
                            }
                        }
                        break;
                    case 3: //cellname
                        cellNames[cellNameCount] = readString();
                        cellNameCount++;
                        break;
                    case 4: //cellname
                        modal.s = readString();
                        i = readUnsignedInteger();
                        cellNames[i] = modal.s;
                        break;
                    case 5: //textname
                        textNames[textNameCount] = readString();
                        textNameCount++;
                        break;
                    case 6: //textname
                        modal.s = readString();
                        i = readUnsignedInteger();
                        textNames[i] = modal.s;
                        break;
                    case 7: //property
                        modal.s = readString();
                        break;
                    case 8: //property
                        modal.s = readString();
                        i = readUnsignedInteger();
                        break;
                    case 9: //property string
                        modal.s = readString();
                        break;
                    case 10: //property string
                        modal.s = readString();
                        i = readUnsignedInteger();
                        break;
                    case 11: //layername textlayername 
                    case 12:
                        modal.s = readString();
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
                        }
                        layerNames.Add("L" + l + "D" + i, modal.s);
                        break;
                    case 13: // cellrecord
                        cell_ = drawing_.addCell();
                        i = readUnsignedInteger();
                        cellNames[i] = cellNames[i] switch
                        {
                            "" => "layout#cell~" + i,
                            _ => cellNames[i] != null ? cellNames[i] : modal.s
                        };
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
                                cellNames[i] = cellNames[i] switch
                                {
                                    "" => "layout#cell~" + i,
                                    _ => cellNames[i]
                                };
                                modal.placement_cell = cellNames[i];
                            }
                            else
                            {
                                modal.placement_cell = readString();
                            }
                        }
                        if ((info_byte & 1) != 0)
                        {
                            modal.mirror_x = true;
                        }

                        modal.angle = (info_byte & 6) switch
                        {
                            6 => 270,
                            _ => (info_byte & 6) switch
                            {
                                4 => 180,
                                _ => (info_byte & 6) switch
                                {
                                    2 => 90,
                                    _ => modal.angle
                                }
                            }
                        };
                        switch (modal.absoluteMode)
                        {
                            case true:
                            {
                                if ((info_byte & 32) != 0)
                                {
                                    modal.placement_x = readSignedInteger();
                                }
                                if ((info_byte & 16) != 0)
                                {
                                    modal.placement_y = readSignedInteger();
                                }

                                break;
                            }
                            default:
                            {
                                if ((info_byte & 32) != 0)
                                {
                                    modal.placement_x += readSignedInteger();
                                }
                                if ((info_byte & 16) != 0)
                                {
                                    modal.placement_y += readSignedInteger();
                                }

                                break;
                            }
                        }
                        if ((info_byte & 8) != 0)
                        {
                            readRepetition();
                            processRepetition(elementType.cellrefElement);
                        }
                        else
                        {
                            addCellref();
                        }

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
                                cellNames[i] = cellNames[i] switch
                                {
                                    "" => "layout#cell~" + i,
                                    _ => cellNames[i]
                                };
                                modal.placement_cell = cellNames[i];
                            }
                            else
                            {
                                modal.placement_cell = readString();
                            }
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
                        switch (modal.absoluteMode)
                        {
                            case true:
                            {
                                if ((info_byte & 32) != 0)
                                {
                                    modal.placement_x = readSignedInteger();
                                }
                                if ((info_byte & 16) != 0)
                                {
                                    modal.placement_y = readSignedInteger();
                                }

                                break;
                            }
                            default:
                            {
                                if ((info_byte & 32) != 0)
                                {
                                    modal.placement_x += readSignedInteger();
                                }
                                if ((info_byte & 16) != 0)
                                {
                                    modal.placement_y += readSignedInteger();
                                }

                                break;
                            }
                        }
                        if ((info_byte & 8) != 0)
                        {
                            readRepetition();
                            processRepetition(elementType.cellrefElement);
                        }
                        else
                        {
                            addCellref();
                        }

                        break;
                    case 19: //text
                        info_byte = readRaw();
                        if ((info_byte & 64) != 0)
                        {
                            if ((info_byte & 32) != 0)
                            {
                                i = readUnsignedInteger();
                                textNames[i] = textNames[i] switch
                                {
                                    "" => "layout#text~" + i,
                                    _ => textNames[i]
                                };
                                modal.text_string = textNames[i];
                            }
                            else
                            {
                                modal.text_string = readString();
                            }
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
                        switch (modal.absoluteMode)
                        {
                            case true:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.text_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.text_y = readSignedInteger();
                                }

                                break;
                            }
                            default:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.text_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.text_y += readSignedInteger();
                                }

                                break;
                            }
                        }
                        if ((info_byte & 4) != 0)
                        {
                            readRepetition();
                            processRepetition(elementType.textElement);
                        }
                        else
                        {
                            addText();
                        }

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
                        switch (modal.absoluteMode)
                        {
                            case true:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }

                                break;
                            }
                            default:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }

                                break;
                            }
                        }
                        if ((info_byte & 4) != 0)
                        {
                            readRepetition();
                            processRepetition(elementType.boxElement);
                        }
                        else
                        {
                            addBox();
                        }

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
                        switch (modal.absoluteMode)
                        {
                            case true:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }

                                break;
                            }
                            default:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }

                                break;
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
                            switch (modal.geometry_w)
                            {
                                // Avoid zero width path.
                                case 0:
                                    modal.geometry_w = 10;
                                    break;
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
                        switch (modal.absoluteMode)
                        {
                            case true:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }

                                break;
                            }
                            default:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }

                                break;
                            }
                        }
                        if ((info_byte & 4) != 0)
                        {
                            readRepetition();
                            processRepetition(elementType.pathElement);
                        }
                        else
                        {
                            addPath();
                        }

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
                        switch (modal.absoluteMode)
                        {
                            case true:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }

                                break;
                            }
                            default:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }

                                break;
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
                        switch (record)
                        {
                            case 24:
                                modal.trapezoid_delta_a = read1Delta(false).X;
                                break;
                            default:
                                modal.trapezoid_delta_b = read1Delta(false).X;
                                break;
                        }
                        switch (modal.absoluteMode)
                        {
                            case true:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }

                                break;
                            }
                            default:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }

                                break;
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
                        switch (modal.absoluteMode)
                        {
                            case true:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }

                                break;
                            }
                            default:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }

                                break;
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
                        switch (modal.absoluteMode)
                        {
                            case true:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }

                                break;
                            }
                            default:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }

                                break;
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
                        if ((info_byte & 4) != 0 && (info_byte & 2) != 0)
                        {
                            readUnsignedInteger();
                        }
                        if ((info_byte & 4) != 0 && (info_byte & 2) == 0)
                        {
                            readString();
                        }
                        switch (info_byte & 8)
                        {
                            case 0:
                            {
                                int count = info_byte >> 4;
                                count = count switch
                                {
                                    15 => readUnsignedInteger(),
                                    _ => count
                                };
                                for (int j = 0; j < count; j++)
                                {
                                    readProperty();
                                }

                                break;
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
                        switch (modal.absoluteMode)
                        {
                            case true:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x = readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y = readSignedInteger();
                                }

                                break;
                            }
                            default:
                            {
                                if ((info_byte & 16) != 0)
                                {
                                    modal.geometry_x += readSignedInteger();
                                }
                                if ((info_byte & 8) != 0)
                                {
                                    modal.geometry_y += readSignedInteger();
                                }

                                break;
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
            foreach (GCCell t in drawing_.cellList.Where(t => t != null))
            {
                s1 = t.cellName;
                if (s1.Left(12) != "layout#cell~")
                {
                    continue;
                }

                s1 = s1.Substring(12, s1.Length - 12);
                t.cellName = cellNames[Convert.ToInt32(s1)];
            }
            foreach (GCElement t1 in drawing_.cellList.Where(t => t != null).SelectMany(t => t.elementList))
            {
                if (t1.isCellref() || t1.isCellrefArray())
                {
                    if (t1.depend() == null)
                    {
                        s1 = t1.getName();
                        if (s1 != null && s1.Left(12) == "layout#cell~")
                        {
                            s1 = s1.Substring(12, s1.Length - 12);
                            t1.setName(cellNames[Convert.ToInt32(s1)]);
                            t1.setCellRef(drawing_.findCell(cellNames[Convert.ToInt32(s1)]));
                        }
                    }
                }

                if (!t1.isText())
                {
                    continue;
                }

                s1 = t1.getName();
                if (s1.Left(12) != "layout#text~")
                {
                    continue;
                }

                s1 = s1.Substring(12, s1.Length - 12);
                t1.setName(textNames[Convert.ToInt32(s1)]);
            }

            try
            {
                drawing_.active_cell = drawing.findCellIndex(cell_.cellName);
                // drawing_.resize(1000.0 / drawing_.databaseunits);
            }
            catch (Exception)
            {
                const string err = "Unable to find any cells. This library only supports Oasis saved in strict mode.";
                error_msgs.Add(err);
                throw new Exception(err);
            }

            statusUpdateUI?.Invoke("Done");
            progressUpdateUI?.Invoke(1.0f);
        }
        catch (Exception e)
        {
            valid = false;
            switch (error_msgs.IndexOf(e.Message))
            {
                case -1:
                    error_msgs.Add(e.Message);
                    break;
            }
        }
        return valid;
    }
}