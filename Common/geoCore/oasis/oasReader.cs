using geoCoreLib;
using MiscUtil.Conversion;
using MiscUtil.IO;
using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Clipper2Lib;

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
    public struct Modals
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
        public Path64 polygon_point_list { get; set; }
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
        public Repetition repetition { get; set; }
        public int trapezoid_delta_a { get; set; }
        public int trapezoid_delta_b { get; set; }
        public bool trapezoid_orientation { get; set; }
        
        public string s { get; set; }
    }

    private GCCell cell_;

    // This is used to allow for implicit (late) definitions
    // to be applied to the layout data after everything is loaded.
    class CellData
    {
        public string cellName;
        public int textIndex = 0; // this tracks the insertion point for text in textNames.
        public int pp_textIndex = 0; // used in the implicit process to track how many text elements have been processed in the cell.
        public string[] textNames = new string[1024]; // arbitrary limit
    }

    private CellData[] cellData = new CellData[1024000];
//    private string[] cellNames = new string[1024000]; // maxLayers
//    private string[] textNames = new string[1024000]; // maxLayers
    public Dictionary<string, string> layerNames { get; set; }
    private Modals modal;

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
        modal.repetition = new ();
        modal.polygon_point_list = new ();
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
        modal.repetition = new();
        modal.polygon_point_list = new ();
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

            for (int cd = 0; cd < cellData.Length; cd++)
            {
                cellData[cd] = new();
            }
            int record;
            bool tableAtEnd = false;
            cell_ = null;
            int cellIndex = 0;
            do
            {
                if (cellIndex > cellData.Length)
                {
                    string err = "More cells (" + cellIndex + ") than are able to be supported (" + cellData.Length + ")!";
                    error_msgs.Add(err);
                    throw new Exception(err);
                }
                if (cellData[cellIndex].textIndex > cellData[cellIndex].textNames.Length)
                {
                    string err = "More text names (" + cellData[cellIndex].textIndex + ") than are able to be supported (" + cellData[cellIndex].textNames.Length + ")!";
                    error_msgs.Add(err);
                    throw new Exception(err);
                }
                record = readUnsignedInteger();
                byte info_byte;
                switch (record)
                {
                    case oasValues.PAD:
                        break;
                    case oasValues.START:
                        string test = readString();
                        if (test != "1.0")
                        {
                            string err2 = "Unknown/unsupported version of OASIS: " + test;
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
                        else
                        {
                            tableAtEnd = true;
                        }

                        break;
                    case oasValues.END:
                        if (tableAtEnd)
                            for (i = 0; i < 12; i++)
                            {
                                readUnsignedInteger();
                            }

                        break;
                    case oasValues.CELLNAME_IMPLICIT:
                        cellData[cellIndex].cellName = readString();
                        cellIndex++;
                        break;
                    case oasValues.CELLNAME:
                        modal.s = readString();
                        i = readUnsignedInteger();
                        cellData[i].cellName = modal.s;
                        break;
                    case oasValues.TEXTSTRING_IMPLICIT:
                        cellData[cellIndex - 1].textNames[cellData[cellIndex - 1].textIndex] = readString();
                        cellData[cellIndex - 1].textIndex++;
                        break;
                    case oasValues.TEXTSTRING:
                        modal.s = readString();
                        i = readUnsignedInteger();
                        cellData[cellIndex - 1].textNames[i] = modal.s;
                        break;
                    case oasValues.PROPNAME_IMPLICIT:
                        modal.s = readString();
                        break;
                    case oasValues.PROPNAME:
                        modal.s = readString();
                        i = readUnsignedInteger();
                        break;
                    case oasValues.PROPSTRING_IMPLICIT:
                        modal.s = readString();
                        break;
                    case oasValues.PROPSTRING:
                        modal.s = readString();
                        i = readUnsignedInteger();
                        break;
                    case oasValues.LAYERNAME_DATA: 
                    case oasValues.LAYERNAME_TEXT:
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

                        if (!layerNames.ContainsKey("L" + l + "D" + i))
                        {
                            layerNames.Add("L" + l + "D" + i, modal.s);
                        }
                        break;
                    case oasValues.CELL_REF_NUM:
                        cell_ = drawing_.addCell();
                        i = readUnsignedInteger();
                        if (cellData[i].cellName == "")
                        {
                            cellData[i].cellName = "layout#cell~" + i;
                        }
                        else
                        {
                            if (cellData[i].cellName == null)
                            {
                                cellData[i].cellName = modal.s;
                            }
                        }

                        cell_.cellName = cellData[i].cellName;
                        resetModal();
                        break;
                    case oasValues.CELL:
                        cell_ = drawing_.addCell();
                        cell_.cellName = readString();
                        resetModal();
                        break;
                    case oasValues.XYABSOLUTE:
                        modal.absoluteMode = true;
                        break;
                    case oasValues.XYRELATIVE:
                        modal.absoluteMode = false;
                        break;
                    case oasValues.PLACEMENT:
                        modal.mag = 1;
                        modal.angle = 0;
                        modal.mirror_x = false;
                        info_byte = readRaw();
                        if ((info_byte & 128) != 0)
                        {
                            if ((info_byte & 64) != 0)
                            {
                                i = readUnsignedInteger();
                                if (cellData[i].cellName == "")
                                {
                                    cellData[i].cellName = "layout#cell~" + i;
                                }

                                modal.placement_cell = cellData[i].cellName;
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

                        if ((info_byte & 6) == 6)
                        {
                            modal.angle = 270;
                        }
                        else if ((info_byte & 6) == 4)
                        {
                            modal.angle = 180;
                        }
                        else
                        {
                            if ((info_byte & 6) == 2)
                            {
                                modal.angle = 90;
                            }
                            else
                            {
                                modal.angle = modal.angle;
                            }
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
                        else
                        {
                            addCellref();
                        }

                        break;
                    case oasValues.PLACEMENT_TRANSFORM:
                        modal.mag = 1;
                        modal.angle = 0;
                        modal.mirror_x = false;
                        info_byte = readRaw();
                        if ((info_byte & 128) != 0)
                        {
                            if ((info_byte & 64) != 0)
                            {
                                i = readUnsignedInteger();
                                if (cellData[i].cellName == "")
                                {
                                    cellData[i].cellName = "layout#cell~" + i;
                                }

                                modal.placement_cell = cellData[i].cellName;
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
                        else
                        {
                            addCellref();
                        }

                        break;
                    case oasValues.TEXT:
                        info_byte = readRaw();
                        if ((info_byte & 64) != 0)
                        {
                            if ((info_byte & 32) != 0)
                            {
                                i = readUnsignedInteger();
                                if (cellData[cellIndex - 1].textNames[i] == "")
                                {
                                    cellData[cellIndex - 1].textNames[i] = "layout#text~" + i;
                                }
                                else
                                {
                                    if (cellData[cellIndex - 1].textNames[i] == null)
                                    {
                                        cellData[cellIndex - 1].textNames[i] = modal.s;
                                    }
                                }

                                modal.text_string = cellData[cellIndex - 1].textNames[i];
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
                        else
                        {
                            if (cellData[cellIndex - 1].textNames[cellData[cellIndex - 1].textIndex] == "")
                            {
                                cellData[cellIndex - 1].textNames[cellData[cellIndex - 1].textIndex] = "layout#text~" + cellData[cellIndex - 1].textIndex;
                                cellData[cellIndex - 1].textIndex++;
                            }
                            else
                            {
                                if (cellData[cellIndex - 1].textNames[i] == null)
                                {
                                    cellData[cellIndex - 1].textNames[i] = modal.text_string;
                                    cellData[cellIndex - 1].textIndex++;
                                }
                            }
                            addText();
                        }

                        break;
                    case oasValues.RECTANGLE:
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
                        else
                        {
                            addBox();
                        }

                        break;
                    case oasValues.POLYGON:
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
                    case oasValues.PATH:
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
                        else
                        {
                            addPath();
                        }

                        break;
                    case oasValues.TRAPEZOID_AB:
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
                        modal.trapezoid_delta_a = (int)read1Delta(false).X;
                        modal.trapezoid_delta_b = (int)read1Delta(false).X;
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
                    case oasValues.TRAPEZOID_A:
                    case oasValues.TRAPEZOID_B:
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
                            modal.trapezoid_delta_a = (int)read1Delta(false).X;
                        }
                        else
                        {
                            modal.trapezoid_delta_b = (int)read1Delta(false).X;
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
                    case oasValues.CTRAPEZOID:
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
                    case oasValues.CIRCLE:
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
                    case oasValues.PROPERTY: //property
                    case oasValues.LAST_PROPERTY:
                        info_byte = readRaw();
                        if ((info_byte & 4) != 0 && (info_byte & 2) != 0)
                        {
                            readUnsignedInteger();
                        }
                        if ((info_byte & 4) != 0 && (info_byte & 2) == 0)
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
                            else
                            {
                                count = count;
                            }

                            for (int j = 0; j < count; j++)
                            {
                                readProperty();
                            }
                        }

                        break;
                    case oasValues.XNAME_IMPLICIT:
                        readUnsignedInteger();
                        readString();
                        readUnsignedInteger();
                        break;
                    case oasValues.XNAME:
                        readUnsignedInteger();
                        readString();
                        break;
                    case oasValues.XELEMENT:
                        readUnsignedInteger();
                        readString();
                        break;
                    case oasValues.XGEOMETRY:
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
                    case oasValues.CBLOCK:
                        uint comp_type = (uint)readUnsignedInteger();
                        if (comp_type != 0)
                        {
                            string err2 = "Unsupported compression. " + record;
                            error_msgs.Add(err2);
                            throw new Exception(err2);
                        }
                        uint after = (uint)readUnsignedInteger();
                        uint before = (uint)readUnsignedInteger();
                        zLibInit(before, after);
                        break;
                    default:
                        string err = "Unknown/unsupported Record." + record;
                        error_msgs.Add(err);
                        // throw new Exception(err);
                        break;
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
                t.cellName = cellData[Convert.ToInt32(s1)].cellName;
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
                            t1.setName(cellData[Convert.ToInt32(s1)].cellName);
                            t1.setCellRef(drawing_.findCell(cellData[Convert.ToInt32(s1)].cellName));
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
                t1.setName(cellData[cellIndex - 1].textNames[Convert.ToInt32(s1)]);
            }

            // Deal with implicit values.
            // These seem to be reverse-ordered compared to the list order, so use cdIndex to transform the index.
            // Here we walk our array of cell data and then each instance of text data
            // to set any deferred values.
            int cellCount = drawing_.cellList.Count;
            // for (int cell = 0; cell < cellIndex; cell++)
            Parallel.For(0, cellCount, (cell) =>
            {
                int cdIndex = (cellCount - 1) - cell;
                drawing_.cellList[cell].cellName = cellData[cdIndex].cellName;
                int elementCount = drawing_.cellList[cell].elementList.Count;
                // for (int element = 0; element < drawing_.cellList[cell].elementList.Count; element++)
                Parallel.For(0, elementCount, (element) =>
                {
                    // Not sure if index will be weird here like above. Need a test case.
                    if (drawing_.cellList[cell].elementList[element].isText())
                    {
                        int txtIndex = (cellData[cdIndex].textIndex - 1) - cellData[cdIndex].pp_textIndex;
                        cellData[cdIndex].pp_textIndex++;
                        drawing_.cellList[cell].elementList[element].setName(cellData[cdIndex].textNames[txtIndex]);
                    }
                });
            });
            
            try
            {
                drawing_.active_cell = drawing.findCellIndex(cell_.cellName);
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