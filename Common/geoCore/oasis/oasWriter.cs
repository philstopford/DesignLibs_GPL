using geoCoreLib;
using MiscUtil.Conversion;
using MiscUtil.IO;
using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using Clipper2Lib;

namespace oasis;

public partial class oasWriter
{
    public delegate void StatusUpdateUI(string text);
    public StatusUpdateUI statusUpdateUI { get; set; }
    public delegate void ProgressUpdateUI(double progress);
    public ProgressUpdateUI progressUpdateUI { get; set; }

    private EndianBinaryWriter bw;
    private GCDrawingfield drawing_;
    private Dictionary<string, string> namedLayers;
    private string filename_;
    public string error_message;

    public Modals modal;

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
        public int texttype { get; set; }
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
        public int path_end_extension{ get; set; }
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
        public Path64 repArray { get; set; }
        public double trapezonid_delta_a { get; set; }
        public double trapezonid_delta_b { get; set; }
        public bool trapezonid_orientation { get; set; }
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
        modal.polygon_point_list = new ();
        modal.repArray = new ();
    }

    public oasWriter(GeoCore gc, string filename)
    {
        pOASWriter(gc, filename);
    }

    private void pOASWriter(GeoCore gc, string filename)
    {
        if (!gc.isValid())
        {
            throw new("Provided GeoCore instance is not marked as valid");
        }
        drawing_ = gc.getDrawing();
        filename_ = filename;
        namedLayers = gc.getLayerNames();
    }

    public bool save()
    {
        return pSave_setup();
    }

    private bool pSave_setup()
    {

        bool compressed = filename_.ToLower().EndsWith(".gz");

        statusUpdateUI?.Invoke("Saving Oasis");
        progressUpdateUI?.Invoke(0);

        FileStream s = File.OpenWrite(filename_);

        bool ret = false;

        if (compressed)
        {
            using GZipStream gzs = new(s, CompressionMode.Compress);
            bw = new EndianBinaryWriter(EndianBitConverter.Little, gzs);
        }
        else
        {
            bw = new EndianBinaryWriter(EndianBitConverter.Little, s);
        }

        try
        {
            pSave_write();
            ret = true;
        }
        catch (Exception e)
        {
            error_message = e.ToString();
        }

        s.Close();
        s.Dispose();

        return ret;
    }

    private void pSave_write()
    {
        bw.Write("%SEMI-OASIS".ToCharArray());
        bw.Write((byte)13);
        bw.Write((byte)10);

        // start record
        bw.Write((byte)1);
        writeString("1.0");
        double tmp = drawing_.getDrawingScale();
        writeReal(tmp * 1E3);

        //  offset table:
        for (int i = 0; i < 13; ++i)
        {
            writeRaw(0);
        }

        writeUnsignedInteger(28);
        writeUnsignedInteger(21);
        writeString("S_MAX_SIGNED_INTEGER_WIDTH");
        writeUnsignedInteger(8);
        writeUnsignedInteger(4);
        writeUnsignedInteger(28);
        writeUnsignedInteger(13);
        writeString("S_MAX_UNSIGNED_INTEGER_WIDTH");
        writeUnsignedInteger(28);
        writeUnsignedInteger(21);
        writeString("S_TOP_CELL");

        int cellCount = 0;
        foreach (GCCell t in drawing_.cellList)
        {
            t.saved = false;
            cellCount++;
        }

        List<string> cellNames_beingWritten = new();
        for (int i = 0; i < drawing_.cellList.Count; i++)
        {
            switch (drawing_.cellList[i].elementList)
            {
                case null:
                    continue;
            }

            switch (drawing_.cellList[i].saved)
            {
                case false:
                {
                    writeUnsignedInteger(12);
                    cellNames_beingWritten.Add(drawing_.cellList[i].cellName);
                    writeString(drawing_.cellList[i].cellName);
                    if (i != drawing_.cellList.Count - 1)
                    {
                        writeUnsignedInteger(28);
                        writeUnsignedInteger(17);
                    }
                    drawing_.cellList[i].saved = true;
                    break;
                }
            }
        }

        cellNames_beingWritten.Reverse();

        foreach (string t in cellNames_beingWritten)
        {
            bw.Write((byte)3);
            writeString(t);
        }

        // Write out the layer names.
        // Layer name mapping
        foreach (KeyValuePair<string, string> entry in namedLayers)
        {
            // split the key.
            string key = entry.Key; // LxxxDyyy
            string[] tokens = key.Split('D');
            int datatype = Convert.ToInt32(tokens[1]);
            int layer = Convert.ToInt32(tokens[0].Split('L')[1]);

            string name = entry.Value;

            bw.Write((byte)11);
            writeString(name);
            bw.Write((byte)3);
            writeUnsignedInteger((uint)layer);
            bw.Write((byte)3);
            writeUnsignedInteger((uint)datatype);
            bw.Write((byte)12);
            writeString(name);
            bw.Write((byte)3);
            writeUnsignedInteger((uint)layer);
            bw.Write((byte)3);
            writeUnsignedInteger((uint)datatype);
        }

        foreach (GCCell t in drawing_.cellList)
        {
            t.saved = false;
        }

        int cellCounter = cellNames_beingWritten.Count;
        int cc = 0;
        int updateInterval = cellCount / 100;
        updateInterval = updateInterval switch
        {
            0 => 1,
            _ => updateInterval
        };
        double progress = 0;
        foreach (GCCell t in drawing_.cellList)
        {
            switch (cc % updateInterval)
            {
                case 0:
                    statusUpdateUI?.Invoke(t.cellName);
                    progressUpdateUI?.Invoke(progress);
                    progress += 0.01;
                    break;
            }
            switch (t.elementList)
            {
                case null:
                    continue;
            }
            cellCounter--;
            switch (t.saved)
            {
                case false:
                {
                    if (!t.dependNotSaved())
                    {
                        resetModal();
                        bw.Write((byte)13);
                        writeUnsignedInteger((uint)cellCounter);
                        t.saveOASIS(this);
                        t.saved = true;
                    }

                    break;
                }
            }
            cc++;
        }

        writeUnsignedInteger(2);
        for (int i = 0; i < 253; i++)
        {
            bw.Write((byte)128);
        }
        writeUnsignedInteger(0);
        writeUnsignedInteger(0);

        bw.Close();
        bw.Dispose();
    }
}