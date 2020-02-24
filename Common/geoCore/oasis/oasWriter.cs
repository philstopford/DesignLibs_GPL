using geoCoreLib;
using geoLib;
using MiscUtil.Conversion;
using MiscUtil.IO;
using System;
using System.Collections.Generic;
using System.IO;

namespace oasis
{
    public partial class oasWriter
    {
        public delegate void StatusUpdateUI(string text);
        public StatusUpdateUI statusUpdateUI { get; set; }
        public delegate void ProgressUpdateUI(double progress);
        public ProgressUpdateUI progressUpdateUI { get; set; }

        EndianBinaryWriter bw;
        GCDrawingfield drawing_;
        Dictionary<string, string> namedLayers;
        string filename_;

        public modals modal;

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
            public Int32 texttype { get; set; }
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
            public Int32 path_end_extension{ get; set; }
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
            public double trapezonid_delta_a { get; set; }
            public double trapezonid_delta_b { get; set; }
            public bool trapezonid_orientation { get; set; }
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

        public oasWriter(GeoCore gc, String filename)
        {
            pOASWriter(gc, filename);
        }

        void pOASWriter(GeoCore gc, String filename)
        {
            drawing_ = gc.getDrawing();
            filename_ = filename;
            namedLayers = gc.getLayerNames();
        }

        public bool save()
        {
            return pSave();
        }

        bool pSave()
        {
            bool ok = true;

            statusUpdateUI?.Invoke("Saving Oasis");
            progressUpdateUI?.Invoke(0);

            Stream s = File.OpenWrite(filename_);
            bw = new EndianBinaryWriter(EndianBitConverter.Little, s);

            bw.Write("%SEMI-OASIS".ToCharArray());
            bw.Write((byte)13);
            bw.Write((byte)10);

            // start record
            bw.Write((byte)1);
            writeString("1.0");
            writeReal(drawing_.databaseunits);

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
            for (int i = 0; i < drawing_.cellList.Count; i++)
            {
                drawing_.cellList[i].saved = false;
                cellCount++;
            }

            List<string> cellNames_beingWritten = new List<string>();
            for (int i = 0; i < drawing_.cellList.Count; i++)
            {
                if (drawing_.cellList[i].elementList == null)
                {
                    continue;
                }

                if (!drawing_.cellList[i].saved)
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
                }
            }

            cellNames_beingWritten.Reverse();

            for (int i = 0; i < cellNames_beingWritten.Count; i++)
            {
                bw.Write((byte)3);
                writeString(cellNames_beingWritten[i]);
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

            for (int i = 0; i < drawing_.cellList.Count; i++)
            {
                drawing_.cellList[i].saved = false;
            }

            int cellCounter = cellNames_beingWritten.Count;
            int cc = 0;
            int updateInterval = cellCount / 100;
            if (updateInterval == 0)
            {
                updateInterval = 1;
            }
            double progress = 0;
            for (int i = 0; i < drawing_.cellList.Count; i++)
            {
                if (cc % updateInterval == 0)
                {
                    statusUpdateUI?.Invoke(drawing_.cellList[i].cellName);
                    progressUpdateUI?.Invoke(progress);
                    progress += 0.01;
                }
                if (drawing_.cellList[i].elementList == null)
                {
                    continue;
                }
                cellCounter--;
                if (drawing_.cellList[i].saved == false)
                {
                    if (!drawing_.cellList[i].dependNotSaved())
                    {
                        resetModal();
                        bw.Write((byte)13);
                        writeUnsignedInteger((uint)cellCounter);
                        drawing_.cellList[i].saveOASIS(this);
                        drawing_.cellList[i].saved = true;
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
            s.Close();
            s.Dispose();
            return ok;
        }
    }
}
