using geoCoreLib;
using MiscUtil.Conversion;
using MiscUtil.IO;
using System;
using System.IO;
using System.IO.Compression;

namespace gds
{
    public partial class gdsWriter
    {
        public delegate void StatusUpdateUI(string text);
        public StatusUpdateUI statusUpdateUI { get; set; }
        public delegate void ProgressUpdateUI(double progress);
        public ProgressUpdateUI progressUpdateUI { get; set; }

        public EndianBinaryWriter bw { get; set; }
        public GCDrawingfield drawing_ { get; set; }
        string filename_;
        bool ok;
        bool noTerminate = false; // debug to allow file comparison during write.

        public gdsWriter(GeoCore gc, String filename)
        {
            pGDSWriter(gc, filename);
        }

        void pGDSWriter(GeoCore gc, String filename)
        {
            drawing_ = gc.getDrawing();
            filename_ = filename;
        }

        public bool save()
        {
            return pSave();
        }

        bool pSave()
        {
            ok = true;

            bool compressed = filename_.ToLower().EndsWith(".gz");

            Stream s = File.Create(filename_);

            statusUpdateUI?.Invoke("Saving GDS");
            progressUpdateUI?.Invoke(0);

            if (compressed)
            {
                using (GZipStream gzs = new GZipStream(s, CompressionMode.Compress))
                {
                    bw = new EndianBinaryWriter(EndianBitConverter.Big, gzs);
                }
            }
            else
            {
                bw = new EndianBinaryWriter(EndianBitConverter.Big, s);
            }

            bw.Write((UInt16)6);
            bw.Write((UInt16)2);
            bw.Write((UInt16)600);
            // bgnlib
            bw.Write((UInt16)28);
            bw.Write((byte)1);
            bw.Write((byte)2);

            // Get date and time.

            // Modification
            bw.Write((UInt16)(drawing_.modyear));
            bw.Write((UInt16)(drawing_.modmonth));
            bw.Write((UInt16)(drawing_.modday));
            bw.Write((UInt16)(drawing_.modhour));
            bw.Write((UInt16)(drawing_.modmin));
            bw.Write((UInt16)(drawing_.modsec));

            // Access
            bw.Write((UInt16)(drawing_.accyear));
            bw.Write((UInt16)(drawing_.accmonth));
            bw.Write((UInt16)(drawing_.accday));
            bw.Write((UInt16)(drawing_.acchour));
            bw.Write((UInt16)(drawing_.accmin));
            bw.Write((UInt16)(drawing_.accsec));

            writeString(drawing_.libname, 2);

            //units
            bw.Write((UInt16)20);
            bw.Write((byte)3);
            bw.Write((byte)5);
            write8ByteReal(drawing_.userunits);
            write8ByteReal(1E-6 / drawing_.databaseunits);

            int cellCount = 0;
            for (int i = 0; i < drawing_.cellList.Count; i++)
            {
                drawing_.cellList[i].saved = false;
                cellCount++;
            }

            bool saved = false;
            int cc = 0;
            int updateInterval = cellCount / 100;
            if (updateInterval == 0)
            {
                updateInterval = 1;
            }
            double progress = 0;
            while (!saved)
            {
                saved = true;
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
                    if (drawing_.cellList[i].saved == false)
                    {
                        if (!drawing_.cellList[i].dependNotSaved())
                        {
                            drawing_.cellList[i].saveGDS(this);
                        }
                        else
                        {
                            saved = false;
                        }
                    };
                    cc++;
                }
            }
            //endlib
            if (!noTerminate)
            {
                bw.Write((UInt16)4);
                bw.Write((byte)4);
                bw.Write((byte)0);
            }

            bw.Close();
            bw.Dispose();
            s.Close();
            s.Dispose();

            return ok;
        }
    }
}