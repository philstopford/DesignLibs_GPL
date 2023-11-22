using geoCoreLib;
using MiscUtil.Conversion;
using MiscUtil.IO;
using System;
using System.IO;
using System.IO.Compression;

namespace gds;

public partial class gdsWriter
{
    public delegate void StatusUpdateUI(string text);
    public StatusUpdateUI statusUpdateUI { get; set; }
    public delegate void ProgressUpdateUI(double progress);
    public ProgressUpdateUI progressUpdateUI { get; set; }

    public EndianBinaryWriter bw { get; set; }
    public GCDrawingfield drawing_ { get; set; }
    private string filename_;
    private const bool noTerminate = false; // debug to allow file comparison during write.

    public gdsWriter(GeoCore gc, string filename)
    {
        pGDSWriter(gc, filename);
    }

    private void pGDSWriter(GeoCore gc, string filename)
    {
        if (!gc.isValid())
        {
            throw new("Provided GeoCore instance is not marked as valid");
        }
        drawing_ = gc.getDrawing();
        filename_ = filename;
    }

    public bool save()
    {
        return pSave_setup();
    }

    private bool pSave_setup()
    {
        bool compressed = filename_.ToLower().EndsWith(".gz");

        Stream s = File.Create(filename_);

        statusUpdateUI?.Invoke("Saving GDS");
        progressUpdateUI?.Invoke(0);

        bool ret = false;

        switch (compressed)
        {
            case true:
            {
                using GZipStream gzs = new(s, CompressionMode.Compress);
                bw = new EndianBinaryWriter(EndianBitConverter.Big, gzs);
                try
                {
                    pSave_write();
                    ret = true;
                }
                catch (Exception)
                {

                }

                break;
            }
            default:
                bw = new EndianBinaryWriter(EndianBitConverter.Big, s);
                try
                {
                    pSave_write();
                    ret = true;
                }
                catch (Exception)
                {

                }

                break;
        }

        s.Close();
        s.Dispose();

        return ret;
    }

    private void pSave_write()
    {
        bw.Write((ushort)6);
        bw.Write((ushort)2);
        bw.Write((ushort)600);
        // bgnlib
        bw.Write((ushort)28);
        bw.Write((byte)1);
        bw.Write((byte)2);

        // Get date and time.

        // Modification
        bw.Write((ushort)drawing_.modyear);
        bw.Write((ushort)drawing_.modmonth);
        bw.Write((ushort)drawing_.modday);
        bw.Write((ushort)drawing_.modhour);
        bw.Write((ushort)drawing_.modmin);
        bw.Write((ushort)drawing_.modsec);

        // Access
        bw.Write((ushort)drawing_.accyear);
        bw.Write((ushort)drawing_.accmonth);
        bw.Write((ushort)drawing_.accday);
        bw.Write((ushort)drawing_.acchour);
        bw.Write((ushort)drawing_.accmin);
        bw.Write((ushort)drawing_.accsec);

        writeString(drawing_.libname, 2);

        //units
        bw.Write((ushort)20);
        bw.Write((byte)3);
        bw.Write((byte)5);
        write8ByteReal(drawing_.userunits);
        write8ByteReal(1E-6 / drawing_.databaseunits);

        int cellCount = 0;
        foreach (GCCell t in drawing_.cellList)
        {
            t.saved = false;
            cellCount++;
        }

        bool saved = false;
        int cc = 0;
        int updateInterval = cellCount / 100;
        updateInterval = updateInterval switch
        {
            0 => 1,
            _ => updateInterval
        };
        double progress = 0;
        while (!saved)
        {
            saved = true;
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
                switch (t.saved)
                {
                    case false when !t.dependNotSaved():
                        t.saveGDS(this);
                        break;
                    case false:
                        saved = false;
                        break;
                }
                cc++;
            }
        }

        switch (noTerminate)
        {
            //endlib
            case false:
                bw.Write((ushort)4);
                bw.Write((byte)4);
                bw.Write((byte)0);
                break;
        }

        bw.Close();
        bw.Dispose();
    }
}