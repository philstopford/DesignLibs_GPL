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
    public StatusUpdateUI statusUpdateUI { get; init; }
    public delegate void ProgressUpdateUI(double progress);
    public ProgressUpdateUI progressUpdateUI { get; init; }

    public EndianBinaryWriter bw { get; set; }
    public GCDrawingfield drawing_ { get; set; }
    private string filename_;
    private const bool noTerminate = false; // debug to allow file comparison during write.
    public string error_message;

    public gdsWriter(GeoCore gc, string filename)
    {
        pGDSWriter(gc, filename);
    }

    private void pGDSWriter(GeoCore gc, string filename)
    {
        if (!gc.isValid())
        {
            throw new Exception("Provided GeoCore instance is not marked as valid");
        }
        drawing_ = gc.getDrawing();
        
        // Calculate scaling needed to convert to standardized 1E-11 database units
        double targetDbUnits = 1E-11; // 0.01 nm in meters
        double scaleFactor = drawing_.databaseunits / targetDbUnits;
        
        // Apply scaling to coordinates if needed
        if (Math.Abs(scaleFactor - 1.0) > 1E-15) // Only scale if significantly different
        {
            drawing_.resize(scaleFactor);
        }
        
        // Set standardized units
        drawing_.userunits = 1E-3;   // 1 mm in meters  
        drawing_.databaseunits = 1E-11; // 0.01 nm in meters
        
        filename_ = filename;
    }

    public bool save()
    {
        return pSave_setup();
    }

    private bool pSave_setup()
    {
        bool compressed = filename_.ToLower().EndsWith(".gz");

        FileStream s = File.Create(filename_);

        statusUpdateUI?.Invoke("Saving GDS");
        progressUpdateUI?.Invoke(0);

        bool ret = false;
        error_message = "";

        if (compressed)
        {
            using GZipStream gzs = new(s, CompressionMode.Compress);
            bw = new EndianBinaryWriter(EndianBitConverter.Big, gzs);
        }
        else
        {
            bw = new EndianBinaryWriter(EndianBitConverter.Big, s);
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
        bw.Write((ushort)6);
        bw.Write(gdsValues.sHEADER);
        bw.Write((ushort)600);
        // bgnlib
        bw.Write((ushort)28);
        bw.Write(gdsValues.sBGNLIB);

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
        bw.Write(gdsValues.sUNITS);
        // Standardize to 0.01 nm (1E-11 m) database units for consistency and 3rd party compatibility
        write8ByteReal(1E-3);  // userunits: 1 mm in meters
        write8ByteReal(1E-11); // databaseunits: 0.01 nm in meters

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
                if (cc % updateInterval == 0)
                {
                    statusUpdateUI?.Invoke(t.cellName);
                    progressUpdateUI?.Invoke(progress);
                    progress += 0.01;
                }

                if (t.elementList == null)
                {
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
                bw.Write(gdsValues.sENDLIB);
                break;
        }

        bw.Close();
        bw.Dispose();
    }
}