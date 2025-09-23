using System;
using System.Text;
using Clipper2Lib;
using geoCoreLib;
using geoWrangler;

namespace oasis;

internal partial class oasReader
{
    private int readUnsignedInteger()
    {
        byte help;
        int result = 0;
        int pos = 0;
        
        do
        {
            help = readRaw();
            int h = help & 127; // Remove continuation bit
            result += h << pos;
            pos += 7;
            
            // Check for overflow (32-bit unsigned integer limit)
            if (pos > 32)
            {
                string errorMsg = $"Unsigned integer overflow: more than 32 bits of data (pos={pos})";
                error_msgs.Add(errorMsg);
                throw new OverflowException(errorMsg);
            }
        }
        while (help >= 128); // Continue while continuation bit is set

        return result;
    }

    private byte readRaw()
    {
        switch (zLibUsed)
        {
            case true:
                return zLibReadRaw();
            default:
            {
                byte ret = br.ReadByte();
                return ret;
            }
        }
    }

    private void readProperty()
    {
        int valueType = readUnsignedInteger();
        
        try
        {
            switch (valueType)
            {
                case 0: // Real (positive)
                    uint realValue0 = (uint)readUnsignedInteger();
                    // Store property value for potential later processing
                    break;
                case 1: // Real (negative)
                    uint realValue1 = (uint)readUnsignedInteger();
                    // Store property value for potential later processing
                    break;
                case 2: // Real (positive reciprocal)
                    uint realValue2 = (uint)readUnsignedInteger();
                    // Store property value for potential later processing
                    break;
                case 3: // Real (negative reciprocal) 
                    uint realValue3 = (uint)readUnsignedInteger();
                    // Store property value for potential later processing
                    break;
                case 4: // Real (positive ratio)
                    uint numerator4 = (uint)readUnsignedInteger();
                    uint denominator4 = (uint)readUnsignedInteger();
                    // Store property value for potential later processing
                    break;
                case 5: // Real (negative ratio)
                    uint numerator5 = (uint)readUnsignedInteger();
                    uint denominator5 = (uint)readUnsignedInteger();
                    // Store property value for potential later processing
                    break;
                case 6: // Single-precision float
                {
                    byte[] a = new byte[4];
                    for (uint z = 0; z < 4; z++)
                    {
                        a[z] = readRaw();
                    }
                    float floatValue = BitConverter.ToSingle(a, 0);
                    // Store property value for potential later processing
                }
                    break;
                case 7: // Double-precision float
                {
                    byte[] a = new byte[8];
                    for (uint z = 0; z < 8; z++)
                    {
                        a[z] = readRaw();
                    }
                    double doubleValue = BitConverter.ToDouble(a, 0);
                    // Store property value for potential later processing
                }
                    break;
                case 8: // Unsigned integer
                    uint uintValue = (uint)readUnsignedInteger();
                    // Store property value for potential later processing
                    break;
                case 9: // Signed integer
                    int sintValue = readSignedInteger();
                    // Store property value for potential later processing
                    break;
                default: // String or reference
                    string stringValue = readString();
                    // Store property value for potential later processing
                    break;
            }
        }
        catch (Exception ex)
        {
            string errorMsg = $"Error reading property value type {valueType}: {ex.Message}";
            error_msgs.Add(errorMsg);
            throw new Exception(errorMsg, ex);
        }
    }

    private double readReal()
    {
        int formatType = readUnsignedInteger();
        
        try
        {
            switch (formatType)
            {
                case 0: // Positive integer
                    return readUnsignedInteger();
                case 1: // Negative integer
                    return -readUnsignedInteger();
                case 2: // Positive reciprocal
                {
                    uint denominator = (uint)readUnsignedInteger();
                    if (denominator == 0)
                    {
                        throw new DivideByZeroException("Reciprocal with zero denominator");
                    }
                    return 1.0 / denominator;
                }
                case 3: // Negative reciprocal
                {
                    uint denominator = (uint)readUnsignedInteger();
                    if (denominator == 0)
                    {
                        throw new DivideByZeroException("Reciprocal with zero denominator");
                    }
                    return -1.0 / denominator;
                }
                case 4: // Positive ratio
                {
                    uint numerator = (uint)readUnsignedInteger();
                    uint denominator = (uint)readUnsignedInteger();
                    if (denominator == 0)
                    {
                        throw new DivideByZeroException("Ratio with zero denominator");
                    }
                    return (double)numerator / denominator;
                }
                case 5: // Negative ratio
                {
                    uint numerator = (uint)readUnsignedInteger();
                    uint denominator = (uint)readUnsignedInteger();
                    if (denominator == 0)
                    {
                        throw new DivideByZeroException("Ratio with zero denominator");
                    }
                    return -(double)numerator / denominator;
                }
                case 6: // Single-precision IEEE float
                {
                    byte[] a = new byte[4];
                    for (uint z = 0; z < 4; z++)
                    {
                        a[z] = readRaw();
                    }
                    return BitConverter.ToSingle(a, 0);
                }
                case 7: // Double-precision IEEE float
                {
                    byte[] a = new byte[8];
                    for (uint z = 0; z < 8; z++)
                    {
                        a[z] = readRaw();
                    }
                    return BitConverter.ToDouble(a, 0);
                }
                default:
                    string errorMsg = $"Unknown real format type: {formatType}";
                    error_msgs.Add(errorMsg);
                    throw new Exception(errorMsg);
            }
        }
        catch (Exception ex) when (!(ex is DivideByZeroException))
        {
            string errorMsg = $"Error reading real value with format type {formatType}: {ex.Message}";
            error_msgs.Add(errorMsg);
            throw new Exception(errorMsg, ex);
        }
    }

    private string readString()
    {
        int length = readUnsignedInteger();
        
        if (length == 0)
        {
            return string.Empty;
        }
        
        if (length < 0)
        {
            string errorMsg = $"Invalid string length: {length}";
            error_msgs.Add(errorMsg);
            throw new Exception(errorMsg);
        }
        
        try
        {
            byte[] buffer = new byte[length];
            for (int i = 0; i < length; i++)
            {
                buffer[i] = readRaw();
            }
            
            // OASIS strings are null-terminated, so we need to handle that properly
            int actualLength = Array.IndexOf(buffer, (byte)0);
            if (actualLength >= 0)
            {
                byte[] trimmedBuffer = new byte[actualLength];
                Array.Copy(buffer, trimmedBuffer, actualLength);
                return Encoding.UTF8.GetString(trimmedBuffer);
            }
            
            return Encoding.UTF8.GetString(buffer);
        }
        catch (Exception ex)
        {
            string errorMsg = $"Error reading string of length {length}: {ex.Message}";
            error_msgs.Add(errorMsg);
            throw new Exception(errorMsg, ex);
        }
    }

    private int readSignedInteger()
    {
        int pos = 0;
        byte help = readRaw();
        bool negative = false;
        int h = help & 127;
        
        // Check sign bit (LSB)
        negative = (h & 1) == 1;
        int result = h >> 1; // Remove sign bit
        pos += 6; // 7 bits total, minus 1 sign bit
        
        while (help >= 128) // Continue bit set
        {
            help = readRaw();
            h = help & 127;
            result += h << pos;
            pos += 7;
        }
        
        if (pos > 31) // Check for overflow (32-bit signed integer limit)
        {
            string errorMsg = $"Signed integer overflow: more than 31 bits of data (pos={pos})";
            error_msgs.Add(errorMsg);
            throw new OverflowException(errorMsg);
        }

        return negative ? -result : result;
    }

    private Point64 read1Delta(bool dir)
    {
        // dir true-> vertical
        int i = readUnsignedInteger();

        switch (i % 2)
        {
            case 0:
                i >>= 1;
                break;
            case 1:
                i = -(i >> 1);
                break;
        }

        return dir switch
        {
            true => new Point64(0, i),
            _ => new Point64(i, 0)
        };
    }

    private Point64 readGDelta()
    {
        int i = readUnsignedInteger();
        int k;
        switch (i % 2)
        {
            case 0: //form 0
                i >>= 1;
                k = i % 8;
                i >>= 3;
                switch (k)
                {
                    case 0:
                        return new Point64(i, 0);
                    case 1:
                        return new Point64(0, i);
                    case 2:
                        return new Point64(-i, 0);
                    case 3:
                        return new Point64(0, -i);
                    case 4:
                        return new Point64(i, i);
                    case 5:
                        return new Point64(-i, i);
                    case 6:
                        return new Point64(-i, -i);
                    case 7:
                        return new Point64(i, -i);
                }
                break;
            case 1: //form 1
                i >>= 1;
                switch (i % 2)
                {
                    case 0:
                        i >>= 1;
                        break;
                    default:
                        i = -(i >> 1);
                        break;
                }
                k = readUnsignedInteger();
                switch (k % 2)
                {
                    case 0:
                        k >>= 1;
                        break;
                    default:
                        k = -(k >> 1);
                        break;
                }
                return new Point64(i, k);
        }
        return new Point64(0, 0);
    }

    private Point64 read2Delta()
    {
        int i = readUnsignedInteger();
        switch (i % 4)
        {
            case 0:
                i >>= 2;
                return new Point64(i, 0);
            case 1:
                i >>= 2;
                return new Point64(0, i);
            case 2:
                i = -(i >> 2);
                return new Point64(i, 0);
            case 3:
                i = -(i >> 2);
                return new Point64(0, i);
        }
        return new Point64(0, 0);
    }

    private Point64 read3Delta()
    {
        int i = readUnsignedInteger();
        switch (i % 8)
        {
            case 0:
                i >>= 3;
                return new Point64(i, 0);
            case 1:
                i >>= 3;
                return new Point64(0, i);
            case 2:
                i >>= 3;
                return new Point64(-i, 0);
            case 3:
                i >>= 3;
                return new Point64(0, -i);
            case 4:
                i >>= 3;
                return new Point64(i, i);
            case 5:
                i >>= 3;
                return new Point64(-i, i);
            case 6:
                i >>= 3;
                return new Point64(-i, -i);
            case 7:
                i >>= 3;
                return new Point64(i, -i);
        }
        return new Point64(0, 0);
    }

    private void readRepetition()
    {
        modal.repetition = new Repetition();
        int repetitionType = readUnsignedInteger();
        
        try
        {
            switch (repetitionType)
            {
                case 0: // No repetition
                    modal.repetition.type = Repetition.RepetitionType.None;
                    break;
                    
                case 1: // Rectangular array, regular spacing
                    modal.repetition.type = Repetition.RepetitionType.Rectangular;
                    modal.repetition.columns = readUnsignedInteger() + 2;
                    modal.repetition.rows = readUnsignedInteger() + 2;
                    modal.repetition.spacing.X = readUnsignedInteger();
                    modal.repetition.spacing.Y = readUnsignedInteger();
                    
                    if (modal.repetition.columns < 2 || modal.repetition.rows < 2)
                    {
                        throw new Exception($"Invalid repetition dimensions: {modal.repetition.columns}x{modal.repetition.rows}");
                    }
                    break;
                    
                case 2: // One-dimensional array, X-axis spacing
                    modal.repetition.type = Repetition.RepetitionType.Rectangular;
                    modal.repetition.columns = readUnsignedInteger() + 2;
                    modal.repetition.rows = 1;
                    modal.repetition.spacing.X = readUnsignedInteger();
                    modal.repetition.spacing.Y = 0;
                    
                    if (modal.repetition.columns < 2)
                    {
                        throw new Exception($"Invalid repetition columns: {modal.repetition.columns}");
                    }
                    break;
                    
                case 3: // One-dimensional array, Y-axis spacing
                    modal.repetition.type = Repetition.RepetitionType.Rectangular;
                    modal.repetition.columns = 1;
                    modal.repetition.rows = readUnsignedInteger() + 2;
                    modal.repetition.spacing.X = 0;
                    modal.repetition.spacing.Y = readUnsignedInteger();
                    
                    if (modal.repetition.rows < 2)
                    {
                        throw new Exception($"Invalid repetition rows: {modal.repetition.rows}");
                    }
                    break;
                    
                case 4: // One-dimensional array, explicit X coordinates
                case 5: // One-dimensional array, explicit X coordinates with grid
                {
                    modal.repetition.type = Repetition.RepetitionType.ExplicitX;
                    int count = readUnsignedInteger() + 1;
                    
                    if (count <= 0 || count > 100000) // Sanity check
                    {
                        throw new Exception($"Invalid repetition count: {count}");
                    }
                    
                    double grid_factor = 1.0;
                    if (repetitionType == 5)
                    {
                        grid_factor = readUnsignedInteger();
                        if (grid_factor <= 0)
                        {
                            throw new Exception($"Invalid grid factor: {grid_factor}");
                        }
                    }

                    modal.repetition.coords.Clear();
                    for (int i = 0; i < count; i++)
                    {
                        double coord = grid_factor * readUnsignedInteger();
                        modal.repetition.coords.Add(coord);
                    }
                    break;
                }
                
                case 6: // One-dimensional array, explicit Y coordinates
                case 7: // One-dimensional array, explicit Y coordinates with grid
                {
                    modal.repetition.type = Repetition.RepetitionType.ExplicitY;
                    int count = readUnsignedInteger() + 1;
                    
                    if (count <= 0 || count > 100000) // Sanity check
                    {
                        throw new Exception($"Invalid repetition count: {count}");
                    }
                    
                    double grid_factor = 1.0;
                    if (repetitionType == 7)
                    {
                        grid_factor = readUnsignedInteger();
                        if (grid_factor <= 0)
                        {
                            throw new Exception($"Invalid grid factor: {grid_factor}");
                        }
                    }

                    modal.repetition.coords.Clear();
                    for (int i = 0; i < count; i++)
                    {
                        double coord = grid_factor * readUnsignedInteger();
                        modal.repetition.coords.Add(coord);
                    }
                    break;
                }
                
                case 8: // Two-dimensional array, arbitrary vectors
                    modal.repetition.type = Repetition.RepetitionType.Regular;
                    modal.repetition.columns = readUnsignedInteger() + 2;
                    modal.repetition.rows = readUnsignedInteger() + 2;
                    modal.repetition.colVector = readGDelta();
                    modal.repetition.rowVector = readGDelta();
                    
                    if (modal.repetition.columns < 2 || modal.repetition.rows < 2)
                    {
                        throw new Exception($"Invalid repetition dimensions: {modal.repetition.columns}x{modal.repetition.rows}");
                    }
                    break;
                    
                case 9: // One-dimensional array with vector
                    modal.repetition.type = Repetition.RepetitionType.Regular;
                    modal.repetition.columns = readUnsignedInteger() + 2;
                    modal.repetition.rows = 1;
                    modal.repetition.colVector = readGDelta();
                    modal.repetition.rowVector = new Point64(modal.repetition.colVector.Y, modal.repetition.colVector.X);
                    
                    if (modal.repetition.columns < 2)
                    {
                        throw new Exception($"Invalid repetition columns: {modal.repetition.columns}");
                    }
                    break;
                    
                case 10: // Arbitrary array, explicit coordinate pairs
                case 11: // Arbitrary array, explicit coordinate pairs with grid
                {
                    modal.repetition.type = Repetition.RepetitionType.Explicit;
                    int count = readUnsignedInteger() + 1;
                    
                    if (count <= 0 || count > 100000) // Sanity check
                    {
                        throw new Exception($"Invalid repetition count: {count}");
                    }
                    
                    int grid_factor = 1;
                    if (repetitionType == 11)
                    {
                        grid_factor = readUnsignedInteger();
                        if (grid_factor <= 0)
                        {
                            throw new Exception($"Invalid grid factor: {grid_factor}");
                        }
                    }

                    modal.repetition.offsets.Clear();
                    for (int i = 0; i < count; i++)
                    {
                        Point64 p = readGDelta();
                        modal.repetition.offsets.Add(new Point64(p.X * grid_factor, p.Y * grid_factor));
                    }
                    break;
                }
                
                default:
                    string errorMsg = $"Unknown repetition type: {repetitionType}";
                    error_msgs.Add(errorMsg);
                    throw new Exception(errorMsg);
            }
        }
        catch (Exception ex)
        {
            string errorMsg = $"Error reading repetition type {repetitionType}: {ex.Message}";
            error_msgs.Add(errorMsg);
            throw new Exception(errorMsg, ex);
        }
    }

    private void readExtension()
    {
        int i = readUnsignedInteger();
        switch (i & 3)
        {
            case 3:
                modal.path_start_extension_value = readSignedInteger();
                modal.path_start_extension = 4;
                break;
            case 2:
                modal.path_start_extension = 2;
                break;
            case 1:
                modal.path_start_extension = 0;
                break;
            default:
                modal.path_start_extension = modal.path_start_extension;
                break;
        }

        switch (i & 12)
        {
            case 12:
                modal.path_end_extension_value = readSignedInteger();
                modal.path_end_extension = 4;
                break;
            case 8:
                modal.path_end_extension = 2;
                break;
            case 4:
                modal.path_end_extension = 0;
                break;
            default:
                modal.path_end_extension = modal.path_end_extension;
                break;
        }
    }

    private void readPointList(bool addImplecid)
    {
        int type = readUnsignedInteger();
        int count = readUnsignedInteger();
        Path64 pointlist = Helper.initedPath64(count + 1);
        pointlist[0] = new Point64(0, 0);
        Point64 last = new(0, 0);
        Point64 oPt;
        int i;
        bool dir = true;
        switch (type)
        {
            case 0:
                dir = false;
                for (i = 1; i <= count; i++)
                {
                    oPt = read1Delta(dir);
                    last = GeoWrangler.move(last, oPt.X, oPt.Y);
                    pointlist[i] = new Point64(last);
                    dir = !dir;
                }

                if (addImplecid)
                {
                    pointlist.Add(new Point64(pointlist[^1].X, 0));
                    pointlist.Add(new Point64(pointlist[0]));
                }

                break;
            case 1:
                for (i = 1; i <= count; i++)
                {
                    oPt = read1Delta(dir);
                    last = GeoWrangler.move(last, oPt.X, oPt.Y);
                    pointlist[i] = new Point64(last);
                    dir = !dir;
                }

                if (addImplecid)
                {
                    pointlist.Add(new Point64(pointlist[^1].X, 0));
                    pointlist.Add(new Point64(pointlist[0]));
                }

                break;
            case 2:
                for (i = 1; i <= count; i++)
                {
                    oPt = read2Delta();
                    last = GeoWrangler.move(last, oPt.X, oPt.Y);
                    pointlist[i] = new Point64(last);
                }

                if (addImplecid)
                {
                    pointlist.Add(new Point64(pointlist[0]));
                }

                break;
            case 3:
                for (i = 1; i <= count; i++)
                {
                    oPt = read3Delta();
                    last = GeoWrangler.move(last, oPt.X, oPt.Y);
                    pointlist[i] = new Point64(last);
                }

                if (addImplecid)
                {
                    pointlist.Add(new Point64(pointlist[0]));
                }

                break;
            case 4:
                for (i = 1; i <= count; i++)
                {
                    oPt = readGDelta();
                    last = GeoWrangler.move(last, oPt.X, oPt.Y);
                    pointlist[i] = new Point64(last);
                }

                if (addImplecid)
                {
                    pointlist.Add(new Point64(pointlist[0]));
                }

                break;
            case 5:
            {
                Point64 l = new(0, 0);
                for (i = 1; i <= count; i++)
                {
                    oPt = readGDelta();
                    last = GeoWrangler.move(last, oPt.X, oPt.Y);
                    l = GeoWrangler.move(l, last.X, last.Y);
                    pointlist[i] = new Point64(l);
                }

                if (addImplecid)
                {
                    pointlist.Add(new Point64(pointlist[0]));
                }
            }
                break;
        }
        modal.polygon_point_list = new Path64(pointlist);
    }

    private void zLibInit(uint before, uint after)
    {
        if (before == 0)
        {
            string errorMsg = "Invalid compressed data size: 0 bytes";
            error_msgs.Add(errorMsg);
            throw new Exception(errorMsg);
        }
        
        if (after == 0)
        {
            string errorMsg = "Invalid uncompressed data size: 0 bytes";
            error_msgs.Add(errorMsg);
            throw new Exception(errorMsg);
        }
        
        try
        {
            zLibUsed = true;
            byte[] compressedData = br.ReadBytes((int)before);
            
            if (compressedData.Length != before)
            {
                string errorMsg = $"Could not read expected {before} bytes of compressed data, got {compressedData.Length} bytes";
                error_msgs.Add(errorMsg);
                throw new Exception(errorMsg);
            }
            
            zLibOut = utility.Utils.decompress(compressedData);
            
            if (zLibOut.Length != after)
            {
                string errorMsg = $"Decompressed data size mismatch: expected {after} bytes, got {zLibOut.Length} bytes";
                error_msgs.Add(errorMsg);
                throw new Exception(errorMsg);
            }
            
            zlibOutPos = 0;
        }
        catch (Exception ex)
        {
            zLibUsed = false;
            string errorMsg = $"Failed to initialize zlib decompression (before={before}, after={after}): {ex.Message}";
            error_msgs.Add(errorMsg);
            throw new Exception(errorMsg, ex);
        }
    }

    private byte zLibReadRaw()
    {
        if (!zLibUsed)
        {
            throw new InvalidOperationException("Attempting to read from zlib when not initialized");
        }
        
        if (zlibOutPos >= zLibOut.Length)
        {
            throw new InvalidOperationException("Attempting to read beyond end of zlib data");
        }
        
        byte result = zLibOut[zlibOutPos];
        zlibOutPos++;
        
        if (zlibOutPos >= zLibOut.Length)
        {
            zLibUsed = false;
        }
        
        return result;
    }
}