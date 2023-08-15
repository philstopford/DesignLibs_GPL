using System.Text;

namespace Veldrid.SPIRV
{
    internal static class Util
    {
        public static Encoding UTF8 { get; } = new UTF8Encoding(encoderShouldEmitUTF8Identifier: false);

        internal static unsafe string GetString(byte* data, uint length)
        {
            if (data == null) { return null; }

            return UTF8.GetString(data, (int)length);
        }

        internal static bool HasSpirvHeader(byte[] bytes)
        {
            return bytes.Length >= 4
                && bytes[0] == 0x03
                && bytes[1] == 0x02
                && bytes[2] == 0x23
                && bytes[3] == 0x07;
        }
    }
}
