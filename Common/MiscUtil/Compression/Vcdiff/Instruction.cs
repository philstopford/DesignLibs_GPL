
namespace MiscUtil.Compression.Vcdiff;

/// <summary>
/// Contains the information for a single instruction
/// </summary>
internal struct Instruction
{
    internal InstructionType Type { get; }

    internal byte Size { get; }

    internal byte Mode { get; }

    internal Instruction(InstructionType type, byte size, byte mode)
    {
        this.Type = type;
        this.Size = size;
        this.Mode = mode;
    }


}