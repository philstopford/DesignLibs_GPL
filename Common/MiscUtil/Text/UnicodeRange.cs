using MiscUtil.Collections;
using System.Collections.Generic;
using System.Linq;

namespace MiscUtil.Text;

/// <summary>
/// Utility class providing a number of singleton instances of
/// Range&lt;char&gt; to indicate the various ranges of unicode characters,
/// as documented at http://msdn.microsoft.com/en-us/library/20bw873z.aspx.
/// Note that this does not indicate the Unicode category of a character,
/// merely which range it's in.
/// TODO: Work out how to include names. Can't derive from Range[char].
/// </summary>
public static class UnicodeRange
{
    private static readonly List<Range<char>> allRanges = [];

    private static Range<char> CreateRange(char from, char to)
    {
        // TODO: Check for overlaps
        Range<char> ret = new(from, to);
        allRanges.Add(ret);
        return ret;
    }

#pragma warning disable 1591
    public static Range<char> BasicLatin { get; } = CreateRange('\u0000', '\u007f');

    public static Range<char> Latin1Supplement { get; } = CreateRange('\u0080', '\u00ff');

    public static Range<char> LatinExtendedA { get; } = CreateRange('\u0100', '\u017f');

    public static Range<char> LatinExtendedB { get; } = CreateRange('\u0180', '\u024f');

    public static Range<char> IpaExtensions { get; } = CreateRange('\u0250', '\u02af');

    public static Range<char> SpacingModifierLetters { get; } = CreateRange('\u02b0', '\u02ff');

    public static Range<char> CombiningDiacriticalMarks { get; } = CreateRange('\u0300', '\u036f');

    public static Range<char> GreekAndCoptic { get; } = CreateRange('\u0370', '\u03ff');

    public static Range<char> Cyrillic { get; } = CreateRange('\u0400', '\u04ff');

    public static Range<char> CyrillicSupplement { get; } = CreateRange('\u0500', '\u052f');

    public static Range<char> Armenian { get; } = CreateRange('\u0530', '\u058f');

    public static Range<char> Hebrew { get; } = CreateRange('\u0590', '\u05FF');

    public static Range<char> Arabic { get; } = CreateRange('\u0600', '\u06ff');

    public static Range<char> Syriac { get; } = CreateRange('\u0700', '\u074f');

    public static Range<char> Thaana { get; } = CreateRange('\u0780', '\u07bf');

    public static Range<char> Devangari { get; } = CreateRange('\u0900', '\u097f');

    public static Range<char> Bengali { get; } = CreateRange('\u0980', '\u09ff');

    public static Range<char> Gurmukhi { get; } = CreateRange('\u0a00', '\u0a7f');

    public static Range<char> Gujarati { get; } = CreateRange('\u0a80', '\u0aff');

    public static Range<char> Oriya { get; } = CreateRange('\u0b00', '\u0b7f');

    public static Range<char> Tamil { get; } = CreateRange('\u0b80', '\u0bff');

    public static Range<char> Telugu { get; } = CreateRange('\u0c00', '\u0c7f');

    public static Range<char> Kannada { get; } = CreateRange('\u0c80', '\u0cff');

    public static Range<char> Malayalam { get; } = CreateRange('\u0d00', '\u0d7f');

    public static Range<char> Sinhala { get; } = CreateRange('\u0d80', '\u0dff');

    public static Range<char> Thai { get; } = CreateRange('\u0e00', '\u0e7f');

    public static Range<char> Lao { get; } = CreateRange('\u0e80', '\u0eff');

    public static Range<char> Tibetan { get; } = CreateRange('\u0f00', '\u0fff');

    public static Range<char> Myanmar { get; } = CreateRange('\u1000', '\u109f');

    public static Range<char> Georgian { get; } = CreateRange('\u10a0', '\u10ff');

    public static Range<char> HangulJamo { get; } = CreateRange('\u1100', '\u11ff');

    public static Range<char> Ethiopic { get; } = CreateRange('\u1200', '\u137f');

    public static Range<char> Cherokee { get; } = CreateRange('\u13a0', '\u13ff');

    public static Range<char> UnifiedCanadianAboriginalSyllabics { get; } = CreateRange('\u1400', '\u167f');

    public static Range<char> Ogham { get; } = CreateRange('\u1680', '\u169f');

    public static Range<char> Runic { get; } = CreateRange('\u16a0', '\u16ff');

    public static Range<char> Tagalog { get; } = CreateRange('\u1700', '\u171f');

    public static Range<char> Hanunoo { get; } = CreateRange('\u1720', '\u173f');

    public static Range<char> Buhid { get; } = CreateRange('\u1740', '\u175f');

    public static Range<char> Tagbanwa { get; } = CreateRange('\u1760', '\u177f');

    public static Range<char> Khmer { get; } = CreateRange('\u1780', '\u17ff');

    public static Range<char> Mongolian { get; } = CreateRange('\u1800', '\u18af');

    public static Range<char> Limbu { get; } = CreateRange('\u1900', '\u194f');

    public static Range<char> TaiLe { get; } = CreateRange('\u1950', '\u197f');

    public static Range<char> KhmerSymbols { get; } = CreateRange('\u19e0', '\u19ff');

    public static Range<char> PhoneticExtensions { get; } = CreateRange('\u1d00', '\u1d7f');

    public static Range<char> LatinExtendedAdditional { get; } = CreateRange('\u1e00', '\u1eff');

    public static Range<char> GreekExtended { get; } = CreateRange('\u1f00', '\u1fff');

    public static Range<char> GeneralPunctuation { get; } = CreateRange('\u2000', '\u206f');

    public static Range<char> SuperscriptsandSubscripts { get; } = CreateRange('\u2070', '\u209f');

    public static Range<char> CurrencySymbols { get; } = CreateRange('\u20a0', '\u20cf');

    public static Range<char> CombiningDiacriticalMarksforSymbols { get; } = CreateRange('\u20d0', '\u20ff');

    public static Range<char> LetterlikeSymbols { get; } = CreateRange('\u2100', '\u214f');

    public static Range<char> NumberForms { get; } = CreateRange('\u2150', '\u218f');

    public static Range<char> Arrows { get; } = CreateRange('\u2190', '\u21ff');

    public static Range<char> MathematicalOperators { get; } = CreateRange('\u2200', '\u22ff');

    public static Range<char> MiscellaneousTechnical { get; } = CreateRange('\u2300', '\u23ff');

    public static Range<char> ControlPictures { get; } = CreateRange('\u2400', '\u243f');

    public static Range<char> OpticalCharacterRecognition { get; } = CreateRange('\u2440', '\u245f');

    public static Range<char> EnclosedAlphanumerics { get; } = CreateRange('\u2460', '\u24ff');

    public static Range<char> BoxDrawing { get; } = CreateRange('\u2500', '\u257f');

    public static Range<char> BlockElements { get; } = CreateRange('\u2580', '\u259f');

    public static Range<char> GeometricShapes { get; } = CreateRange('\u25a0', '\u25ff');

    public static Range<char> MiscellaneousSymbols { get; } = CreateRange('\u2600', '\u26ff');

    public static Range<char> Dingbats { get; } = CreateRange('\u2700', '\u27bf');

    public static Range<char> MiscellaneousMathematicalSymbolsA { get; } = CreateRange('\u27c0', '\u27ef');

    public static Range<char> SupplementalArrowsA { get; } = CreateRange('\u27f0', '\u27ff');

    public static Range<char> BraillePatterns { get; } = CreateRange('\u2800', '\u28ff');

    public static Range<char> SupplementalArrowsB { get; } = CreateRange('\u2900', '\u297f');

    public static Range<char> MiscellaneousMathematicalSymbolsB { get; } = CreateRange('\u2980', '\u29ff');

    public static Range<char> SupplementalMathematicalOperators { get; } = CreateRange('\u2a00', '\u2aff');

    public static Range<char> MiscellaneousSymbolsandArrows { get; } = CreateRange('\u2b00', '\u2bff');

    public static Range<char> CjkRadicalsSupplement { get; } = CreateRange('\u2e80', '\u2eff');

    public static Range<char> KangxiRadicals { get; } = CreateRange('\u2f00', '\u2fdf');

    public static Range<char> IdeographicDescriptionCharacters { get; } = CreateRange('\u2ff0', '\u2fff');

    public static Range<char> CjkSymbolsandPunctuation { get; } = CreateRange('\u3000', '\u303f');

    public static Range<char> Hiragana { get; } = CreateRange('\u3040', '\u309f');

    public static Range<char> Katakana { get; } = CreateRange('\u30a0', '\u30ff');

    public static Range<char> Bopomofo { get; } = CreateRange('\u3100', '\u312f');

    public static Range<char> HangulCompatibilityJamo { get; } = CreateRange('\u3130', '\u318f');

    public static Range<char> Kanbun { get; } = CreateRange('\u3190', '\u319f');

    public static Range<char> BopomofoExtended { get; } = CreateRange('\u31a0', '\u31bf');

    public static Range<char> KatakanaPhoneticExtensions { get; } = CreateRange('\u31f0', '\u31ff');

    public static Range<char> EnclosedCjkLettersandMonths { get; } = CreateRange('\u3200', '\u32ff');

    public static Range<char> CjkCompatibility { get; } = CreateRange('\u3300', '\u33ff');

    public static Range<char> CjkUnifiedIdeographsExtensionA { get; } = CreateRange('\u3400', '\u4dbf');

    public static Range<char> YijingHexagramSymbols { get; } = CreateRange('\u4dc0', '\u4dff');

    public static Range<char> CjkUnifiedIdeographs { get; } = CreateRange('\u4e00', '\u9fff');

    public static Range<char> YiSyllables { get; } = CreateRange('\ua000', '\ua48f');

    public static Range<char> YiRadicals { get; } = CreateRange('\ua490', '\ua4cf');

    public static Range<char> HangulSyllables { get; } = CreateRange('\uac00', '\ud7af');

    public static Range<char> HighSurrogates { get; } = CreateRange('\ud800', '\udb7f');

    public static Range<char> HighPrivateUseSurrogates { get; } = CreateRange('\udb80', '\udbff');

    public static Range<char> LowSurrogates { get; } = CreateRange('\udc00', '\udfff');

    public static Range<char> PrivateUse { get; } = CreateRange('\ue000', '\uf8ff');

    public static Range<char> PrivateUseArea { get; } = CreateRange('\uf900', '\ufaff');

    public static Range<char> CjkCompatibilityIdeographs { get; } = CreateRange('\ufb00', '\ufb4f');

    public static Range<char> AlphabeticPresentationForms { get; } = CreateRange('\ufb50', '\ufdff');

    public static Range<char> ArabicPresentationFormsA { get; } = CreateRange('\ufe00', '\ufe0f');

    public static Range<char> VariationSelectors { get; } = CreateRange('\ufe20', '\ufe2f');

    public static Range<char> CombiningHalfMarks { get; } = CreateRange('\ufe30', '\ufe4f');

    public static Range<char> CjkCompatibilityForms { get; } = CreateRange('\ufe50', '\ufe6f');

    public static Range<char> SmallFormVariants { get; } = CreateRange('\ufe70', '\ufeff');

    public static Range<char> ArabicPresentationFormsB { get; } = CreateRange('\uff00', '\uffef');

    public static Range<char> HalfwidthandFullwidthForms { get; } = CreateRange('\ufff0', '\uffff');
#pragma warning restore 1591

    /// <summary>
    /// Returns the unicode range containing the specified character.
    /// </summary>
    /// <param name="c">Character to look for</param>
    /// <returns>The unicode range containing the specified character, or null if the character
    /// is not in a unicode range.</returns>
    public static Range<char> GetRange(char c)
    {
        // TODO: Make this efficient. SortedList should do it with a binary search, but it
        // doesn't give us quite what we want
        return allRanges.FirstOrDefault(range => range.Contains(c));
    }
}