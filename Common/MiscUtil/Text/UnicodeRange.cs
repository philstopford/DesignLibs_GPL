using MiscUtil.Collections;
using System.Collections.Generic;

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
    private static readonly List<Range<char>> allRanges = new();

    private static Range<char> CreateRange(char from, char to)
    {
        // TODO: Check for overlaps
        Range<char> ret = new(from, to);
        allRanges.Add(ret);
        return ret;
    }

    private static readonly Range<char> basicLatin = CreateRange('\u0000', '\u007f');
    private static readonly Range<char> latin1Supplement = CreateRange('\u0080', '\u00ff');
    private static readonly Range<char> latinExtendedA = CreateRange('\u0100', '\u017f');
    private static readonly Range<char> latinExtendedB = CreateRange('\u0180', '\u024f');
    private static readonly Range<char> ipaExtensions = CreateRange('\u0250', '\u02af');
    private static readonly Range<char> spacingModifierLetters = CreateRange('\u02b0', '\u02ff');
    private static readonly Range<char> combiningDiacriticalMarks = CreateRange('\u0300', '\u036f');
    private static readonly Range<char> greekAndCoptic = CreateRange('\u0370', '\u03ff');
    private static readonly Range<char> cyrillic = CreateRange('\u0400', '\u04ff');
    private static readonly Range<char> cyrillicSupplement = CreateRange('\u0500', '\u052f');
    private static readonly Range<char> armenian = CreateRange('\u0530', '\u058f');
    private static readonly Range<char> hebrew = CreateRange('\u0590', '\u05FF');
    private static readonly Range<char> arabic = CreateRange('\u0600', '\u06ff');
    private static readonly Range<char> syriac = CreateRange('\u0700', '\u074f');
    private static readonly Range<char> thaana = CreateRange('\u0780', '\u07bf');
    private static readonly Range<char> devangari = CreateRange('\u0900', '\u097f');
    private static readonly Range<char> bengali = CreateRange('\u0980', '\u09ff');
    private static readonly Range<char> gurmukhi = CreateRange('\u0a00', '\u0a7f');
    private static readonly Range<char> gujarati = CreateRange('\u0a80', '\u0aff');
    private static readonly Range<char> oriya = CreateRange('\u0b00', '\u0b7f');
    private static readonly Range<char> tamil = CreateRange('\u0b80', '\u0bff');
    private static readonly Range<char> telugu = CreateRange('\u0c00', '\u0c7f');
    private static readonly Range<char> kannada = CreateRange('\u0c80', '\u0cff');
    private static readonly Range<char> malayalam = CreateRange('\u0d00', '\u0d7f');
    private static readonly Range<char> sinhala = CreateRange('\u0d80', '\u0dff');
    private static readonly Range<char> thai = CreateRange('\u0e00', '\u0e7f');
    private static readonly Range<char> lao = CreateRange('\u0e80', '\u0eff');
    private static readonly Range<char> tibetan = CreateRange('\u0f00', '\u0fff');
    private static readonly Range<char> myanmar = CreateRange('\u1000', '\u109f');
    private static readonly Range<char> georgian = CreateRange('\u10a0', '\u10ff');
    private static readonly Range<char> hangulJamo = CreateRange('\u1100', '\u11ff');
    private static readonly Range<char> ethiopic = CreateRange('\u1200', '\u137f');
    private static readonly Range<char> cherokee = CreateRange('\u13a0', '\u13ff');
    private static readonly Range<char> unifiedCanadianAboriginalSyllabics = CreateRange('\u1400', '\u167f');
    private static readonly Range<char> ogham = CreateRange('\u1680', '\u169f');
    private static readonly Range<char> runic = CreateRange('\u16a0', '\u16ff');
    private static readonly Range<char> tagalog = CreateRange('\u1700', '\u171f');
    private static readonly Range<char> hanunoo = CreateRange('\u1720', '\u173f');
    private static readonly Range<char> buhid = CreateRange('\u1740', '\u175f');
    private static readonly Range<char> tagbanwa = CreateRange('\u1760', '\u177f');
    private static readonly Range<char> khmer = CreateRange('\u1780', '\u17ff');
    private static readonly Range<char> mongolian = CreateRange('\u1800', '\u18af');
    private static readonly Range<char> limbu = CreateRange('\u1900', '\u194f');
    private static readonly Range<char> taiLe = CreateRange('\u1950', '\u197f');
    private static readonly Range<char> khmerSymbols = CreateRange('\u19e0', '\u19ff');
    private static readonly Range<char> phoneticExtensions = CreateRange('\u1d00', '\u1d7f');
    private static readonly Range<char> latinExtendedAdditional = CreateRange('\u1e00', '\u1eff');
    private static readonly Range<char> greekExtended = CreateRange('\u1f00', '\u1fff');
    private static readonly Range<char> generalPunctuation = CreateRange('\u2000', '\u206f');
    private static readonly Range<char> superscriptsandSubscripts = CreateRange('\u2070', '\u209f');
    private static readonly Range<char> currencySymbols = CreateRange('\u20a0', '\u20cf');
    private static readonly Range<char> combiningDiacriticalMarksforSymbols = CreateRange('\u20d0', '\u20ff');
    private static readonly Range<char> letterlikeSymbols = CreateRange('\u2100', '\u214f');
    private static readonly Range<char> numberForms = CreateRange('\u2150', '\u218f');
    private static readonly Range<char> arrows = CreateRange('\u2190', '\u21ff');
    private static readonly Range<char> mathematicalOperators = CreateRange('\u2200', '\u22ff');
    private static readonly Range<char> miscellaneousTechnical = CreateRange('\u2300', '\u23ff');
    private static readonly Range<char> controlPictures = CreateRange('\u2400', '\u243f');
    private static readonly Range<char> opticalCharacterRecognition = CreateRange('\u2440', '\u245f');
    private static readonly Range<char> enclosedAlphanumerics = CreateRange('\u2460', '\u24ff');
    private static readonly Range<char> boxDrawing = CreateRange('\u2500', '\u257f');
    private static readonly Range<char> blockElements = CreateRange('\u2580', '\u259f');
    private static readonly Range<char> geometricShapes = CreateRange('\u25a0', '\u25ff');
    private static readonly Range<char> miscellaneousSymbols = CreateRange('\u2600', '\u26ff');
    private static readonly Range<char> dingbats = CreateRange('\u2700', '\u27bf');
    private static readonly Range<char> miscellaneousMathematicalSymbolsA = CreateRange('\u27c0', '\u27ef');
    private static readonly Range<char> supplementalArrowsA = CreateRange('\u27f0', '\u27ff');
    private static readonly Range<char> braillePatterns = CreateRange('\u2800', '\u28ff');
    private static readonly Range<char> supplementalArrowsB = CreateRange('\u2900', '\u297f');
    private static readonly Range<char> miscellaneousMathematicalSymbolsB = CreateRange('\u2980', '\u29ff');
    private static readonly Range<char> supplementalMathematicalOperators = CreateRange('\u2a00', '\u2aff');
    private static readonly Range<char> miscellaneousSymbolsandArrows = CreateRange('\u2b00', '\u2bff');
    private static readonly Range<char> cjkRadicalsSupplement = CreateRange('\u2e80', '\u2eff');
    private static readonly Range<char> kangxiRadicals = CreateRange('\u2f00', '\u2fdf');
    private static readonly Range<char> ideographicDescriptionCharacters = CreateRange('\u2ff0', '\u2fff');
    private static readonly Range<char> cjkSymbolsandPunctuation = CreateRange('\u3000', '\u303f');
    private static readonly Range<char> hiragana = CreateRange('\u3040', '\u309f');
    private static readonly Range<char> katakana = CreateRange('\u30a0', '\u30ff');
    private static readonly Range<char> bopomofo = CreateRange('\u3100', '\u312f');
    private static readonly Range<char> hangulCompatibilityJamo = CreateRange('\u3130', '\u318f');
    private static readonly Range<char> kanbun = CreateRange('\u3190', '\u319f');
    private static readonly Range<char> bopomofoExtended = CreateRange('\u31a0', '\u31bf');
    private static readonly Range<char> katakanaPhoneticExtensions = CreateRange('\u31f0', '\u31ff');
    private static readonly Range<char> enclosedCjkLettersandMonths = CreateRange('\u3200', '\u32ff');
    private static readonly Range<char> cjkCompatibility = CreateRange('\u3300', '\u33ff');
    private static readonly Range<char> cjkUnifiedIdeographsExtensionA = CreateRange('\u3400', '\u4dbf');
    private static readonly Range<char> yijingHexagramSymbols = CreateRange('\u4dc0', '\u4dff');
    private static readonly Range<char> cjkUnifiedIdeographs = CreateRange('\u4e00', '\u9fff');
    private static readonly Range<char> yiSyllables = CreateRange('\ua000', '\ua48f');
    private static readonly Range<char> yiRadicals = CreateRange('\ua490', '\ua4cf');
    private static readonly Range<char> hangulSyllables = CreateRange('\uac00', '\ud7af');
    private static readonly Range<char> highSurrogates = CreateRange('\ud800', '\udb7f');
    private static readonly Range<char> highPrivateUseSurrogates = CreateRange('\udb80', '\udbff');
    private static readonly Range<char> lowSurrogates = CreateRange('\udc00', '\udfff');
    private static readonly Range<char> privateUse = CreateRange('\ue000', '\uf8ff');
    private static readonly Range<char> privateUseArea = CreateRange('\uf900', '\ufaff');
    private static readonly Range<char> cjkCompatibilityIdeographs = CreateRange('\ufb00', '\ufb4f');
    private static readonly Range<char> alphabeticPresentationForms = CreateRange('\ufb50', '\ufdff');
    private static readonly Range<char> arabicPresentationFormsA = CreateRange('\ufe00', '\ufe0f');
    private static readonly Range<char> variationSelectors = CreateRange('\ufe20', '\ufe2f');
    private static readonly Range<char> combiningHalfMarks = CreateRange('\ufe30', '\ufe4f');
    private static readonly Range<char> cjkCompatibilityForms = CreateRange('\ufe50', '\ufe6f');
    private static readonly Range<char> smallFormVariants = CreateRange('\ufe70', '\ufeff');
    private static readonly Range<char> arabicPresentationFormsB = CreateRange('\uff00', '\uffef');
    private static readonly Range<char> halfwidthandFullwidthForms = CreateRange('\ufff0', '\uffff');

#pragma warning disable 1591
    public static Range<char> BasicLatin => basicLatin;
    public static Range<char> Latin1Supplement => latin1Supplement;
    public static Range<char> LatinExtendedA => latinExtendedA;
    public static Range<char> LatinExtendedB => latinExtendedB;
    public static Range<char> IpaExtensions => ipaExtensions;
    public static Range<char> SpacingModifierLetters => spacingModifierLetters;
    public static Range<char> CombiningDiacriticalMarks => combiningDiacriticalMarks;
    public static Range<char> GreekAndCoptic => greekAndCoptic;
    public static Range<char> Cyrillic => cyrillic;
    public static Range<char> CyrillicSupplement => cyrillicSupplement;
    public static Range<char> Armenian => armenian;
    public static Range<char> Hebrew => hebrew;
    public static Range<char> Arabic => arabic;
    public static Range<char> Syriac => syriac;
    public static Range<char> Thaana => thaana;
    public static Range<char> Devangari => devangari;
    public static Range<char> Bengali => bengali;
    public static Range<char> Gurmukhi => gurmukhi;
    public static Range<char> Gujarati => gujarati;
    public static Range<char> Oriya => oriya;
    public static Range<char> Tamil => tamil;
    public static Range<char> Telugu => telugu;
    public static Range<char> Kannada => kannada;
    public static Range<char> Malayalam => malayalam;
    public static Range<char> Sinhala => sinhala;
    public static Range<char> Thai => thai;
    public static Range<char> Lao => lao;
    public static Range<char> Tibetan => tibetan;
    public static Range<char> Myanmar => myanmar;
    public static Range<char> Georgian => georgian;
    public static Range<char> HangulJamo => hangulJamo;
    public static Range<char> Ethiopic => ethiopic;
    public static Range<char> Cherokee => cherokee;
    public static Range<char> UnifiedCanadianAboriginalSyllabics => unifiedCanadianAboriginalSyllabics;
    public static Range<char> Ogham => ogham;
    public static Range<char> Runic => runic;
    public static Range<char> Tagalog => tagalog;
    public static Range<char> Hanunoo => hanunoo;
    public static Range<char> Buhid => buhid;
    public static Range<char> Tagbanwa => tagbanwa;
    public static Range<char> Khmer => khmer;
    public static Range<char> Mongolian => mongolian;
    public static Range<char> Limbu => limbu;
    public static Range<char> TaiLe => taiLe;
    public static Range<char> KhmerSymbols => khmerSymbols;
    public static Range<char> PhoneticExtensions => phoneticExtensions;
    public static Range<char> LatinExtendedAdditional => latinExtendedAdditional;
    public static Range<char> GreekExtended => greekExtended;
    public static Range<char> GeneralPunctuation => generalPunctuation;
    public static Range<char> SuperscriptsandSubscripts => superscriptsandSubscripts;
    public static Range<char> CurrencySymbols => currencySymbols;
    public static Range<char> CombiningDiacriticalMarksforSymbols => combiningDiacriticalMarksforSymbols;
    public static Range<char> LetterlikeSymbols => letterlikeSymbols;
    public static Range<char> NumberForms => numberForms;
    public static Range<char> Arrows => arrows;
    public static Range<char> MathematicalOperators => mathematicalOperators;
    public static Range<char> MiscellaneousTechnical => miscellaneousTechnical;
    public static Range<char> ControlPictures => controlPictures;
    public static Range<char> OpticalCharacterRecognition => opticalCharacterRecognition;
    public static Range<char> EnclosedAlphanumerics => enclosedAlphanumerics;
    public static Range<char> BoxDrawing => boxDrawing;
    public static Range<char> BlockElements => blockElements;
    public static Range<char> GeometricShapes => geometricShapes;
    public static Range<char> MiscellaneousSymbols => miscellaneousSymbols;
    public static Range<char> Dingbats => dingbats;
    public static Range<char> MiscellaneousMathematicalSymbolsA => miscellaneousMathematicalSymbolsA;
    public static Range<char> SupplementalArrowsA => supplementalArrowsA;
    public static Range<char> BraillePatterns => braillePatterns;
    public static Range<char> SupplementalArrowsB => supplementalArrowsB;
    public static Range<char> MiscellaneousMathematicalSymbolsB => miscellaneousMathematicalSymbolsB;
    public static Range<char> SupplementalMathematicalOperators => supplementalMathematicalOperators;
    public static Range<char> MiscellaneousSymbolsandArrows => miscellaneousSymbolsandArrows;
    public static Range<char> CjkRadicalsSupplement => cjkRadicalsSupplement;
    public static Range<char> KangxiRadicals => kangxiRadicals;
    public static Range<char> IdeographicDescriptionCharacters => ideographicDescriptionCharacters;
    public static Range<char> CjkSymbolsandPunctuation => cjkSymbolsandPunctuation;
    public static Range<char> Hiragana => hiragana;
    public static Range<char> Katakana => katakana;
    public static Range<char> Bopomofo => bopomofo;
    public static Range<char> HangulCompatibilityJamo => hangulCompatibilityJamo;
    public static Range<char> Kanbun => kanbun;
    public static Range<char> BopomofoExtended => bopomofoExtended;
    public static Range<char> KatakanaPhoneticExtensions => katakanaPhoneticExtensions;
    public static Range<char> EnclosedCjkLettersandMonths => enclosedCjkLettersandMonths;
    public static Range<char> CjkCompatibility => cjkCompatibility;
    public static Range<char> CjkUnifiedIdeographsExtensionA => cjkUnifiedIdeographsExtensionA;
    public static Range<char> YijingHexagramSymbols => yijingHexagramSymbols;
    public static Range<char> CjkUnifiedIdeographs => cjkUnifiedIdeographs;
    public static Range<char> YiSyllables => yiSyllables;
    public static Range<char> YiRadicals => yiRadicals;
    public static Range<char> HangulSyllables => hangulSyllables;
    public static Range<char> HighSurrogates => highSurrogates;
    public static Range<char> HighPrivateUseSurrogates => highPrivateUseSurrogates;
    public static Range<char> LowSurrogates => lowSurrogates;
    public static Range<char> PrivateUse => privateUse;
    public static Range<char> PrivateUseArea => privateUseArea;
    public static Range<char> CjkCompatibilityIdeographs => cjkCompatibilityIdeographs;
    public static Range<char> AlphabeticPresentationForms => alphabeticPresentationForms;
    public static Range<char> ArabicPresentationFormsA => arabicPresentationFormsA;
    public static Range<char> VariationSelectors => variationSelectors;
    public static Range<char> CombiningHalfMarks => combiningHalfMarks;
    public static Range<char> CjkCompatibilityForms => cjkCompatibilityForms;
    public static Range<char> SmallFormVariants => smallFormVariants;
    public static Range<char> ArabicPresentationFormsB => arabicPresentationFormsB;
    public static Range<char> HalfwidthandFullwidthForms => halfwidthandFullwidthForms;
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
        foreach (Range<char> range in allRanges)
        {
            if (range.Contains(c))
            {
                return range;
            }
        }
        return null;
    }
}