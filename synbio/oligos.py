"""Functions for oligos. Tm calc."""

import math

from Bio.Seq import Seq


DNA_NN = {
    "init": (0.2, -5.7),
    "init_A/T": (2.2, 6.9),
    "init_G/C": (0, 0),
    "init_oneG/C": (0, 0),
    "init_allA/T": (0, 0),
    "init_5T/A": (0, 0),
    "sym": (0, -1.4),
    "AA/TT": (-7.6, -21.3),
    "AT/TA": (-7.2, -20.4),
    "TA/AT": (-7.2, -20.4),
    "CA/GT": (-8.5, -22.7),
    "GT/CA": (-8.4, -22.4),
    "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2),
    "CG/GC": (-10.6, -27.2),
    "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.0),
}
"""
SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440

from Bio.SeqUtils
"""

DNA_IMM = {
    "AG/TT": (1.0, 0.9),
    "AT/TG": (-2.5, -8.3),
    "CG/GT": (-4.1, -11.7),
    "CT/GG": (-2.8, -8.0),
    "GG/CT": (3.3, 10.4),
    "GG/TT": (5.8, 16.3),
    "GT/CG": (-4.4, -12.3),
    "GT/TG": (4.1, 9.5),
    "TG/AT": (-0.1, -1.7),
    "TG/GT": (-1.4, -6.2),
    "TT/AG": (-1.3, -5.3),
    "AA/TG": (-0.6, -2.3),
    "AG/TA": (-0.7, -2.3),
    "CA/GG": (-0.7, -2.3),
    "CG/GA": (-4.0, -13.2),
    "GA/CG": (-0.6, -1.0),
    "GG/CA": (0.5, 3.2),
    "TA/AG": (0.7, 0.7),
    "TG/AA": (3.0, 7.4),
    "AC/TT": (0.7, 0.2),
    "AT/TC": (-1.2, -6.2),
    "CC/GT": (-0.8, -4.5),
    "CT/GC": (-1.5, -6.1),
    "GC/CT": (2.3, 5.4),
    "GT/CC": (5.2, 13.5),
    "TC/AT": (1.2, 0.7),
    "TT/AC": (1.0, 0.7),
    "AA/TC": (2.3, 4.6),
    "AC/TA": (5.3, 14.6),
    "CA/GC": (1.9, 3.7),
    "CC/GA": (0.6, -0.6),
    "GA/CC": (5.2, 14.2),
    "GC/CA": (-0.7, -3.8),
    "TA/AC": (3.4, 8.0),
    "TC/AA": (7.6, 20.2),
    "AA/TA": (1.2, 1.7),
    "CA/GA": (-0.9, -4.2),
    "GA/CA": (-2.9, -9.8),
    "TA/AA": (4.7, 12.9),
    "AC/TC": (0.0, -4.4),
    "CC/GC": (-1.5, -7.2),
    "GC/CC": (3.6, 8.9),
    "TC/AC": (6.1, 16.4),
    "AG/TG": (-3.1, -9.5),
    "CG/GG": (-4.9, -15.3),
    "GG/CG": (-6.0, -15.8),
    "TG/AG": (1.6, 3.6),
    "AT/TT": (-2.7, -10.8),
    "CT/GT": (-5.0, -15.8),
    "GT/CT": (-2.2, -8.4),
    "TT/AT": (0.2, -1.5),
    "AI/TC": (-8.9, -25.5),
    "TI/AC": (-5.9, -17.4),
    "AC/TI": (-8.8, -25.4),
    "TC/AI": (-4.9, -13.9),
    "CI/GC": (-5.4, -13.7),
    "GI/CC": (-6.8, -19.1),
    "CC/GI": (-8.3, -23.8),
    "GC/CI": (-5.0, -12.6),
    "AI/TA": (-8.3, -25.0),
    "TI/AA": (-3.4, -11.2),
    "AA/TI": (-0.7, -2.6),
    "TA/AI": (-1.3, -4.6),
    "CI/GA": (2.6, 8.9),
    "GI/CA": (-7.8, -21.1),
    "CA/GI": (-7.0, -20.0),
    "GA/CI": (-7.6, -20.2),
    "AI/TT": (0.49, -0.7),
    "TI/AT": (-6.5, -22.0),
    "AT/TI": (-5.6, -18.7),
    "TT/AI": (-0.8, -4.3),
    "CI/GT": (-1.0, -2.4),
    "GI/CT": (-3.5, -10.6),
    "CT/GI": (0.1, -1.0),
    "GT/CI": (-4.3, -12.1),
    "AI/TG": (-4.9, -15.8),
    "TI/AG": (-1.9, -8.5),
    "AG/TI": (0.1, -1.8),
    "TG/AI": (1.0, 1.0),
    "CI/GG": (7.1, 21.3),
    "GI/CG": (-1.1, -3.2),
    "CG/GI": (5.8, 16.9),
    "GG/CI": (-7.6, -22.0),
    "AI/TI": (-3.3, -11.9),
    "TI/AI": (0.1, -2.3),
    "CI/GI": (1.3, 3.0),
    "GI/CI": (-0.5, -1.3),
}
"""
Internal mismatch and inosine table (DNA)
Allawi & SantaLucia (1997), Biochemistry 36: 10581-10594
Allawi & SantaLucia (1998), Biochemistry 37: 9435-9444
Allawi & SantaLucia (1998), Biochemistry 37: 2170-2179
Allawi & SantaLucia (1998), Nucl Acids Res 26: 2694-2701
Peyret et al. (1999), Biochemistry 38: 3468-3477
Watkins & SantaLucia (2005), Nucl Acids Res 33: 6258-6267

from Bio.SeqUtils
"""

DNA_TMM = {
    "AA/TA": (-3.1, -7.8),
    "TA/AA": (-2.5, -6.3),
    "CA/GA": (-4.3, -10.7),
    "GA/CA": (-8.0, -22.5),
    "AC/TC": (-0.1, 0.5),
    "TC/AC": (-0.7, -1.3),
    "CC/GC": (-2.1, -5.1),
    "GC/CC": (-3.9, -10.6),
    "AG/TG": (-1.1, -2.1),
    "TG/AG": (-1.1, -2.7),
    "CG/GG": (-3.8, -9.5),
    "GG/CG": (-0.7, -19.2),
    "AT/TT": (-2.4, -6.5),
    "TT/AT": (-3.2, -8.9),
    "CT/GT": (-6.1, -16.9),
    "GT/CT": (-7.4, -21.2),
    "AA/TC": (-1.6, -4.0),
    "AC/TA": (-1.8, -3.8),
    "CA/GC": (-2.6, -5.9),
    "CC/GA": (-2.7, -6.0),
    "GA/CC": (-5.0, -13.8),
    "GC/CA": (-3.2, -7.1),
    "TA/AC": (-2.3, -5.9),
    "TC/AA": (-2.7, -7.0),
    "AC/TT": (-0.9, -1.7),
    "AT/TC": (-2.3, -6.3),
    "CC/GT": (-3.2, -8.0),
    "CT/GC": (-3.9, -10.6),
    "GC/CT": (-4.9, -13.5),
    "GT/CC": (-3.0, -7.8),
    "TC/AT": (-2.5, -6.3),
    "TT/AC": (-0.7, -1.2),
    "AA/TG": (-1.9, -4.4),
    "AG/TA": (-2.5, -5.9),
    "CA/GG": (-3.9, -9.6),
    "CG/GA": (-6.0, -15.5),
    "GA/CG": (-4.3, -11.1),
    "GG/CA": (-4.6, -11.4),
    "TA/AG": (-2.0, -4.7),
    "TG/AA": (-2.4, -5.8),
    "AG/TT": (-3.2, -8.7),
    "AT/TG": (-3.5, -9.4),
    "CG/GT": (-3.8, -9.0),
    "CT/GG": (-6.6, -18.7),
    "GG/CT": (-5.7, -15.9),
    "GT/CG": (-5.9, -16.1),
    "TG/AT": (-3.9, -10.5),
    "TT/AG": (-3.6, -9.8),
}
"""
Terminal mismatch table (DNA)
SantaLucia & Peyret (2001) Patent Application WO 01/94611

from Bio.SeqUtils
"""


def calc_tm(seq1: str, seq2: str = "", pcr: bool = True) -> float:
    """Calculate the annealing temperature between seq1 and seq2.

    If seq2 is not provided, it's exact complement is used.
    In otherwords, it's assumed to be an exact match. This tm
    calculate does not account for pseudoknots or anything other
    than an exact, unpadded alignment between seq1 and seq2.

    This is largely a copy-paste from Bio.SeqUtils.MeltingTemp with
    some different defaults. Here, the reaction mixture is assumed to
    be PCR and concentrations for Mg, Tris, K, and dNTPs are included
    that match a typical PCR reaction according to Thermo and NEB. Additionally,
    the salt correction formula from IDT's Owczarzy et al. (2008) is used.

    NEB: https://www.neb.com/tools-and-resources/usage-guidelines/guidelines-for-pcr-optimization-with-taq-dna-polymerase
    ThermoFisher: https://www.thermofisher.com/order/catalog/product/18067017?SID=srch-srp-18067017

    NOTE: Sequences are assumed not to be symmetrical. Oligo not binding to self.

    Args:
        seq1: The seq whose tm is calculated

    Keyword Args:
        seq2: The seq that seq1 anneals to in 3' -> 5' direction
        pcr: Whether tm is being calculated for the oligo is in a
            PCR reaction mixture. If so, ion and Tris concentrations
            that match a typical NEB/ThermoFisher PCR mixture are used

    Returns:
        The estimated tm as a float
    """

    if not seq2:
        seq2 = str(Seq(seq1).complement())

    if len(seq1) != len(seq2):
        raise ValueError(
            f"length mismatch between seq1 {len(seq1)} and seq2 {len(seq2)}"
        )

    # sum enthalpy and entropy. Enthaply is first value of each tuple and
    # entropy is the second value of each tuple in:
    # SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440

    # start with initiation enthalpy and entropy
    delta_h, delta_s = DNA_NN["init"]

    # add in initial A/T and initial G/Cs
    init = seq1[0] + seq1[-1]
    init_at = init.count("A") + init.count("T")
    init_gc = init.count("G") + init.count("C")
    init_at_h, init_at_s = DNA_NN["init_A/T"]
    init_gc_h, init_gc_s = DNA_NN["init_G/C"]
    delta_h += init_at * init_at_h + init_gc * init_gc_h
    delta_s += init_at * init_at_s + init_gc * init_gc_s

    # work through each nearest neighbor pair
    for i in range(len(seq1) - 1):
        pair1 = seq1[i] + seq1[i + 1]
        pair2 = seq2[i] + seq2[i + 1]
        pair = pair1 + "/" + pair2
        pair_rev = pair[::-1]

        # assuming internal neighbor pair
        pair_delta_h, pair_delta_s = 0.0, 0.0
        if pair in DNA_NN:
            pair_delta_h, pair_delta_s = DNA_NN[pair]
        elif pair_rev in DNA_NN:
            pair_delta_h, pair_delta_s = DNA_NN[pair_rev]
        elif pair in DNA_IMM:
            pair_delta_h, pair_delta_s = DNA_IMM[pair]
        elif pair_rev in DNA_IMM:
            pair_delta_h, pair_delta_s = DNA_IMM[pair_rev]

        # overwrite if it's a terminal pair
        if i in (0, len(seq1) - 2):
            if pair in DNA_TMM:
                pair_delta_h, pair_delta_s = DNA_TMM[pair]
            elif pair_rev in DNA_TMM:
                pair_delta_h, pair_delta_s = DNA_TMM[pair_rev]

        delta_h += pair_delta_h
        delta_s += pair_delta_s

    # adjust salt based on mode
    if pcr:
        seq1_conc = 250.0
        seq2_conc = 0.0
        Na = 0
        K = 50
        Tris = 2  # see Thermo
        Mg = 1.5  # see NEB
        dNTPs = 0.2  # see NEB
    else:
        seq1_conc = 25.0
        seq2_conc = 25.0
        Na = 50
        K = 0
        Tris = 0
        Mg = 0
        dNTPs = 0

    # salt correction for deltaS
    # copied-pasted from Bio.SeqUtils' use of a decision tree by:
    # Owczarzy et al. (2008), Biochemistry 4 7: 5336-5353
    Mon = Na + K + Tris / 2.0  # monovalent ions
    mg = Mg * 1e-3  # Lowercase ions (mg, mon, dntps) are molar
    mon = Mon * 1e-3

    # coefficients to a multi-variate from the paper
    a, b, c, d = 3.92, -0.911, 6.26, 1.42
    e, f, g = -48.2, 52.5, 8.31

    if dNTPs > 0:
        dntps = dNTPs * 1e-3
        ka = 3e4  # Dissociation constant for Mg:dNTP
        # Free Mg2+ calculation:
        mg = (
            -(ka * dntps - ka * mg + 1.0)
            + math.sqrt((ka * dntps - ka * mg + 1.0) ** 2 + 4.0 * ka * mg)
        ) / (2.0 * ka)
    if Mon > 0:
        R = math.sqrt(mg) / mon
        if R < 0.22:
            corr = (4.29 * _gc(seq1) / 100 - 3.95) * 1e-5 * math.log(
                mon
            ) + 9.40e-6 * math.log(mon) ** 2
            return corr
        elif R < 6.0:
            a = 3.92 * (0.843 - 0.352 * math.sqrt(mon) * math.log(mon))
            d = 1.42 * (1.279 - 4.03e-3 * math.log(mon) - 8.03e-3 * math.log(mon) ** 2)
            g = 8.31 * (0.486 - 0.258 * math.log(mon) + 5.25e-3 * math.log(mon) ** 3)
    corr = (
        a
        + b * math.log(mg)
        + (_gc(seq1) / 100) * (c + d * math.log(mg))
        + (1 / (2.0 * (len(seq1) - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
    ) * 1e-5

    # tm with concentration consideration
    k = (seq1_conc - (seq2_conc / 2.0)) * 1e-9
    R = 1.9872
    tm = (delta_h * 1000.0) / (delta_s + R * math.log(k)) - 273.15

    # add in salt correction
    tm = 1 / (1 / (tm + 273.15) + corr) - 273.15

    return tm


def _gc(seq: str) -> float:
    """Return the GC ratio of a sequence."""

    seq = seq.upper()

    return float(seq.count("G") + seq.count("C")) / float(len(seq))

