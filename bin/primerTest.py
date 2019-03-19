from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import primer3

def calculateTm(primer):
    return mt.Tm_NN(Seq(primer), Na=50, Tris=10, Mg=2, dNTPs=0.6,saltcorr=7)

# print(calculateTm('GCTACGGCTTCGACAATTCT'))
# print(primer3.calcTm('GCTACGGCTTCGACAATTCT'))
# print(primer3.calcHairpin('GCTACGGCTTCGACAATTCT'))
# d = {'foo': 100, 'bar': 200}
#
#
#
# print(2224 + len('IKCFCSCSKLPFALSTFGKITKLGLHPDVVTFNTLLHGLCVEDR')*3)
#
# print(Seq('GGTGAAGCTGTATATATCACA').reverse_complement())

def check(seqPr):
    print(calculateTm(seqPr))
    print(primer3.calcHairpin(seqPr))
    print("Rev_com")
    print(Seq(seqPr).reverse_complement())

