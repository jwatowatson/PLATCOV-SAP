import argparse
import pdb
from pango_aliasor.aliasor import Aliasor
import pandas as pd
import sys

def Get_Args():
    parser = argparse.ArgumentParser(description="Convert a PANGO lineage into PLATCOV classification")
    parser.add_argument("--input", help="Input file - one variant per line without headers", required=True)
    parser.add_argument("--output", help="Output file", required=True)
    args = parser.parse_args()

    return args
    
def Classify_Lineage(original):
    #current classes of variants:
    #B.1.1.7 is Alpha
    #B.1.351 is Beta
    #B.1.617.2* is Delta (inc all AY lineages)
    #P.1 is Gamma
    #BA.1 is BA.1.*
    #BA.2 is BA.2* EXCLUDING BA.2.75*
    #BA.4 is BA.4
    #BA.5 is BA.5
    #BQ.1* is BQ.1
    '''Take lineage names, expand, assign to WHO lineage'''
    try:
        uncompressed = aliasor.uncompress(original)
    except AttributeError: #catch 'nan' values which cannot be parsed
        return "none"
    if uncompressed.startswith('B.1.1.7'):
        lineage = "Alpha"
    elif uncompressed.startswith('B.1.351'):
        lineage = "Beta"
    elif uncompressed.startswith('B.1.617.2'):
        lineage = "Delta"
    elif uncompressed.startswith('P.1'):
        lineage = "Gamma"
    #ordering is important - here we extract the various BA.2 and BA.5 sublineages first
    #unwieldy but the best way to pull out all the BQ.1 descendants
    elif uncompressed.startswith('B.1.1.529.5.3.1.1.1.1.1'):
        lineage = "BQ.1"
    #also includes intra-Omicron recombinants with the same spike as BA.2.75    
    elif uncompressed.startswith("B.1.1.529.2.75") or uncompressed.startswith('XBF') or uncompressed.startswith('XBK'):
        lineage = "BA.2.75"
    elif uncompressed.startswith("B.1.1.529.2"):
        lineage = "BA.2"
    elif uncompressed.startswith("B.1.1.529.5"):
        lineage = "BA.5"
    elif uncompressed.startswith("B.1.1.529.1"):
        lineage = "BA.1"
    elif uncompressed.startswith("B.1.1.529.4"):
        lineage = "BA.4"
    elif uncompressed.startswith("XBB.1.5") or uncompressed.startswith("XBB.1.9") or uncompressed.startswith("XBB.1.16") or uncompressed.startswith("XBL") or uncompressed.startswith("XBB.2.3"):
        lineage = "XBB.1.5-like"
    elif uncompressed.startswith('XBB'):
        lineage = 'XBB'
     
    else:
        lineage = 'Other'
    return lineage
    
if __name__=="__main__":
    args = Get_Args()
    #relies on PANGO Aliasor: https://github.com/corneliusroemer/pango_aliasor  
    aliasor = Aliasor()
    
    #try:
    inputdata = pd.read_csv(args.input, names=['Original'])
    #except:
    #    print("Error: input file {0} does not exists".format(args.input))
    #    sys.exit(1)
    
    #add a column with the VariantClass
    inputdata['VariantClass'] = inputdata.apply(lambda x: "{0}".format(Classify_Lineage(x.Original)), axis=1)

    #write CSV output
    inputdata.to_csv(args.output, index=False)    
