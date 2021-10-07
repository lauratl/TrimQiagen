#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 09:31:14 2021

@author: laura
"""

import regex
import sys
import argparse



def outputRead( fq1out, fq2out, umi):
    fq1out.write(read_info["r1"]["name"].split(" ")[0] + "_"  + umi + " " + read_info["r1"]["name"].split(" ")[1] + "\n" + read_info["r1"]["seq"] + "\n+\n" + read_info["r1"]["qual"] + "\n")
    fq2out.write(read_info["r2"]["name"].split(" ")[0] + "_"  + umi + " " + read_info["r2"]["name"].split(" ")[1] + "\n" + read_info["r2"]["seq"] + "\n+\n" + read_info["r2"]["qual"] + "\n")

    
def checkQual(qual, max_low_qual, adapter_start):
    # print("Adapter start position: " + str(adapter_start))
    # print("Position to ignore: " + args.adapter_r2_ignorequalpos)
    # print("In the string: " + str([int(x) - adapter_start  for x in args.adapter_r2_ignorequalpos.split(",")]))
    qual_list = [char for char in qual]
    high_quals = 0
    for i in range(len(qual_list)):
        if i in [int(x) - adapter_start  for x in args.adapter_r2_ignorequalpos.split(",")]:
       # if i == 6 : # Skip qual validation in pos 18
            continue
        if ord(qual_list[i]) - 33 > max_low_qual:
            high_quals+=1
    
    if high_quals > args.adapter_r2_highqualsallowed:
        return(False)
    return(True)


def longEnough(seq, minlen):
    if(len(seq)) >= minlen:
        return(True)
    return(False)



def findRequiredAdapter():
    
    adapter_r2_start = args.adapter_r2_start
    adapter_r2_end
    adapter_seq = adapters["a2"]
 
    if regex.search("(" + adapter_seq + "){s<=" +str(args.adapter_r2_mismatches) + "}", read_info["r2"]["seq"][(adapter_r2_start-1):adapter_r2_end]) is not None:
        adapter_counter["r2"]["a2"]["n"]["seq"]+=1
        
    elif checkQual(read_info["r2"]["qual"][(adapter_r2_start-1):adapter_r2_end]  , args.adapter_r2_maxlowqual, adapter_r2_start ):
        adapter_counter["r2"]["a2"]["n"]["qual"]+=1
    else:
        return( False)
    
    umi = read_info["r2"]["seq"][:args.umi_end]
    read_info["r2"]["seq"] = read_info["r2"]["seq"][adapter_r2_end:]
    read_info["r2"]["qual"] = read_info["r2"]["qual"][adapter_r2_end:]
    return(umi)
    
def findAdapter(adapter, seq, minmatches, position):
    
    full_match = adapter
    end_match = adapter[:minmatches] + "".join(["("+nuc for nuc in adapter[minmatches:]]) + "".join([")?" for nuc in adapter[minmatches:]]) + "$"
    start_match = "^" + "".join(["(" for nuc in adapter[:-minmatches]]) + "".join([nuc + ")?" for nuc in adapter[:-minmatches]]) + adapter[-minmatches:] 
    
    if position == "any":
        match = "(" + "|".join([full_match, end_match, start_match]) + ")"
    if position == "end":
        match = "(" + "|".join([full_match, end_match]) + ")"
    if position == "start":
        match = "(" + "|".join([full_match, start_match]) + ")"
        
    return(regex.search(match, seq))
        


def trimRead(read, adapter_name, form, position):
    
    adapter = adapters[adapter_name]
    
    if form == "rc":
        adapter = revComp(adapter)
        
    m = findAdapter(adapter, read_info[read]["seq"], minmatches = 1, position = position)
    
    
    if m is not None:
        
        if args.verbose:
            print("Adapter " + adapter_name + " in the " + form + " form (" + adapter + ") trimmed at the " + position + " of read " + read)
        if m.span()[1] - m.span()[0] >= args.min_matches_trimmed_stats:
            adapter_counter[read][adapter_name][form] += 1
        if position == "end": # If the adapter is at the end, we trim the sequence from the adapter start
                    
            read_info[read]["seq"] = read_info[read]["seq"][:m.span()[0]] 
            read_info[read]["qual"] = read_info[read]["qual"][:m.span()[0]] 
            
        elif position =="start": # if the adapter is at the start, we keep the sequence from the adapter end
                       
            read_info[read]["seq"] = read_info[read]["seq"][m.span()[1]:] 
            read_info[read]["qual"] =read_info[read]["qual"][m.span()[1]:] 
          
        

def revComp(seq):
    
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    return("".join([complement[nuc] for nuc in seq])[::-1])


def findWeirdAdapters(adapter, seq, minmatches):
    
    
    if findAdapter(adapter, seq, minmatches, position = "any") is not None:
        
        return(True)
    

def noWeirdAdapters():
    
    noWeird = True
    for read in ["r1", "r2"]:
        if read =="r1":
            form = "n"
            mm = args.min_matches_a1
        elif read == "r2":
            form = "rc"
            mm = args.min_matches_a2
        for adapter_name in ["a1", "a2"]:
            adapter = adapters[adapter_name]
            if form =="rc":
                adapter = revComp(adapter)
            if findWeirdAdapters(adapter, read_info[read]["seq"], mm):
                adapter_counter[read][adapter_name][form]+=1
                if args.verbose:
                    print("Weird adapter " + adapter_name + " in the " + form + " form found in read " + read)
                noWeird = False
    
    
    return(noWeird)
        

def createParser():
    
    
    parser = argparse.ArgumentParser(description = "Trim Qiagen reads. Remove R1 and R2 adapters and append UMI to the read name. Reads are only kept if they meet one of these criteria:\n a) They have the R2 adapter sequence in the expected position [adapter_r2_start, adapter_r2_end], allowing adapter_r2_mismatches mismatches.\n b)They show a decrease in the base qualities at [adapter_r2_start, adapter_r2_end], being lower than adapter_r2_maxlowqual in all positions but adapter_r2_ignorequalpos ")
    
    parser.add_argument("--r1", help = "Input file for R1")
    parser.add_argument("--r2", help = "Input file for R2")
    parser.add_argument("--o1", help = "Output file for R1 (Default: R1.fastq.trimmed)")
    parser.add_argument("--o2", help = "Output file for R2 (Default: R2.fastq.trimmed)")
    
    parser.add_argument("--umi_start", help = "1-based position for UMI start (Default: 1)", default = 1, type = int)
    parser.add_argument("--umi_end", help = "1-based position for UMI end (Default: 12)",   default = 12, type = int)
    
    parser.add_argument("--adapter_r2_start", help = "1-based position for R2 adapter start (Default: 13)", default = 13, type = int)
    parser.add_argument("--adapter_r2_end", help = "1-based position for R2 adapter end (Default: adapter_r2_start + len(adapter_r2_seq) -1)", type = int)
    parser.add_argument("--adapter_r2_seq", help = "R2 adapter sequence (Default: ATTGGAGTCCT)", default = "ATTGGAGTCCT")
    parser.add_argument("--adapter_r2_mismatches", help = "Max number of mismatches allowed for R2 adapter (Default: 1)", default = 1, type = int)
    parser.add_argument("--adapter_r2_maxlowqual", help = "Maximum base quality score to be considered low quality base in R2 adapter (Default: 25)", default = 25,  type = int)
    parser.add_argument("--adapter_r2_ignorequalpos", help = "Comma-separated positions to ignore the low quality test in R2 adapter (Default: 18)", default = "18")
    parser.add_argument("--adapter_r2_highqualsallowed", help = "Maximum positions allowed to have higher qual than adapter_r2_maxlowqual and still be considered adapter (Default: 2)", type = int,  default = 2)

    parser.add_argument("--adapter_r1_seq", help = "R1 adapter sequence. It appears as reverse complement in R1. (Default: CAAAACGCAATACTGTACATT)", default = "CAAAACGCAATACTGTACATT")
    
    parser.add_argument("--min_read_length", help = "Min read length to keep the read pair (Default: 70)", default = 70,  type = int)
    
    
    parser.add_argument("--explain", help="Explain the tool behaviour", action="store_true")
    parser.add_argument("--verbose", help="Show the process read by read", action="store_true")

    parser.add_argument("--test_parameters", help="Show which paremeters will be used", action="store_true")
    
    parser.add_argument("--min_matches_a1", help="Minimum length of initial/terminal matches to consider R1 adapter as present in R1 or revComp R1 adapter as present in R2, and thus remove the read pair (Default: 5)",  default = 5, type = int)
    parser.add_argument("--min_matches_a2", help="Minimum length of initial/terminal matches to consider R2 adapter as present in R1 or revComp R2 adapter as present in R2, and thus remove the read pair (Default: 5)",  default = 5, type = int)
    parser.add_argument("--min_matches_trimmed_stats", help="Minimum length of proper adapters matches when trimming to count it for the stats (Default: 4)",  default = 4, type = int)

    return(parser)


def showDoc():
    print(sys.argv[0] + " assumes that the QIAseq targeted DNA panels library looks like this:\n")
    print("    UMI      R2 Adapter           Sequence to keep                 R1 adapter")      
    print("┏----------┒┏---------┒┏-------------------------------------┒┏-------------------┒")
    print("UUUUUUUUUUUUAAAAAAAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXAAAAAAAAAAAAAAAAAAAAA\n")
    print("                                 <-------------------------R1----------------------")
    print("-------------------------R2---------------------->                                 ")

    print("- UMIs are random sequences (default length 12) used downstream in the pipeline to remove duplicates. In order to post-process the FASTQs, UMIs must be removed from the read sequence and appended to the read name.  ")
    print("- Adapters must be removed before mapping to the reference genome, to avoid mismatches")

    print("    - R2 adapter is by default the sequence ATTGGAGTCCT, or AGGACTCCAAT when it appears as the reverse complement.")
    print("    - R1 adapter is by default the sequence CAAAACGCAATACTGTACATT, or AATGTACAGTATTGCGTTTTG when it appears as the reverse complement.")
    


    print("\nR2 reads would  have the following structure:\n")
    print("    UMI      R2 Adapter           Sequence to keep                 [R1 adapter]")
    print("┏----------┒┏---------┒┏-------------------------------------┒┏-------------------┒")
    print("UUUUUUUUUUUUAAAAAAAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXAAAAAAAAAAAAAAAAAAAAA\n")
    print("- UMIs must be present at the beginning of the read (or at the position specified by the user)")
    print("- R2 adapter must be present at position 13 (or at the position specified by the user). It should be removed before mapping to the reference genome to avoid mismatches. ")
    print("- R1 adapter may be (partially or totally) present at the end of the read")

    print("\nR1 reads would have the following structure:\n")
    print(" [R1 Adapter revComp]        Sequence to keep          [R2 Adapter revComp] [UMI]")
    print("┏-------------------┒┏-------------------------------------┒┏---------┒┏----------┒")
    print("AAAAAAAAAAAAAAAAAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXAAAAAAAAAAAUUUUUUUUUUUU\n")
    print("- R1 and R2 adapters can be totally or partially present at the beginning or at the end of the read, respectively. They should be in the reverse complement form \n\n")
    print("Apart from extracting the UMIs and removing the adapters, we want to discard those read pairs that don't present the expected structure. Since UMIs are random sequences, we have to rely mainly on the adapters for this. We remove read pairs that:")
    print("- Don't have the R2 Adapter at the expected position in R2. But sequencing machines behave in two different ways: ")
    print("    a) Some machines have no problem sequencing the adapter. In this case, we just look for the adapter sequence, optionally allowing for mismatches")
    print("    b) Other machines have problems reading the same sequencing at the same position across all the reads, leading to a decrease in the base qualities. In this case we cannot search for the sequence itself, but we can detect the decrease in the base quality scores in its expected position.")
    print("  If we find the adapter following one of these approaches, the read pair is considered correct and it will be written to the output files.")
    print("- Have any of the adapters in the normal form in R1 (they should be in the reverse complement form)")
    print("- Have any of the adapters in the reverse complement form in R2 (they should be in the normal form)")
    print("- Have at least one of the reads (R1 or R2) shorter than 70 (or the value specified by the user)")


    print("\n\nIn the output fastqs, adapters are removed and UMIs are appended to the read name to match downstream analysis requirements. UMI and adapter positions are removed from both the sequence and the base quality string. For example, a read like this with "+'\x1b[1;35;40m'+"UMI"+'\x1b[0m'+" and "+'\x1b[0;36;40m'+"adapter"+'\x1b[0m'+":\n")
    print('\x1b[1;30;40m' + "@NB500939:97:HTNM3AFX2:1:11101:24584:1043 2:N:0:25" )
    print('\x1b[1;35;40m' + "CCTATCGTGGAA" + '\x1b[1;36;40m'+ "GACTAGTACTC"+'\x1b[1;30;40m'+"CTAGTTCTAGAGCAGTAACCAGTACAACTATTCCAACTCTGACTATTTCTTCTGATGAACCAGAGACCACAACTTCATTGNTNNNNCATTCTGAGGCAAAGATGATTTCAGCCATTCCAACTTTAG" )
    print("+" )
    print('\x1b[1;35;40m' + "AAAAAEEEEEEE"+'\x1b[0;36;40m' +"///////////"+'\x1b[1;30;40m'+"///E/<6<6E/EE6EA//EE6E/////////</<//A/E/E/AAA/AAA//<AE<<E//<E/EEE/E//<//<//</A/E#<####<///A/A/EEE///E/<E//<6A/A<////6///AE/E/E" + '\x1b[0m')

    print("\nwill be outputed like this:\n")

    print('\x1b[1;30;40m' + "@NB500939:97:HTNM3AFX2:1:11101:24584:1043_"+ '\x1b[1;35;40m'+ "CCTATCGTGGAA"+'\x1b[1;30;40m' +" 2:N:0:25" )
    print( "CTAGTTCTAGAGCAGTAACCAGTACAACTATTCCAACTCTGACTATTTCTTCTGATGAACCAGAGACCACAACTTCATTGNTNNNNCATTCTGAGGCAAAGATGATTTCAGCCATTCCAACTTTAG" )
    print("+" )
    print("///E/<6<6E/EE6EA//EE6E/////////</<//A/E/E/AAA/AAA//<AE<<E//<E/EEE/E//<//<//</A/E#<####<///A/A/EEE///E/<E//<6A/A<////6///AE/E/E" + '\x1b[0m' )
    print("\n\n")
   
    print("Examples of R1 trimming:\n")
    print("- If we find the R2 adapter reverse complement (AGGACTCCAAT) at the end, we trim from there:")
    print("    XXXXXXXXXXXXXXXXAGGACTCCAAT -> XXXXXXXXXXXXXXXX   # Complete")
    print("    XXXXXXXXXXXXXXXXAGGACT      -> XXXXXXXXXXXXXXXX   # Partial (start)")
    print("    XXXXXXXXXXXXXXXXA           -> XXXXXXXXXXXXXXXX   # Partial (just the first A)\n")
    print("- If we find the R1 adapter reverse complement (AATGTACAGTATTGCGTTTTG ) at the beginning, we trim until there ")
    print("    CAAAACGCAATACTGTACATTXXXXXX -> XXXXXX             # Complete")
    print("    AATACTGTACATTXXXXXXXXXXXXXX -> XXXXXXXXXXXXXX     # Partial (end)\n")
    print("- If we find the normal forms of R1 or R2 adapte (CAAAACGCAATACTGTACATT or ATTGGAGTCCT), we remove the read pair")
    print("    XXXATTGGAGTCCTXXXXXXXXXXXXXXXXXXXXXXXX -> .       # Complete")
    print("    XXXCAAAACGCAATACTGTACATTXXXXXXXXXXXXXX -> .       # Complete")
    print("    ATTGGAGTCCTXXXXXXXXXXXXXXXXXXXXXXXXXXX -> .       # Complete")
    print("    ACTGTACATTXXXXXXXXXXXXXXXXXXXXXXXXXXXX -> .       # Partial")
    print("")
   
    
    
parser = createParser()
args = parser.parse_args()  

if args.explain:
    
    showDoc()
    sys.exit()


if not args.r1 or not args.r2:
    print("Usage: " + sys.argv[0] + " --r1 <fastq_R1> --r2 <fastq_R2> [--o1 <output_fastq_R1> --o2 <output_fastq_R2>] ")
    sys.exit()

if args.o1:
    o1 = args.o1
else:
    o1 = args.r1 + ".trimmed"
    
if args.o2:
    o2 = args.o2
else:
    o2 = args.r2 + ".trimmed"
    
if args.adapter_r2_end:
    if 1 + args.adapter_r2_end - args.adapter_r2_start != len(args.adapter_r2_seq):
        print("R2 adapter start and end positions don't match adapter sequence length")
        sys.exit()
    else:
        adapter_r2_end = args.adapter_r2_end
else:
    adapter_r2_end = args.adapter_r2_start + len(args.adapter_r2_seq) -1


print("\nTrimming " + args.r1 + " and " + args.r2 + " into " + o1 + " and " + o2 )
print("Extracting UMI from position " + str(args.umi_start) +  " to position " + str(args.umi_end))
print("Removing R2 adapter " + args.adapter_r2_seq + " from position " + str(args.adapter_r2_start) +  " to position " + str(adapter_r2_end) + " allowing up to " +str(args.adapter_r2_mismatches)+ " mismatches")
print("If the R2 adapter sequence is not found, we will look for a decrease in base quality scores between " +str(args.adapter_r2_start) + " and " + str(adapter_r2_end) +". Quality must be below "  + str(args.adapter_r2_maxlowqual) + " except for position(s) " + args.adapter_r2_ignorequalpos)

if args.test_parameters:
    sys.exit()

fq1 = open(args.r1, 'r')
fq2 = open(args.r2, 'r')
fq1out = open(o1, 'w')
fq2out = open(o2, 'w')


line_counter = 0
line_counter = 0



read_info = {"r1":{}, "r2":{}}

# r1_info = {}
# r2_info = {}

adapter_counter = {"r1" : {"a1" : {"n": 0, "rc": 0}, "a2" : {"n": 0, "rc": 0}},
           "r2" : {"a1" : {"n": 0, "rc": 0}, "a2" : {"n" : {"seq": 0, "qual":0}, "rc": 0}}}

stats = {"no_umi":0, "weird_adapter": 0, "short_r1": 0, "short_r2":0, "ok": 0}

adapters = {"a1": args.adapter_r1_seq, "a2": args.adapter_r2_seq}

for r1 in fq1:
    if r1[-1]=="\n":
        r1 = r1[:-1]
    if line_counter == 0:
        read_info["r1"]["name"] = r1
    elif line_counter == 1:
        read_info["r1"]["seq"] = r1
    # elif line_counter == 2:
    #     r1_info["plus"] = r1
    elif line_counter == 3:
        read_info["r1"]["qual"] = r1
        
        
    
    for r2 in fq2:
        if r2[-1]=="\n":
            r2 = r2[:-1]
        if line_counter == 0:
            read_info["r2"]["name"] = r2
        elif line_counter == 1:
            read_info["r2"]["seq"] = r2
        # elif line_counter == 2:
        #     r2_info["plus"] = r2
        elif line_counter == 3:
            read_info["r2"]["qual"] = r2
        
        
        break
    
    line_counter += 1
    
    
    if line_counter == 4:
        
        if args.verbose:
            print("R1:")
            print(read_info["r1"])
            print("R2:")
            print(read_info["r2"])
        
        
        
        # Check for the presence of weird adapters
        
        weirdAdaptersNotFound = noWeirdAdapters()  # Detect here to avoid small trimming or proper adapters to interfer with a longer matching of wrong adapters. Don't continue the loop here to keep the whole count of adapters
        
        # Find required adapter in R2
        
        umi = findRequiredAdapter()
    
        
         
        # Trim R1 end
        
        trimRead(read = "r1", adapter_name = "a2", form = "rc" , position = "end")

        # Trim R1 start
        
        trimRead(read = "r1", adapter_name = "a1", form = "rc" , position = "start")
        
        # Trim R2 end
        
        trimRead(read = "r2", adapter_name = "a1", form = "n" , position = "end")
        
        


        if not umi:
            if args.verbose:
                print("R2 adapter is not in the right place -> Discard read pair\n")
            stats["no_umi"]+=1
        elif not weirdAdaptersNotFound:
            if args.verbose:
                print("There were some unexpected adapters -> Discard read pair\n")
            stats["weird_adapter"]+=1
        elif not longEnough(read_info["r1"]["seq"], args.min_read_length):
            if args.verbose:
                print("R1 is too short -> Discard read pair\n")
            stats["short_r1"]+=1
        elif not longEnough(read_info["r2"]["seq"], args.min_read_length):
            if args.verbose:
                print("R2 is too short -> Discard read pair\n")         
            stats["short_r2"]+=1
       
        else:
            if args.verbose:
                print("Read pair should be in the output\n")
                print("UMI is " + umi + "\n")
            stats["ok"]+=1
            outputRead( fq1out, fq2out, umi)
        
        
           
        
   

        line_counter = 0
        
fq1.close()
fq2.close()
fq1out.close()
fq2out.close()

# print(adapter_counter)
print("\n\n*****************\nSummary:\n")
print(stats)

print("\nRead\tA1 normal\tA1 rc\tA2 normal\tA2 rc ")
for read in adapter_counter.keys():
    print(f'{read}\t{adapter_counter[read]["a1"]["n"]}\t{adapter_counter[read]["a1"]["rc"]}\t{adapter_counter[read]["a2"]["n"]}\t{adapter_counter[read]["a2"]["rc"]}')