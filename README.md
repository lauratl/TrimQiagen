# TrimQiagen

## Background

TrimQiagen  assumes that the QIAseq targeted DNA panels library looks like this:

        UMI      R2 Adapter           Sequence to keep                 R1 adapter      
    ┏----------┒┏---------┒┏-------------------------------------┒┏-------------------┒
    UUUUUUUUUUUUAAAAAAAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXAAAAAAAAAAAAAAAAAAAAA
                                     <-------------------------R1----------------------
    -------------------------R2---------------------->                                 

- UMIs are random sequences default length 12 used downstream in the pipeline to remove duplicates. In order to post-process the FASTQs, UMIs must be removed from the read sequence and appended to the read name.  
- Adapters must be removed before mapping to the reference genome, to avoid mismatches
	- R2 adapter is by default the sequence ATTGGAGTCCT, or AGGACTCCAAT when it appears as the reverse complement.
	- R1 adapter is by default the sequence CAAAACGCAATACTGTACATT, or AATGTACAGTATTGCGTTTTG when it appears as the reverse complement.
    

R2 reads would  have the following structure:

        UMI      R2 Adapter           Sequence to keep                 [R1 adapter]
    ┏----------┒┏---------┒┏-------------------------------------┒┏-------------------┒
    UUUUUUUUUUUUAAAAAAAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXAAAAAAAAAAAAAAAAAAAAA

    - UMIs must be present at the beginning of the read or at the position specified by the user
    - R2 adapter must be present at position 13 or at the position specified by the user. It should be removed before mapping to the reference genome to avoid mismatches. 
    - R1 adapter may be partially or totally present at the end of the read

R1 reads would have the following structure:

     [R1 Adapter revComp]        Sequence to keep          [R2 Adapter revComp] [UMI]
    ┏-------------------┒┏-------------------------------------┒┏---------┒┏----------┒
    AAAAAAAAAAAAAAAAAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXAAAAAAAAAAAUUUUUUUUUUUU

    - R1 and R2 adapters can be totally or partially present at the beginning or at the end of the read, respectively. They should be in the reverse complement form 



Apart from extracting the UMIs and removing the adapters, we want to discard those read pairs that don't present the expected structure. Since UMIs are random sequences, we have to rely mainly on the adapters for this. We remove read pairs that:
- Don't have the R2 Adapter at the expected position in R2. But sequencing machines behave in two different ways:
	- Some machines have no problem sequencing the adapter. In this case, we just look for the adapter sequence, optionally allowing for mismatches
	- Other machines have problems reading the same sequencing at the same position across all the reads, leading to a decrease in the base qualities. In this case we cannot search for the sequence itself, but we can detect the decrease in the base quality scores in its expected position.
     If we find the adapter following one of these approaches, the read pair is considered correct and it will be written to the output files.
- Have any of the adapters in the normal form in R1 they should be in the reverse complement form
- Have any of the adapters in the reverse complement form in R2 they should be in the normal form
- Have at least one of the reads R1 or R2 shorter than 70 or the value specified by the user


In the output fastqs, adapters are removed and UMIs are appended to the read name to match downstream analysis requirements. UMI and adapter positions are removed from both the sequence and the base quality string. 
    
Examples of R1 trimming:

    - If we find the R2 adapter reverse complement AGGACTCCAAT at the end, we trim from there:
        XXXXXXXXXXXXXXXXAGGACTCCAAT -> XXXXXXXXXXXXXXXX   # Complete
        XXXXXXXXXXXXXXXXAGGACT      -> XXXXXXXXXXXXXXXX   # Partial start
        XXXXXXXXXXXXXXXXA           -> XXXXXXXXXXXXXXXX   # Partial just the first A

    - If we find the R1 adapter reverse complement AATGTACAGTATTGCGTTTTG  at the beginning, we trim until there 
        CAAAACGCAATACTGTACATTXXXXXX -> XXXXXX             # Complete
        AATACTGTACATTXXXXXXXXXXXXXX -> XXXXXXXXXXXXXX     # Partial end

    - If we find the normal forms of R1 or R2 adapte CAAAACGCAATACTGTACATT or ATTGGAGTCCT, we remove the read pair
        XXXATTGGAGTCCTXXXXXXXXXXXXXXXXXXXXXXXX -> .       # Complete
        XXXCAAAACGCAATACTGTACATTXXXXXXXXXXXXXX -> .       # Complete
        ATTGGAGTCCTXXXXXXXXXXXXXXXXXXXXXXXXXXX -> .       # Complete
        ACTGTACATTXXXXXXXXXXXXXXXXXXXXXXXXXXXX -> .       # Partial
    


## Usage

TrimQiagen works with Python3. Default parameters are our preferred ones for QIAseq gene panels. 


```
$ python src/TrimQiagen.py -h
usage: TrimQiagen.py [-h] [--r1 R1] [--r2 R2] [--o1 O1] [--o2 O2] [--umi_start UMI_START] [--umi_end UMI_END]
                     [--adapter_r2_start ADAPTER_R2_START] [--adapter_r2_end ADAPTER_R2_END]
                     [--adapter_r2_seq ADAPTER_R2_SEQ] [--adapter_r2_mismatches ADAPTER_R2_MISMATCHES]
                     [--adapter_r2_maxlowqual ADAPTER_R2_MAXLOWQUAL]
                     [--adapter_r2_ignorequalpos ADAPTER_R2_IGNOREQUALPOS]
                     [--adapter_r2_highqualsallowed ADAPTER_R2_HIGHQUALSALLOWED] [--adapter_r1_seq ADAPTER_R1_SEQ]
                     [--min_read_length MIN_READ_LENGTH] [--explain] [--verbose] [--test_parameters]
                     [--min_matches_a1 MIN_MATCHES_A1] [--min_matches_a2 MIN_MATCHES_A2]
                     [--min_matches_trimmed_stats MIN_MATCHES_TRIMMED_STATS]

Trim Qiagen reads. Remove R1 and R2 adapters and append UMI to the read name. Reads are only kept if they meet one
of these criteria: a) They have the R2 adapter sequence in the expected position [adapter_r2_start, adapter_r2_end],
allowing adapter_r2_mismatches mismatches. b)They show a decrease in the base qualities at [adapter_r2_start,
adapter_r2_end], being lower than adapter_r2_maxlowqual in all positions but adapter_r2_ignorequalpos

optional arguments:
  -h, --help            show this help message and exit
  --r1 R1               Input file for R1
  --r2 R2               Input file for R2
  --o1 O1               Output file for R1 (Default: R1.fastq.trimmed)
  --o2 O2               Output file for R2 (Default: R2.fastq.trimmed)
  --umi_start UMI_START
                        1-based position for UMI start (Default: 1)
  --umi_end UMI_END     1-based position for UMI end (Default: 12)
  --adapter_r2_start ADAPTER_R2_START
                        1-based position for R2 adapter start (Default: 13)
  --adapter_r2_end ADAPTER_R2_END
                        1-based position for R2 adapter end (Default: adapter_r2_start + len(adapter_r2_seq) -1)
  --adapter_r2_seq ADAPTER_R2_SEQ
                        R2 adapter sequence (Default: ATTGGAGTCCT)
  --adapter_r2_mismatches ADAPTER_R2_MISMATCHES
                        Max number of mismatches allowed for R2 adapter (Default: 1)
  --adapter_r2_maxlowqual ADAPTER_R2_MAXLOWQUAL
                        Maximum base quality score to be considered low quality base in R2 adapter (Default: 25)
  --adapter_r2_ignorequalpos ADAPTER_R2_IGNOREQUALPOS
                        Comma-separated positions to ignore the low quality test in R2 adapter (Default: 18)
  --adapter_r2_highqualsallowed ADAPTER_R2_HIGHQUALSALLOWED
                        Maximum positions allowed to have higher qual than adapter_r2_maxlowqual and still be
                        considered adapter (Default: 2)
  --adapter_r1_seq ADAPTER_R1_SEQ
                        R1 adapter sequence. It appears as reverse complement in R1. (Default:
                        CAAAACGCAATACTGTACATT)
  --min_read_length MIN_READ_LENGTH
                        Min read length to keep the read pair (Default: 70)
  --explain             Explain the tool behaviour
  --verbose             Show the process read by read
  --test_parameters     Show which paremeters will be used
  --min_matches_a1 MIN_MATCHES_A1
                        Minimum length of initial/terminal matches to consider R1 adapter as present in R1 or
                        revComp R1 adapter as present in R2, and thus remove the read pair (Default: 5)
  --min_matches_a2 MIN_MATCHES_A2
                        Minimum length of initial/terminal matches to consider R2 adapter as present in R1 or
                        revComp R2 adapter as present in R2, and thus remove the read pair (Default: 5)
  --min_matches_trimmed_stats MIN_MATCHES_TRIMMED_STATS
                        Minimum length of proper adapters matches when trimming to count it for the stats (Default:
                        4) 

```  
