# TrimQiagen

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
    
   
