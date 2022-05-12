# KINN

KINN is a Novel Alignment-Free Method Based on the Inner Distance Distribution of Kmer Pairs,  
The program is written in python

# Usage

  python kinn_main.py protein --seqs hrv.fasta --medatacsv hrv.csv --k 2 --savefolder $savefolder
  
  ## Parameters   
  
     dna or protein  (must input)  
    --seqs  sequences file of species  
    --savefolder   position of saving the distance  
  
     --k     kmer(optional)  
    --medatacsv Species specific label (optional)  
   
   # Output 
   
      output is the mega style, you can directly pass it into mega X to obtain the phylogenetic tree.
      
   
