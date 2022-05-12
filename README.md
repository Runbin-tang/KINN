# KINN

KINN is a Novel Alignment-Free Method Based on the Inner Distance Distribution of Kmer Pairs,  
The program is written in python (in Linux)

# Usage  

    python kinn_main.py protein --seqs hrv.fasta --medatacsv hrv.csv --k 2 --savefolder $savefolder
  
  ## Parameters   
  
     dna or protein  (must input)   
     
    --seqs  sequences file of species  
    --savefolder   position of saving the distance  
  
     --k     kmer(optional)  
    --medatacsv Species specific label (optional)   
    --metric    distnace function, default='cosine'
    --seqformat'  default='fasta  
    
   # Output 
   
      output is the mega style, you can directly pass it into mega X to obtain the phylogenetic tree.
      
   # Example  
      
       python kinn_main.py --seqtype protein --k 2 --seqs hrv.fasta --medatacsv hrv.csv --savefolder dis
      
      
   
