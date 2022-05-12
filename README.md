# KINN

KINN is a Novel Alignment-Free Method Based on the Inner Distance Distribution of Kmer Pairs,  
The program is written in python (in Linux)

# Usage  
    
    python kinn_main.py protein --seqs hrv.fasta --medatacsv hrv.csv --k 2 --savefolder $savefolder
    
   ## Help  
     python kinn_main.py -h  
     
  ## Parameters   
  
     dna or protein  (must input)   
     
    --seqs  sequences file of species  
    --savefolder   position of saving the distance  
  
     --k     kmer(optional)  
    --medatacsv Species specific label (optional)   
    --metric    distnace function, default='cosine'
    --seqformat'  default='fasta  
    
   # Output 
   
      Output is the mega style, you can directly pass it into mega X to obtain the phylogenetic tree.
      
   # Example  
      
       python kinn_main.py --seqtype protein --k 2 --seqs example.fasta --medatacsv example.csv --savefolder dis
      
   # Citation
    Runbin Tang, Zuguo Yu, Jinyan Li, Inner distance distributions of k-mer pairs in biological sequences for alignment-free construction of accurate phylogeny trees
   
