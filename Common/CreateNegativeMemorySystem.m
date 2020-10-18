function [Pre, Post, c, X0] = CreateNegativeMemorySystem()


         %geneOff geneI geneOn gR A I Dummy
  	X0 = [1       0      0     0  0 0 1 ]';

    GeneOn = 0.006;
    GeneOff = 0.005;
    Transcription = 1;
    mRNADegradation = 0.02;
    Translation = 0.004;
    ProteinDegradation = 0.0004;
    
    dummy = 0.1;
    
    c = [
         GeneOn
         GeneOff
         
         GeneOn
         GeneOff
         
         0
         0

         0
         0
         
         dummy
         ];

           %geneOff geneOn gR A    I dummy
    Pre = [1 0       0     0  0    0 0
           0 1       0     0  0    0 0
        
           0 1       0     0  1    0 0   
           0 0       1     0  0    0 0  
           
           1 0       0     0  0    1 0 
           0 0       0     1  0    0 0 
           
           0 1       0     0  0    0 0
           1 0       0     0  1    0 0

            
           0 0       0     0  0    0 1
           ];
        
    Post = [0 1       0     0  0    0 0
            1 0       0     0  0    0 0
            
            0 0       1     0  1    0 0  
            0 1       0     0  0    0 0  
            
            0 0       0     1  0    1 0 
            1 0       0     0  0    0 0

            0 1       0     0  1    0 0
            1 0       0     0  0    1 0
            
            0 0       0     0  0    0 1
             
            ];


end