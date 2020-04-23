function [Pre, Post, c, X0] = CreateThreeStateExpressionSystemLeaky()


         %geneOff geneI geneOn  mRNA      P  Dummy
  	X0 = [1       0      0       0        0  1]';

    GeneOn = 0.006;
    GeneOff = 0.005;
    LeakyTranscription = 0.05;
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
         
    
         LeakyTranscription
         Transcription
         0*mRNADegradation

         Translation
         ProteinDegradation

         dummy
         ];

    Pre = [1 0       0       0     0    0
           0 1       0       0     0    0
        
           0 1       0       0     0    0   
           0 0       1       0     0    0    
            
           0 1       0       0     0    0  
           0 0       1       0     0    0   
           0 0       0       1     0    0   
                 
           0 0       0       1     0    0
           0 0       0       0     1    0
            
           0 0       0       0     0    1
           ];
        
    Post = [0 1       0       0     0    0
            1 0       0       0     0    0
            
            0 0       1       0      0    0   
            0 1       0       0      0    0   
        
            0 1       0       1      0    0
            0 0       1       1      0    0 
            0 0       0       0      0    0
            
            0 0       0       1     1    0
            0 0       0       0     0    0
            
            
            0 0       0       0     0    1
             
            ];


end