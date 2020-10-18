function [Pre, Post, c, X0] = CreateSystemForFullMoments()


         %geneOff geneI geneOn  mRNA      P  Z
  	X0 = [1       0      0       0        0  1]';

    GeneOn = 0.006;
    GeneOff = 0.005;
    LeakyTranscription = 0.05;
    Transcription = 1;
    mRNADegradation = 0.02;
    Translation = 0.004;
    ProteinDegradation = 0.0004;
    
    
    c = [
         GeneOn
         GeneOff
         GeneOn
         GeneOff
         
    
         LeakyTranscription
         Transcription
         mRNADegradation

         Translation
         ProteinDegradation

         
         ];

    Pre = [1 0       0       0     0    0
           0 1       0       0     0    0
        
           0 1       0       0     0    0   
           0 0       1       0     0    0    
            
           0 1       0       0     0    0  
           0 0       1       0     0    0   
           0 0       0       1     0    0   
                 
           0 0       0       1     0    1
           0 0       0       0     1    0

           ];
        
    Post = [0 1       0       0     0    0
            1 0       0       0     0    0
            
            0 0       1       0      0    0   
            0 1       0       0      0    0   
        
            0 1       0       1      0    0
            0 0       1       1      0    0 
            0 0       0       0      0    0
            
            0 0       0       1     1    1
            0 0       0       0     0    0
            

             
            ];


end