function [Pre, Post, c, X0] = CreateMSN2PromoterSystem()


         %geneOff geneI geneOn  mRNA      P  Dummy
  	X0 = [1       0      0       0         0  1]';

    GeneOn = 0.006;
    GeneOff = 0.005;
    LeakyTranscription1 = 0.001;
    LeakyTranscription2 = 0.05;
    Transcription = 0.1;
    mRNADegradation = 0.02;
    Translation = 0.004;
    ProteinDegradation = 0.0004;

    
    c = [
         GeneOn
         GeneOff
         GeneOn
         GeneOff
    
         LeakyTranscription1
         LeakyTranscription2
         Transcription
         mRNADegradation

         Translation
         ProteinDegradation

         0.1
         ];

           %geneOff geneOn  mRNA       dummy
    Pre = [1 0       0       0     0    0
           0 1       0       0     0    0
        
           0 1       0       0     0    0   
           0 0       1       0     0    0    
            
           1 0       0       0     0    0 
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
        
            1 0       0       1     0    0 
            0 1       0       1      0    0
            0 0       1       1      0    0 
            0 0       0       0      0    0
            
            0 0       0       1     1    0
            0 0       0       0     0    0
            

            0 0       0       0     0    1
             
            ];


end