function [Pre, Post, c, X0] = CreateTwoStageExpressionSystemZ()


         %mRNA      P     Z
  	X0 = [1         0     1]';

    Transcription = 0.1;
    mRNADegradation = 0.03;
    Translation = 0.01;
    ProteinDegradation = 0.001;
    
    
    c = [
         Transcription
         mRNADegradation
         Translation
         ProteinDegradation
         
         ];

           %mRNA        P   
    Pre = [ 
            0           0   0
            1           0   0
            
            1           0   1
            0           1   0
           ];
        
    Post = [
            1           0   0
            0           0   0
            
            1           1   1
            0           0   0 
            ];


end