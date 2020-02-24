function [Pre, Post, c, X0] = CreateComplexPromoterModel()


         %geneOff geneI geneOn A Dummy
  	X0 = [1       0      0     0 1 ]';

    GeneOn = 0.006;
    GeneOff = 0.01;
    GeneAct = 0.0001;
    GeneDeact = 0.001;
    prodA = 0.01;
    degA = 0.0001;

    
    dummy = 0.1;
    
    c = [
         GeneOn
         GeneOff
         
         prodA
         degA
         
         GeneAct
         GeneDeact

         
         dummy
         ];

           %geneOff geneOn A    dummy
    Pre = [1 0       0     0    0
           0 1       0     0    0
        
           0 1       0     0    0
           0 0       0     1    0
           
           0 1       0     0    0
           0 0       1     0    0
           
           0 0       0     0    1

           ];
        
    Post = [0 1       0     0    0
            1 0       0     0    0
            
            0 1       0     1    0
            0 0       0     0    0
            
            0 0       1     0    0
            1 0       0     0    0
            
            0 0       0     0    1
            ];


end