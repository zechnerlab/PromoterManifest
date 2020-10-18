function [Pre, Post, c, X0] = CreateBirthDeathSystem()


         %A         
  	X0 = [10]';

    Birth = 10;
    Death = 0.01;

    c = [Birth
         Death
         ];

           %A     
    Pre = [ 0
            1
           ];
        
    Post = [1
            0
            ];


end