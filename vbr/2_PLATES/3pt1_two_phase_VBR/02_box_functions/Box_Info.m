%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defining the Box_Info as a class for fun/practice....
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Box_Info
    properties
        var1name
        var1val
        var1range
        var1units
        var2name
        var2val
        var2range
        var2units
    end
    
    methods        
     function obj = Box_Info(settings,ivar1,ivar2)            
         
       % organize the box
           obj.var1name =  settings.Box.var1name;
           obj.var1val = settings.Box.var1range(ivar1);
           obj.var1range = settings.Box.var1range(:);
           obj.var1units = settings.Box.var1units;
           
           if settings.Box.nvar2 > 1
               obj.var2name = settings.Box.var2name;
               obj.var2val = settings.Box.var2range(ivar2);
               obj.var2range = settings.Box.var2range(:);
               obj.var2units = settings.Box.var2units;
           end
       
     end
    end
end



