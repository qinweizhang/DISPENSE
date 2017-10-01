classdef MCFopClass < handle
    properties
        sens
        ncoils
        adjoint
        imsize
        phase
        dimsize
    end
    methods
        function set_phase(obj,phasemap)
            obj.phase=phasemap;
            
        end
        
        function set_MCFop_Params(obj,sens,imsize,dimsize)
            obj.sens=sens;
            obj.ncoils=size(sens,3);
            obj.adjoint = 0;
            obj.imsize=imsize;
            obj.dimsize=dimsize;
            
        end
        
        function set_MCFop_adjoint(obj,a)
            obj.adjoint = a;
            
        end
        
    end
    
end
