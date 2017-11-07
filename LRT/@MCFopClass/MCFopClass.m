classdef MCFopClass < handle
    properties
        sens
        ncoils
        adjoint
        imsize
        phase
        dimsize
        NUFFT_op
        trj_length
        nav_dim
        nav_mask
    end
    methods
        function set_phase(obj,phasemap)
            obj.phase=phasemap;
            
        end
        
        function set_MCFop_Params(obj,sens,imsize,dimsize, varargin)
            obj.sens=sens;
            obj.ncoils=size(sens,3);
            obj.adjoint = 0;
            obj.imsize=imsize;
            obj.dimsize=dimsize;
            if(nargin >= 5)
                obj.NUFFT_op = varargin{1};
                obj.trj_length = varargin{2};
                obj.nav_dim = varargin{3};
                obj.nav_mask = varargin{4};
            else
                obj.NUFFT_op = [];
                obj.trj_length = [];
                obj.nav_dim = [];
            end
                
            
        end
        
        function set_MCFop_adjoint(obj,a)
            obj.adjoint = a;
            
        end
        
    end
    
end
