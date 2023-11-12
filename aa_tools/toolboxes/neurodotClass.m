% aa toolbox interface for Neurodot toolbox

classdef neurodotClass < toolboxClass
    
    properties (Access = protected)
        hGUI = []% GUI handles
    end
    
    methods
        
        function obj = neurodotClass(path,varargin)
            
            defaultAddToPath = false;
            
            argParse = inputParser;
            argParse.addRequired('path', @ischar);
            argParse.addParameter('name','', @ischar);
            argParse.addParameter('doAddToPath', defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.parse(path,varargin{:});
            
            obj = obj@toolboxClass(argParse.Results.name,argParse.Results.path,argParse.Results.doAddToPath,{});
        end
        
        function load(obj)
            addpath(genpath(obj.toolPath));
            load@toolboxClass(obj)
        end
        
    end
    
end
