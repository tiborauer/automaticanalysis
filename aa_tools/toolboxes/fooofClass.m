% This requires FOOOF to be installed and added as an aa toolbox.
% Also, FOOOF requires python 3.6-3.11.
%
% The easiest way to install using conda (https://fooof-tools.github.io/fooof/#installation)
%
% add the correspodning entry to your parameterset
%   <toolbox desc='Toolbox with implemented interface in extrafunctions/toolboxes' ui='custom'>
%       <name desc='Name corresponding to the name of the interface without the "Class" suffix' ui='text'>fooof</name>
%           <dir ui='dir'>/users/abcd1234/tools/conda/envs/fooof</dir>
%   </toolbox>
%
classdef fooofClass < toolboxClass
    properties (Access = protected)
        hGUI = []% GUI handles
    end
    
    properties (SetAccess = protected)
        condaEnvironment
        script = 'run_fooof.py' % wrapper script inside fooof_mods
    end

    methods
        function obj = fooofClass(path,varargin)
            defaultAddToPath = false;
            
            argParse = inputParser;
            argParse.addRequired('path',@ischar);
            argParse.addParameter('name','',@ischar);
            argParse.addParameter('doAddToPath',defaultAddToPath,@(x) islogical(x) || isnumeric(x));
            argParse.addParameter('condaEnvironment','fooof',@ischar);
            argParse.parse(path,varargin{:});
            
            obj = obj@toolboxClass(argParse.Results.name,argParse.Results.path,argParse.Results.doAddToPath,{});

            obj.condaEnvironment = argParse.Results.condaEnvironment;
        end
        
        function load(obj)
            load@toolboxClass(obj)

        end
    end
end
