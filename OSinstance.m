% Matlab class for Optimization Services instance
classdef OSinstance < handle
    properties
        document
        instanceHeader
        instanceData
    end
    methods
        function instance = OSinstance
            % default constructor
            instance.document = ...
                com.mathworks.xml.XMLUtils.createDocument('osil');
            docRootNode = instance.document.getDocumentElement;
            docRootNode.setAttribute('xmlns', ...
                'os.optimizationservices.org');
            docRootNode.setAttribute('xmlns:xsi', ...
                'http://www.w3.org/2001/XMLSchema-instance');
            docRootNode.setAttribute('xsi:schemaLocation', ...
                ['os.optimizationservices.org ' ...
                'http://www.optimizationservices.org/schemas/2.0/OSiL.xsd']);
            
            instance.instanceHeader.element = ...
                instance.document.createElement('instanceHeader');
            instance.instanceData.element = ...
                instance.document.createElement('instanceData');
            
            docRootNode.appendChild(instance.instanceHeader.element);
            docRootNode.appendChild(instance.instanceData.element);
            
            variables = instance.document.createElement('variables');
            variables.setAttribute('numberOfVariables', '0');
            instance.instanceData.variables.element = variables;
            instance.instanceData.variables.numberOfVariables = 0;
            instance.instanceData.element.appendChild(variables);
            
            % only implementing single-objective problems for now
            objectives = instance.document.createElement('objectives');
            obj = instance.document.createElement('obj');
            obj.setAttribute('maxOrMin', 'min');
            obj.setAttribute('numberOfObjCoef', '0');
            objectives.appendChild(obj);
            instance.instanceData.objectives.element = objectives;
            instance.instanceData.objectives.obj.element = obj;
            instance.instanceData.element.appendChild(objectives);
            
            constraints = instance.document.createElement('constraints');
            constraints.setAttribute('numberOfConstraints', '0');
            instance.instanceData.constraints.element = constraints;
            instance.instanceData.constraints.numberOfConstraints = 0;
            instance.instanceData.element.appendChild(constraints);
            
            linearConstraintCoefficients = ...
                instance.document.createElement('linearConstraintCoefficients');
            linearConstraintCoefficients.setAttribute('numberOfValues', '0');
            start = instance.document.createElement('start');
            el = instance.document.createElement('el');
            txt = instance.document.createTextNode('0');
            el.appendChild(txt);
            start.appendChild(el);
            colIdx = instance.document.createElement('colIdx');
            value = instance.document.createElement('value');
            linearConstraintCoefficients.appendChild(start);
            linearConstraintCoefficients.appendChild(colIdx);
            linearConstraintCoefficients.appendChild(value);
            instance.instanceData.linearConstraintCoefficients = ...
                struct('element', linearConstraintCoefficients, ...
                'start', start, 'colIdx', colIdx, 'value', value, ...
                'numberOfValues', 0);
            
            nonlinearExpressions = instance.document.createElement( ...
                'nonlinearExpressions');
            nonlinearExpressions.setAttribute( ...
                'numberOfNonlinearExpressions', '1');
            objectiveNL = instance.document.createElement('nl');
            objectiveNL.setAttribute('idx', '-1');
            dummyobjective = OSnonlinear(objectiveNL, 0);
            objectiveNL.appendChild(dummyobjective.element);
            nonlinearExpressions.appendChild(objectiveNL);
            instance.instanceData.nonlinearExpressions = ...
                struct('element', nonlinearExpressions, 'objectiveNL', ...
                objectiveNL, 'numberOfNonlinearExpressions', 1);
            instance.instanceData.element.appendChild(nonlinearExpressions);
        end
        
        function set.instanceHeader(instance, value)
            % property set function for instanceHeader.name,
            % instanceHeader.source, and instanceHeader.description
            if ~isstruct(value)
                error('instanceHeader must be a structure')
            end
            fnames = fieldnames(value);
            for i=1:length(fnames)
                switch fnames{i}
                    case 'element'
                        instance.instanceHeader.element = value.element;
                    case {'name', 'source', 'description'}
                        instance.instanceHeader.(fnames{i}) = ...
                            value.(fnames{i});
                        if isfield(instance.instanceHeader, ...
                                [fnames{i} 'Text'])
                            instance.instanceHeader.([fnames{i} ...
                                'Text']).setTextContent(value.(fnames{i}));
                        else
                            docu = instance.instanceHeader.element.getOwnerDocument;
                            newelem = docu.createElement(fnames{i});
                            instance.instanceHeader.([fnames{i} 'Text']) = ...
                                docu.createTextNode(value.(fnames{i}));
                            newelem.appendChild( ...
                                instance.instanceHeader.([fnames{i} 'Text']));
                            instance.instanceHeader.element.appendChild(newelem);
                        end
                    case {'nameText', 'sourceText', 'descriptionText'}
                        % ignore these
                    otherwise
                        error('invalid field for instanceHeader')
                end
            end
        end
        
        function variable = createVariable(instance, name, type, lb, ub, mult)
            variables = instance.instanceData.variables;
            numvars_prev = variables.numberOfVariables;
            
            % first new element goes in instanceData.variables
            elem1 = instance.document.createElement('var');
            if nargin > 1 && ~isempty(name)
                elem1.setAttribute('name', name);
            end
            if nargin > 2 && ~isempty(type)
                elem1.setAttribute('type', type);
            end
            if nargin > 3 && ~isempty(lb)
                elem1.setAttribute('lb', sprintf('%.17g', lb));
            else
                % either set lb to -INF by default, or make very clear in docu?
                %warning('OSinstance:lb0','note that lb = 0 by default in OSiL')
            end
            if nargin > 4 && ~isempty(ub)
                elem1.setAttribute('ub', sprintf('%.17g', ub));
            end
            if nargin > 5 && ~isempty(mult)
                if mult ~= 1
                    error('mult ~= 1 not currently supported')
                end
                elem1.setAttribute('mult', sprintf('%d', mult));
            end
            variables.element.appendChild(elem1);
            instance.instanceData.variables.numberOfVariables = ...
                numvars_prev + 1;
            variables.element.setAttribute('numberOfVariables', ...
                sprintf('%d', numvars_prev + 1));
            
            % second new element goes in output nonlinear expression
            elem2 = instance.document.createElement('variable');
            elem2.setAttribute('idx', sprintf('%d', numvars_prev));
            elem2.setAttribute('coef', '1.0');
            variable = OSnonlinear(elem2);
            % ideally want a way of setting variable properties after
            % creation too, probably through the OSnonlinear object
        end
        
        function setObjective(instance, expression, name, maxOrMin)
            % numberOfObjCoef? constant? weight? mult?
            % clear attributes or set to default if not specified?
            obj = instance.instanceData.objectives.obj.element;
            if nargin > 2 && ~isempty(name)
                obj.setAttribute('name', name);
            end
            if nargin > 3 && ~isempty(maxOrMin)
                obj.setAttribute('maxOrMin', maxOrMin);
            end
            objectiveNL = instance.instanceData.nonlinearExpressions.objectiveNL;
            objectiveNL.replaceChild(expression.element, ...
                objectiveNL.getFirstChild);
        end
        
        function createConstraint(instance, expression, ...
                name, constant, lb, ub, mult)
            constraints = instance.instanceData.constraints;
            numconstr_prev = constraints.numberOfConstraints;
            linearConstraintCoefficients = ...
                instance.instanceData.linearConstraintCoefficients;
            numlinvals_prev = linearConstraintCoefficients.numberOfValues;
            numlinvals_new = numlinvals_prev;
            
            % check if new constraint is linear
            [idx, val] = linearCoefficients(expression);
            if ~isempty(idx)
                offset = sum(val(idx == -1));
                % use sparse to combine coefficients for duplicate indices
                coefvec = sparse(idx(idx > 0), 1, val(idx > 0));
                numlinvals_new = numlinvals_new + nnz(coefvec);
                instance.instanceData.linearConstraintCoefficients.numberOfValues = ...
                    numlinvals_new;
                linearConstraintCoefficients.element.setAttribute( ...
                    'numberOfValues', sprintf('%d', numlinvals_new));
                
                [idx, ~, val] = find(coefvec);
                for i=1:nnz(coefvec)
                    % add new entry to linearConstraintCoefficients.colIdx
                    el = instance.document.createElement('el');
                    txt = instance.document.createTextNode( ...
                        sprintf('%d', idx(i) - 1)); % zero-based colIdx
                    el.appendChild(txt);
                    linearConstraintCoefficients.colIdx.appendChild(el);
                    % add new entry to linearConstraintCoefficients.value
                    el = instance.document.createElement('el');
                    txt = instance.document.createTextNode( ...
                        sprintf('%.17g', val(i)));
                    el.appendChild(txt);
                    linearConstraintCoefficients.value.appendChild(el);
                end
                
                if numlinvals_prev == 0
                    % add linearConstraintCoefficients section to
                    % instanceData if it was previously empty
                    instance.instanceData.element.insertBefore( ...
                        linearConstraintCoefficients.element, ...
                        instance.instanceData.nonlinearExpressions.element);
                end                
                % adjust constant to account for offset
                if nargin > 3 && ~isempty(constant)
                    constant = constant + offset;
                elseif offset ~= 0
                    constant = offset;
                end
            else
                % constraint is nonlinear, so add new element
                % in instanceData.nonlinearExpressions
                conNL = instance.document.createElement('nl');
                conNL.setAttribute('idx', sprintf('%d', numconstr_prev));
                conNL.appendChild(expression.element.cloneNode(true));
                nonlinearExpressions = instance.instanceData.nonlinearExpressions;
                numNL_prev = nonlinearExpressions.numberOfNonlinearExpressions;
                nonlinearExpressions.element.appendChild(conNL);
                instance.instanceData.nonlinearExpressions.numberOfNonlinearExpressions = ...
                    numNL_prev + 1;
                nonlinearExpressions.element.setAttribute( ...
                    'numberOfNonlinearExpressions', ...
                    sprintf('%d', numNL_prev + 1));
            end
            
            % add new element to instanceData.constraints
            con = instance.document.createElement('con');
            if nargin > 2 && ~isempty(name)
                con.setAttribute('name', name);
            end
            if exist('constant', 'var') && ~isempty(constant)
                con.setAttribute('constant', sprintf('%.17g', constant));
            end
            if nargin > 4 && ~isempty(lb)
                con.setAttribute('lb', sprintf('%.17g', lb));
            end
            if nargin > 5 && ~isempty(ub)
                con.setAttribute('ub', sprintf('%.17g', ub));
            end
            if nargin > 6 && ~isempty(mult)
                if mult ~= 1
                    error('mult ~= 1 not currently supported')
                end
                con.setAttribute('mult', sprintf('%d', mult));
            end
            constraints.element.appendChild(con);
            instance.instanceData.constraints.numberOfConstraints = ...
                numconstr_prev + 1;
            constraints.element.setAttribute('numberOfConstraints', ...
                sprintf('%d', numconstr_prev + 1));
            
            % add new entry to linearConstraintCoefficients.start
            el = instance.document.createElement('el');
            txt = instance.document.createTextNode( ...
                sprintf('%d', numlinvals_new));
            el.appendChild(txt);
            linearConstraintCoefficients.start.appendChild(el);
        end
        
        function variable = substitutionVariable(instance, expression, varargin)
            % create a new variable to replace a nonlinear expression
            % and add a corresponding equality constraint
            variable = instance.createVariable(varargin{:});
            instance.createConstraint(-variable + expression, [], [], 0, 0);
        end
        
        function matrix = linearConstraintMatrix(instance)
            % assemble the linear constraint matrix for Matlab use
            lcc = instance.instanceData.linearConstraintCoefficients;
            colIdx = zeros(lcc.numberOfValues, 1); % preallocate
            value = zeros(lcc.numberOfValues, 1); % preallocate
            for i=1:lcc.numberOfValues
                colIdx(i) = str2double(lcc.colIdx.item(i-1).getTextContent);
                value(i) = str2double(lcc.value.item(i-1).getTextContent);
            end
            start = zeros(lcc.start.getLength, 1); % preallocate
            rowIdx = zeros(lcc.numberOfValues, 1); % preallocate
            start(1) = str2double(lcc.start.item(0).getTextContent);
            for i=1:length(start)-1
                start(i+1) = str2double(lcc.start.item(i).getTextContent);
                rowIdx(start(i)+1 : start(i+1)) = i;
            end
            matrix = sparse(rowIdx, colIdx + 1, value, length(start)-1, ...
                instance.instanceData.variables.numberOfVariables);
        end
    end
end
