% Matlab class for Optimization Services nonlinear expression
classdef OSnonlinear < handle
    properties
        element % xml element
        fcn % Matlab function representation
    end
    methods
        function nonlinear = OSnonlinear(element, number)
            if nargin == 1
                nonlinear.element = element;
                if element.getNodeName.equals('variable')
                    idx = str2double(element.getAttribute('idx'));
                    coef = str2double(element.getAttribute('coef'));
                    nonlinear.fcn = @(x) coef*x(idx+1); % 1-based indexing
                end
            elseif nargin > 1
                % create nonlinear expression from a number
                % in parent document of first input
                if length(number) > 1
                    error('only scalars currently supported')
                end
                if isa(element, 'OSinstance')
                    document = element.document;
                else
                    document = element.getOwnerDocument;
                end
                nonlinear.element = document.createElement('number');
                nonlinear.element.setAttribute('value', ...
                    sprintf('%.17g', number));
                nonlinear.fcn = @(x) number;
            end
        end
        
        function [idx, val] = linearCoefficients(in)
            % determine whether an OSnonlinear object is linear
            % if nonlinear, return empty
            % if linear, return vector of (1-based) indices
            % and vector of coefficient values
            % constant terms are assigned idx = -1
            % currently only handles sum, minus, variable, or number elements
            idx = [];
            val = [];
            if in.element.getNodeName.equals('sum')
                imax = in.element.getLength;
                elem = in.element.getFirstChild;
                idx = zeros(imax, 1); % preallocate
                val = zeros(imax, 1); % preallocate
            elseif in.element.getNodeName.equals('minus')
                % handle subtraction recursively on each child
                [idx, val] = linearCoefficients( ...
                    OSnonlinear(in.element.getFirstChild));
                if isempty(idx) || isempty(val)
                    % first child is nonlinear
                    return
                end
                [idx2, val2] = linearCoefficients( ...
                    OSnonlinear(in.element.getLastChild));
                if isempty(idx2) || isempty(val2)
                    % second child is nonlinear
                    idx = idx2;
                    val = val2;
                else
                    % both children are linear
                    idx = [idx; idx2];
                    val = [val; -val2];
                end
                return
            else
                imax = 1;
                elem = in.element;
            end
            for i=1:imax
                if elem.getNodeName.equals('variable')
                    idx(i) = str2double(elem.getAttribute('idx')) + 1;
                    val(i) = str2double(elem.getAttribute('coef'));
                elseif elem.getNodeName.equals('number')
                    idx(i) = -1;
                    val(i) = str2double(elem.getAttribute('value'));
                else
                    idx = [];
                    val = [];
                    return
                end
                elem = elem.getNextSibling;
            end
        end
        
        % should check for non-scalar inputs, add vector functions like
        % norm, dot, sum, prod, concatenation?
        
        function out = plus(in1, in2)
            if isnumeric(in1)
                if isequal(in1, 0)
                    % special case for 0
                    out = in2;
                    return
                end
                in1 = OSnonlinear(in2.element, in1);
                elem1 = in1.element;
                elem2 = in2.element.cloneNode(true);
            elseif isnumeric(in2)
                if isequal(in2, 0)
                    % special case for 0
                    out = in1;
                    return
                end
                in2 = OSnonlinear(in1.element, in2);
                elem2 = in2.element;
                elem1 = in1.element.cloneNode(true);
            else
                elem1 = in1.element.cloneNode(true);
                elem2 = in2.element.cloneNode(true);
            end
            document = elem1.getOwnerDocument;
            if ~isequal(document, elem2.getOwnerDocument)
                error('inputs must be expressions from same OSinstance')
            end
            if elem1.getNodeName.equals('sum') && ...
                    elem2.getNodeName.equals('sum')
                outelem = elem1;
                for i=1:elem2.getLength
                    outelem.appendChild(elem2.getFirstChild);
                end
            elseif elem1.getNodeName.equals('sum')
                outelem = elem1;
                outelem.appendChild(elem2);
            elseif elem2.getNodeName.equals('sum')
                outelem = elem2;
                outelem.insertBefore(elem1, outelem.getFirstChild);
            else
                outelem = document.createElement('sum');
                outelem.appendChild(elem1);
                outelem.appendChild(elem2);
            end
            out = OSnonlinear(outelem);
            out.fcn = @(x) in1.fcn(x) + in2.fcn(x);
        end
        
        function out = minus(in1, in2)
            if isnumeric(in2) || in2.element.getNodeName.equals('variable')
                % subtract constant by adding opposite value,
                % subtract variable by swapping sign of coef (see uminus)
                out = in1 + (-in2);
                return
            elseif isnumeric(in1)
                if isequal(in1, 0)
                    % special case for 0
                    out = -in2;
                    return
                end
                in1 = OSnonlinear(in2.element, in1);
                elem1 = in1.element;
                elem2 = in2.element.cloneNode(true);
            else
                elem1 = in1.element.cloneNode(true);
                elem2 = in2.element.cloneNode(true);
            end
            document = elem1.getOwnerDocument;
            if ~isequal(document, elem2.getOwnerDocument)
                error('inputs must be expressions from same OSinstance')
            end
            outelem = document.createElement('minus');
            outelem.appendChild(elem1);
            outelem.appendChild(elem2);
            out = OSnonlinear(outelem);
            out.fcn = @(x) in1.fcn(x) - in2.fcn(x);
        end
        
        function out = mtimes(in1, in2)
            if isnumeric(in1)
                elem2 = in2.element.cloneNode(true);
                if elem2.getNodeName.equals('variable')
                    % constant times variable, fold into coef
                    outelem = elem2;
                    coef = str2double(outelem.getAttribute('coef'));
                    outelem.setAttribute('coef', sprintf('%.17g', in1*coef));
                    out = OSnonlinear(outelem);
                    return
                end
                in1 = OSnonlinear(in2.element, in1);
                elem1 = in1.element;
            elseif isnumeric(in2)
                elem1 = in1.element.cloneNode(true);
                if elem1.getNodeName.equals('variable')
                    % variable times constant, fold into coef
                    outelem = elem1;
                    coef = str2double(outelem.getAttribute('coef'));
                    outelem.setAttribute('coef', sprintf('%.17g', coef*in2));
                    out = OSnonlinear(outelem);
                    return
                end
                in2 = OSnonlinear(in1.element, in2);
                elem2 = in2.element;
            else
                elem1 = in1.element.cloneNode(true);
                elem2 = in2.element.cloneNode(true);
            end
            document = elem1.getOwnerDocument;
            if ~isequal(document, elem2.getOwnerDocument)
                error('inputs must be expressions from same OSinstance')
            end
            if elem1.getNodeName.equals('product') && ...
                    elem2.getNodeName.equals('product')
                outelem = elem1;
                for i=1:elem2.getLength
                    outelem.appendChild(elem2.getFirstChild);
                end
            elseif elem1.getNodeName.equals('product')
                outelem = elem1;
                outelem.appendChild(elem2);
            elseif elem2.getNodeName.equals('product')
                outelem = elem2;
                outelem.insertBefore(elem1, outelem.getFirstChild);
            else
                outelem = document.createElement('product');
                outelem.appendChild(elem1);
                outelem.appendChild(elem2);
            end
            out = OSnonlinear(outelem);
            out.fcn = @(x) in1.fcn(x) * in2.fcn(x);
        end
        
        function out = mrdivide(in1, in2)
            if isnumeric(in1)
                in1 = OSnonlinear(in2.element, in1);
                elem1 = in1.element;
                elem2 = in2.element.cloneNode(true);
            elseif isnumeric(in2)
                elem1 = in1.element.cloneNode(true);
                if elem1.getNodeName.equals('variable')
                    % variable divided by constant, fold into coef
                    outelem = elem1;
                    coef = str2double(outelem.getAttribute('coef'));
                    outelem.setAttribute('coef', sprintf('%.17g', coef/in2));
                    out = OSnonlinear(outelem);
                    return
                end
                in2 = OSnonlinear(in1.element, in2);
                elem2 = in2.element;
            else
                elem1 = in1.element.cloneNode(true);
                elem2 = in2.element.cloneNode(true);
            end
            document = elem1.getOwnerDocument;
            if ~isequal(document, elem2.getOwnerDocument)
                error('inputs must be expressions from same OSinstance')
            end
            outelem = document.createElement('divide');
            outelem.appendChild(elem1);
            outelem.appendChild(elem2);
            out = OSnonlinear(outelem);
            out.fcn = @(x) in1.fcn(x) / in2.fcn(x);
        end
        
        function out = mpower(in1, in2)
            if isnumeric(in1)
                in1 = OSnonlinear(in2.element, in1);
                elem1 = in1.element;
                elem2 = in2.element.cloneNode(true);
            elseif isnumeric(in2)
                if isequal(in2, 2)
                    % special case for square
                    out = univariate_fcn(in1, 'square');
                    out.fcn = @(x) in1.fcn(x)^2;
                    return
                end
                in2 = OSnonlinear(in1.element, in2);
                elem2 = in2.element;
                elem1 = in1.element.cloneNode(true);
            else
                elem1 = in1.element.cloneNode(true);
                elem2 = in2.element.cloneNode(true);
            end
            document = elem1.getOwnerDocument;
            if ~isequal(document, elem2.getOwnerDocument)
                error('inputs must be expressions from same OSinstance')
            end
            outelem = document.createElement('power');
            outelem.appendChild(elem1);
            outelem.appendChild(elem2);
            out = OSnonlinear(outelem);
            out.fcn = @(x) in1.fcn(x) ^ in2.fcn(x);
        end
        
        function out = uplus(in1) % y = +x;
            out = in1;
        end
        
        function out = univariate_fcn(in1, fcn)
            elem1 = in1.element.cloneNode(true);
            document = elem1.getOwnerDocument;
            outelem = document.createElement(fcn);
            outelem.appendChild(elem1);
            out = OSnonlinear(outelem);
        end
        
        function out = uminus(in1)
            if in1.element.getNodeName.equals('variable')
                % unary minus of a variable, swap sign of coef
                outelem = in1.element.cloneNode(true);
                coefstr = outelem.getAttribute('coef').toCharArray;
                if coefstr(1) == '-'
                    outelem.setAttribute('coef', coefstr(2:end));
                else
                    outelem.setAttribute('coef', ['-'; coefstr]);
                end
                out = OSnonlinear(outelem);
            else
                out = univariate_fcn(in1, 'negate');
                out.fcn = @(x) -in1.fcn(x);
            end
        end
        
        function out = sqrt(in1)
            out = univariate_fcn(in1, 'squareRoot');
            out.fcn = @(x) sqrt(in1.fcn(x));
        end
        
        function out = log(in1)
            out = univariate_fcn(in1, 'ln');
            out.fcn = @(x) log(in1.fcn(x));
        end
        
        function out = log10(in1)
            out = univariate_fcn(in1, 'log10');
            out.fcn = @(x) log10(in1.fcn(x));
        end
        
        function out = exp(in1)
            out = univariate_fcn(in1, 'exp');
            out.fcn = @(x) exp(in1.fcn(x));
        end
        
        function out = sin(in1)
            out = univariate_fcn(in1, 'sin');
            out.fcn = @(x) sin(in1.fcn(x));
        end
        
        function out = cos(in1)
            out = univariate_fcn(in1, 'cos');
            out.fcn = @(x) cos(in1.fcn(x));
        end
        
        function out = tan(in1)
            out = univariate_fcn(in1, 'tan');
            out.fcn = @(x) tan(in1.fcn(x));
        end
        
        function out = asin(in1)
            out = univariate_fcn(in1, 'arcsin');
            out.fcn = @(x) asin(in1.fcn(x));
        end
        
        function out = acos(in1)
            out = univariate_fcn(in1, 'arccos');
            out.fcn = @(x) acos(in1.fcn(x));
        end
        
        function out = atan(in1)
            out = univariate_fcn(in1, 'arctan');
            out.fcn = @(x) atan(in1.fcn(x));
        end
        
        function out = sinh(in1)
            out = univariate_fcn(in1, 'sinh');
            out.fcn = @(x) sinh(in1.fcn(x));
        end
        
        function out = cosh(in1)
            out = univariate_fcn(in1, 'cosh');
            out.fcn = @(x) cosh(in1.fcn(x));
        end
        
        function out = tanh(in1)
            out = univariate_fcn(in1, 'tanh');
            out.fcn = @(x) tanh(in1.fcn(x));
        end
        
        function out = asinh(in1)
            out = univariate_fcn(in1, 'arcsinh');
            out.fcn = @(x) asinh(in1.fcn(x));
        end
        
        function out = acosh(in1)
            out = univariate_fcn(in1, 'arccosh');
            out.fcn = @(x) acosh(in1.fcn(x));
        end
        
        function out = atanh(in1)
            out = univariate_fcn(in1, 'arctanh');
            out.fcn = @(x) atanh(in1.fcn(x));
        end
        
        function out = cot(in1)
            out = univariate_fcn(in1, 'cot');
            out.fcn = @(x) cot(in1.fcn(x));
        end
        
        function out = coth(in1)
            out = univariate_fcn(in1, 'coth');
            out.fcn = @(x) coth(in1.fcn(x));
        end
        
        function out = acot(in1)
            out = univariate_fcn(in1, 'arccot');
            out.fcn = @(x) acot(in1.fcn(x));
        end
        
        function out = acoth(in1)
            out = univariate_fcn(in1, 'arccoth');
            out.fcn = @(x) acoth(in1.fcn(x));
        end
        
        function out = sec(in1)
            out = univariate_fcn(in1, 'sec');
            out.fcn = @(x) sec(in1.fcn(x));
        end
        
        function out = sech(in1)
            out = univariate_fcn(in1, 'sech');
            out.fcn = @(x) sech(in1.fcn(x));
        end
        
        function out = asec(in1)
            out = univariate_fcn(in1, 'arcsec');
            out.fcn = @(x) asec(in1.fcn(x));
        end
        
        function out = asech(in1)
            out = univariate_fcn(in1, 'arcsech');
            out.fcn = @(x) asech(in1.fcn(x));
        end
        
        function out = csc(in1)
            out = univariate_fcn(in1, 'csc');
            out.fcn = @(x) csc(in1.fcn(x));
        end
        
        function out = csch(in1)
            out = univariate_fcn(in1, 'csch');
            out.fcn = @(x) csch(in1.fcn(x));
        end
        
        function out = acsc(in1)
            out = univariate_fcn(in1, 'arccsc');
            out.fcn = @(x) acsc(in1.fcn(x));
        end
        
        function out = acsch(in1)
            out = univariate_fcn(in1, 'arccsch');
            out.fcn = @(x) acsch(in1.fcn(x));
        end
        
        function out = erf(in1)
            out = univariate_fcn(in1, 'erf');
            out.fcn = @(x) erf(in1.fcn(x));
        end
    end
end
