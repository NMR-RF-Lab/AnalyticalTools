function pkResults = amplitudeRatio(pkResults, fitResults, refPeak)
% Inputs
%   - pkResults: custom OXSA output results structure
%   - fitResults: default OXSA output results structure
%   - refPeak: peak to calculate ratios in reference to 

% Ensure cell string format
if ~iscell(refPeak)
    refPeak = cellstr(refPeak);
end

% Determine which name for 'pk'
if isfield(fitResults, 'pk')
    pkFieldName = 'pk';
elseif isfield(fitResults, 'pkWithLinLsq')
    pkFieldName = 'pkWithLinLsq';
else
    error('cannot find pk structure!')
end

% Grab data parameters
pp = numel(refPeak);
[ii,jj,kk] = size(pkResults); % size
pkNames = [pkResults(:,1,1).peakName];
pkIndices = dictionary(pkNames, 1:ii);

for p = 1:pp
    for j = 1:jj
        for k = 1:kk
            % Special case of "total" or "all"
            if strcmpi(refPeak{p},"total") || strcmpi(refPeak{p},"all") % separate checks necessary in case refPeak{p} is nonscalar
                PKS = pkNames;
                [num, mask] = convToSym(pkResults,pkIndices,j,k,false,PKS);
            else
                PKS = refPeak{p};
                [num, mask] = convToSym(pkResults,pkIndices,j,k,false,pkNames(~strcmp(pkNames,PKS))); % doesn't work when refPeak{p} is nonscalar, need to fix
            end
            denom = convToSym(pkResults,pkIndices,j,k,true,PKS);
            

            % Combine & convert to char
            derivedStr = cellfun(@(x) x + " / " + denom, num, 'UniformOutput', false);
            derivedStr = cellfun(@char, derivedStr, 'UniformOutput', false);

            % Calculate ratio
            [rVal, rCRB] = AMARES.estimateDerivedParamAndCRB(fitResults(j,k).(pkFieldName), ...
                                                             fitResults(j,k).xFit, ...
                                                             fitResults(j,k).constraintsCellArray, ...
                                                             fitResults(j,k).covariance, ...
                                                             derivedStr);

            % Reorganize by peak
            Dx = 0;
            for i = 1:ii
                if mask(i)
                    Dx = Dx + 1;
                    pkResults(i,j,k).ratioDenom{p} = refPeak{p};
                    pkResults(i,j,k).ratio(p) = rVal(Dx);
                    pkResults(i,j,k).ratioCRB(p) = rCRB(Dx);
                else
                    pkResults(i,j,k).ratioDenom{p} = NaN;
                    pkResults(i,j,k).ratio(p) = NaN;
                    pkResults(i,j,k).ratioCRB(p) = NaN;
                end
            end
        end
    end
end
    function [pk, mask] = convToSym(pkResults,pkIndices,J,K,FLAG,refPeak)
        % Initialize pk indices & logical mask
        I = 1:size(pkResults,1);
        mask = false(size(I));

        % Select peak indices
        mask(pkIndices(refPeak)) = true;

        % Ensure none are empty
        mask(cellfun(@(x) isempty(x) || (isscalar(x) && strcmp(x,"")), {pkResults(:,J,K).subPeakNames})) = false;
        
        I = I(mask);

        % Determine numerator & denominator
        pk = {pkResults(I,J,K).subPeakNames};

        % Grab total number of peaks
        pks = cellfun(@(x) cellfun(@numel,x), {pkResults(I,J,K).multDx}, 'UniformOutput',false); % embedded cellfun to account for grouped peaks (Pitot, ATP)
        pk = cellfun(@(x,y) arrayfun(@(a,b) repmat(a,1,b), x, y, 'UniformOutput',false), pk, pks,'UniformOutput',false); % duplicate based on number of multiplets
        for n = 1:numel(pks)
            for subn = 1:numel(pks{n})
                if pks{n}(subn) > 1
                    pk{n}{subn} = pk{n}{subn} + string(1:pks{n}(subn)); % add numbering to multiplets
                end
            end
        end

        % Extract from the added bottom cell layer
        pk = cellfun(@(x) horzcat(x{:}), pk, 'UniformOutput', false);

        % Add amplitude symbol
        pk = cellfun(@(x) x + "_am", pk, 'UniformOutput', false);

        % Add summation symbol, where applicable
        for n = 1:numel(pk)
            if ~isscalar(pk{n})
                pk{n} = join(pk{n}," + ");
            end
        end

        % If desired, combine all
        if FLAG && ~isscalar(pk)
            pk = {join([pk{:}]," + ")};
        end

        % Add parenthesis
        pk = cellfun(@(x) "(" + x + ")", pk, 'UniformOutput', false);
    end
end