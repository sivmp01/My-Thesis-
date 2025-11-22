load('JanFeb_Prices.mat');
whos
if exist('TT_full','var')
    disp(class(TT_full));
    if istimetable(TT_full)
        % Case A: TT_full is a timetable → use row times
        t_all = TT_full.Properties.RowTimes;
    else
        % Case B: TT_full is a table → find first datetime column
        dtMask = varfun(@(x) isa(x,'datetime'), TT_full, 'OutputFormat','uniform');
        if any(dtMask)
            dtVar = TT_full.Properties.VariableNames{find(dtMask,1,'first')};
            t_all = TT_full.(dtVar);
        else
            error('No datetime column found in TT_full. Available vars: %s', strjoin(TT_full.Properties.VariableNames,', '));
        end
    end
elseif exist('DateTime','var')
    % Case C: datetime is a standalone variable
    t_all = DateTime;
else
    error('No datetime vector found. Check variables in JanFeb_Prices.mat');
end
