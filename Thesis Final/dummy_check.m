% --- ultra-small model: 1 variable, 0 constraints ---

model.A          = sparse(0,1);   % 0-by-1 empty matrix
model.obj        = 0;             % objective coeff for that lone var
model.modelsense = 'min';
model.rhs        = [];            % no rows → empty RHS
model.sense      = '';            % no rows → empty sense string

params.OutputFlag = 0;            % silence solver log

info = gurobi(model, params);     % call succeeds
disp(info.vers)                   % prints the Gurobi version
