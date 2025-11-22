% %% --- Gurobi one-line diagnostics ---
% 
% % 1) show what licence variables MATLAB sees
% disp('--- ENV ---')
% disp(getenv('GRB_WLSACCESSID'))
% disp(getenv('GRB_WLSSECRET'))
% disp(getenv('GRB_LICENSEID'))
% disp(getenv('GRB_LICENSE_FILE'))
% 
% % 2) show what mex file is being called
% disp('--- WHICH GUROBI ---')
% which gurobi -all
% 
% % 3) run the tiniest possible model (no params)
% disp('--- SOLVE TEST ---')
% model.A          = sparse(1,3,1);   % 1×3 matrix with one "1"
% model.obj        = [1 0 0];
% model.modelsense = 'min';
% model.rhs        = 1;
% model.sense      = '=';
% try
%     result = gurobi(model);
%     disp(result.objval)
% catch ME
%     disp('*** ERROR CAUGHT ***')
%     disp(ME.message)
% end
% OLD (causes error)
% gurobi('version')

% NEW: either query version like this …
v = gurobi([]);           % returns a struct with version info
disp(v.vers)

% … or, much simpler, just solve a tiny LP to confirm everything
model.A          = sparse(1,3,1);
model.obj        = [1 0 0];
model.modelsense = 'min';
model.rhs        = 1;
model.sense      = '=';
result = gurobi(model);
fprintf('Solver OK, obj = %.4f\n', result.objval);
