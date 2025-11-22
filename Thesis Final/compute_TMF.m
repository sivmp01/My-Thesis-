%function [TMF_gnd,TMF_middle,TMF_roof] = compute_TMF(num_f,wall_A,window_A,f_A,roof_area,C_air,d_wall,d_floor,d_window,d_roof,rho_wall,rho_window,rho_floor,rho_roof,cp_wall,cp_roof,cp_window,cp_floor)
function [TMF_gnd, TMF_middle, TMF_roof] = compute_TMF(num_f, wall_A, window_A, f_A, roof_area, C_air, ...
    d_wall, rho_wall, cp_wall, ...
    d_window, rho_window, cp_window, ...
    d_floor, rho_floor, cp_floor, ...
    d_roof, rho_roof, cp_roof, ...
    d_ceiling, rho_ceiling, cp_ceiling)
%COMPUTE_TMF Summary of this function goes here
C_air_f= C_air/num_f;
floor_num=zeros(num_f,1);
floor_num(1)=1;
    for i=1:num_f-1
    floor_num(i+1)=floor_num(i)+1;
    end
%COMPUTE_TMF Summary of this function goes here
%%Capacitance of each element
% Capacitance of each multilayer component
p=0.3; % participation factor
C_wall   = wall_A   * sum(d_wall   .* rho_wall   .* cp_wall)*p   ;
C_window = window_A * sum(d_window .* rho_window .* cp_window)*p ;
C_floor  = f_A      * sum(d_floor  .* rho_floor  .* cp_floor) *p ;
C_roof   = roof_area* sum(d_roof   .* rho_roof   .* cp_roof)  *p ;
C_ceiling= f_A      * sum(d_ceiling.* rho_ceiling.* cp_ceiling)*p ;

    if num_f > 2
    % Calculate TMF for middle floor only if >2 floors
    C_total_middle = C_wall + C_window + C_floor + C_ceiling;
    TMF_middle = 1 + (C_total_middle / C_air_f);
    else
    TMF_middle = NaN;
    end
% TMF for ground floor (no ceiling, no roof)
C_total_gnd = C_wall + C_window + C_floor;
TMF_gnd = 1 + (C_total_gnd / C_air_f); 
% TMF for top floor (has roof and floor)
C_total_roof = C_wall + C_window + C_floor + C_roof;
TMF_roof = 1 + (C_total_roof / C_air_f);

disp([TMF_gnd,TMF_middle,TMF_roof])


end

