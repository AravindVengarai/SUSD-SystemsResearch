global exponential
close all
clear all
Rand_Seeds=load('input_data/Rand_Seeds.txt');

num_runs = 10;
num_funcs = 10;
func_range = 1:10;

dimension = 10;
num_agents = floor(dimension*2);
num_swarms = floor(dimension*5);

max_FEs = getMaxFEs(dimension); %change this if dim ~= 5,10,15,20
%max_FEs = 50000;
max_iterations = ceil(max_FEs/(num_agents*num_swarms))+1;
exponential = false; %MAKE SURE TO CHANGE!
use_ms = true;
use_variable_gains = false; 

graph_bounds = 100;

% LINEAR PARAMETERS
l_susd_gain = 8; %changing these do nothing if use_variable_gains is true
l_form_gain = 1e-3;
l_form_dist = 10; %changing this will change the starting form dist
stopping_condition = 10^-100;

% EXPONENTIAL PARAMETERS
e_susd_gain = 3; %changing these do nothing if use_variable_gains is true
e_form_gain = 0;
e_form_dist = 2;
epsilon = 0.01;
window = 25000;


%variable tuning arrays to store data
global num_same_iterations_before_shrink num_same_iterations prev_min_pos ...
    dist_scale form_gain_scale susd_gain_scale num_wrong_dir prev_min_func_val
num_same_iterations_before_shrink = 10;
num_same_iterations = ones(1, num_swarms);
prev_min_pos = zeros(dimension, num_swarms);
num_wrong_dir = ones(1, num_swarms);
prev_min_func_val = ones(1, num_swarms);
dist_scale = ones(1, num_swarms);
form_gain_scale = ones(1, num_swarms);
susd_gain_scale = ones(1, num_swarms);

optimum= [100, 1100 ,700 ,1900 ,1700 ,1600 ,2100 ,2200 ,2400 ,2500];

C = [
    0 0 0;
    0 0 1;
    0 1 0;
    0 1 1;
    1 0 0;
    1 0 1;
    1 1 0;
    1 1 1;
];
for config=1:size(C,1)
    Alg_Name=[ 'SUSD_(' num2str(C(config,1)) num2str(C(config,2)) num2str(C(config,3)) ')'];
    func_list = getFunc(C(config,:));
    all_errors = zeros(16,num_runs, num_funcs);
    for func_num = func_range
        fprintf('func number: %i\n', func_num)
        func = @(x) func_list(x,func_num);
        func_min = optimum(func_num);
        for run = 1:num_runs
            fprintf('run number: %i\n', run)
            run_seed=Rand_Seeds(mod(dimension*func_num*num_runs+run-num_runs,length(Rand_Seeds))+1);
             x0 = generateStartingPositions(dimension, -graph_bounds,...
                     graph_bounds, num_agents, run_seed, num_swarms, true, l_form_dist);
            if exponential
                if use_variable_gains
                    [xx, xmins, fmins, iter, s_conditions] = susd_exp_search(x0,...
                                                                           func,...
                                                                           epsilon,...
                                                                           window,...
                                                                           max_iterations,...
                                                                           max_FEs,...
                                                                           e_susd_gain,...
                                                                           e_form_gain,...
                                                                           e_form_dist,...
                                                                           @getGains);
                else
                [xx, xmins, fmins, iter, s_conditions] = susd_exp_search(x0,...
                                                                       func,...
                                                                       epsilon,...
                                                                       window,...
                                                                       max_iterations,...
                                                                       max_FEs,...
                                                                       e_susd_gain,...
                                                                       e_form_gain,...
                                                                       e_form_dist);
                end
            else
                if use_variable_gains
                    if use_ms
                        [xx, xmins, fmins, iter, s_conditions] = susd_search_ms(x0,...
                                                                           func,...
                                                                           stopping_condition,...
                                                                           max_iterations,...
                                                                           max_FEs,...
                                                                           l_susd_gain,...
                                                                           l_form_gain,...
                                                                           l_form_dist,...
                                                                           @getGains);
                    else
                        [xx, xmins, fmins, iter, s_conditions] = susd_search(x0,...
                                                       func,...
                                                       stopping_condition,...
                                                       max_iterations,...
                                                       max_FEs,...
                                                       l_susd_gain,...
                                                       l_form_gain,...
                                                       l_form_dist,...
                                                       @getGains);
                    end
                else
                    if use_ms
                        [xx, xmins, fmins, iter, s_conditions] = susd_search_ms(x0,...
                                                                           func,...
                                                                           stopping_condition,...
                                                                           max_iterations,...
                                                                           max_FEs,...
                                                                           l_susd_gain,...
                                                                           l_form_gain,...
                                                                           l_form_dist);
                    else
                        [xx, xmins, fmins, iter, s_conditions] = susd_search(x0,...
                                                       func,...
                                                       stopping_condition,...
                                                       max_iterations,...
                                                       max_FEs,...
                                                       l_susd_gain,...
                                                       l_form_gain,...
                                                       l_form_dist);
                    end
                end
            end

            %Calculating error matrix for cec
            errors = zeros(1, 16);
            for i=0:15
                ii=floor(dimension^((i/5)-3) * max_FEs/(num_agents*num_swarms));
                if ii==0
                    ii=1;
                end
               errors(i+1) = fmins(ii);
               if C(1) == 1
                  errors(i+1) =  errors(i+1) - func_min; 
               end
            end 
            all_errors(:,run, func_num) = errors';


            if dimension == 2 && run == 5
                graph2Dsusd(func, xx, [-graph_bounds graph_bounds...
                                     -graph_bounds graph_bounds]);
            end
        end
        all_errors_this_func = all_errors(:,:,func_num);

        file_name=sprintf('Results\\%s_%s_%s.txt',Alg_Name, int2str(func_num),int2str(dimension));
        save(file_name, 'all_errors_this_func', '-ascii');
    end
    %generate table II
    table2 = zeros(10,6);
    for j = 1:10
       table2(j,1) = j;
       table2(j,2) = min(all_errors(end,:,j));
       table2(j,3) = max(all_errors(end,:,j));
       table2(j,4) = median(all_errors(end,:,j));
       table2(j,5) = mean(all_errors(end,:,j));
       table2(j,6) = std(all_errors(end,:,j));
    end
    file_name=sprintf('Results\\%s_%s_%s.txt',Alg_Name, "all_funcs",int2str(dimension));
    save(file_name, 'table2', '-ascii');
end

%%%%%%%%%%%%%%FUNCTIONS USED FOR VARIABLE GAINS%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [form_dist, form_gain, susd_gain] = getGains(avg_feval, min_pos, swarm_num)
    global exponential num_same_iterations_before_shrink ...
        num_same_iterations prev_min_pos dist_scale susd_gain_scale ...
        form_gain_scale prev_min_func_val num_wrong_dir

    D = size(min_pos,2);
    
    % scale by some num if the min_pos doesn't change for some iterations
    tolerance = min_pos * .01;
    if abs((min_pos - tolerance)) < abs(prev_min_pos) & abs(prev_min_pos) < abs((min_pos + tolerance))
        num_same_iterations(swarm_num) = num_same_iterations(swarm_num) + 1;
    else
        if prev_min_func_val < avg_feval
           num_wrong_dir = num_wrong_dir + 1; 
        else
           num_wrong_dir = 0;
        end
        num_same_iterations(swarm_num) = 1;
        prev_min_pos(:,swarm_num) = min_pos;
    end
    if num_wrong_dir > 1
       susd_gain_scale(swarm_num) = susd_gain_scale(swarm_num)*0.5e-2; 
       form_gain_scale(swarm_num) = form_gain_scale(swarm_num)*0.5e-2;
       dist_scale(swarm_num) = dist_scale(swarm_num)*1.2; 
       num_wrong_dir = 0;
    else
       susd_gain_scale(swarm_num) = susd_gain_scale(swarm_num)*1.015;
       form_gain_scale(swarm_num) = form_gain_scale(swarm_num)*1.0001;
       dist_scale(swarm_num) = dist_scale(swarm_num)*0.9; 
    end
    if num_same_iterations(swarm_num) >= num_same_iterations_before_shrink
       dist_scale(swarm_num) = dist_scale(swarm_num)*0.95;
       form_gain_scale(swarm_num) = form_gain_scale(swarm_num)*1.2;
       susd_gain_scale(swarm_num) = susd_gain_scale(swarm_num)*0.99;
       num_same_iterations(swarm_num) = 1;
    end
        
    prev_min_func_val = avg_feval;
    
    % form dist
    if exponential
        form_dist = min(50*(dist_scale(swarm_num)+1e-14),100);
    else
        form_dist = min(5*(dist_scale(swarm_num)+1e-14), 50);
    end
    
    %form gain
    if exponential
        form_gain = 1.9e-2  * (form_gain_scale(swarm_num)+1e-10);
    else
        form_gain = 1e-3 ;
    end
    
    % susd gain
    if exponential
        susd_gain = 7000 * (susd_gain_scale(swarm_num)+1e-10);
    else
        susd_gain = 2;
    end
end