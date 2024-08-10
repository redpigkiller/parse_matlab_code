addpath(genpath('src'));

[answer, IsCancelled] = request_user_input();

if IsCancelled == 1
    return
end

remove_dir(answer.sim_src, answer.sim_exp);

%% functions
function [answer, IsCancelled] = request_user_input()
    % %%%%%%%%%%%%%%%%%%%%
    % get dump
    % %%%%%%%%%%%%%%%%%%%%
    dir_dmp_now = get_subfolders_r('_dump', 2);

    if dir_dmp_now.count == 0
		answer = struct;
		IsCancelled = 1;
		fprintf(1, "Can not find any simulation.\n");
		return;
    end

    % %%%%%%%%%%%%%%%%%%%%
    % verify the results
    % %%%%%%%%%%%%%%%%%%%%
    dir_now = struct;

    count = 0;
    name = {};
    array = [];
    for i_arr = 1 : dir_dmp_now.count
        exp_struct = dir_dmp_now.array(i_arr);
        if exp_struct.count ~= 0
            count = count + 1;
            name{end+1} = dir_dmp_now.name{i_arr}; %#ok<AGROW> 
            array = [array; exp_struct]; %#ok<AGROW> 
        end
    end
    
    dir_now.count = count;
    dir_now.name = name;
    dir_now.array = array;

    % %%%%%%%%%%%%%%%%%%%%
    % default answer
    % %%%%%%%%%%%%%%%%%%%%
	defAns = struct;
	defAns.sim_src = dir_now.name(1);
	defAns.sim_exp = dir_now.array(1).name(1);

	% ask for user input
	prompt = {};
	title = 'Remove simulation data';
	formats = {};
	options = struct;

	prompt(1,:) = {'Simulation source:', 'sim_src', []};
	formats(1,1).type = 'list';
	formats(1,1).style = 'listbox';
	formats(1,1).size = [0 200];
	formats(1,1).format = 'text';
	formats(1,1).items = dir_now.name;
	formats(1,1).callback = @(~,~,h,k)(set(h(k+1), ...
		'String', dir_now.array(get(h(k),'Value')).name, ...
		'Value', 1));
	
	prompt(2,:) = {'Simulation experiment:', 'sim_exp',[]};
	formats(2,1).type = 'list';
	formats(2,1).style = 'listbox';
	formats(2,1).size = [0 200];
	formats(2,1).format = 'text';
	formats(2,1).items = dir_now.array(1).name;
    formats(2,1).limits = [0 2]; % multi-select
	
	options.Resize = 'on';
	options.AlignControls = 'on';
    options.ButtonNames = {'Continue for sure?', 'Cancel'};
	
	[answer, IsCancelled] = inputsdlg(prompt, title, formats, defAns, options);
end

function remove_dir(sim_src, sim_exp)
    answer = questdlg('remove for sure?', 'double check', 'yes', 'no', 'no');
    
    switch answer
    case 'yes'

        sim_name = sim_src{1};

        dump_dir = fullfile('_dump', sim_name);
        save_dir = fullfile('_sim_result_est', sim_name);
        res_dir = fullfile('results', sim_name);
        fig_dir = fullfile('figure', sim_name);
        for idx = 1 : length(sim_exp)
            full_path = fullfile(dump_dir, sim_exp{idx});
            status = rmdir(full_path, 's');

            full_path = fullfile(save_dir, sim_exp{idx});
            status = rmdir(full_path, 's');

            full_path = fullfile(res_dir, sim_exp{idx});
            status = rmdir(full_path, 's');

            full_path = fullfile(fig_dir, sim_exp{idx});
            status = rmdir(full_path, 's');
        end
        
    case 'no'
        % do nothing
    end
end
