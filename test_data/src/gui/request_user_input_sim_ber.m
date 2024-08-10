function [answer, IsCancelled] = request_user_input_sim_ber(default_Mc, default_exp_name)
    % %%%%%%%%%%%%%%%%%%%%
    % get gen_cases
    % %%%%%%%%%%%%%%%%%%%%
    dir_gen_now = get_subfolders_r('_gendata', 2);

	if dir_gen_now.count == 0
		answer      = struct;
		IsCancelled = 1;
		fprintf(1, "Can not find any simulation source.\n" + ...
            "Please execute gen_cases.m first!\n");
		return;
	end

    % %%%%%%%%%%%%%%%%%%%%
    % get simulation result for estimated DFO and Channel
    % %%%%%%%%%%%%%%%%%%%%
    dir_sim_now = get_subfolders_r('_sim_result_est', 3);

    if dir_sim_now.count == 0
		answer      = struct;
		IsCancelled = 1;
        fprintf(1, "Can not find any simulation results.\n" + ...
            "Please execute sim_estimate_channel.m first!\n");
		return;
    end

    % %%%%%%%%%%%%%%%%%%%%
    % combine two source and merge them (also delete invalid folder)
    % %%%%%%%%%%%%%%%%%%%%
    dir_now = struct;

    tmp_names       = intersect(dir_gen_now.name, dir_sim_now.name);
    tmp_name_idx    = zeros(length(tmp_names), 1);
    for i_name = 1 : length(tmp_names)
        name                    = tmp_names(i_name);
        tmp_name_idx(i_name)    = find(strcmp(dir_sim_now.name, name));
    end
    tmp_name_idx(tmp_name_idx == 0) = [];
    tmp_name_idx                    = sort(tmp_name_idx);
    
    count   = 0;
    name    = {};
    array   = [];
    for i_arr = 1 : length(tmp_name_idx)
        exp_struct = dir_sim_now.array(tmp_name_idx(i_arr));
        if exp_struct.count ~= 0
            count       = count + 1;
            name{end+1} = dir_sim_now.name{tmp_name_idx(i_arr)}; %#ok<AGROW> 
            array       = [array; exp_struct]; %#ok<AGROW> 
        end
    end
    
    dir_now.count   = count;
    dir_now.name    = name;
    dir_now.array   = array;

    % %%%%%%%%%%%%%%%%%%%%
    % default answer 1
    % %%%%%%%%%%%%%%%%%%%%
    defInd              = struct;
    defInd.i_sim_src    = 1;
    defInd.i_exp_src    = 1;

	defAns                      = struct;
	defAns.Mc                   = default_Mc;
	defAns.exp_name             = default_exp_name;
	defAns.flag_collect_mode    = false;
	defAns.flag_skip_simulated  = true;
	defAns.flag_use_parfor      = true;
	defAns.flag_use_gpu         = false;

    % %%%%%%%%%%%%%%%%%%%%
    % get MD5 of dir_now
    % %%%%%%%%%%%%%%%%%%%%
    md5 = GetMD5(dir_now, 'Array');

	defVal_path = fullfile('_tmp', 'defaultInput_sim_ber.mat');
    if isfile(defVal_path)
		defVal      = load(defVal_path).defVal;
        md5_prev    = defVal.md5;

        if strcmp(md5, md5_prev)
            defInd = defVal.defInd;
			defAns = defVal.defAns;
		else
			defAns.exp_name             = defVal.defAns.exp_name;
	        defAns.flag_collect_mode    = defVal.defAns.flag_collect_mode;
	        defAns.flag_skip_simulated  = defVal.defAns.flag_skip_simulated;
			defAns.flag_use_parfor      = defVal.defAns.flag_use_parfor;
			defAns.flag_use_gpu         = defVal.defAns.flag_use_gpu;
        end
    end

    % %%%%%%%%%%%%%%%%%%%%
    % default answer 2
    % %%%%%%%%%%%%%%%%%%%%
	defAns.sim_src_name = dir_now.name(defInd.i_sim_src);
	defAns.exp_src_name = dir_now.array(defInd.i_sim_src).name(defInd.i_exp_src);
	
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%% ask for user input %%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	prompt  = {};
	title   = 'Simulation setting';
	formats = {};
	options = struct;
	
	prompt(1,:)             = {'Enter simulation data source:', 'sim_src_name', []};
	formats(1,1).type       = 'list';
	formats(1,1).style      = 'popupmenu';
	formats(1,1).size       = -1;
	formats(1,1).format     = 'text';
	formats(1,1).items      = dir_now.name;
	formats(1,1).callback   = @(~,~,h,k)(set(h(k+1), ...
		'String', dir_now.array(get(h(k),'Value')).name, ...
		'Value', 1));
	formats(1,1).span       = [1 4];
	
	prompt(2,:)         = {'Enter experiment source:', 'exp_src_name',[]};
	formats(2,1).type   = 'list';
	formats(2,1).style  = 'popupmenu';
	formats(2,1).size   = -1;
	formats(2,1).format = 'text';
	formats(2,1).items  = dir_now.array(defInd.i_sim_src).name;
	formats(2,1).span   = [1 4];

	prompt(3,:)         = {'number of simulation trials', 'Mc',[]};
	formats(3,1).type   = 'edit';
	formats(3,1).size   = -1;
	formats(3,1).format = 'integer';
	formats(3,1).span   = [1 4];

	prompt(4,:)         = {'experiment', 'exp_name',[]};
	formats(4,1).type   = 'edit';
	formats(4,1).size   = -1;
	formats(4,1).format = 'text';
	formats(4,1).span   = [1 4];

	prompt(5,:)         = {'collect mode', 'flag_collect_mode',[]};
	formats(5,1).type   = 'check';

	prompt(6,:)         = {'skip simulated', 'flag_skip_simulated',[]};
	formats(5,2).type   = 'check';

	prompt(7,:)         = {'use parfor', 'flag_use_parfor',[]};
	formats(5,3).type   = 'check';

	prompt(8,:)         = {'use GPU', 'flag_use_gpu',[]};
	formats(5,4).type   = 'check';

	options.Resize          = 'off';
	options.AlignControls   = 'on';
	
	[answer, IsCancelled] = inputsdlg(prompt, title, formats, defAns, options);
	if IsCancelled == 0
		answer.sim_src_name = answer.sim_src_name{1};
		answer.exp_src_name = answer.exp_src_name{1};

        % record defAns
        defVal      = struct;
        defVal.md5  = md5;
        % compute ind
        ind             = struct;
        ind.i_sim_src   = find(ismember(dir_now.name, answer.sim_src_name));
        ind.i_exp_src   = find(ismember(dir_now.array(ind.i_sim_src).name, answer.exp_src_name));
        % save to defVal
        defVal.defInd = ind;
        defVal.defAns = answer;
        if ~exist('_tmp', 'dir');	mkdir('_tmp'); end
        save(defVal_path, 'defVal');
	end
end
