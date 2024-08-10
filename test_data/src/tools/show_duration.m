function oupt_str = show_duration(time_in_sec)
%SHOW_DURATION Simply convert the time in seconds to day, hour, minute, ...
%and return the string.
%

    time_in_ms = ceil(1000 * time_in_sec);

    % milliseconds
    t_ms = mod(time_in_ms, 1000);

    % seconds
    time_in_ms = fix(time_in_ms / 1000);
    t_sec = mod(time_in_ms, 60);

    % minutes
    time_in_ms = fix(time_in_ms / 60);
    t_min = mod(time_in_ms, 60);

    % hours
    time_in_ms = fix(time_in_ms / 60);
    t_hr = mod(time_in_ms, 24);

    % days
    time_in_ms = fix(time_in_ms / 24);
    t_day = mod(time_in_ms, 24);

    % generate output string
    oupt_str = "";

    flag_out = false;
    if flag_out || t_day ~= 0
        flag_out = true;
        oupt_str = oupt_str + num2str(t_day) + " d, ";
    end

    if flag_out || t_hr ~= 0
        flag_out = true;
        oupt_str = oupt_str + num2str(t_hr) + " hr, ";
    end

    if flag_out || t_min ~= 0
        flag_out = true;
        oupt_str = oupt_str + num2str(t_min) + " min, ";
    end

    if flag_out || t_sec ~= 0
        flag_out = true;
        oupt_str = oupt_str + num2str(t_sec) + " sec, ";
    end

    if flag_out || t_ms ~= 0
        oupt_str = oupt_str + num2str(t_ms) + " ms";
    end

end