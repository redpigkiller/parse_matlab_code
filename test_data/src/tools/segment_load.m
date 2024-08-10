function [load_idx, sub_idx] = segment_load(load_Mc_idx, seg_size)
%SEGMENT_LOAD Given the index of Mc (you want to load) and the size of
% segment, return the index of filename and the index of the data in that
% file.

	load_idx = floor((load_Mc_idx - 1) / seg_size) + 1;
	sub_idx = mod(load_Mc_idx - 1, seg_size) + 1;
end

