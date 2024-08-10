function result = iif(cond, t, f)
%IIF Evaluates a condition and returns one of two parameters based on the
%result. Emulates the ?: functionality of C and C++.
% 
% Usage:
% iif(cond, t, f)
% Returns t if cond == true, and f otherwise.

    if cond
        result = t;
    else
        result = f;
    end
end