function [a_str,font_size] = get_asterisks_from_pval(p,noval)

if nargin < 2
    noval = 1;
end

font_size = 30;

if p > 0.05
    if noval
        a_str = '';
    else
        a_str = 'n.s.';
    end
    font_size = 24;
end

if p < 0.05
    a_str = '*';
end

if p < 0.01
    a_str = '**';
end

if p < 0.001
    a_str = '***';
end

if p < 0.0001
    a_str = '****';
end
