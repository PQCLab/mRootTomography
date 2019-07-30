function h = rt_fprintreplace(str,h)

if nargin > 1
    fprintf(repmat('\b',1,h.nbytes));
end

h.str = str;
h.nbytes = fprintf(str);

end

