function [results, h2_array_avg]=ins_countlinks(H,thr_h2, onchan, type)

% MODIFIED BY SARA 22/03/21 (normalization of degrees)

% H: H2 from Anywave mat file
% threshold: threshold on graph (=0: strength, >0: degrees)
% links{1}: OUT, {2}: IN, {3}: TOT (strength)
% [links,tit,times]=ins_countlinks(H,thr_h2,thr_lag, verbose)


% links{1}: OUT, {2} IN, {3} OUT/IN (degrees) or OUT-INT (strength)
%
numchan=size(H.aw_h2,1);
if nargin < 3
    onchan = 0;
end
if nargin <4 && onchan==1
    type = 'mean';
end
% pick highest h2 across directions
[h2_array,~]=ins_findmaxh2(H);
results = struct();
if onchan == 1
    if strcmpi(type, 'mean')
        h2_array =mean(h2_array, 3);
    elseif strcmpi(meth, 'max')
        h2_array =max(h2_array, [], 3);
    end
    h2_array_avg =h2_array;
    totdegree = (h2_array >= thr_h2);
    results.indegree = triu(totdegree)';
    results.outdegree =  tril(totdegree);
    results.totdegree = results.indegree + results.outdegree;
    totstrength = h2_array.*(h2_array >= thr_h2);
    results.instrength = triu(totstrength)';
    results.outstrength =  tril(totstrength);
    results.totstrength = results.instrength + results.outstrength;
else
    %%Calcul In, Out and tot 
    %degree
    h2_array_sup = (h2_array >= thr_h2);
    results.indegree = squeeze(sum(h2_array_sup, 1));
    results.outdegree = squeeze(sum(h2_array_sup, 2));
    results.totdegree = results.indegree + results.outdegree;

    %strength
    results.instrength = squeeze(sum(h2_array, 1));
    results.outstrength = squeeze(sum(h2_array, 2));
    results.totstrength = results.instrength + results.outstrength;

    %strength normalized
    results.instrength_norm = results.instrength/(numchan-1);
    results.outstrength_norm = results.outstrength/(numchan-1);
    results.totstrength_norm = results.instrength_norm + results.outstrength_norm;
    
    h2_array_avg =mean(h2_array, 3);
end

end