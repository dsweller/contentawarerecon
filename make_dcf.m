function w = make_dcf(ks,extrapolate)
%
% assume ks are on scale of -FOV/2 to FOV/2 (not -pi to pi)
% each column represents a dimension (2D is two columns, etc.)
%
% This code builds on the voronoi-based density compensation in Fessler's
% image reconstruction toolbox, but with additional capabilities for
% dealing with the outermost sampled points of k-space.

if ~exist('extrapolate','var') || isempty(extrapolate), extrapolate = 'max'; end

% get unique points within tolerance
ks_round = round(ks.*1e4).*1e-4;
[~,ia,ic] = unique(ks_round,'rows','stable');
ks = ks(ia,:);

nrep = accumarray(ic,1,[size(ks,1),1]); % # of repetitions

[v,c] = voronoin(ks);

% get convex hull of ks
inds_hull = convhulln(ks);

% get voronoi points outside of hull
nullsp = arrayfun(@(ii) null(bsxfun(@minus,ks(inds_hull(ii,2:end),:),ks(inds_hull(ii,1),:))),1:size(inds_hull,1),'UniformOutput',false); % normal vector for each facet of hull
if any(cellfun(@(nullsp1) size(nullsp1,2) ~= 1,nullsp))
    warning('Degenerate boundary; assuming all Voronoi points are good.');
    c_in = true(size(c));
else
    nullsp = cellfun(@(nullsp1) nullsp1(:,1).',nullsp,'UniformOutput',false);
    nullsp = cat(1,nullsp{:});
    
    ks_hull = ks(unique(inds_hull(:)),:);
    center_hull = mean(ks_hull,1);
    flip_nullsp = sum(bsxfun(@minus,ks(inds_hull(:,1),:),center_hull).*nullsp,2) < 0;
    nullsp(flip_nullsp,:) = -nullsp(flip_nullsp,:); % ensure normal vector pointing outwards

    v_in = arrayfun(@(ii) bsxfun(@minus,v,ks(inds_hull(ii,1),:))*(nullsp(ii,:).') <= 0,1:size(inds_hull,1),'UniformOutput',false);
    v_in = all(isfinite(v),2) & all(cat(2,v_in{:}),2);
    c_in = cellfun(@(c1) all(v_in(c1)),c);
end

[~,w_in] = cellfun(@(c1) convhulln(v(c1,:)),c(c_in),'UniformOutput',false);
w_in = cat(1,w_in{:});
w = zeros(size(ks,1),1);
w(c_in) = w_in;
if isnumeric(extrapolate)
    w(~c_in) = extrapolate;
else
    switch lower(extrapolate)
        case 'zero'
        case 'nan'
            w(~c_in) = NaN;
        case 'min'
            w(~c_in) = min(w_in);
        case 'max'
            w(~c_in) = max(w_in);
        case 'prev'
            inds_out = cumsum(c_in);
            inds_out(inds_out < 1) = 1;
            w(~c_in) = w(inds_out(~c_in));
        case 'next'
            inds_out = length(w_in) - cumsum(c_in(end:-1:1)) + 1;
            inds_out(inds_out > length(w_in)) = length(w_in);
            w(~c_in) = w(inds_out(~c_in));
        otherwise
            error('Unrecognized extrapolation type');
    end
end

w = w./nrep;
w = w(ic);

end
