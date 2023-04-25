function lines = connect_points(pts3d, dist_thresh)

pdists = squareform(pdist(pts3d'));
nPts = length(pts3d);

lines = {};
% pick a starting point, it might as well be 1.

curr_pts = 1:nPts;
all_pts = 1:nPts;

linecounter = 1;

while ~isempty(curr_pts)

    rand_pt = curr_pts(randi(length(curr_pts)));
    curr_pts = setdiff(curr_pts,[rand_pt]);

    this_line = [rand_pt];
    vol_break=0;

    while vol_break==0
        curr_pt = this_line(end);
        dists_to_curr_pt = pdists(curr_pt,:);
        dists_to_curr_pt(curr_pt)=9999;
        dists_to_curr_pt(setdiff(all_pts,curr_pts))=9999;
        
        [dist_to_next,closest_next] = min(dists_to_curr_pt);

        if dist_to_next> dist_thresh
            vol_break = 1; % exit this line and start over
        else
            this_line = [this_line;closest_next];
            curr_pts = setdiff(curr_pts,[closest_next]);
        end


    end
    lines{linecounter} = this_line;
    linecounter = linecounter+1;
end

end