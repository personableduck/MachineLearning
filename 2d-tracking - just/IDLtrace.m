% re-generate the tracks in 3D tracking data with IDLtrack.m

clear;
[list_mat, path_mat] = uigetfile({'*.mat'},...
    'Select the data file to CONVERT','MultiSelect','on');
if path_mat == 0    % restore path and fname if cancelled
    return;
end
if iscell(list_mat) == 0    % single file
    list_mat = {list_mat};
end
n_mat = length(list_mat);

for i_mat = 1:n_mat
    fprintf('\nLoading ROI %g/%g...\n\n',i_mat,n_mat)
    fname_mat = list_mat{i_mat}; 

    load([path_mat,fname_mat]);

    n_frame = length(xyz_frame);
    xyzs = [];
    for i_f = 1:n_frame
        if isempty(xyz_frame{i_f})==0
            xyzs_frame = xyz_frame{i_f};
            xyzs_frame(:,4) = i_frame_all(i_f);      % use frame as the time unit
            xyzs = [xyzs;xyzs_frame];
        end
    end
    xyzs = xyzs(isnan(xyzs(:,1))==0,:);

    commandwindow;
    maxdisp = speed_max*dt/1000;
    param = struct('mem',10,'good',0,'dim',3,'quiet',0);
    xyzst = IDLtrack(xyzs,maxdisp,param);

    % convert xyzst to tracksSimple
    n_track = xyzst(end,5);
    tracksSimple = repmat(struct('coor',[],'start_frame',[],'end_frame',[]),[n_track,1]);
    for i_track = 1:n_track
        xyzs_track = xyzst(xyzst(:,5)==i_track,1:4);
        f_track = xyzs_track(:,4);
        tracksSimple(i_track).start_frame = min(f_track);
        tracksSimple(i_track).end_frame = max(f_track);
        tracksSimple(i_track).coor = nan(max(f_track)-min(f_track)+1,3);
        for i_f = 1:max(f_track)-min(f_track)+1
            coor = xyzs_track(f_track==(min(f_track)+i_f-1),1:3);
            if ~isempty(coor)
                tracksSimple(i_track).coor(i_f,1:3) = coor;
            end
        end
    end

    save([path_mat,fname_mat],'tracksSimple','n_track','-append');

    if exist('tracksFinal','var')
        tracksFinal=[];
        save([path_mat,fname_mat],'tracksFinal','-append');
    end
end

fprintf('Done\n');