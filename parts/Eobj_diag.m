classdef Eobj_diag < handle
methods
%% RAW DIAGNOSTICS
    function [] = run_diagnostics(obj)
        % TODAY
        obj.check_data_dir_orphan();
        obj.check_block_parity_all();
        %obj.check_data_all();
    end
%% ORHPAN
    function [] = check_data_dir_orphan(obj)
        [oDir,oDB]=obj.get_data_dir_orphan();
        if ~isempty(oDir)
            disp('FILES NOT IN DB:')
            disp(oDir)
        else
            disp('NO ORPHANED FILES IN DIRECTORY')
        end
        if ~isempty(oDB)
            disp('DB LABELS WITHOUT FILE')
            disp(oDB)
        else
            disp('NO ORPHANED DB LABELS')
        end
    end
    function [oDir,oDB] = get_data_dir_orphan(obj)
        ddir=obj.get_all_matching_raw_files_dir();
        ddb =obj.get_raw_files_all(1);
        oDir=setdiff(ddir,ddb); % orphan files
        oDB=setdiff(ddb,ddir);   % orhpan db string
    end
    %%
    function [] = check_block_parity_all(obj)
        [noflag,nodata]=obj.get_block_data_parity_all();
        if ~isempty(noflag)
            disp('DATA FILE EXISTS, BUT BLOCK INDICATES INCOMPLETE')
            disp(noflag);
        end
        if ~isempty(nodata)
            disp('BLOCK INDICATES COMPLETE, BUT NO DATA FILE')
            disp(nodata);
        end
        if isempty(nodata) && isempty(noflag)
            disp('ALL COMPLETION BLOCKS MATCH FILES')
        end
    end
%% PARITY
    function [noflag,nodata] = get_block_data_parity_all(obj)
        IND=obj.get_all_data_combs();
        bad=zeros(size(IND,1),1);
        for i = 1:size(IND,1)
            ind=IND(i,:);
            bad(i)=obj.check_block_data_parity(ind{:});
        end
        noflag=IND(bad==2,:);
        nodata=IND(bad==3,:);
    end
    function out=check_block_data_parity(obj,subj,method,std,blk)
        bData=~isempty(obj.rawData.(subj).(method).(std).data{blk});
        bFlag=logical(obj.rawData.(subj).(method).(std).blk(blk));
        if ((bData && bFlag) || ~(bData && bFlag))
            out=1;
        elseif (bData && ~bFlag)
            out=2;
        elseif (~bData && bFlag)
            out=3;
        end
    end
    function files = get_all_matching_raw_files_dir(obj)
        rex=[obj.name '.*\.mat'];
        files=matchingFilesInDir(obj.dir.RAW,rex);
        files=cellfun(@filepartsfun,files,'UniformOutput',0);
        function out=filepartsfun(files)
            [~,out]=filePartsSane(files);
        end
    end
%% TODO SORT
    function [out,fname]=check_block(obj,subj,method,std,blk,bSkipCheck)
        if ~exist('bSkipCheck','var') || isempty(bSkipCheck)
            bSkipCheck=1;
        end
        if bSkipCheck
            exitflag=obj.check_raw_dir();
            if exitflag
                out=0;
                return
            end
        end
        fname=obj.get_block_fname(subj,method,std,blk);
        out=exist(fname,'file')==2;
    end
    function out=check_block_all(obj,bSkipCheck)
        if ~exist('bSkipCheck','var') || isempty(bSkipCheck)
            bSkipCheck=0;
        end
        obj.init_raw_empty_all(obj);
        IND=obj.get_all_data_combs();
        if ~bSkipSheck
            out=obj.checkDataDir();
        else
            out=1;
        end
        if out==0
            error('DataDir does not exist')
        end

        for i = 1:length(IND)
            ind=IND(i,:);
            out=obj.check_block(ind{:},1);
        end
    end
%% CHECK
    function out=check_prjDir(obj)
        % XXX  new px function
    end
    function out = check_redo(obj)
       % TODO
    end
end
end
