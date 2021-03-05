classdef Eobj_fix < handle
methods(Hidden=true)
    function obj=fix_single_cell_fnames(obj)
        dflds=fieldnames(obj.fnames);
        for d = 1:length(dflds)
            dfld=dflds{d};
            flds=fieldnames(obj.fnames.(dfld))
            for f = 1:length(flds)
                fld=flds{f};
                val=obj.fnames.(dfld).(fld);
                if iscell(val) && numel(val) ==1
                    obj.fnames.(dfld).(fld)=val{1};
                end

            end
        end
    end
    function obj=fix_redo_from_file_all(obj)
        modes=obj.modeflds;
        stds=1:obj.nStd;
        blks=1:obj.nBlk;
        subjs=1:obj.nSubj;
        IND=distribute(subjs,modes,stds,blks);
        for i = 1:size(IND,1)
            ind=IND(i,:);
            obj.fix_redo_from_file(ind{:});
        end
    end
    function obj=fix_redo_from_file(obj,subj,mode,std,blk)
        subj=obj.auto_subj(subj);
        mode=obj.auto_mode(mode);
        std=obj.auto_std_fld(std);
        blk=obj.auto_blk_num(blk);

        fname=obj.rawData.(subj).(mode).(std).data{blk};
        if isempty(fname)
            obj.rawData.(subj).(mode).(std).redo(blk)=0;
            return
        end
        [r1,~]=obj.get_r1_r2_from_fname(fname);
        obj.rawData.(subj).(mode).(std).redo(blk)=r1;

    end
%% NAMES
    function obj=fix_exp_names_all(obj)
        mthds=obj.modeflds;
        stds=1:obj.nStd;
        blks=1:obj.nBlk;
        IND=distribute(mthds,stds,blks);
        fnames=cell(size(IND,1),1);
        for i = 1:size(IND,1)
            ind=IND(i,:);
            obj.fix_exp_name(ind{:});
        end
        fnames(cellfun(@isempty,fnames))=[];
    end
    function obj=fix_exp_name(obj,mode,std,blk)
        mode=obj.auto_mode(mode);
        std=obj.auto_std_fld(std);
        blk=obj.auto_blk_num(blk);
        name=obj.expData.(mode).(std){blk};
        if isempty(name)
            return
        end

        % NOTE
        [~,name]=fileparts(name);
        obj.expData.(mode).(std){blk}=name;

    end
    function obj=fix_raw_names_all(obj)
        mthds=obj.modeflds;
        subjs=obj.subjs;
        stds=1:obj.nStd;
        blks=1:obj.nBlk;
        IND=distribute(subjs,mthds,stds,blks);
        fnames=cell(size(IND,1),1);
        for i = 1:size(IND,1)
            ind=IND(i,:);
            obj.fix_raw_name(ind{:});
        end
        fnames(cellfun(@isempty,fnames))=[];
    end
    function obj=fix_raw_name(obj,subj,method,std,blk)
        subj=obj.auto_subj(subj);
        method=obj.auto_mode(method);
        std=obj.auto_std_fld(std);
        blk=obj.auto_blk_num(blk);
        name=obj.rawData.(subj).(method).(std).data{blk};
        if isempty(name)
            return
        end

        % NOTE
        name=strrep(name,'.mat','');
        if startsWith(name,'RAW_')
            error('write code')
        elseif ~startsWith(name,'raw_')
            name=['raw_' name];
        end

        name
        obj.rawData.(subj).(method).(std).data{blk}=name;
    end
%%
end
end
