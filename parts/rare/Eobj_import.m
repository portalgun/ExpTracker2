classdef Eobj_import < handle
methods
    function obj=import_exp_in_all(obj)
        % ATTEMPTS TO ADOPT EXP DATA (stim) IF IT EXISTS
        INDS=distribute(obj.modeflds,1:obj.nStd,1:obj.nBlk);
        for i = 1:size(INDS,1)
            ind=INDS(i,:);
            obj.import_exp_in(ind{:});
        end
    end
    function obj=import_exp_in(obj,mode,std,blk)
        % ATTEMPTS TO LOAD IN EXP DATA
        fname=obj.gen_fname_exp(mode,std,blk,0);
        if ~exist(fname,'file')
            fname=obj.gen_fname_exp(mode,std,blk,1);
        end
        if ~exist(fname,'file')
            return
        end
        name=obj.gen_name_exp(mode,std,blk);
        mode=obj.auto_mode(mode);
        std=obj.auto_std_fld(std);
        obj.expData.(mode).(std){blk}=name;
    end
%%
    function merge_exp(obj)
        % ???
        %fldsrm={'stdXunq','cmpXunq'};
        mthds={'test','train','pilot'};
        stds=obj.get_std_flds();
        blks=1:obj.nBlk;

        %%% across subjects
        INDS=distribute(mthds,stds,blks);

        %% across subjects and standards
        %INDS=distribute(mthds,blks);
        %INDS=[INDS(:,1) [repmat({[]},size(INDS,1),1)] INDS(:,2)]

        for i = 1:size(INDS,1)
            ind=INDS(i,:);
            names=obj.get_raw_names([],ind{:});
            names=cellfun(@(x) [obj.dir.EXP strrep(x,'raw_','exp_') '.mat'],names,'UniformOutput',false);
            if any(cellfun(@isempty,names))
                continue
            end

            SS=cell(length(names),1);
            for j=1:length(names)
                continueflag=0;
                try
                    S=load(names{j});
                catch
                    disp(names{j})
                    continueflag=1;
                    break
                end
                if isstruct(S) && nflds(S)==1
                    flds=fieldnames(S);
                    S=S.(flds{1});
                end

                SS{j}=S;
            end
            if continueflag==1
                continue
            end

            if numel(names) == 1
                b=1;
            else
                bIndFldsAll=cellStructCmp(SS);
                b=all(bIndFldsAll);
            end
            if b
                name=names{1};
                parts=strsplit(name,'_');
                subj=[und parts{end-2}];
                blk=[und parts{end}];
                [~,nblk]=strtok(blk,'-');
                name=strrep(name,subj,'');

                expName=[strrep(name,nblk,'') '.mat'];
                exp=SS{1};
                save(expName,'exp','-mat');
                for k=1:length(names)
                    delete(names{k});
                end
            end
        end
    end
end
end
