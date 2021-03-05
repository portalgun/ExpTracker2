classdef Eobj_wrangle < handle
%% WRANGLE
    % COMBINE
    % CLEAN
    % MERGE
    % PRUNE
    % TABLE
methods

    function [obj,DT]=wrangle_all(obj,bSave,bNested,bOverwrite)
        if ~exist('bOverwrite','var')
            bOverwrite=0;
        end
       % XXX ADD MODE COMPATibility
        % XXX ADD FORCE option
        % XXX ADD LOADING IF FILE EXISTS
        if ~exist('bNested','var') || isempty(bNested)
            bNested=0;
        end
        INDS=distribute(1:obj.nSubj,1:obj.nStd);
        n=size(INDS,1);

        CM=obj.block_combine_all(bNested);
        Cl=obj.block_clean(CM,bSave,bOverwrite);    clear CM;
        assignin('base','Cl',Cl);

        for i = 1:n
            break % XXX
            subj=INDS(i,1);
            std=INDS(i,2);
            if nflds(Cl{std,subj})==0
                continue
            end
            % TODO MODE HERE
            [Pr,~]=obj.block_prune(Cl{std,subj},bSave,bOverwrite,subj,1,std);
            % BROKEN XXX
            %assignin('base','Pr',Pr);
            %DT=obj.block_table(Pr,bSave,bOverwrite,subj,1,std);
        end

        % MERGE
        MR=obj.block_merge(Cl,bSave,bOverwrite);    clear Cl;
        assignin('base','MR',MR);

        % PRUNE
        [Pr,~]=obj.block_prune(MR,bSave,bOverwrite); clear MR;
        assignin('base','Pr',Pr);

        % TABLE
        % XXX Not working, need to handle struct instead of cell struct
        DT=obj.block_table(Pr,bSave,bOverwrite);
        assignin('base','DT',DT);

        obj.save();
    end
    function [CM]=block_combine_all(obj,bNested)
        if ~exist('bNested','var') || isempty(bNested)
            bNested=0;
        end

        CM=cell(obj.nStd,numel(obj.subjs));
        if ~bNested; p=pr(obj.nSubj,[],'loading blocks','...'); end

        for s=1:obj.nSubj
            if ~bNested; p.u(); end
            subj=obj.subjs{s};
            Cm=obj.block_combine_stds(subj,1);
            CM(:,s)=Cm;
        end
        if ~bNested; p.c(); end
    end
    function [Cm]=block_combine_stds(obj,subj,bNested)
        if ~exist('bNested','var') || isempty(bNested)
            bNested=0;
        end
        Cm=cell(obj.nStd,1);
        if ~bNested; p=pr(obj.nStd,[],'loading blocks','...'); end
        for std = 1:obj.nStd
            if ~bNested; p.u(); end
           Cm{std}=obj.block_combine(subj,std,1);
        end
        if ~bNested; p.c(); end
    end
    function [cm,INDS,Fnames]=block_combine(obj,subj,std,bNested)
        if ~exist('bNested','var') || isempty(bNested)
            bNested=0;
        end
        [Sall,Fnames,INDS]=obj.load_raw_test_blocks(subj,std,[],bNested);
        ind=cellfun(@isempty,Fnames);
        Sall=obj.wrangle_raw_objects(Sall);
               

        % ADD BLOCK NUMBER
        for i = 1:numel(Sall)
            if ind(i)
                continue
            end
            if isfield(Sall{i},'S')
                n=size(Sall{i}.S.stdX);
            else
                n=size(Sall{i}.stdX);
            end
            Sall{i}.block=ones(n)*INDS{i,4};
        end

        Fnames(ind)=[];
        INDS(ind)=[];
        Sall(ind)=[];
        if ~isempty(Sall)
            cm=cellStructMerge(Sall,0,1);
        else
            cm=struct();
        end
    end
    function Sall=wrangle_raw_objects(obj,Sall)
        for i = 1:numel(Sall)
            if isempty(Sall{i})
                continue
            end
            if isobject(Sall{i})
                Sall{i}=obj2struct(Sall{i},1);
                Sall{i}=obj.fix_RSP(Sall{i}); % UNCOMMENT TO FIX RSP % XXX
            end
        end
    end
    function S=fix_RSP(obj,S)
    % NOTE TEMPORARY FIX
    % remove columns that are just zero
        if isfield(S.RSP,'responses')
            S.RSP.responses( :, ~any(S.RSP.responses,1) ) = [];
        end
        disp(num2str(size(S.RSP.responses)))
    end
    function Cl=block_clean(obj,CM,bSave,bOverwrite)
        if ~exist('bOverwrite','var')
            bOverwrite=0;
        end
        if ~exist('CM','var') || isempty(CM)
            % TODO
        end

        if strcmp(obj.method,'2IFC')
            keyFld='cmpIind';
        else
            error('write code');
        end
        assignin('base','CM',CM)
        out=cellStructCleanup(CM,keyFld);
        Cl=out.S;

        if ~bSave
            return
        end

        assignin('base','Cl',Cl);
        obj.save_data(Cl,'OUT','cleaned',[],[],[],bOverwrite);

        for i = 1:obj.nSubj
        for j = 1:obj.nStd
            if nflds(Cl{j,i})==0
                continue
            end
            obj.save_data(Cl{j,i},'OUT','cleaned_ind',i,1,j,bOverwrite);
        end
        end
    end
    function MR=block_merge(obj,Cl,bSave,bOverwrite)
        if ~exist('Cl','var') || isempty(Cl)
        end
        MR=cellStructMerge(Cl,0,0);
        if ~bSave
            return
        end
        obj.save_data(MR,'OUT','merged',[],[],[],bOverwrite);
    end
    function [Pr,Pp]=block_prune(obj,MR,bSave,bOverwrite,subj,mode,std)
        if ~exist('bOverwrite','var')
            bOverwrite=0;
        end
        if ~exist('MR','var') || isempty(MR)
            % TODO
        end
        flds2IFC={'block','R','RcmpChosen','Rcorrect','cmpIind','stdIind','stdXind','cmpXind','cmpIntrvl','stdIntrvl','stdX','cmpX','responses'};

        Pp=struct();
        Pr=MR;
        clear MR

        if isfield(Pr,'RSP')
            Pr=obj.flatten_psycho(Pr);
        end
        flds=fieldnames(Pr);

        for i = 1:length(flds)
            fld=flds{i};

            sz=size(Pr.(fld));
            dims=num2strSane([obj.nStd numel(obj.subjs)]);
            c1=~(ismember(obj.nTrl,sz) | ismember(obj.nTrlPerBlk,sz));
            c2=~contains(num2strSane(sz),dims);
            c3=strcmp(obj.method,'2IFC') && ~ismember(fld,flds2IFC);
            if strcmp(obj.method,'2IFC') && ismember(fld,flds2IFC)
                continue
            elseif c3 || c1 || c2
                Pp.(fld)=Pr.(fld);
                Pr=rmfield(Pr,fld);
            end
        end
        if isfield(Pr,'cmpIind')
            Pr=dataTable.fix_shrink(Pr,'cmpIind');
        end

        if ~bSave
            return
        end
        if nargin > 5
            obj.save_data(Pr,'OUT','pruned_ind',subj,mode,std,bOverwrite);
            obj.save_data(Pp,'OUT','prunedP_ind',subj,mode,std,bOverwrite);
        else
            obj.save_data(Pr,'OUT','pruned',[],[],[],bOverwrite);
            obj.save_data(Pp,'OUT','prunedP',[],[],[],bOverwrite);
        end
    end
    function PrNew=flatten_psycho(obj,Pr)

        PrNew=struct();
        flds=fieldnames(Pr);
        psychoFlds={'RSP','S'};
        for i = 1:length(flds)
            fld=flds{i};
            if ismember(fld,psychoFlds)
                flds2=fieldnames(Pr.(fld));
                for j = 1:length(flds2)
                    fld2=flds2{j};
                    PrNew.(fld2)=Pr.(fld).(fld2);
                end
            end
        end
    end
    function DT =block_table(obj,Pr,bSave,bOverwrite,subj,mode,std)
        if ~exist('bOverwrite','var')
            bOverwrite=0;
        end
        labelVals={obj.subjs,obj.get_std_nums};
        labels={'subj','std'};
        DT=cellStruct2dataTable(Pr,labels,labelVals);
        if ~bSave
            return
        end
        if nargin> 6
            obj.save_data(DT,'OUT','tabled_ind',subj,mode,std,bOverwrite);
        else
            obj.save_data(DT,'OUT','tabled',[],[],[],bOverwrite);
        end
    end
end
end
