classdef EAobj_inds < handle
methods
    function ind=get_ind(obj,subj,mode,std,blk,rule)
        if ~exist('mode','var') || isempty('mode')
            mode=obj.subjStatus.(subj).status;
        end
        if ~exist('std','var')
            std=[];
        end
        if ~exist('blk','var')
            blk=[];
        end
        if ~exist('rule','var') || isempty(rule)
            rule='rand_all';
        end

        ind=obj.(['get_ind_' rule])(subj,mode,std,blk);
    end
%%
    function Iall=get_ind_rand_all(obj,subj,mode,std,blk)

        [Iall,Cinds,~]=obj.get_completion_inds(subj,mode,std,blk);

        ind=randperm(size(Iall,1));
        Cinds=logical(Cinds(ind));
        Iall=Iall(ind,:);
        Iall(Cinds,:)=[];
    end
    function Iall=get_ind_rand_incomplete_all(obj,subj,mode)
        [Iall,Cinds,~]=get_completion_inds(obj,subj,mode)
        Iall(Cinds)=[];
        ind=randperm(1:size(Iall,1));
        Iall=Iall(ind,:)
    end
    function Iall=get_ind_rand_low_pass_first(obj,subj,mode)
        % TODO
    end
%%
    function ind=get_same_pass_ind(obj)
        n=length(obj.Eall);
        names=cell(n,1);
        for i = 1:n
            names{i}=sed('s',obj.Eall{i}.name,'_pass[0-9]+$','');
        end
        [~,~,ind]=unique(names);
    end
    function [Iall,CINDS,Call]=get_completion_inds(obj,subj,mode,std,blk)
    %function [INDS,Iall]=get_completion_inds(obj,subj,mode,sortOrder)
    %   INDS = completion std, blk, pass_group, pass
        %IND = completion
        %INDS = completion, EAllind std, blk, pass_group, pass
        psInds=obj.get_same_pass_ind();
        n=length(obj.Eall);
        IND=cell(n,1);
        for i = 1:n
            [IND{i},inds]=obj.Eall{i}.get_std_block_ind_completion(subj,mode,std,blk);

            I=repmat(i,size(inds,1),1);
            a=repmat(psInds(i),size(inds,1),1);
            b=repmat(obj.Eall{i}.pass,size(inds,1),1);
            INDS{i}=[I inds a b];
        end
        Iall=vertcat(INDS{:});
        CINDS=vertcat(IND{:});
        [C,Call]=obj.get_basic_completions(Iall,CINDS);

    end
    function [C,Call]=get_basic_completions(obj,Iall,CINDS);
        %exp competion
        %std cmpletion
        %blk completion
        %pass completion
        %Iall(CINDS,
        C=struct();
        ALL=ones(size(Iall,1),1);
        C.ExpC=accumarray(Iall(:,1),CINDS)./accumarray(Iall(:,1),ALL);
        C.stdC=accumarray(Iall(:,2),CINDS)./accumarray(Iall(:,2),ALL);
        C.blkC=accumarray(Iall(:,3),CINDS)./accumarray(Iall(:,3),ALL);
        C.AcPassC=accumarray(Iall(:,4),CINDS)./accumarray(Iall(:,4),ALL);
        C.PassC=accumarray(Iall(:,5),CINDS)./accumarray(Iall(:,5),ALL);

        Call=zeros(size(Iall));
        flds=fieldnames(C);
        for c = 1:numel(flds)
            fld=flds{c};
            for r = 1:numel(C.(fld))
                Call(Iall(:,c)==r,c)=C.(fld)(r);
            end
        end

    end
    function [IallNew,labels,usInd]=get_permuted_completion_inds(obj,subj,mode,indPerm)
    % indPerm
    % 1=Eind
    % 2=std
    % 3=blk
    % 4=passgroup
    % 5=pass
        [~,si]=sort(indPerm);
        usInd=sort(si);
        labels={'Eind','std','blk','passGroup','pass'}
        [Inds,Cind]=obj.get_sort_completion_inds(subj,mode)
        indPerm=rowVec(indPerm);
        IallAll=Iall(:,indPerm)
    end
    function [Inds]=get_sorted_completion_inds(obj,subj,mode,rule)
        if ~exist('rule','var') || ~isempty(rule)
            rule='none';
        end
        switch rule
            case 'none'
                indPerm=1:6
        end

        %sort columns
        [Iall,labels,usInd]=obj.get_permuted_completion_inds(subj,mode,indPerm,labels);

        %sort rows
        Inds=sortrows(Iall);

        %unsort columns
        Inds=Inds(:,usInd);


    end
end
end
