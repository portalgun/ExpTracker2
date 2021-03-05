classdef dataTable < handle
properties
    DT
    dt
    curIndex=struct()
    index=struct()
end
properties(Hidden = true)
    colNames

    curColInd
    curColName
    curLabels
    curLabelInds

    keyInd

    nLabelVals
    labels
    labelNames
    labelInds
    S
end
methods
    function obj=dataTable(varargin)
        if numel(varargin) > 1 && ischar(varargin{1})
            obj.select_construct(varargin{:});
        end
    end
    function obj=select_construct(obj,from,varargin)
        switch from
            case 'CELLSTRUCT'
                obj.construct_from_cellStruct(varargin{:});
            otherwise
                error(['Create construct method for ' from ]);
        end
    end
    function obj=construct_from_cellStruct(obj,S,labelNames,labelVals)
        obj.labelNames=['flds' labelNames];
        obj.labels=[labelVals];


        flds=fieldnames(S);
        obj.keyInd=getKeyInd(S,flds);
        nFlds=numel(flds);

        %szs=zeros(length(flds),1);
        %for i = 1:nFlds
        %    fld=flds{i};
        %    szs(i)=size(S.(fld),1);
        %end

        %X=unique(szs);
       
        %if numel(X) == 1
        %    obj.keyInd=1;
        %else
        %    counts=hist(szs,X);
        %    obj.keyInd=szs(counts==max(counts));
        %end

        for i = 1:nFlds
            fld=flds{i};
            if ~isequal(size(S.(fld),1),obj.keyInd)
                obj.index.(fld)=S.(fld);
                S=rmfield(S,fld);
            end
        end
        obj.S=S;
        obj.expand_labels();
        obj.expand_cellStruct();
        % obj.expand_index XXX DONT USE
    end
    function obj=expand_labels(obj)
        if ~isequal(numel(obj.labels),numel(obj.nLabelVals))
            flds=fieldnames(obj.S);
            out=distribute(flds,obj.labels{:});
        else
            out=distribute(obj.labels{:});
        end
        out=distribute2chars(out);
        C=zeros(size(out));
        A=cell(1,size(out,2));
        for i = 1:size(out,2)
            [a,~,c]=unique(out(:,i));
            C(:,i)=c;
            A{i}=a;
        end
        obj.labels=A;
        obj.labelInds=C;
        obj.colNames=join(out);
        obj.nLabelVals=cellfun(@numel,obj.labels);
    end
    function obj=expand_cellStruct(obj)
        nFlds=numel(obj.labels{1});
        nrows=prod(obj.nLabelVals);
        nlbls=prod(obj.nLabelVals(2:end));
        obj.DT=zeros(obj.keyInd,nrows);

        for f=1:nFlds
            ind=(1:nlbls)+(nlbls*(f-1));
            fld=obj.labels{1}{f};
            obj.DT(:,ind)=reshape(obj.S.(fld),obj.keyInd,nlbls);
        end
        %obj.labelVals=[{fieldnames(obj.S)} obj.labelVals]
        % XXX
        %obj.S=[];
    end
    function obj=index_to_fld(obj,fld)
        if ismember(fld,obj.labels{1})
            warning([fld ' is already a field.']);
            return
        end
        rep=obj.keyInd/size(obj.index.(fld),1);
        repsz=[rep ones(1,ndims(obj.index.(fld))-1)];
        new=repmat(obj.index.(fld),repsz);
        obj.S.(fld)=new;
        obj.labels{1}{end+1}=fld;
        obj.expand_labels();
        obj.expand_cellStruct();
    end
    %function obj=expand_index(obj)
    %    nrows=prod(obj.nLabelVals(2:end));

    %    flds=fieldnames(obj.index)
    %    nFlds=numel(flds);
    %    for f =1:nFlds
    %        fld=flds{f};
    %        sz=size(obj.index.(fld),1);
    %        obj.index.(fld)=reshape(obj.index.(fld),sz,nrows);
    %    end
    %end
    function obj=select(obj,varargin)
        val=cell(size(varargin));
        key=cell(size(varargin));
        for i = 1:length(varargin)
            val{i}=varargin{i};
            ind= find(cellfun(@(x) ismember(val{i},x) ,obj.labels));
            key{i}=obj.labelNames{ind};
            % XXX CHECK FOR NAME TOO
            if numel(ind) > 1
                if isnumeric(val{i})
                    num2str(val{i});
                else
                    str=val{i};
                end
                error(['value ' str ' is ambiguous'])
            elseif isempty(key{i})
                if isnumeric(val{i})
                    num2str(val{i});
                else
                    str=val{i};
                end
                error(['value ' str ' does not exist ambiguous'])
            end
        end
        obj.select_by_label(key,val);
    end
    function obj=select_by_label(obj,varargin)
        if mod(numel(varargin),2)~=0
            error('requires keypairs')
        end

        %RESHAPE IF CELLS
        if iscell(varargin{1}) && iscell(varargin{2}) && length(varargin) == 2
            for i = 1:length(varargin{1})
                v=i*2-1;
                in{v}=varargin{1}{i};
                in{v+1}=varargin{2}{i};
            end
        else
            in=varargin;
        end

        COL=ones(size(obj.colNames));
        c=0;
        SUBS=zeros(1,numel(obj.nLabelVals));
        for i = 1:2:length(in)
            c=c+1;
            key=in{i};
            if strcmp(key,'fld')
                key='flds'
            end
            val=in{i+1};

            % XXX handle OR

            ind1=find(ismember(obj.labelNames,key));
            ind2=find(ismember(obj.labels{ind1},val));
            ind=obj.labelInds(:,ind1)==ind2;
            COL=ind & COL;

            if ind1 > 1
                SUBS(ind1)=ind2;
            end

        end

        % SELECT INDEX FLDS
        SUBS=strrep(num2strSane(SUBS),'0',':');
        obj.curIndex=[];
        str=['obj.curIndex.(fld)=obj.index.(fld)(' SUBS ');'];
        flds=fieldnames(obj.index);
        for i = 1:length(flds)
            fld=flds{i};
            eval(str);
        end

        obj.curColInd=COL;
        obj.curColName=obj.colNames(COL);
        obj.curLabelInds=obj.labelInds(COL,:);
        obj.dt=obj.DT(:,COL);

        obj.curLabels=cell(size(obj.labelInds,2),1);
        inds=obj.labelInds(COL,:);
        for i = 1:size(obj.labelInds,2)
            INDS=unique(obj.curLabelInds(:,i));
            LABELS=obj.labels{i};
            obj.curLabels{i}=LABELS(INDS);
        end

    end
    function obj=combine(obj,new)
    end
%% LIST
    function []=ls_flds(obj)
        obj.ls_fields;
    end
    function []=ls_fields(obj)
        for i = 1:length(obj.labels{1})
            disp(obj.labels{1}{i});
        end
    end
    function []=ls_subj(obj)
       ind=find(ismember(obj.labelNames,'subj'));
       if isempty(ind)
           disp('No subjects in data table')
           return
       end
       for i = 1:length(obj.labels{ind})
           disp(obj.labels{ind}{i});
       end
    end
    function []=ls_std(obj)
       ind=find(ismember(obj.labelNames,'std'));
       if isempty(ind)
           disp('No standards in data table')
           return
       end
       for i = 1:length(obj.labels{ind})
           disp(obj.labels{ind}{i});
       end
    end
    function []=ls_cur(obj)
        obj.print_helper(obj.curLabels);
    end
    function []=ls_all(obj)
        obj.print_helper(obj.labels);
    end
    function []=print_helper(obj,labels)
        in=cell(1,length(obj.labelNames));
        for i = 1:length(obj.labelNames)
            in{1,i}=vertcat(obj.labelNames{i},labels{i});
        end
        out=printColumn(in);
        [a,b]=strtok(out,newline);
        disp([newline a])
        disp([b newline])
    end
end
methods(Static=true)
    % DO IN ORDER
    function S=standardize_2IFC(S,bFlip)
        if ~exist('bFlip','var') || isempty(bFlip)
            bFlip=0;
        end
        dict=containers.Map();
        names={'responses','R';
               'Responses','R';
               'Rsp','R';
               'rsp','R';
               'RcmpChosen','RcmpChs';
               'RCmpChosen','RcmpChs';
               'Rchs','RcmpChs';
               'RChs','RcmpChs';
              };
        for i = 1:size(names,1)
            dict(names{i,1})=names{i,2};
        end
        for i = 1:numel(S)
            flds=fieldnames(S{i});
            for j=1:numel(flds)
                fld=flds{j};
                try
                    fix=dict(fld);
                catch
                    continue
                end
                if ~strcmp(fix,fld)
                    S{i}=fldRename(S{i},fld,fix);
                end
            end
            if bFlip
                nanind=isnan(S{i}.R);
                S{i}.R(nanind)=0;
                S{i}.R=double(~S{i}.R);
                S{i}.R(nanind)=nan;
            end
            if ~isfield(S{i},'RcmpChs')
                S{i}.RcmpChs=S{i}.cmpIntrvl==S{i}.R;
            end
        end
    end
    function S=fix_sizes(S)
        %% trl, std, subj, pass
        %% get maximum dimensions of each field
        flds=cellStructGetAllFlds(S);
        nS=numel(S);
        nF=numel(flds);
        ndims=zeros(nS,nF);
        sizes=cell(nS,nF);
        for i = 1:nS
        for j = 1:nF
            fld=flds{j};
            n=numel(size(S{i}.(fld)));
            ndims(i,j)=numel(size(S{i}.(fld)));
            sizes{i,j}=size(S{i}.(fld));
        end
        end
        maxd=max(ndims);

        % get max size of each field
        maxSizes=cell(size(maxd));
        for j = 1:length(maxd)
            ind=[ndims(:,j)==maxd(j)];
            maxSizes{j}=max(vertcat(sizes{ind,j}));
        end

        for i = 1:nS
        for j = 1:nF
            fld=flds{j};
            dd=maxd(j)-ndims(i,j);
            if dd ~=0
                % TODO
            end
            ds=maxSizes{j}-sizes{i,j};
            if any(ds~=0)
                ds(ds==0)=1;
                S{i}.(fld)=[S{i}.(fld); nan(ds)];
            end
        end
        end

    end
    function S=fix_shrink(S,indFld)
        %% for when stdx and cmpx get shrunk to 180
        if iscell(S)
            for i =1:numel(S)
                S{i}=shrink_fun(S{i},indFld);
            end
        else
            S=shrink_fun(S,indFld);
        end
        function S=shrink_fun(S,indFld)
            flds=fieldnames(S);
            N=size(S.(indFld),1);
            for f=1:length(flds)
                fld=flds{f};
                n=size(S.(fld),1);
                if ~isequal(n,N) && mod(N,n)==0
                    k=N/n;
                    S.(fld)=repmat(S.(fld),k,1,1);
                end
            end
        end
    end
    function out=flatten(S,indFld)
        out=dataTable.fix_shrink(S,indFld);
        out=dataTable.sort_by_ind(out,indFld);
        out=flatten_fun(out);

        function Snew=flatten_fun(S)
            Snew=S{1};
            flds=fieldnames(Snew);
            for f=1:length(flds)
                fld=flds{f};
                for i = 2:length(out)
                    Snew.(fld)=cat(4,Snew.(fld),S{i}.(fld));
                end
            end
        end
    end
%% START UNTIL
    function out=pass_sort(S,indFld,bRm)
        % use cmpIind
        % RUN LAST IF POSSIBLE
        if ~exist('bRm','var') || isempty(bRm)
            bRm=0;
        end
        out=dataTable.fix_shrink(S,indFld);
        out=dataTable.sort_by_ind(out,indFld);
        if bRm
            out=rmNanRows(out,indFld);
        end
        function Snew=combine_fun(S)
            flds=fieldnames(S{1});
            Snew=S{1};
            for f= 1:length(flds)
                fld=flds{f};
                for i = 2:numel(S)
                    Snew.(fld)=cat(4,Snew.(fld),S{i}.(fld));
                end
            end
        end
        function S=nanRows(S,indFld)
            nStd=size(S{1}.(indFld),2);
            n=size(S{1}.(indFld),1);
            ind=false(n,nStd);
            for s=1:nStd
            for i=1:numel(S)
                k=any(isnan(S{i}.(indFld)(:,s,:)),3);
                ind(:,s)=ind(:,s) | k;
            end
            end
            for i = 1:numel(S)
                flds=fieldnames(S{i});
                for f=  1:length(flds)
                    fld=flds{f};
                    for s = 1:nStd
                        S{i}.(fld)(ind(:,s),s,:)=nan;
                    end
                end
            end
        end
        function S=rmNanRows(S,indFld)
            inds=cellfun(@(x) isnan(x.(indFld)),S,UO,false);
            ind=horzcat(inds{:});
            ind=~any(ind,2);
            for i = 1:length(S)
                S{i}=structfun(@(x) x(ind),S{i},UO,false);
            end
        end
        function Snew=combineStds(S,indFld)
            nStd=numel(S{1}.(indFld));
            nSubj=size(S{1}.(indFld){1},2);
            nExp=numel(S);

            flds=cellStructGetSameFlds(S);
            Snew=cell(nStd,1);
            for s=1:nStd
                Snew{s}=struct();
                for f=1:length(flds)
                    fld=flds{f};
                    n=size(S{1}.(fld){s},1);
                    Snew{s}.(fld)=zeros(n,nSubj,nExp);
                end
            end

            for s= 1:nStd
            for i=1:nExp
            for f = 1:length(flds)
                fld=flds{f};

                Snew{s}.(fld)(:,:,i)=S{i}.(fld){s};
            end
            end
            end
        end

    end
    function Snew=sort_by_ind(S,indFld)
        % sorts fields by indFld between experiments, filling in nans where indFld value is not present
        nStd=size(S{1}.(indFld),2);
        nSubj=size(S{1}.(indFld),3);

        % get all values
        vals=cell(nStd,1);
        for i = 1:numel(S)
        for s = 1:nStd
            vals{s}=[vals{s}; S{i}.(indFld)(:,s,:)];
        end
        end
        for s=1:nStd
            vals{s}=unique(vals{s});
        end

        counts=cell(nStd);
        new=cell(nStd,numel(S));
        IND=cell(nStd,numel(S));
        INDB=cell(nStd,numel(S));
        for s = 1:nStd
            VAL=vals{s};
            counts{s}=zeros(length(VAL),numel(S));
            for i = 1:numel(S)

                val=S{i}.(indFld)(:,s,:);
                counts{s}(:,i)=histc(val,VAL);
                % TODO expand counts
            %end
                new{s,i}=repmat(VAL ,1,numel(S)) .* (counts{s} > 0);
                new{s,i}(new{s,i}==0)=nan;
                tmp=cell2mat(arrayfun(@(x) find(x==val),VAL,UO,false));
                IND{s,i}=tmp;
                %INDB{s,i}=cell2mat(arrayfun(@(x) find(x==VAL),val,UO,false));
            end
        end

        Snew=cell(nStd,1);
        size(Snew);
        flds=cellStructGetSameFlds(S);
        for s = 1:nStd
            N=nan(size(vals{s},1),nStd,nSubj);
            Snew{s}=struct();
            for f = 1:numel(flds)
                fld=flds{f};
                for i = 1:numel(S)
                    Snew{i,1}.(fld)=N;
                    val=S{i}.(fld)(IND{s,i},:,:);
                    Snew{i,1}.(fld)(IND{s,i},:,:)=val;
                end
            end
        end

    end
%% END UNTIL
    function S=merge_at_pass(S)
    % input cellstruct

        [dims,FLDS]=structGetCatDims(S{:});
        d=4-dims; % TODO?
        dims=dims+d;
        S=structMergeFldsAtDims(FLDS,dims,S{:});

    end
end
end
