classdef dataTable < handle
properties
    M
    X
    cCrit

    fM
    fX

    cmpX
    stdX
    RcmpChs
    color
    shape

    keySubj
    keyExp
    keyPrj
    keyPass
    nSubj
    nPrj
    nExp
    nPass

end

properties(Hidden=true)
    i
    E
    S

    n
    m

    nX
    nM
    MM

    crit
end
methods
    function obj=data_Table(Mall,cCrit,obj)
        % XXX AFTER COMBINING BLOCKS, SORTING, PRUNING SAVING
        obj.keySubj = containers.Map;
        obj.keyExp  = containers.Map;
        obj.keyPrj  = containers.Map;
        obj.keyPass = containers.Map;
        obj.nSubj   =0;
        obj.nPrj    =0;
        obj.nExp    =0;
        obj.nPass   =0;
        obj.Mall=sortrows(Mall);
        obj.expand_cCrit(cCrit);
    end
%% Construct
    function obj=get_mats();
        for i = 1:size(obj.Mall,1)
            obj.i=i;
            get_cols();
        end
        obj.expand_meta_data();
        obj.convert_cCrit();
    end
    function obj=get_cols(obj)
        obj.get_meta_data();
        obj.loadE();
        obj.loadS();
        obj.get_new_data();
        obj.append_meta_data();
        obj.append_data();
    end
%%
    function obj=get_meta_data(obj)
        obj.nMprev=obj.nM;
        obj.nM=obj.Mall(obj.i,:);
    end
    function out=check_load(obj)
        if ~strcmp(obj.nMprev{1},nM{1}) || ~strcmp(obj.nMprev{1},nM{1})
            out=1;
        else
            out=0;
        end
    end
    function obj=load_E(obj)
        if obj.checkload();
            obj.E=loadExpStruct(obj.nM{1},obj.nM{2});
        end
    end
    function obj=load_S(obj)
       obj.E= ;% XXX
    end
    function obj=get_new_data(obj)
        obj.nX=[obj.S.stdX obj.S.cmpX obj.S.RcmpChs];
        n=size(nX,1);
        m=size(nX,2);
        reshape(obj,nX,n,1,m);
    end
    function obj=append_meta_data(obj)
        if ~ismember(obj.nM{1},keys(obj.keySubj))
            obj.nSubj=obj.Subj+1;
            obj.keySubj(obj.nM{3})   = obj.nSubj;
        end
        if ~ismember(obj.nM{2},keys(obj.keyPrj))
            obj.nPrj=obj.nPrj+1;
            obj.keyPrj(obj.nM{3})    = obj.nPrj;
        end
        if ~ismember(obj.nM{3},keys(obj.keyExp))
            obj.nExp=obj.nExp+1;
            obj.keyExp(obj.nM{3})    = obj.nExp;
        end
        if ~ismember(obj.nM{4},keys(obj.keyPass))
            obj.nPass=obj.nPass+1;
            obj.keyPass(obj.nM{4})     = obj.nPass;
        end

        s =obj.keySubj(S.subj);
        p =obj.keyPrj(S.Prj);
        e =obj.keyExp(S.Exp);
        ps=obj.kenPass(S.Pass);

        obj.M(1,end+1,:)=[s e p ps];
    end
    function obj=append_data(obj)
        obj.n=size(obj.X,1);
        obj.m=size(obj.X,2);
        obj.o=size(obj.X,3);

        n=size(obj.nX,1);
        o=size(obj.nX,3);

        % add nan to obj.X
        if n > obj.n
            nn=n-obj.n;
            Na=nan(nn,obj.m,obj.o);
            obj.X=[obj.X; Na];
        end

        % add nan to nX
        if n < obj.n
            nn=obj.n-nN;
            Na=nan(nn,1,m);
            obj.nX=[obj.nX; Na];
        end

        obj.X=[obj.X; obj.nX];
    end
%%
    function obj=expand_meta_data(obj)
        n=size(obj.X,1);
        obj.MM=repmat(obj.M,obj.n,1,1);
    end
    function obj=convert_cCrit(obj)
        obj.crit=cell(size(obj.Crit,1));
        for j = 1:size(obj.cCrit,1)
            c=obj.cCrit(j,:);
            key=c{1};
            val=c{2};
            if ~iscell(val)
                val={val};
            end
            switch key
                case {'subj','SUBJ','subjs','SUBJS','Subj','Subjs'}
                    KEY=obj.keySubj;
                    IND=1;
                case {'exp','EXP','Exp','expID'};
                    KEY=obj.keyExp;
                    IND=2;
                case {'prj','PRJ','proj','PROJ','Prj','Proj','prjCode'}
                    KEY=obj.keyPrj;
                    IND=3;
                case {'pass','PASS','passes','PASSES','Pass','Passes'}
                    KEY=obj.keyPass;
                    IND=4;
            end

            VV=[];
            for k = 1:numel(val)
                if strcmp(val{k},'all')
                    V=values(KEY);
                    V=[V{:}];
                else
                    V=KEY('val');
                end
                VV=[VV, V];
            end
            obj.crit(j,:)=[IND,{VV}];
        end
    end
    function obj=combine(obj)
        % XXX combine more specific conditions
        obj.fX=[];
        Used=zeros(size(obj.X),2);
        for j = 1:size(obj.crit,1)
            crit=obj.crit(j,:);
            ind=crit{1};
            VAL=crit{2};

            %GET COLUMNS TO COMBINE
            IND=zeros(size(obj.MM,1),size(obj.MM,2));
            for k = 1:numel(VAL)
                val=VAL{k};
                IND=IND | obj.MM(:,:,ind)==val;
            end
            IND=repmat(IND,1,1,3);

            % COMBINE COLUMNS
            col=obj.X(IND);
            col=col(:);

            % APPEND COLUMN TO fX
            obj.fX=appendColIsz(obj.fX,col)

            % XXX index
            obj.mF;

            USED=USED | IND(:,1);
        end

        % SELECT REST OF M
        used=repmat(USED,1,1,size(obj.M,3));
        M=obj.M(~used);

        % APPEND REST OF M TO MF
        obj.mF=[obj.mF M];


        % Select REST OF X
        used=repmat(USED,size(obj.X,1),1);
        X=obj.X(~used);

        % APPEND NANS to REST OF X
        nN=size(obj.fX))-size(X,1);
        Na=nan(nN,size(X,2));
        X=[X; Na];

        % APPEND REST OF X TO fX
        obj.fX=[obj.fX X];
    end
end
end

