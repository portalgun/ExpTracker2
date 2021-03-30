classdef exp_details < handle;
properties
    N

    Aname
    AprjCode
    AimgDTB
    AnatORflt
    AimgDim
    Amethod
    AprjInd;
    Asubjs;
    Apass;
    AstdXunqAll;
    AcmpXunqAll;
    AXunits;
    AXname;
    bname;
    bprjCode;
    bimgDTB;
    bnatORflt;
    bimgDim;
    bmethod;
    bprjInd;
    bsubjs;
    bpass;
    bstdXunqAll;
    bcmpXunqAll;
    bXunits;
    bXname;
    bnCmp;
    bnStd;
    AnCmp;
    AnStd;

    prjIndS;
    passS;
    passSfull;
    stdXunqAllS;
    cmpXunqAllS;

    name;
    prjCode;
    imgDTB;
    natORflt;
    imgDim;
    method;
    prjInd;
    subjs;
    pass;
    stdXunqAll;
    cmpXunqAll;
    Xunits;
    Xname;

    nCmp;
    nStd;
end
methods
    function obj=exp_details(varargin);
        if isa(varargin{1},'exp_details');
            EXPS=varargin;
        elseif isa(varargin{1},'EAobj');
            EXPS=obj.combine_from_EAobj(varargin{:});
        else;
            EXPS=varargin;
        end
        obj.combine_details(EXPS{:});
        obj.shrink_all();
        obj.get_name();
    end
    function EXPS=combine_from_EAobj(obj,EA);
        EXPS=cell(numel(EA.Eall),1);
        for i = 1:size(EXPS,1);
            EXPS{i}=EA.select_detail(i);
        end
    end
    function obj=combine_details(obj,varargin);
        EXPS=varargin;
        flds=fieldnames(EXPS{1});
        obj.N=numel(EXPS);
        for i =1:length(flds);
            fld=flds{i};
            if startsWith(fld,'A');
                fldname=fld;
            else;
                fldname=['A' fld ];
            end
            obj.(fldname)=cell(obj.N,1);
            for k = 1:obj.N;
                obj.(fldname){k}=EXPS{k}.(fld);
            end
        end
    end
    function obj=shrink_all(obj);
        flds=fieldnames(obj);
        flds=flds(startsWith(flds,'A'));
        for i = 1:length(flds);
            obj.shrink(flds{i});
        end
    end
    function obj=shrink(obj,fld);
        fl=fld(2:end);
        fls=[fl 'S'];
        flb=['b' fl];
        if isempty(obj.(fld)) || strcmp(fld,'Aname');
            return;
        end
        obj.(flb)=0;
        vals=obj.(fld);
        val=obj.(fld){1};
        if isuniform(vals) && ( startsWith(fld,'An') || strcmp(fld,'Asubjs'));
            obj.(fl)=val;
            obj.(flb)=1;
        elseif isuniform(vals) && endsWith(fld,'stdXunqAll') && iscell(val);
            val=num2cell(distribute(val{:}));
            obj.(fl)=val;
            obj.(fls)=Eobj.stdXunq2string(val);
            obj.(flb)=1;
        elseif isuniform(vals) && endsWith(fld,'unqAll') && iscell(val);
            val=distribute(val{:});
            obj.(fl)=val;
            obj.(fls)=strsplit(num2strSane(val),';');
            obj.(flb)=1;
        elseif isuniform(vals) && endsWith(fld,'unqAll');
            obj.(fl)=val;
            obj.(fls)=strsplit(num2strSane(val),';');
            obj.(flb)=1;
        elseif isuniform(vals) && strcmp(fl,'prjInd');
            obj.(fl)=val;
            obj.(fls)=strrep(num2strSane(val),',','.');
            obj.(flb)=1;
        elseif isuniform(vals) && ischar(val);
            obj.(fl)=val;
            obj.(flb)=1;
        elseif isuniform(vals) && isnumeric(val) && strcmp(fl,'pass');
            obj.pass=val;
            obj.passS=strsplit(num2strSane(val),',');
            obj.passSfull=strrep(num2strSane(val),',','-');
        elseif isuniform(vals) && isnumeric(val);
            obj.(fl)=strrep(num2strSane(val),',','-');
            obj.(flb)=1;
        elseif strcmp(fl,'pass');
            n=horzcat(obj.(fld){:});
            obj.pass=n;
            obj.passS=strsplit(num2strSane(n),',');
            obj.passSfull=strrep(num2strSane(n),',','-');
        elseif strcmp(fl,'prjInd');
            v=cellfun(@join_fun,obj.(fld),UO,false);
            %v=join(obj.(fld),'.');
            %obj.(fls)=v{1};
        elseif isnumeric(val);
            n=horzcat(obj.(fld){:});
            obj.(fls)=strrep(num2strSane(n),',','-');
        elseif ischar(val);
            v=join(obj.(fld),'-');
            obj.(fl)=v{1};
        end
        function out=join_fun(x)
            if isnumeric(x)
                x=num2strSane(x);
            end
            out=strrep(x,',','-');
        end
    end
    function obj=get_name(obj);
        obj.name=[obj.prjCode und obj.imgDTB und obj.natORflt und obj.imgDim und obj.method und obj.prjIndS und 'pass' obj.passSfull];
    end
    function obj=rm_subjs(obj,subjs);
        for i = 1:length(subjs);
            subj=subjs{i};
            obj.rm_subj(subj);
        end
    end
    function obj=rm_subj(obj,subj);
        for i = 1:obj.N;
            ind=ismember(obj.Asubjs{i},subj);
            if sum(ind)>0;
                obj.Asubjs{i}(ind)=[];
            end
        end
        obj.shrink('Asubjs');
        % XXX get_name();
    end
    function STR=get_rc_title(obj,varargin);
        if (strcmp(varargin{1},'std') || strcmp(varargin{1},'stdX')) && (numel(varargin)==1 || (numel(varargin)==2 && isnumeric(numel(varargin{2}))));
            STR=obj.stdXunqAllS;
            if numel(varargin) > 1;
                i=varargin{2};
                STR=STR{i};
            end
            st=repmat('Std ',numel(STR),1);
            STR=char(string(transpose(STR)));
            STR=cellstr([st STR]);
            return;
        elseif (strcmp(varargin{1},'cmp') || strcmp(varargin{1},'cmpX')) && (numel(varargin)==1 || (numel(varargin)==2 && isnumeric(numel(varargin{2}))));
            STR=cellfun(@(x) strsplit(x,',') ,transpose(obj.cmpXunqAllS),UO,false);
            STR=vertcat(STR{:});
            sz=size(STR);
            STR=STR(:);
            if numel(varargin) > 1;
                i=varargin{2};
                STR=STR{i};
            end
            st=repmat('Cmp ',numel(STR),1);
            STR=squeeze(char(string(STR)));
            STR=cellstr([st STR]);
            STR=reshape(STR,sz);
            return;
        elseif (strcmp(varargin{1},'subjs') || strcmp(varargin{1},'subj')) &&  (numel(varargin)==1 || (numel(varargin)==2 && isnumeric(numel(varargin{2}))));
            STR=obj.subjs;
            if numel(varargin) > 1;
                i=varargin{2};
                STR=STR{i};
            end
            return;
        end
        STR=cell(obj.N,1);
        for j = 1:length(varargin);
            for i = 1:obj.N;
                str=obj.get_string_cell(varargin{j});
                if strcmp(varargin{j},'pass');
                    STR{i}=[STR{i} ' pass-' str{i}];
                else;
                    STR{i}=[STR{i} ' ' str{i}];
                end
            end
        end
        for i = 1:obj.N;
            STR{i}=STR{i}(2:end);
            STR{i}(1)=makeUpperCase(STR{i}(1));
        end
    end
    function xtitl=get_xtitl(obj)
        if isempty(obj.Xname)
            Xname='X';
        else
            Xname=obj.Xname;
            Xname(1)=makeUpperCase(Xname(1));
        end
        if isempty(obj.Xunits)
            xtitl=Xname;
        else
            xtitl=[Xname ' (' obj.Xunits ')' ];
        end
    end
    function out=get_string_cell(obj,fl);
        if startsWith(fl,'A');
            fl=fl(2:end);
        end
        flb=['b' fl];
        fls=[fl 'S'];
        fla=['A' fl ];
        bfl=  isprop(obj,fl)  && iscell(obj.(fl))  && numel(obj.(fl))> 1  && ischar(obj.(fl){1});
        bfls= isprop(obj,fls) && iscell(obj.(fls)) && numel(obj.(fls))> 1 && ischar(obj.(fls){1});
        if bfl;
            out=obj.(fl);
        elseif bfls;
            out=obj.(fls);
        else;
            out=obj.(fla);
        end
    end
    function ind=get_ind_on_criteria(obj,dim,pass)
        ind=find(ismember(obj.AimgDim,dim) & transpose(obj.pass)==pass);
    end
    function ind=get_unique_exp_ind(obj);
        names=regexprep(obj.Aname,'_pass.*','');
        [~,~,ind]=unique(names);
    end
    function EXPS=combine_passes(obj);
        ind=obj.get_unique_exp_ind();
        N=max(ind);
        flds=fieldnames(obj);
        flds(~startsWith(flds,'A'))=[];
        EXPS=cell(N,1);
        %passes=cell(N,1);
        for i = 1:N;
            I=find(ind==i,1,'first');
            EXPS{i}=struct();
            for j = 1:length(flds);
                fld=flds{j};
                EXPS{i}.(fld)=obj.(fld){I};
                if strcmp(fld,'Aname') || strcmp(fld,'Apass');
                    continue;
                end
            end
            EXPS{i}.pass=[obj.Apass{ind==i}];
            passes=['pass' strrep(num2strSane(EXPS{i}.pass),',','-')];
            EXPS{i}.Aname=regexprep(obj.Aname{I},'pass.*',passes);
        end
        EXPS=exp_details(EXPS{:});
    end
    function EXP=select(obj,ind);
        EXP=obj.get_EXP(ind);
    end
    function EXP=get_EXP(obj,ind);
        EXP=exp_detail(obj,ind);
    end
    function EXPS=split(obj,inds)
        EXPS=cell(length(inds),1);
        inds=rowVec(inds);
        for ii = 1:length(inds)
            i=inds(ii);
            EXPS{ii}=obj.select(i);
        end
        EXPS=exp_details(EXPS{:});
    end
end
end
