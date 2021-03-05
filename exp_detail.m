classdef exp_detail < handle
properties
    name
    prjCode
    imgDTB
    natORflt
    imgDim
    method
    prjInd
    subjs
    pass
    stdXunqAll
    cmpXunqAll
    nCmp
    nStd
    Xunits
    Xname
end
methods
    function obj=exp_detail(IN,varargin)
        if isa(IN,'Eobj')
            obj.construct_from_Eobj(IN);
        elseif isa(IN,'exp_details')
            obj.construct_from_details(IN,varargin{1});
        end
    end
    function obj=construct_from_Eobj(obj,IN)
        obj.prjCode =IN.prjCode;
        obj.imgDTB  =IN.imgDTB;
        obj.natORflt=IN.natORflt;
        obj.imgDim  =IN.imgDim;
        obj.method  =IN.method;
        obj.prjInd  =IN.prjInd;
        obj.subjs   =IN.subjs;
        obj.pass    =IN.pass;
        obj.name    =IN.name;
        obj.Xunits  =IN.Xunits;
        obj.Xname  = IN.Xname;
        if isfield(IN.methodVars,'stdXunqAll')
            obj.stdXunqAll=IN.methodVars.stdXunqAll;
        end
        if isfield(IN.methodVars,'cmpXunqAll')
            obj.cmpXunqAll=IN.methodVars.cmpXunqAll;
        end
        obj.nCmp=IN.nCmp;
        obj.nStd=IN.nStd;
    end
    function construct_from_details(obj,IN,IND)
        flds=fieldnames(IN);
        ind=startsWith(flds,'A');
        flds=flds(ind);
        for i = 1:length(flds)
            fld=flds{i};
            fl=fld(2:end);
            if ~isempty(IN.(fld))
                obj.(fl)=IN.(fld){IND};
            end
        end
    end
    function stri=std_str(obj)
        stri=num2str(obj.std);
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
    function titl=get_titl(obj)
        titl=strrep(obj.name,'_',' ');
    end
end
end
