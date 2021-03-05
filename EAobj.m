classdef EAobj < handle & EAobj_inds
% TODO
% DB
% check pass status
% run subj
% sort by completion
% sort by status
% locked experiments
% mass mod
% check_all
% load_all
% summary
properties
    Eall
    details
    namesAll
    namePartsAll

    names

    prjCodes
    imgDTB
    natORflt
    imgDim
    method
    passes
    prjInds

    subjs
    modes
    stds
    blocks
end
properties(Hidden=true)
    aliases
    aliasesHash
    dbDir
    dbURL
    bLegacy=0;
end
methods
    function obj=EAobj(Opts)
        if ~exist('Opts','var')
            Opts=struct();
        end

        dire=dbDirs('EXP');
        if obj.bLegacy
            obj.dbDir=[dire 'legacy' filesep];
        else
            obj.dbDir=dire;
        end

        p=inputParser();
        p.addParameter('names',[]);

        p.addParameter('prjCodes',[]);
        p.addParameter('imgDTB',[]);
        p.addParameter('natORflt',[]);
        p.addParameter('imgDim',[]);
        p.addParameter('method',[]);
        p.addParameter('passes',[]);
        p.addParameter('prjInds',[]);

        p.addParameter('subjs',[]);
        p.addParameter('modes',[]);
        p.addParameter('stds',[]);
        p.addParameter('blocks',[]);
        p.addParameter('aliases',[]);

        p=parseStruct(Opts,p);

        flds=fieldnames(p.Results);
        for i = 1:length(flds)
            fld=flds{i};
            obj.(fld)=p.Results.(fld);
        end

        obj.get_names_from_aliases();
        obj.get_all_names();
        if isempty(obj.names)
            obj.names=obj.namesAll;
        else
            obj.check_names();
        end
        obj.get_nameParts();
        obj.load_Eall();
        obj.select_details();
    end
    function keyPairs=get_names_keyPairs(obj)
        keyPairs={...
            ;'prjCodes', 1 ...
            ;'imgDTB'  , 2 ...
            ;'natORflt', 3 ...
            ;'imgDim'  , 4 ...
            ;'method'  , 5 ...
            ;'prjInds' , 7 ...
            ;'passes'  , 6 ...
        };
    end
    function keyPairs=get_narrow_keyPairs(obj)
        keyPairs={...
            ;'subjs'  , 1 ...
            ;'modse'  , 2 ...
            ;'stds'   , 3 ...
            ;'blocks' , 4 ...
        };
    end
    function obj=get_all_names(obj)
        rex= '^(?!_+).*\.mat';
        files=matchingFilesInDir(obj.dbDir,rex);
        obj.namesAll=cellfun(@filepartsfun,files,'UniformOutput',0);
        obj.namePartsAll=cellfun(@splitfun,obj.namesAll,'UniformOutput',0);
        obj.namePartsAll=vertcat(obj.namePartsAll{:});
        obj.namePartsAll(:,end )=cellfun(@(x) str2double(x(end)),obj.namePartsAll(:,end),'UniformOutput',0);
        function out=filepartsfun(files)
            [~,out]=fileparts(files);
        end
        function out=splitfun(names)
            out=strsplit(names,'_');
        end
    end
    function obj=check_names(obj)
        if ischar(obj.names)
            obj.names={obj.names};
        end
        obj.names = transpose({ obj.names{:} }) ;
        ind = ~ismember(obj.names,obj.namesAll);
        if any(ind)
            bad=obj.names(ind);
            for i = 1:length(bad)
                warning(['Not valid project inds: '  bad{i}] )
            end
        end
        obj.names(ind)=[];
    end
    function obj=get_nameParts(obj,nameParts)
        keyPairs=obj.get_names_keyPairs;
        ind=zeros(size(obj.namePartsAll));
        for i = 1:size(keyPairs,1)
            fld=keyPairs{i,1};
            if isempty(obj.(fld))
                continue
            elseif isnumeric(obj.namePartsAll{1,i})
                val=vertcat(obj.namePartsAll{:,i});
            else
                val=obj.namePartsAll(:,i);
            end
            ind(:,i)=ismember(val,obj.(fld));
        end
        ind=logical(prod(ind,2));
        names=obj.namesAll(ind);
        obj.names=unique([obj.names;names]);
    end
    function obj=load_Eall(obj)
        obj.Eall=cell(length(obj.names),1);
        for i = 1:length(obj.names)
            obj.Eall{i}=Eobj.load(obj.names{i});
        end
    end
    function obj=save_Eall(obj)
        for i = 1:length(obj.names)
            obj.Eall{i}.save();
        end
    end
    function obj=get_narrow_from_parts()
        % TODO
    end
    function out=get_name_parts()
        %flds=
        for i =1:length(flds)
            fld=flds{i};
            out=~isempty(obj.(fld));
            if out
                break
            end
        end
    end
    function obj=apply_prop_all(obj,fld,val)
        for i = 1:length(obj.names)
            obj.Eall{i}.(fld)=val
        end
    end
    function VALS=get_prop_all(obj,varargin)
        n=length(obj.Eall);
        VALS=cell(n,1);
        flds=cellfun(@(x) ['(' '''' x '''' ').'],varargin,'UniformOutput',0);
        flds=horzcat(flds{:});
        flds=flds(1:end-1);
        for i = 1:n
            str=['obj.Eall{i}.' flds ';'];
            VALS{i}=eval(str);
        end
        VALS=cell2matSafe(VALS);
    end
    function out=get_narrow_parts(obj)
        % TODO
    end
    function obj=wrangle_Eall(obj,bSave)
        n=length(obj.names);
        p=pr(n,[],'Mass wrangling data','...');
        for i = 1:n
            p.u();
            obj.Eall{i}.wrangle_all(bSave,1);
        end
        p.c();
    end
    function [out,EXPS,fnames]=load_out(obj,fld,subj,mode,std,rmsubjs,rmstds)
        if ~exist('rmsubjs','var')
            rmsubjs=[];
        end
        if ~exist('rmstds','var')
            rmstds=[];
        end
        if ~exist('std','var')
            std=[];
        end
        if ~exist('mode','var')
            mode=[];
        end
        if ~exist('blk','var')
            blk=[];
        end
        n=length(obj.names);
        out=cell(n,1);
        fnames=cell(n,1);
        for i = 1:n
            [out{i},fnames{i}]=obj.Eall{i}.load_out(fld,subj,mode,std);

            if isempty(rmsubjs)
                continue
            end

            if ~isempty(rmsubjs)
                out{i}=obj.Eall{i}.(['rm_subj_' fld])(out{i},rmsubjs);
            end
            if ~isempty(rmstds)
                out{i}=obj.Eall{i}.(['rm_std_' fld])(out{i},rmstds);
            end

        end
        if nargout > 1
            EXPS=obj.select_details();
            EXPS.rm_subjs(rmsubjs);
        end
    end
    function EXP=select_detail(obj,ind)
        EXP=obj.Eall{ind}.select_detail();
    end
    function EXPS=select_details(obj)
        EXPS=exp_details(obj);
    end
%%
    function obj=get_names_from_aliases(obj)
        obj.aliasesHash=EAobj.get_aliases_hash();
        if isempty(obj.aliases)
            return
        end
        obj.names=cell(length(obj.aliases),1);
        for i = 1:length(obj.aliases)
            obj.names{i}=obj.aliasesHash(obj.aliases{i});
        end
    end
end
methods(Static=true)
    function hash=get_aliases_hash()
        hash=containers.Map;
        dirn=autodir(dbDirs('EXP'));
        cel=readOptsFile([dirn 'aliases.txt']);
        for i = 1:size(cel,1)
            key=cel{i,1};
            hash(key)=cel{i,2};
        end
    end
    function obj=load_by_aliases(Aliases)
        Opts=struct();
        Opts.aliases=Aliases;
        obj=EAobj(Opts);
    end
    function obj=run_loop_old(Aliases,subj,mode,std)
       
        if ~exist('Aliases','var')
            Aliases=[];
        end
        if ~exist('mode','var')
            mode=[];
        end
        if ~exist('std','var')
            std=[];
        end
        obj=EAobj.load_by_aliases(Aliases);

        i=0;
        exitflags=zeros(length(Aliases));
        while ~all(exitflags)
            i=i+1;
            if i > length(Aliases)
                i=1;
            end

            if ~exitflags(i)
                [obj.Eall{i},exitflags(i)]=obj.Eall{i}.run(subj,mode,std);
            end
        
        end
    end
    function obj=run_loop_all(subj,Aliases,rule,mode,std,blk)
        bTest=1;
        if ~exist('Aliases','var')
            Aliases=[];
        end
        if ~exist('rule','var') || isempty(rule)
            rule='rand_all';
        end
        if ~exist('mode','var')
            mode=[];
        end
        if ~exist('std','var')
            std=[];
        end
        if ~exist('blk','var')
            blk=[];
        end
        obj=EAobj.load_by_aliases(Aliases);

        if contains(rule,'all')
            bReset=0;
            INDS=obj.get_ind(subj,mode,std,blk,rule);
        else
            bReset=1;
        end

        i=0;
        while true
            if bReset
                i=1;
                INDS=obj.get_ind(subj,mode,std,blk,rule);
            else
                i=i+1;
            end
            ind=INDS(i,1);
            std=INDS(i,2);
            blk=INDS(i,3);
            if bTest==2
                disp(num2str(size(INDS)));
                disp('i  s  b  pg p')
                disp(num2str(INDS(i,:)));
            end

            obj.Eall{ind}.run(subj, mode,std,blk, 1, bTest);
        end
    end
end
end
