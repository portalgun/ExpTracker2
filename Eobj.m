classdef Eobj < handle & Eobj_auto & Eobj_diag & Eobj_init & Eobj_raw & Eobj_data & Eobj_status & Eobj_wrangle &  Eobj_expRunner & Eobj_psy & Eobj_rename & Eobj_fix & Eobj_import %& Eobj_fix 
% TODO
%% in needs db
%%
% fix data
% redo
% -----
% verify indTrl
% verify methodvars
% verify rand
% verify data contents
% run individual trials
% plot training progress
% plot flag and label
%
% subj,exp symbols colors
%
% OTHER OBJECTS
%    psychToolbox
%    exp
%    AMA ?
%    LRSI2DSP2AMA2
%    exp_data
%    plots
%    ROC
properties
    name
    alias % TODO
    % DSP_1 DSP_2 DPSRMS_1 DPA_1-1
    % ALIAS FILE
    prjCode
    imgDTB
    natORflt
    imgDim
    method
    prjInd
    subjs={'TST'};
    pass

    nSubj
    nStd
    nCmp=1;
    nBlk
    nTrl
    nTrlPerBlk
    nTrlPerLvl
    Xname
    Xunits

    flags=struct('bComplete',[],'bDBlocked',[]);
    subjStatus=struct(struct('TST',struct('status','','message','')));

    dirs=struct();
    dir=struct()
    fnames=struct()

    expData=struct();
    rawData=struct();
    indTrl

    methodVars=struct();
        % bUseFeedback
        % nTrialsPerLevel
        % stdXunqAll
        % cmpXunqAll
        % cmpXspacing
        % stdXblkTable
        % cmpXblkTable

    rnd=struct();
    % Intervl also probe placement

    expHost
    pAuthor
    creationDate
    description
    relPapers
    language='matlab';

    LRSI
end
properties(Hidden = true)
    aliases
    prjHier
    expID
    bSkipSyncTest

    dbDir
    
    bLegacy=0;
    legacy
    keyPairs
    burge
    notused
    Eflds

    indBlk
    indStd

    % TESTING
    E
    flds2IFC={'R','RcmpChosen','Rcorrect','cmpIind','stdIind','stdXind','cmpXind','cmpIntrvl','stdIntrvl','stdX','cmpX','responses'};
    LRSIflds={'CP','CL','BI','PRJ','AMA'};
    modeflds={'test','train','pilot'};

    bLoaded=0;
end
methods
    function obj=Eobj(name,prjCode,expID)
        if exist('prjCode','var') && ~isempty(prjCode)
            obj.prjCode=prjCode;
        end
        if exist('expID','var') && ~isempty(expID)
            obj.expID=expID;
        end
        if exist('name','var') && ~isempty(name)
            obj.name=name;
        elseif ~isempty(obj.expID) && ~isempty(obj.prjCode)
            obj.name=[obj.prjCode '_' expID];
        end
    end
    function obj=init(obj,name,prjCode,expID)
        obj.init_dirs();
        obj.dbDir=dbDirs('EXP');

        if isempty(obj.nTrlPerLvl)
            obj.nTrlPerLvl=obj.nTrl/obj.nCmp;
        end
        obj.init_exp_empty_all();
        obj.init_raw_empty_all();
        obj.init_rnd();
        obj.init_fnames();
        obj.nSubj=numel(obj.subjs);

        obj.convert_dirs();
    end
%% META
    function obj=save(obj,bTest)
        if ~exist('bTest','var') || isempty(bTest)
            bTest=0;
        end
        obj.check_data_dir('META',1);
        obj.dir.META=dbDirs('EXP');
        obj.bLegacy=0;

        fname=obj.name;
        dirn=obj.dir.META;
        fnamefull=[dirn fname];

        % obj.check_meta_lock % NOTE
        if ~bTest
            save(fnamefull,'obj');
        end

        % update DB
        % NOTE
    end
    function obj=load_fun(obj)
        obj.aliases=Eobj.getAliases();
        if obj.bLoaded
            out=basicYN('Reload current object?');
            if ~out
                return
            end
        end
        if isempty(obj.name)
            obj.get_fname_by_num();
        end

        out=obj.check_data_dir('META',1);
        if out == 0
            return
        end

        if obj.bLegacy
            EE=obj.load_legacy();
            obj.convert_from_legacy(EE);
            obj.bLoaded=1;
            return
        end

        dirn=obj.dir.META;
        fname=obj.name;

        obj=load([dirn fname]);
        if isstruct(obj)
            obj=obj.obj;
        end
        obj.bLoaded=1;
    end
    function []=get_fname_by_num(obj)
        sizes=cellfun(@numel,obj.aliases);
        maxSz=max(sizes(:,1));
        noAlias=repmat('.',1,maxSz);

        rex= '^(?!_+).*\.mat';
        files=matchingFilesInDir(obj.dir.META,rex);
        files=strrep(files,'.mat','');
        n=length(files);

        for i=1:n
            ind=ismember(obj.aliases(:,2),files{i});
            if ~any(ind)
                alias=noAlias;
                l=0;
            else
                alias=obj.aliases{ind};
                l=length(alias);
            end
            if l > 0 & l < maxSz
                N=maxSz-l;
                alias=[alias repmat(' ',1,N)];
            end
            str=[sprintf('%02d',i) '  ' alias '   ' files{i}];
            disp(str)
        end
        out=basicRange(n,1);
        obj.name=files{out};
    end
% ---
    function obj=reset_exp(obj)
        obj.expData=struct();
        obj.init_exp_empty_all();
    end
    function obj=reset_raw(obj)
        obj.rawData=struct();
        obj.init_raw_empty_all();
    end
    function nobj=copy_as_new_pass(obj,newpassnum)
        if ~isint(newpassnum)
            error('newpassnum must be integers')
        end
        newpassnum=num2str(newpassnum);
        nobj=copyObj(obj);
        nobj.reset_raw();
        nobj=nobj.rename_pass(newpassnum,1);
        nobj.creationDate=date();

        fname=nobj.name;
        dirn=nobj.dir.META;
        fnamefull=[dirn fname];
        if ~exist(fnamefull,'file')
            nobj.save();
            disp('saved.');
        else
            warning('this project exists! not saving. clearing');
            nobj=[];
        end
        
    end
    function nobj=copy_as_new_ind(obj,a,b)
        if ~isint(a) || ~isint(b)
            error('prjinds must be integers')
        end
        prjInd=strrep(num2strSane([a b]),',','-');

        nobj=copyObj(obj);
        nobj.reset_raw();
        nobj.reset_exp();
        nobj=nobj.rename_prjInd(prjInd,1);
        nobj.creationDate=date();

        fname=nobj.name;
        dirn=nobj.dir.META;
        fnamefull=[dirn fname];
        if ~exist(fnamefull,'file')
            nobj.save();
            disp('saved.');
        else
            warning('this project exists! not saving. clearing');
            nobj=[];
        end

    end
end
methods(Static=true)
    function obj=load(name,prjCode,ExpID)
        if ~exist('name','var')
            name=[];
        end
        if ~exist('prjCode','var')
            prjCode=[];
        end
        if ~exist('expID','var')
            expID=[];
        end
        if ~isempty(name)
            alias=Eobj.nameFromAlias(name);
        else
            alias=[];
        end
        if ~isempty(alias)
            name=alias;
        end

        obj=Eobj(name,prjCode,expID);
        obj.init_dirs();
        obj=obj.load_fun();
        obj.init();
        obj.convert_dirs();

    end
    function name=nameFromAlias(alias)
        aliases=Eobj.getAliases;
        ind=ismember(aliases(:,1),alias);
        if sum(ind) > 1
            error('ambiguous alias???');
        elseif sum(ind) == 0
            name=alias;
        else
            name=aliases{ind,2};
        end
    end
    function aliases=getAliases()
        dirn=dbDirs('EXP');
        aliases=readOptsFile([dirn 'aliases.txt']);
    end
    function out=prj_exists(name)
        dir=dbDirs(['EXP'])
        out=exist([dir name '.mat'],'file') > 0
    end
    function obj=new(opts)
        opts=parse_opts(opts);
        obj=Eobj(opts.prjCode,opts.expID);
        if Eobj.prjExists(otps.name)
            error('project already exists');
        end
        flds=fieldnames(opts);
        for i = 1:numel(flds)
            fld=flds{i};
            if strcmp(fld,'expID'); continue; end
            obj.(fld)=opts.(fld)
        end
        obj.creationDate=date();
        obj.init();
        obj.save();
    end
    function opts=parse_opts(opts)
        P={...
              'prjCode',     [],'ischar' ...
              ;'imgDTB',     [],'ischar' ...
              ;'natORflt',   [],'ischar' ...
              ;'imgDim',     [],'ischar' ...
              ;'method',     [],'ischar'...
              ;'prjInd',     [],[]... % optional
              ;'alias',      [],'ischar_e'... % optional
              ...
              ;'nBlk',        [],'isinit1'...
              ;'nTrlPerLvl',  [],'isinit1'...
              ;'nTrlPerBlk',  [],'isinit1'...
              ;'Xname',       [],'ischar'...
              ;'Xunits',      [],'ischar'...
              ;'rndSd',       [],'isint1'...
              ;'expHost',     [],'ischar'...
              ;'pAuthor',     [],'ischar'...
              ;'description', [],'ischar'...
              ...
              ;'rndSd'    [],'isint1'...
              ...
              ;'subjs',   [],'ischarcell_e'... % optional
              ;'pass',    [],'ischarcell_e'... % optional
              ...
        }
        opts=Parse([],opts,P);

        %NATORFLT

        if opts.description < 20
            error('Description to short')
        end
        if opts.pAuthor < 3  || opts.pAuthor > 3
            error('Author name must be 3 char initials')
        end
        name=Eobj.get_E_name(opts.prjCode,...
                        opts.imgDTB,...
                        opts.natORflt,...
                        opts.imgDim,...
                        opts.method,...
                        opts.prjInd,...
                        opts.pass);
        ind=regexp(name,['^' opts.prjCode '_']);
        ind=ind(end);
        opts.expID=name(ind+1:end);
        opts.name=name;
    end
% XXX parsing, no invalid characters s.a. '_'
    function prjCode=parse_prjCode(prjCode)
        if isempty(natORflt); error('prjCode cannot be empty'); end
        if ~ischar(prjCode); error('prjCode must be a string of 3-5 chars (pref 3)'); end
        if length(prjCode) < 5; error('prjCode must be a string of 3-5 chars (pref 3)'); end
        prjCode=makeUpper(prjCode);
    end
    function imgDim=parse_imgDim(imgDim)
        if isempty(imgDim)
            error('imgDim cannot be empty')
        elseif isnum(imgDim)
            imgDim=str2double(imgDim);
        elseif isalphanum(imgDim) && imgDim(end)=='D'
            imgDim=str2double(imgDim(1:end-1));
        elseif ~isint(imgDim)
            error('unhandled imgDim dimensions')
        end
        if imgDim > 4
            error('unhandled imgDim dimensions')
        end
        imgDim=[num2str(imgDim,'%i') 'D'];
    end
    function natORflt=parse_natORflt(natORflt)
        if isempty(natORflt)
            error('natORflt cannot be empty')
        end
        natORflt=makeUpper(opts.natORflt);
        if ~ismember(opts.natORflt,{'NAT','FLT'})
            error('invalid value for natORflt')
        end
    end
    function method=parse_method(method)
        validMethods={'ADJ','MCH','EST'};
        method=makeUpper(method);
        if regExp(method,'^[1-9][AI]FC$' )
            method=strrep(method,'A','I');
        elseif ~ismember(method,validMethods)
            error('invalid method')
        end
    end
    function prjInd=parse_prjInd(prjInd)
        if ischar(prjInd) && ~regExp(prjInd,'[1-9]{1}[0-9]*-[1-9]{1}[0-9]*')
            error('invalid pojectInd')
        elseif ~isint (prjInd)
            error('prjInds must be integers')
        elseif ~all(isint(prjInd)) || numel(prjInd) > 2 || isempty(prjInd)
            error('prjInds must be 1 or 2 integers')
        elseif ischar(prjInd)
            prjInd=str2double(strsplit(prjInd,'-'));
        end
    end
    function pass=parse_pass(pass)
        if isempty(pass)
            pass=1;
        elseif isnumeric(pass) && ~isint(pass)
            error('invalid pass number');
        elseif ischar(pass) & ~regExp(prjInd,'[Pp]{0,1}(ass|ASS){0,1}[1-9]{1}[0-9]*')
            error('invalid pass number');
        elseif ischar(pass)
            ind=regExp(pass,'[0-9]');
            pass=num2str(pass(ind));
        end
    end
    function mode=parse_mode(mode)
        validModes={'test','train','pilot'};
        if ~ischar(mode)
            error('invalid mode');
        end
        mode=makeLower(mode);
        if ~imember(mode,validModes)
            error('invalid mode');
        end
    end
    function blk=parse_blk(blk)
        if numel(blk) > 1 || isempty(blk)
            error('blk must be a single integer')
        elseif isnumeric(blk) && ~isint(blk)
            error('invalid blk');
        elseif ischar && regexp(blk,'^[0-9]+^')
            error('invalid blk');
        end
    end
    function std=parse_std(std)
        if numel(std) > 1 || isempty(std)
            error('std must be 1 element')
        end
    end
    function name=get_base_name(prjCode,imgDTB,natORflt,imgDim,method,prjInd)
        prjCode  = Eobj.parse_prjCode(prjCode);
        imgDTB   = Eobj.parse_imgDTB(imgDTB);
        natORflt = Eobj.parse_natORflt(natORflt);
        imgDim   = Eobj.parse_imgDim(imgDim);
        method   = Eobj.parse_method(method);
        prjInd   = Eobj.parse_prjInd(prjInd);
        prjInd   = strrep(num2strSane(prjInd),',','-')

        parts={ prjCode, imgDTB, natORflt, imgDim, method, prjInd};
        strs=join(parts,'_');
        name=strs{1};
    end
    function name=get_E_name(prjCode,imgDTB,natORflt,imgDim,method,prjInd,pass)
        bname=Eobj.get_base_name(prjCode,imgDTB,natORflt,imgDim,method,prjInd);
        pass=Eobj.parse_pass(pass);
        pass=['pass' num2str(pass)];
        name=[bname und pass];
    end
    function name=get_expData_name(prjCode,imgDTB,natORflt,imgDim,method,prjInd, mode,std,blk)
        bname=Eobj.get_base_name(prjCode,imgDTB,natORflt,imgDim,method,prjInd);
        mode=Eobj.parse_mode(mode);
        std=Eobj.parse_std(std);
        std=Eobj.auto_std_str(std); % XXX
        blk=Eobj.auto_blk_str(blk); % XXX
        name=['exp' bname und std und blk];
    end
end
end
