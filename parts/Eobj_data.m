classdef Eobj_data < handle
methods
% Open
    function []=open_data(obj,dType,fld)
        out=obj.check_data(dType,fld,1);
        if out == 1
            open(obj.fnames.(dType).(fld));
        end

    end
    function []=open_def(obj,fld)
        obj.open_data('DEF',fld);
    end
    function []=open_code(obj,fld)
        obj.open_data('CODE',fld);
    end

% GEN NAME
    function name=gen_name_def(obj,fld)
        % exp
        fld=makeUpperCase(fld);
        if ~isfield(obj.fnames.DEF,fld)
            disp(['not a valid field ' fld])
        end
        fld=makeLowerCase(fld);
        name=['def' und fld und regexprep(obj.name,'_pass[1-9]+','')];
        name=strrep(name,'-','d');
    end
    function name=gen_name_exp(obj,mode,std,blk,bD)
        if ~exist('bD','var') || isempty(bD)
            bD=0;
        end
        mode=obj.auto_mode(mode);
        std=obj.auto_std_str(std);
        blk=obj.auto_blk_str(blk);
        name=['exp' und regexprep(obj.name,'_pass[1-9]+','') und mode und std und blk];
        if bD
            m=['_' obj.auto_ind_str() '_'];
            r=['_' obj.auto_indd_str() '_'];
            name=sed('s',name,m,r);
        end
    end
    function name=gen_name_in(obj,fld)
        FLD=makeUpperCase(fld);
        if ~isfield(obj.fnames.IN,FLD)
            disp(['not a valid field ' FLD])
            return
        end
        name=[fld und obj.name];
        name=regexprep(name,'_pass[1-9]+','');
    end
    function fname=gen_name_out(obj,fld, subj,mode,std)
        if ~exist('subj','var')
            subj=[];
        end
        if ~exist('mode','var')
            mode=[];
        end
        if ~exist('std','var')
            std=[];
        end
        fld=obj.auto_out_fld(fld);
        code=check_data_dim(obj,'OUT',fld);
        if code==1
            fld=makeLowerCase(fld);
            fld=obj.auto_out_fld(fld);
            fname=[fld und obj.name ];
            strrep(fld,'prunedp','prunedP');
        elseif code==2 && ~isempty(subj) && ~isempty(mode) && ~isempty(std)
            subj=obj.auto_subj(subj);
            mode=obj.auto_mode(mode);
            std=num2str(obj.auto_std_num(std));
            fname=[fld und obj.name und und subj und mode und std];
        elseif code ==2
            fname = [fld und obj.name];
        elseif code == 3
            % XXX
            error('Write Code')
        end
    end
    function fname=gen_fname_def(obj,fld)
        % exp
        name=obj.gen_name_def(fld);
        fname=[obj.dir.DEF makeLowerCase(fld) filesep name];
    end
    function fname=gen_fname_exp(obj,mode,std,blk,bD)
        if ~exist('bD','var') || isempty(bD)
            bD=0;
        end
        name=obj.gen_name_exp(mode,std,blk,bD);
        fname=[obj.dir.EXP name '.mat'];
    end
    function fname=gen_fname_in(obj,fld)
        % exp
        name=obj.gen_name_in(fld);
        fname=[obj.dir.IN name '.mat'];
    end
    function fname=gen_fname_out(obj,fld, subj,mode,std)
        if ~exist('subj','var')
            subj=[];
        end
        if ~exist('mode','var')
            mode=[];
        end
        if ~exist('std','var')
            std=[];
        end
        fld=obj.auto_out_fld(fld);
        name=obj.gen_name_out(fld,subj,mode,std);
        dirn=obj.dir.OUT;
        fname=[dirn name '.mat'];

    end
%% GET NAME
    function name=get_name_in(obj,fld)
        name=obj.get_name_data('IN',fld);
    end
    function name=get_name_data(obj,dType,fld,subj,mode,std)
        if ~exist('subj','var') || isempty(subj)
            subj=[];
        end
        if ~exist('mode','var') || isempty(mode)
            mode=[];
        end
        if ~exist('std','var') || isempty(std)
            std=[];
        end
        dType=makeUpperCase(dType);
        if strcmp(dType,'OUT')
            fld=obj.auto_out_fld(fld);
        else
            fld=makeUpperCase(fld);
        end

        code=obj.check_data_dim(dType,fld);

        %if code==2  && ~isempty(subj) && ~isempty(mode) && ~isempty(std)% inds
        %    subj=obj.auto_subj_num(subj);
        %    std=obj.auto_std_ind(std);
        %    mode=obj.auto_mode_num(mode);
        %    name=obj.fnames.(dType).(fld){subj,mode,std};
        if code ==2
            name=obj.fnames.(dType).(fld){subj,mode,std};
        elseif code == 1
            name=obj.fnames.(dType).(fld);
        elseif code == -1
            name=obj.fnames.(dType).(fld){subj};
        end
        if iscell(name)
            name=name{1};
        end
        [~,name]=filePartsSane(name);
    end
%% GET FNAME
    function fname=get_fname_exp(obj,mode,std,blk)
        if ~exist('mode','var') || isempty(mode)
            mode=[];
        end
        if ~exist('std','var') || isempty(std)
            std=[];
        end
        name=obj.get_name_exp(mode,std,blk);
        if isempty(name)
            fname=[];
            return
        end
        fname=[obj.dir.EXP name '.mat'];
    end
    function fname=get_fname_def(obj,fld)
        % XXX WORKED OFF PC
        name=obj.get_name_data('DEF',fld);
        dirn=obj.dir.DEF;
        ext='.m';
        sdir=makeUpperCase(fld);
        fname=[dirn sdir filesep name ext];
    end
    function fname=get_fname_in(obj,fld)
        name=obj.get_name_data('IN',fld);
        ext='.mat';
        dirn=obj.dir.IN;
        fname=[dirn name ext];
    end
    function name=get_name_exp(obj,mode,std,blk)
        std=obj.auto_std_fld(std);
        mode=obj.auto_mode(mode);
        blk=obj.auto_blk_num(blk);
        name=obj.expData.(mode).(std){blk};
    end
    function fname=get_fname_data(obj,dType,fld,subj,mode,std)
        if ~exist('subj','var') || isempty(subj)
            subj=[];
        end
        if ~exist('mode','var') || isempty(mode)
            mode=[];
        end
        if ~exist('std','var') || isempty(std)
            std=[];
        end
        dType=makeUpperCase(dType);
        fld=makeUpperCase(fld);
        name=obj.get_name_data(dType,fld,subj,mode,std);
        dirn=obj.dir.(dType);
        if strcmp(dType,'CODE')
            ext='.m';
        else
            ext='.mat';
        end
        fname=[dirn name ext];
    end
%%  LOAD
    function Def=load_def(obj,fld)
        if obj.bBlk & obj.bPtchs
            fname=[dbDirs('loc') obj.fnames.DEF.(fld)];
        else
            fname=obj.get_fname_def(fld);
        end
        cur=pwd;
        [path,name]=filePartsSane(fname);
        cd(path);
        run(name);
        cd(cur);
        clear cur fname path name

        vars=whos();
        vars=transpose(vertcat({vars.name}));
        Def=struct();

        contFlds={'obj','Def','vars','var','fld'};
        for i = 1:length(vars)
            var=vars{i};
            if ismember(var,contFlds)
                continue
            end
            eval(['Def.(var)=' var ';']);
        end
    end
    function out=load_in(obj,fld,mthd)
        if isempty(mthd) && contains(fld,und)
            fld=[fld und mthd];
        end

        % CHECK
        out=obj.check_data('IN',fld,1,[],[],[]);
        if out < 1
            return
        end

        % GET DIR
        dirn=obj.dir.IN;


        % GET FNAME
        fname=obj.get_fname_data('IN',fld,subj,mode,std);

        % LOAD
        out=load([dirn fname]);

        while true
            flds=fieldnames(out);
            if isstruct(out) && numel(flds)==1 && isstruct(out.(flds{1}))
                out=out.(flds{1});
            else
                return
            end
        end
    end
    function out=load_exp(obj,mode,std,blk)

        % XXX checking
        %out=obj.check_data('IN','EXP',1,std,mode,blk);
        %
        %if out < 1
        %    return
        %end
        %
        if obj.bPtchs && obj.bBlk
            base=dbDirs('base');
            fname=[base obj.fnames.IN.(['DMP_' makeUpperCase(obj.auto_mode(mode))])];
            if ~endsWith(fname,'.mat')
                fname=[fname '.mat'];
            end
            load(fname);
            P.exp_init(obj.alias,obj.auto_mode_num(mode),obj.auto_std_ind(std),blk);
            out=P;

        else
            fname=obj.get_fname_exp(mode,std,blk);
            out=load(fname);
        end

        while true
            flds=fieldnames(out);
            if isstruct(out) && numel(flds)==1 && isstruct(out.(flds{1}))
                out=out.(flds{1});
            else
                return
            end
        end
    end
    function [out,fnames]=load_out(obj,fld,subj,mode,std,blk)
        %% whether single or blk format
        % 1 reg
        % -1 irreg
        % 2 empty
        if ~exist('subj','var')
            subj=[];
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
        code=check_data_dim(obj,'OUT',fld);
        if code == -1
            [out,fnames]=obj.load_out_irreg(fld,subj,mode,std,blk);
        elseif code == 1
            [out,fnames]=load_out_fun(obj,fld,subj,mode,std,blk);
        end

    end
    function [out,fnames]=load_out_irreg(obj,fld,subj,mode,std,blk)
        if ~exist('subj','var') || isempty(subj)
            subj=1:obj.nSubj;
        end
        if ~exist('std','var') || isempty(std)
            std=1:obj.nStd;
        end
        if ~exist('mode','var') || isempty(mode)
            mode=1;
        end
        if ~exist('blk','var') || isempty(blk) && contains(fld,'clean')
            blk=1:obj.nBlk;
        else
            blk=0;
        end

        bMult= iscell(fld) || ...
        ((iscell(mode) || (all(isint(mode)) && numel(mode) > 1))) || ...
        ((iscell(std) || (all(isnumeric(std)) && numel(std) > 1)));

        if ~bMult
            [out,fnames]=obj.load_out_fun(fld,subj,mode,std,blk);
            return
        end
        subj=obj.auto_subj_num(subj);
        mode=obj.auto_mode_num(mode);
        std=obj.auto_std_ind(std);

        INDS=distribute(subj,std,blk);
        n=size(INDS,1);
        out=cell(n,1);
        fnames=cell(n,1);
        for i = 1:n
            [out{i},fnames{i}]=obj.load_out_fun(fld,INDS(i,1),mode,INDS(i,2),INDS(i,3));
        end
    end
    function [out,fname]=load_out_fun(obj,fld,subj,mode,std,blk)
        if ~exist('subj','var')
            subj=[];
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

        % CHECK
        [out,fname]=obj.check_data_load('OUT',fld,subj,mode,std);

        % LOAD
        if out==0
            error(['No such file: ' newline '    ' fname]);
        else
            out=load(fname);
        end

        while true
            flds=fieldnames(out);
            if isstruct(out) && numel(flds)==1 && isstruct(out.(flds{1}))
                out=out.(flds{1});
            else
                return
            end
        end
    end
    function [out,fname]=check_data_load(obj,dType,fld,subj,mode,std,blk)
        if ~exist('subj','var')
            subj=[];
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

        [codeDB,~,fname,~]=obj.check_data(1,'OUT',fld,subj,mode,std,blk);
        if ismember('W',codeDB)
            strrep(codeDB,'W','');
        end
        if isempty(codeDB)
            out=1;
        else
            out=0;
        end
    end
    function SS=load_out_all(obj,dType,fld,mode)
        SS=cell(numel(obj.subjs),abs(obj.nStd));
        for subj = 1:numel(obj.subjs)
        for std  = 1:obj.nStd
            obj.load_out(SS{subj,std},dType,fld,subj,mode,std);
        end
        end
    end
%% SAVE OUT
    function obj=save_data(obj,S,dType,fld,subj,mode,std,bOverwrite)
        if ~exist('bOverwrite','var')
            bOverwrite=0;
        end
        if isempty(S)
            error('Data not included')
        end
        if ~exist('subj','var')
            subj=[];
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

        % CHECK
        [out,fname]=obj.check_data_save(dType,fld,subj,mode,std,blk,bOverwrite);
        if ~out
            return
        end


        % GET DIR
        dirn=obj.dir.(dType);

        % GET  NAME
        [~,name]=filePartsSane(fname);
        % SAVE
        save(fname,'-v7.3','S');
        disp(['Saving ' fname ' successful.']);


        if ~isempty(subj) && ~isempty(mode) && ~isempty(std)
            s=obj.auto_subj_num(subj);
            m=obj.auto_mode_num(mode);
            d=obj.auto_std_ind(std);
            obj.fnames.(dType).(fld){s,m,d}=name;
        else
            obj.fnames.(dType).(fld)=name;
        end

        % UPDATE META
        obj.save();
    end
    function obj=save_data_all(obj,SS,dType,fld,mode)
        for subj = 1:numel(obj.subjs)
        for std  = 1:obj.nStd
            obj.save_data(SS{subj,std},dType,fld,subj,mode,std);
        end
        end
    end
%% CHECKS
    function [OUT,CODE] = check_data_all(obj,bPrint)
        dTypes=fieldnames(obj.fnames);
        nD=numel(dTypes);
        OUT=cell(nD,1);
        CODE=cell(nD,1);
        if ~exist('bPrint','var')
            bPrint=1;
        end
        for d = 1:nD
            dType=dTypes{d};
            flds=fieldnames(obj.(dType));
            nF=numel(flds);
            out=cell(nF,1);
            code=cell(nF,1);
            for f = 1:nF
                fld=flds{f};
                code=obj.check_data_dim(dType,fld);
                if code==1
                    [out{f},code{f}]=obj.check_data(dType,fld,bPrint);
                elseif code == 2
                    [out{f},code{f}]=obj.check_data_inds(dType,fld,bPrint);
                elseif code == -1
                    [out{f},code{f}]=obj.check_data_irreg(dType,fld,bPrint);
                end
            end
            CODE{d}=code;
            OUT{d}=out;
        end
    end
    function [out,fname]=check_data_save(obj,dType,fld,subj,mode,std,blk,bOverwrite)
        if ~exist('bOverwrite','var') || isempty(bOverwrite)
            bOverwrite=0;
        end
        if ~exist('bPrint','var') || isempty(bPrint)
            bPrint=0;
        end
        if ~exist('subj','var')
            subj=[];
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
        [codeDB,codeFN,fnameDB,fnameGen]=obj.check_data(0,dType,fld,subj,mode,std,blk);
        if ~strcmp(fnameDB,fnameGen)
            disp(['save name: ' fnameDB]);
            disp(['gen  name: ' fnameGen]);
            out1=basicYN('These two do not match. Continue saving under gen name?');
            if ~out1
                out=0;
                return
            end
        end
        if   ismember('E',codeFN) && ismember('-',codeDB)
            out=1;
            fname=fnameGen;
        elseif ~ismember('E',codeFN) && ~ismember('E',codeDB)
            if bOverwrite
                out=1;
            else
                out=basicYN(['File exists, and in database. Overwrite ' fnameGen ' ?']);
            end
            fname=fnameGen;
        elseif ismember('E',codeFN)  &&  ismember('E',codeDB)
            if bOverwrite
                out=1;
            else
                out=basicYN(['DB says ' fnameDB ' exists, but file cannot be found. Continue saving?']);
            end
            fname=fnameDB;
        elseif ismember('E',codeFN)  &&  ~ismember('E',codeDB)
            dk
        elseif ~ismember('E',codeFN)  &&  ismember('-',codeDB)
            if bOverwrite
                out=1;
            else
                out=basicYN(['File exists, but not in db. Overwrite ' fnameGen ' ?']);
            end
            fname=fnameGen;
        else
            out=0;
            %chkAllPrint(codeFN);
        end
    end
    function [codeDB,codeFN,fnameDB,fnameGen]=check_data(obj,bPrint,dType,fld,subj,mode,std,blk)
        if ~exist('bPrint','var') || isempty(bPrint)
            bPrint=0;
        end
        if ~exist('subj','var')
            subj=[];
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
        [codeDB,fnameDB]=obj.check_data_db(dType,fld,subj,mode,std);
        [codeFN,fnameGen]=obj.check_data_gen(dType,fld,subj,mode,std,blk);

        if endsWith(fnameDB,'.mat')
            codeDB=strrep(codeDB,'X','');
        end
        if endsWith(fnameGen,'.mat')
            codeGen=strrep(codeFN,'X','');
        end

        %if bPrint
        %    disp('codeDB:')
        %    str=chkAllPrint(codeDB,fulldir,highestdir,fname);
        %    str=indent(str,4);
        %    disp(str)

        %    str=chkAllPrint(codeFN,fulldir,highestdir,fname);
        %    str=indent(str,4);
        %    disp(str);
        %end
    end
    function [code,fname]=check_data_db(obj,dType,fld,subj,mode,std,bSkipDir)
        if ~exist('bSkipDir','var') || isempty(bSkipDir)
            bSkipDir=0;
        end
        if ~exist('subj','var')
            subj=[];
        end
        if ~exist('mode','var')
            mode=[];
        end
        if ~exist('std','var')
            std=[];
        end
        fname=obj.get_fname_data(dType,fld,subj,mode,std);
        if isempty(fname)
            code='-';
            return
        end
        [~,code]=chkFileAll(fname,1,bSkipDir);
    end
    function [code,fname]=check_data_gen(obj,dType,fld,subj,mode,std,blk,bSkipDir)
        if ~exist('subj','var')
            subj=[];
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
        dType=makeUpperCase(dType);
        fld=makeUpperCase(fld);

        if ~exist('bSkipDir','var') || isempty(bSkipDir)
            bSkipDir=0;
        end

        switch dType
        case 'CODE'
            error('Code fnames cannot be generated. Must be set manually.');
        case 'DEF'
            fname=obj.gen_fname_def(fld);
        case 'IN'
            fname=obj.gen_fname_in(fld);
        case 'EXP'
            fname=obj.gen_fname_exp(mode,std,blk,bD);
        case 'OUT'
            fname=obj.gen_fname_out(fld,subj,mode,std);
        case 'MAP'
            error('Write code');
        otherwise
            error(['Unhandled type ' dType]);
        end


        [~,code]=chkFileAll(fname,1,bSkipDir);
    end
    function [code] = check_data_dir(obj,dType,bPrint)
        dirn=obj.dir.(dType);
        code=chkDirAll(dirn,1);
    end
    function code=check_data_dim(obj,dType,fld)
        if strcmp(dType,'OUT')
            fld=obj.auto_out_fld(fld);
        end
        %% whether single or blk format
        % 1 reg
        % -1 irreg
        % 2 empty
        FLD=obj.fnames.(dType).(fld);
        if iscell(FLD) && numel(FLD)==1
            FLD=FLD{1};
        end
        sz=size(FLD);
        SZ=[numel(obj.nSubj) 3 abs(numel(obj.nStd))];
        bInd=isempty(FLD) || ischar(FLD);
        bIrreg=iscell(obj.fnames.(dType).(fld)) && isequal(sz,SZ);
        if bInd && ~bIrreg
            code=1;
        elseif ~bInd && ~bIrreg
            code=2;
        elseif bIrreg
            code=-1;
        else
            error(['Fname for ' dtype ' ' fld ' is bad']);
        end
    end
    function [out,code]=check_data_inds(obj,dType,fld,bPrint)
        if strcmp(dType,'OUT')
            fld=obj.auto_out_fld(fld);
        end
        if ~exist('bPrint','var')
            bPrint=1;
        end
        nSubj=numel(obj.subjs);
        code=cell(nSubj,3,obj.nStd);
        out=zeros(nSubj,3,obj.nStd);
        for s = 1:numel(obj.subjs)
        for m = 1:3
        for i = 1:obj.nStd
            [out(s,m,i),code{s,m,i}]=obj.check_data(dType,fld,bPrint,s,m,i);
        end
        end
        end
    end
    function [out,code]=check_data_irreg(dType,fld,bPrint)
        if strcmp(dType,'OUT')
            fld=obj.auto_out_fld(fld);
        end
        if ~exist('bPrint','var')
            bPrint=1;
        end
        n=numel(obj.(dType),(fld));
        out=zeros(n,1);
        code=cell(n,1);
        for i = 1:n
            [out(i),code{i}]=obj.check_data(dType,fld,bPrint,i);
        end
    end
% DELETE
    function obj=delete_data(obj,dType,fld)
        % TODO
    end
    function obj=delete_data_all(obj)
        % TODO
    end
% CLEAR
    function obj=clear_data(obj,dType,fld)
        % TODO
    end
    function obj=clear_data_all(obj)
        % TODO
    end
end
end

    %% subj has different role for irregular XXX revisit
    %% OUT CODES
    %% -1 empty name
    %% 0  file doesn't exist
    %% 1  file exists
    %% 2  direcotory doesnt exist
    %% 2  read but not write
