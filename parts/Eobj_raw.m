classdef Eobj_raw < handle
methods
    function [flds,vals]=get_fnames_cell(obj)
        [flds,vals]=struct2Cell(obj.fnames);
    end
    function [flds,vals]=get_raw_fnames_cell(obj)
        [flds,vals]=struct2Cell(obj.rawData);
        inds=cellfun(@(x) strcmp(x{end},'data'),flds);
        flds=flds(inds);
        vals=vals(inds);
    end
    function name=gen_fname_raw_new(obj,subj,mode,std,blk)
        name=obj.gen_fname_block_new(subj,mode,std,blk);
    end
    function name=gen_fname_raw(obj,subj,mode,std,blk,r1,r2)
        name=obj.gen_fname_block(subj,mode,std,blk,r1,r2);
    end
    function name=gen_name_raw_new(obj,subj,method,std,blk)
        name=obj.gen_name_block_new(subj,method,std,blk);
    end
    function name=gen_name_raw(obj,subj,mode,std,blk,r1,r2)
        name=obj.gen_name_block(subj,mode,std,blk,r1,r2);
    end
    function [r1,r2]=get_r1_r2_from_fname(obj,fname)
        [~,fname,~]=filePartsSane(fname);
        spl=strsplit(fname,'_');
        spl=strsplit(spl{end},'-');
        r1=str2double(spl{2});
        r2=str2double(spl{3});
    end
    function fname=gen_fname_block(obj,subj,mode,std,blk,r1,r2)
        name=obj.gen_name_block(subj,mode,std,blk,r1,r2);
        fname=[autodir(obj.dir.RAW) name];
    end
    function fname=gen_fname_block_new(obj,subj,mode,std,blk)
        name=obj.gen_name_block_new(subj,mode,std,blk);
        fname=[autodir(obj.dir.RAW) name];
    end
    function fname=gen_fname_block_cur(obj,subj,mode,std,blk)
        name=obj.gen_name_block_cur(subj,mode,std,blk);
        fname=[autodir(obj.dir.RAW) name];
    end
    function name=gen_name_block(obj,subj,mode,std,blk,r1,r2)
        subj=obj.auto_subj(subj);
        mode=obj.auto_mode(mode);
        std=obj.auto_std_str(std);
        blkstr=obj.auto_blk_str(blk);
        r1=obj.auto_int_str(r1);
        r2=obj.auto_int_str(r2);
        name=['raw_' obj.name und mode und subj und std und blkstr dsh r1 dsh r2 '.mat'];
    end
    function name=gen_name_block_cur(obj,subj,method,std,blk)
        if obj.isempty_raw(subj,method,std,blk)
            r1=0;
        else
            r1=obj.get_block_redo(subj,method,std,blk);
        end
        if ~isfield(obj.rawData,'redo')
            r2=0;
        else
            r2=obj.rawData.redo.nRedo();
        end
            
        name=obj.gen_name_block(subj,method,std,blk,r1,r2);
    end
    function name=gen_name_block_new(obj,subj,method,std,blk)
        % accounts for subj redo or empty
        if obj.isempty_raw(subj,method,std,blk)
            r1=0;
        else
            r1=obj.get_block_redo(subj,method,std,blk)+1;
        end
        r2=obj.rawData.redo.nRedo();
        name=obj.gen_name_block(subj,method,std,blk,r1,r2);
    end
%% LOAD
    function [files,INDS] = get_raw_files(obj,subj,mode,std,blk)
        if ~exist('subj','var') || isempty(subj)
            subj=obj.subjs;
        elseif ischar(subj)
            subj={subj};
        elseif isint(subj)
            subj=obj.subj(subj);
        end
        if ~exist('std','var') || isempty(std)
            std=obj.get_std_flds();
        elseif ischar(std)
            std={auto_std_num(std)};
        elseif isint(std)
            std=obj.get_std_fld(std);
        end
        if ~exist('mode','var') || isempty(mode)
            mode={'test','train','pilot'};
        end
        if ~exist('blk','var') || isempty(blk)
            blk=1:obj.nBlk;
        end
        INDS=distribute(subj,mode,std,blk);
        files=cell(size(INDS,1),1);
        for i = 1:size(INDS,1)
            files{i,1}=obj.get_raw_file(INDS{i,:});
        end
    end
    function [Sall, files, INDS]=load_raw_test_blocks(obj,subjs,stds,blks,bNested)
        if ~exist('bNested','var') || isempty(bNested)
            bNested=0;
        end
        if ~exist('subjs','var') || isempty(subjs)
            subjs=obj.subjs;
        elseif ischar(subjs)
            subjs={subjs};
        elseif isint(subjs)
            subjs=obj.subjs(subjs);
        end
        if ~exist('stds','var') || isempty(stds)
            stds=obj.get_std_flds();
        elseif ischar(stds)
            stds={obj.auto_std_num(stds)};
        elseif isint(stds)
            stds=obj.get_std_fld(stds);
        end
        if ~exist('blks','var') || isempty(blks)
            blks=1:obj.nBlk;
        end
        [files,INDS]=obj.get_raw_test_files(subjs,stds,blks);
        Sall=cell(size(files,1),1);
        %p=pr(length(files));
        if ~bNested; p=pr(length(files),[],'loading blocks','...'); end
        for i = 1:length(files)
            if ~bNested
                p.u();
            end
            if isempty(files{i})
                continue
            end
            if endsWith(files{i},'.mat')
                Sall{i}=load([obj.dir.RAW files{i}]);
            else
                Sall{i}=load([obj.dir.RAW files{i} '.mat']);
            end
            flds=fieldnames(Sall{i});
            if numel(flds)==1
                Sall{i}=Sall{i}.(flds{1});
            end

        end
        if ~bNested; p.c(); end
    end
    function [files,INDS] = get_raw_test_files(obj,subj,std,blk)
        if ~exist('subj','var') || isempty(subj)
            subj=obj.subjs;
        elseif ischar(subj)
            subj={subj};
        elseif isint(subj)
            subj=obj.subj(subj);
        end
        if ~exist('std','var') || isempty(std)
            std=obj.get_std_flds();
        elseif ischar(std)
            std={obj.auto_std_num(std)};
        elseif isint(std)
            std=obj.get_std_fld(std);
        end
        if ~exist('blk','var') || isempty(blk)
            blk=1:obj.nBlk;
        end
        [files,INDS]=obj.get_raw_files(subj,'test',std,blk);
    end
    function [file] = get_raw_file(obj,subj,mode,std,blk)
        [subj,mode,std,blk]=obj.auto_fld(subj,mode,std,blk);
        file=obj.rawData.(subj).(mode).(std).data{blk};
    end
    function [files,IND,ind] = get_raw_files_all(obj,bSqueeze)
        if ~exist('bSqueeze','var') || isempty(bSqueeze)
            bSqueeze=0;
        end
        IND=obj.get_all_data_combs();
        files=cell(size(IND,1),1);
        for i = 1:size(IND,1)
            ind=IND(i,:);
            fname=obj.get_block_fname(ind{:});
            if ~isempty(fname)
                [~,fname,~]=filePartsSane(fname);
                files{i,1}=fname;
            end
        end
        if nargout < 3 && ~bSqueeze
            return
        end
        ind=cellfun(@isempty,files);
        if bSqueeze()
            files(ind)=[];
            IND(ind,:)=[];
        end

    end
    function raw=load_block(obj,subj,method,std,blk,bSkipCheck)
        if ~exist('bSkipCheck','var') || isempty(bSkipCheck)
            bSkipCheck=1;
        end
        out=obj.check_block(subj,method,std,blk,bSkipCheck);
        if ~out
            fname=obj.get_block_fname(subj,method,std,blk);
            if isempty(fname)
                subj
                method
                std
                blk
            end

            raw=load([fname '.mat']);
            if isfield(raw,'raw')
                raw=raw.raw;
            end
        end
    end
%% RAW FNAME
    function out=isempty_raw(obj,subj,mode,std,blk)
        subj=obj.auto_subj(subj);
        mode=obj.auto_mode(mode);
        std=obj.auto_std_fld(std);
        blk=obj.auto_blk_num(blk);
        out=isempty(obj.rawData.(subj).(mode).(std).data{blk});
    end
    function [r1,r2]=get_block_redo(obj,subj,mode,std,blk)
        subj=obj.auto_subj(subj);
        mode=obj.auto_mode(mode);
        std=obj.auto_std_fld(std);
        blk=obj.auto_blk_num(blk);
        r1=obj.rawData.(subj).(mode).(std).redo(blk);
        if ~isfield(obj.rawData,'redo')
            r2=0;
        else
            dk
        end
    end
    function num=get_block_flag(obj,subj,mode,std,blk)
        subj=obj.auto_subj(subj);
        mode=obj.auto_mode(mode);
        std=obj.auto_std_fld(std);
        blk=obj.auto_blk_num(blk);
        num=obj.rawData.(subj).(mode).(std).flag(blk);
    end
    function names=get_raw_names(obj,subjs,mthds,stds,blks)
        if ~exist('subjs','var')
            subjs=[];
        end
        if ~exist('mthds','var')
            mthds=[];
        end
        if ~exist('stds','var')
            stds=[];
        end
        if ~exist('blks','var')
            blks=[];
        end
        names=get_block_names(obj,subjs,mthds,stds,blks);
    end
    function fnames=get_raw_fnames(obj,subjs,mthds,stds,blks)
        if ~exist('subjs','var')
            subjs=[];
        end
        if ~exist('mthds','var')
            mthds=[];
        end
        if ~exist('stds','var')
            stds=[];
        end
        if ~exist('blks','var')
            blks=[];
        end
        fnames=get_block_fnames(obj,subjs,mthds,stds,blks);
    end
    function fnames=get_block_fnames(obj,subjs,mthds,stds,blks)
        if ~exist('subjs','var') || isempty(subjs)
            subjs=obj.subjs;
        elseif ischar(subjs)
            subjs={subjs};
        end
        if ~exist('mthds','var') || isempty(mthds)
            mthds=obj.modeflds;
        elseif ischar(mthds)
            mthds={mthds};
        end
        if ~exist('stds','var') || isempty(stds)
            stds=1:obj.nStd;
        end
        if ~exist('blks','var') || isempty(blks)
            blks=1:obj.nBlk;
        end

        IND=distribute(subjs,mthds,stds,blks);
        fnames=cell(size(IND,1),1);
        for i = 1:size(IND,1)
            ind=IND(i,:);
            fnames{i}=obj.get_block_fname(ind{:});
        end
        fnames(cellfun(@isempty,fnames))=[];
    end
    function names=get_block_names(obj,subjs,mthds,stds,blks)
        if ~exist('subjs','var') || isempty(subjs)
            subjs=obj.subjs;
        elseif ischar(subjs)
            subjs={subjs};
        end
        if ~exist('mthds','var') || isempty(mthds)
            mthds=obj.modeflds;
        elseif ischar(mthds)
            mthds={mthds};
        end
        if ~exist('stds','var') || isempty(stds)
            stds=1:obj.nStd;
        end
        if ~exist('blks','var') || isempty(blks)
            blks=1:obj.nBlk;
        end

        IND=distribute(subjs,mthds,stds,blks);
        names=cell(size(IND,1),1);
        for i = 1:size(IND,1)
            ind=IND(i,:);
            names{i}=obj.get_block_name(ind{:});
        end
        names(cellfun(@isempty,names))=[];
    end
    function fname=get_block_fname(obj,subj,method,std,blk)
        subj=obj.auto_subj(subj);
        method=obj.auto_mode(method);
        std=obj.auto_std_fld(std);
        blk=obj.auto_blk_num(blk);
        fname=obj.rawData.(subj).(method).(std).data{blk};
        if isempty(fname)
            return
        end
        fname=[obj.dir.RAW fname];
    end
    function name=get_block_name(obj,subj,method,std,blk)
        subj=obj.auto_subj(subj);
        method=obj.auto_mode(method);
        std=obj.auto_std_fld(std);
        blk=obj.auto_blk_num(blk);
        name=obj.rawData.(subj).(method).(std).data{blk};
        if isempty(name)
            return
        end
        name=strrep(name,'.mat','');
    end
    function val=get_block_val(obj,subj,method,std,blk)
        subj=obj.auto_subj(subj);
        method=obj.auto_mode(method);
        std=obj.auto_std_fld(std);
        blk=obj.auto_blk_num(blk);
        val=obj.rawData.(subj).(method).(std).blk(blk);
    end
    function out =check_raw_save(obj,fname,subj,mode,std,blk)
        [codeDB,codeFN,fnameDB,fnameGen]=obj.check_raw(0,subj,mode,std,blk);
        if strcmp(fnameDB,fnameGen)
            disp(['save name: ' fnameDB])
            disp(['gen  name: ' fnameGen])
            out1=basicYN('These two do not match. Continue saving under gen name?');
            if ~out1
                out=0;
                return
            end
        end
        if ismember('-',codeDB) &&  ismember('E',codeFN)
            out=1;
        else
            error('write code')
        end
        %elseif ~ismember('E',codeFN) && ~ismember('E',codeDB)
        %    out=basicYN(['File exists, and in database. Overwrite ' fname ' ?']);
        %elseif ismember('E',codeFN)  &&  ~ismember('E',codeDB)
        %    out=basicYN(['DB says ' fname ' exists, but file cannot be found. Continue saving?']);
        %elseif ~ismember('E',codeFN)  &&  ismember('E',codeDB)
        %    out=basicYN(['File exists, but not in db. Overwrite ' fname ' ?']);
        %else
        %    out=0;
        %    chkAlPrint(codeFN);
        %end
    end
    function [codeDB,codeFN,fnameDB,fnameGen]=check_raw(obj,bPrint,subj,mode,std,blk)
        % NOTE STANDARD FOR DATA
        % MM = mismatche
        if ~exist('bPrint','var') || isempty(bPrint)
            bPrint=0;
        end

        % check db file
        [codeDB,fnameDB]=obj.check_raw_db(subj,mode,std,blk);
        [codeFN,fnameGen]=obj.check_raw_gen(subj,mode,std,blk);


        if bPrint
            disp('codeDB:')
            str=chkAllPrint(codeDB,fulldir,highestdir,fname);
            str=indent(str,4);
            disp(str)

            str=chkAllPrint(codeFN,fulldir,highestdir,fname);
            str=indent(str,4);
            disp(str);
        end

    end
    function [code,fname]=check_raw_db(obj,subj,mode,std,blk,bSkipDir)
        if ~exist('bSkipDir','var') || isempty(bSkipDir)
            bSkipDir=0;
        end
        fname=obj.get_block_fname(subj,mode,std,blk);
        if isempty(fname)
            code='-';
            return
        end
        [~,code]=chkFileAll(fname,0,bSkipDir);
    end
    function [code,fname]=check_raw_gen(obj,subj,mode,std,blk,r1,r2,bSkipDir)
        if ~exist('r1','var')
            r1=[];
        end
        if ~exist('r2','var')
            r2=[];
        end
        if ~exist('bSkipDir','var') || isempty(bSkipDir)
            bSkipDir=0;
        end
        %% 0  file doesn't exist
        %% 1  file exists
        if xor(isempty(r1),isempty(r2))
            error('both r1 and r2 must be specifies')
        elseif ~isempty(r1) && ~isempty(r2)
            fname=obj.gen_fname_block(subj,mode,std,blk,r1,r2);
        else
            fname=obj.gen_fname_block_cur(subj,mode,std,blk);
        end
        if isempty(fname)
            code='-';
        else
            [~,code]=chkFileAll(fname,1,bSkipDir);
        end
    end

    function code=check_raw_dir(obj)
        dirn=obj.dir.RAW;
        code=chkDirAll(dirn,1);
    end
%%
end
end
