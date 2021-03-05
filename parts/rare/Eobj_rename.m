classdef Eobj_rename < handle
methods
    function obj=rename_pass(obj,in,bCopy)
        if ~exist('bCopy','var') || isempty(bCopy)
            bCopy=0;
        end
        re= '[0-9]+';
        if exist('in','var') && ~isempty(in)
            if ~regExp(in,re)
                error('Invalid rename string')
            else
                out=in;
            end
        else
            while true
                out=input('input pass num [0-9]: ','s');
                if regExp(out,re)
                    break
                end
                disp('Invalid input')
            end
        end
        val=str2double(out);
        new=['_pass' out];

        re='_pass[0-9]+';
        obj.rename_fun_p(new,re,val,'pass',bCopy);
    end
    function obj=rename_prjInd(obj,in,bCopy)
        if ~exist('bCopy','var') || isempty(bCopy)
            bCopy=0;
        end
        re= '[0-9]+(-|d)[0-9]+';
        if exist('in','var') && ~isempty(in)
            if ~regExp(in,re)
                error('Invalid rename string')
            else
                out=in;
            end
        else
            while true
                out=input('input prjInd str [X-X]: ','s');
                if regExp(out,re)
                    break
                end
                disp('Invalid input')
            end
        end
        val=str2double(split(out,'-'))';

        obj.rename_fun_p(out,re,val,'prjInd',bCopy);
    end
%end
%methods(Access=private)
    function obj=rename_fun_p(obj,new,re,val,NAME,bCopy)
        fun=@(x) sed('s',x,re,new);
        obj.rename_fnames_fun_p(fun,new,bCopy);
        obj.rename_raw_fun_p(fun,new,bCopy);
        obj.rename_db_fun_p(re,new,val,NAME,bCopy);
    end
    function obj=rename_db_fun_p(obj,re,new,val,NAME,bCopy)
        % names
        oldname=obj.name;
        newname=sed('s',oldname,re,new);

        % fnames
        dire=obj.dir.META;
        oldbfile=[dire oldname '.mat'];
        newbfile=[dire newname '.mat'];
        [dire,name,ext]=filePartsSane(oldbfile);
        backupfile=[dire filesep 'backup' filesep name ext];

        % check exist
        if ~bCopy && ~exist(oldbfile,'file')
            error(['db file doesn''t exist ?????'])
        end

        %rename db
        obj.name=newname;


        %rename files
        if ~bCopy
            movefile(oldbfile,backupfile);
            if ~exist(backupfile,'file')
                error('backup db file not created!')
            end
        end

        obj.(NAME)=val;
        obj.save();

        % validate files moved
        if ~exist(newbfile,'file')
            error('new db file not saved !???')
        end

    end
    function obj=rename_fnames_fun_p(obj,fun,new,bCopy)
        % names
        [flds,names]=obj.get_fnames_cell;

        % ----------
        % MANY (names that are blk specific)
        bmany=cellfun(@iscell,names);
        fldsMany=flds(bmany,:);
        namesMany=names(bmany);
        obj.rename_fun_blk_p(fun,new,fldsMany,namesMany,bCopy);

        % ----------
        % single
        bsingle=cellfun(@numel,names)>1 & ~bmany;
        namesSing=cellfun(@(x) {x} ,names(bsingle),UO,false);
        fldsSing=flds(bsingle,:);
        obj.rename_fun_blk_p(fun,new,fldsSing,namesSing,bCopy);

    end
    function obj=rename_raw_fun_p(obj,fun,new,bCopy)
        ext='.mat';
        dire=obj.dir.RAW;
        [flds,oldNames]=get_raw_fnames_cell(obj);
        obj.rename_fun_blk_p(fun,new,flds,oldNames,bCopy);
    end
    function obj=rename_fun_blk_p(obj,fun,new,flds,oldNames,bCopy)
        % names
        newNames=oldNames;
        for i = 1:length(flds)
            newNames{i}=cellfun(fun,oldNames{i},UO,false);
        end
        % fnames
        braw=0;
        bSuccess=0;
        %FLDS=cell(0);
        for i = 1:length(oldNames)
            if strcmp(flds{i}{1},'DEF')
                dire=[obj.dir.(flds{i}{end-1}) flds{i}{end}];
            elseif strcmp(flds{i}{end},'data')
                braw=1;
                dire=obj.dir.RAW;
            elseif strcmp(flds{i}{1},'CODE')
                dire='null';
            else
                dire=obj.dir.(flds{i}{end-1});
            end
            if dire(end) ~= '/'
                dire=[dire '/'];
            end

            if  ismember(flds{i}{end-1},{'DEF','CODE'})
                ext='.m';
            else
                ext='.mat';
            end
            n=length(oldNames{i});
            f=num2cell(repmat(flds{i},n,1),2);
            %FLDS=[FLDS; f];
            of=cell(n,1);
            nf=cell(n,1);
            for j = 1:n
                if ~isempty(oldNames{i}{j})
                    [~,oldNames{i}{j}]=filePartsSane(oldNames{i}{j});
                end
                of{j}=Eobj.extfun(dire,oldNames{i}{j},ext);
                nf{j}=Eobj.extfun(dire,newNames{i}{j},ext);
                if ~exist(nf{j},'file') && ~exist(of{j},'file')
                    nf{j}=of{j};
                    newNames{i}{j}=oldNames{i}{j};
                else
                    nf{j}=Eobj.extfun(dire,newNames{i}{j},ext);
                end
                if strcmp(ext,'.m')
                    newNames{i}{j}=strrep(newNames{i}{j},'-','d');
                    nf{j}=strrep(nf{j},'-','d');
                end
                oldFnames{i,1}=of;
                newFnames{i,1}=nf;
                bSuccess=1;
            end
        end
        if ~bSuccess
            return
        end
        if braw
            fld='rawData';
        elseif isfield(obj.fnames,'DEF')
            fld='fnames';
        else
            error('Where do I put the data ?')
        end

        newFnames=vertcat(newFnames{:});
        oldFnames=vertcat(oldFnames{:});
        %newNames=vertcat(newNames{:});
        %oldNames=vertcat(oldNames{:});

        % SET DB
        % NOTE newNames AND flds DO NOT NEED TO BE EXPANDED!
        fn=obj.(fld);
        for i = 1:length(newNames)
            fn=structSet(fn,flds{i},newNames{i});
        end
        obj.(fld)=fn;

        % mv files
        for i = 1:length(newFnames)

            if isequal(oldFnames{i},newFnames{i})
                continue
            end

            % TEST
            %disp(oldFnames{i})
            %disp(newFnames{i})
            %
            if ~exist(oldFnames{i},'file');
                continue
            end

            if bCopy
                copyfile(oldFnames{i},newFnames{i});
            else
                movefile(oldFnames{i},newFnames{i});
            end
            e=exist(newFnames{i},'file');
            if ~e
                error('new file not created ???')
            end
        end


    end
end
methods(Static=true,Access=private)
    function out=extfun(dire,x,ext)
        out=x;
        if isempty(x)
            return
        end
        dr=filePartsSane(x);
        out=strrep(out,dr,'');
        out=strrep(out,'/','');

        if ~startsWith(x,dire)
            out=[dire out];
        end
        if ~endsWith(x,ext)
            out=[out ext];
        end
    end
end
end
