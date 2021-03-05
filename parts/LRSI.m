
    function codes=get_LRSI_data_codes(obj)
        flds={'CP','CL','BI','PRJ','AMA'};
        codes=zeros(length(flds),1);
        for i = 1:length(flds)
            [~,codes(i)]=check_data('IN',fld);
        end
    end
    function status=check_LRSI_data(obj)
        %% status
        % 0 - listed, missing data
        % 1 - all good
        % 2 - implied missing data
        codes=obj.get_LRSI_data_codes();
        status=codess==1 | codes==-1;
        if status==0
            return
        end

        emptyInd=find(codes==-1);
        existInd=find(codes==1);

        % implied missing data
        if isempty(existInd) || isempty(emptyInd)
            missing=zeros(size(status));
        else
            missing = (emptyInd < min(existInd));
        end

        if sum(missing) == 0 && status==1
            return
        elseif   sum(missing) > 0 && status==1
            status=2;
        end
    end
    function next=get_next_LRSI_mod(obj)
        codes=obj.get_LRSI_data_codes();
        next=find(codes==0,1,'first');
    end
    function ind=get_overwrite(obj,mods)
        codes=get_LRSI_data_codes();
        ind=ismember(1:5,mods) & codes==1;
    end
    function out=run_LRSIgen(obj,mods,bForceOverwrite,bTest)
        if ~exist('var','bTest') || isempty(bTest)
            bTest=0;
        end

        if ~exist('var','bForceOverwrite') || isempty(bForceOverwrite)
            bForceOverwrite=0;
        end
        %% run data
        %% initialize LRSI
        %% determine next to be ran
        %% promppt continue or restart
        %% run


        % CHECK DIRECTORY
        out=obj.check_data_dir('IN',1);
        if out~=1
            return
        end

        % CHECK DATA
        status=obj.check_LRSI_data();
        if status == 0
            disp('Missing data. Fix then try again.')
        elseif status == 2
            disp('Implied missing data.')
            resp=basicYN('Continue?');
            if resp == 0
                return
            end
        end


        % make sure mod values are valid
        if ~isempty(mods)
        end

        % check overwrite
        if ~isempty(mods)
            next=get_next_LRSI_mod;
            mod=next:4;
        else
            ind=obj.get_overwrite(mods);
            if sum(ind)~=0 && ~bforceOverwrite
                modStr=join(obj.LRSIflds{ind},' ');
                modStr=modStr{1};
                disp('Data exists for the following mods:');
                disip(['     ' modStr ]);
                resp=basicYN('Overwrite? If no, existing data will be used as necessary.');
                if resp==0
                    mods(ind)=[];
                end
            end
        end


        % check if generated don't exist
        % XXX

        % CHECK DEFS
        % XXX

        obj.LRSI_main(1,mods);
        obj.LRSI_main(0,mods);
    end
    function obj=LRSI_main(obj,bTest,mods)

        % package_def
        Opts=package_def();

        % INIT
        obj.init_LRSI(Opts);

        mods=mods(:)';
        for mod = mods
            obj.run_LRSI(mod,name,bTest);
            obj.save_LRSI(mod,bTest);
            obj.update_LRSI_META(mod,bTest);
            if ~bTest
                %% NOTE NOT RAN DURING TEST
                obj.save(bTest);
            end
        end
    end
    function Opts=package_def(obj,mod)
        Esess=struct();
        Esess.mod=mod;
        clear mod
        Esses.modStr=obj.LRSIflds{mod};
        Esess.fname=[obj.dir.DEF obj.OUT.DEF.(Esess.modStr)];
        run(Esess.fname);
        vars=who;
        vars(ismember(vars,{'Esess','obj'}))=[];
        Opts=struct();
        for i = 1:length(vars)
            var=vars{i};
            Opts.(var)=eval([var ';']);
        end
    end
    function obj=init_LRSI(obj,Opts)
        obj.LRSI=LRSIgen(Opts);
    end
    function obj=run_LRSI_mod(obj,mod,fname,bTest)
        if ~exist('bTest','var') || isempty(bTest)
            bTest=0;
        end
        modStr=obj.LRSIflds{mod};
        obj.LRSI.fnames.(modStr)=fname;

        fname=obj.get_LRSI_base_fname(mod);
        if bTest
            obj.LRSI.(modStr)=struct();
        else
            %% NOTE NOT RAN DURING TEST
            obj.LRSI.(['mod' num2str(mod)]);
        end
    end
    function name=save_LRSI(obj,mod,bTest)
        if ~exist('bTest','var') || isempty(bTest)
            bTest=0;
        end
        modStr=obj.LRSIflds{mod};
        fname=obj_get_LRSI_fname(obj,mod);

        eval([modStr '=obj.LRSI.(modStr));']);

        eval([' save(fname,''' modStr ''');']);
        if bTest
            delete([fname '.mat']);
        end
    end
    function obj=update_LRSI_meta(obj,mod,bTest)
        if ~exist('bTest','var') || isempty(bTest)
            bTest=0;
        end
        old=obj.fnames.IN.(modStr);
        modStr==obj.LRSIflds{mod};
        name=obj.get_LRSI_base_fname(mod);
        obj.fnames.IN.(modStr)=name;
        if bTest
            obj.fnames.IN.(modStr)=old;
        end
    end
    function name=get_LRSI_base_fname(obj,mod)
        modStr=obj.LRSIflds{mod};
        error('fix this')
        name=obj.gen_fname_data('IN',modStr);
    end
    function name=get_LRSI_fname(obj,mod)
        name=obj.get_LRSI_base_fname(mod);
        fname=[obj.dir.IN name];
    end
%% DEF
    function codes=get_LRSI_def_codes(obj) %% -1 empty name
        %% 0  file doesn't exist
        %% 1  file exists
        %% 2  direcotory doesnt exist
        %% 2  read but not write
        flds={'CP','CL','BI','PRJ','AMA'};
        codes=zeros(length(flds),1);
        for i = 1:length(flds)
            [~,codes(i)]=check_data('DEF',fld);
        end
    end
    function out=check_consistent_def_fnames(obj)
        flds={'CP','CL','BI','PRJ','AMA'};
        fnames=cell(size(flds,1),1);
        SEEN={};
        for i = 1:length(flds)
            fnames{i}=obj.fnames.DEF.(fld);
        end
        [unq,~,ind]=unique(fnames);
        % XXX
    end
    function check_LRSI_defs(obj)
    end
