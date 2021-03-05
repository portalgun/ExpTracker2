classdef Eobj_init < handle
methods
%% INIT
    function obj=init_subj_status(obj)
        % TODO
    end
    function obj=init_flags(obj)
        % TODO
    end
    function obj=init_exp_empty_all(obj)
        INDS=distribute(obj.modeflds,1:obj.nStd);
        for i = 1:size(INDS,1)
            ind=INDS(i,:);
            obj.init_exp_empty(ind{:});
        end
    end
    function obj=init_exp_empty(obj,mode,std)
        mode=obj.auto_mode(mode);
        std=obj.auto_std_fld(std);
        if ~isfield(obj.expData,mode)
            obj.expData.(mode)=struct();
        end
        if ~isfield(obj.expData.(mode),std)
            obj.expData.(mode).(std)=cell(obj.nBlk,1);
        end
    end
    function obj=init_raw_empty_all(obj)
        IND=obj.get_all_data_combs();
        for i = 1:length(IND)
            ind=IND(i,:);
            obj.init_raw_empty(ind{:});
        end
    end
    function obj=init_raw_empty(obj,subj,method,std,blk)
        if ~isfield(obj.rawData,subj)
            obj.rawData.(subj)=struct();
        end
        if ~isfield(obj.rawData.(subj),method)
            obj.rawData.(subj).(method)=struct();
        end
        if ~isfield(obj.rawData.(subj).(method),std)
            obj.rawData.(subj).(method).(std)=struct();
        end
        if ~isfield(obj.rawData.(subj).(method).(std),'blk')
            obj.rawData.(subj).(method).(std).blk=zeros(1,obj.nBlk);
        end
        if ~isfield(obj.rawData.(subj).(method).(std),'rerunCount')
            obj.rawData.(subj).(method).(std).rerunCount=zeros(1,obj.nBlk);
        end
        if ~isfield(obj.rawData.(subj).(method).(std),'flag')
            obj.rawData.(subj).(method).(std).flag=zeros(1,obj.nBlk);
        end
        if ~isfield(obj.rawData.(subj).(method).(std),'redo')
            obj.rawData.(subj).(method).(std).redo=zeros(1,obj.nBlk);
        end
        if ~isfield(obj.rawData.(subj).(method).(std),'data') || isempty(obj.rawData.(subj).(method).(std).data)
            obj.rawData.(subj).(method).(std).data=cell(obj.nBlk,1);
        elseif obj.rawData.(subj).(method).(std).redo(blk)==0 && ~isempty(obj.rawData.(subj).(method).(std).data{blk})
            fname=obj.rawData.(subj).(method).(std).data{blk};
            [obj.rawData.(subj).(method).(std).redo(blk),~]=obj.get_r1_r2_from_fname(fname);
        end
    end
%%
    function obj=init_rnd(obj)
        if ~isfield(obj.rnd,'master') || isempty(obj.rnd.master)
            obj.gen_rnd_master();
        end
        if ~isfield(obj.rnd,'blk') || isempty(obj.rnd.blk)
            obj.gen_rnd_blk();
        end

        if ~isfield(obj.rnd,'trl') || isempty(obj.rnd.trl)
            obj.gen_rnd_trl();
        end
        if ~isfield(obj.rnd,'cmp') || isempty(obj.rnd.cmp)
            obj.gen_rnd_cmp();
        end
        if ~isfield(obj.rnd,'intrvl') || isempty(obj.rnd.intrvl)
            obj.gen_rnd_intrvl();
        end
        obj.apply_rnd();
    end
%%
    function obj=init_dirs(obj)
        obj.dir.META=BLdirs('EXP');
        if obj.bLegacy
            obj.dir.META=[obj.dir.META 'legacy' filesep];
        end
        % TODO: make into config file or dbDirs
        %% TODO MAKE DIRS IF DON't EXIST
        root='/Volumes/Data/.daveDB/';
        obj.dirs.BONICE=struct();
        obj.dirs.BONICE.OUT=[root 'out/'];
        obj.dirs.BONICE.RAW=[root 'raw/'];
        obj.dirs.BONICE.IN =[root 'in/'  obj.imgDTB '/'];
        obj.dirs.BONICE.MAP=[root 'img/' obj.imgDTB '/'];
        obj.dirs.BONICE.EXP=[root 'exp/' obj.imgDTB '/'];
        obj.dirs.BONICE.DEF='/home/dambam/Code/mat/projects/_expVariables/';
        obj.dirs.BONICE.TMP='/home/dambam/tmp/mat/';
        obj.dirs.jburge_wheatstone=obj.dirs.BONICE;
        obj.dirs.jburge_wheatstone.DEF='/homes/davwhite/Code/mat/projects/_expVariables/';
        obj.dirs.jburge_wheatstone.TMP='/homes/davwhite/tmp/mat/';

        root='Z:\Data\.daveDB\';
        obj.dirs.BONICH=struct();
        obj.dirs.BONICH.OUT=[root 'out\'];
        obj.dirs.BONICH.RAW=[root 'raw\'];
        obj.dirs.BONICH.IN =[root 'in\'  obj.imgDTB '\'];
        obj.dirs.BONICH.MAP=[root 'img\' obj.imgDTB '\'];
        obj.dirs.BONICH.EXP=[root 'exp\' obj.imgDTB '\'];
        obj.dirs.BONICH.DEF='Y:\mat\toolboxes\_expVariables\';
        obj.dirs.BONICH.TMP='Y:\tmp\mat';

        obj.dir=obj.dirs.(strrep(hostname,'-','_'));
        obj.dir.META=dbDirs('EXP');
        obj.dir.prjDir=[obj.prjCode filesep];
        % obj.dir.CODE handled by px
    end
    function obj=convert_dirs(obj)
        flds=fieldnames(obj.dir);
        for i = 1:length(flds)
            fld=flds{i};
            obj.dir.(fld)=autodir(obj.dir.(fld));
        end
    end
%%
    function obj=init_fnames(obj)
        obj.init_fname_out();
        obj.init_fname_in();
        obj.init_fname_def();
        obj.init_fname_code();
        obj.init_fname_map();
    end
    function obj=init_fname_out(obj)
        flds={'cleaned','merged','pruned','prunedP','tabled','cleaned_ind','pruned_ind','prunedP_ind','tabled_ind'};

        if ~isfield(obj.fnames,'OUT')
            obj.fnames.OUT=struct();
        end
        for i = 1:length(flds)
            fld=flds{i};
            if ~isfield(obj.fnames.OUT,fld) && endsWith(fld,'_ind')
                obj.fnames.OUT.(fld)=cell(numel(obj.subjs),3,obj.nStd); %3 for modes
            elseif ~isfield(obj.fnames.OUT,fld)
                obj.fnames.OUT.(fld)='';
            end
        end
    end
    function obj=init_fname_map(obj)
        flds={'VET','GEN','BIN','SMP','PCH'};
        if ~isfield(obj.fnames,'MAP')
            obj.fnames.MAP=struct();
        end
        for i = 1:length(flds)
            fld=flds{i};
            if ~isfield(obj.fnames.MAP,fld)
                obj.fnames.MAP.(fld)='';
            end
        end
    end
    function obj=init_fname_in(obj)
        %flds={'IN','CP','CL','BI','PRJ_TEST','PRJ_TRAIN','PRJ_PILOT','EXP','EXPi'};
        %obj.fnames.IN.BI=cell(0);
        %obj.fnames.IN.CL=cell(0);
        %obj.fnames.IN.PRJ_EXPi=cell(0);
        flds={'SEL_PILOT','SEL_TRAIN','SEL_TEST','DMP_PILOT','DMP_TRAIN','DMP_TEST','EXP_PILOT','EXP_TRAIN','EXP_TEST',};

        if ~isfield(obj.fnames,'IN')
            obj.fnames.IN=struct();
        end
        for i = 1:length(flds)
            fld=flds{i};
            if ~isfield(obj.fnames.IN,fld)
                obj.fnames.IN.(fld)='';
            end
        end
        fldscur=fieldnames(obj.fnames.IN);
        for i = 1:length(fldscur)
            fld=fldscur{i};
            if isempty(obj.fnames.IN.(fld)) && ~ismember(fld,flds)
                obj.fnames.IN=rmfield(obj.fnames.IN,fld);
            end
        end
    end
    function obj=init_fname_def(obj)
        if ~isfield(obj.fnames,'DEF')
            obj.fnames.DEF=struct();
        end
        %flds={'CP','CL','BI','PRJ','AMA','EXP',};
        flds={'IMAP','MODEL','EXP'};
        % imap model exp
        for i = 1:length(flds)
            fld=flds{i};
            if ~isfield(obj.fnames.DEF,fld)
                obj.fnames.DEF.(fld)='';
            end
        end
        fldscur=fieldnames(obj.fnames.DEF);
        for i = 1:length(fldscur)
            fld=fldscur{i};
            if isempty(obj.fnames.DEF.(fld)) && ~ismember(fld,flds)
                obj.fnames.DEF=rmfield(obj.fnames.DEF,fld);
            end
        end

    end
    function obj=init_fname_code(obj)
        if ~isfield(obj.fnames,'CODE')
            obj.fnames.CODE=struct();
        end
        if ~isfield(obj.fnames.CODE,'EXP')
            obj.fnames.CODE.EXP='';
        end
        if ~isfield(obj.fnames.CODE,'AMA')
            obj.fnames.CODE.AMA='';
        end
    end
% RAND RND
    function obj=gen_rnd_new(obj)
        out=basicYN('Are you sure you want to regen master rnd?');
        if out == 0
            return
        end
        out=basicYN('ARE YOU REALLY SURE YOU WANT TO REGEN MASTER RND?');
        if out == 0
            return
        end
        obj.gen_rnd_master();
        obj.gen_rnd_blk();
        obj.gen_rnd_trl();
        obj.gen_rnd_blk();
        obj.apply_rnd();
    end
    function obj=gen_rnd_master(obj)
        rng('shuffle');
        obj.rnd.master=rng(2^32-1);
    end
    function obj=gen_rnd_blk(obj)
        % nblk val to randomize block order
        rng(obj.rnd.master,'twister');
        obj.rnd.blk=randi(2^32-1,obj.nBlk,1);
    end
    function obj=gen_rnd_trl(obj)
        % 1 val to randomize all trials
        rng(obj.rnd.master,'twister');
        val=randi(2^32-1,obj.nBlk+1,1);
        obj.rnd.trl=val(end);
    end
    function obj=gen_rnd_cmp(obj)
        rng(obj.rnd.master,'twister');
        val=randi(2^32-1,obj.nBlk+2,1);
        obj.rnd.cmp=val(end);
    end
    function obj=gen_rnd_intrvl(obj)
        % nStd x nBlk to randomize intervals
        rng(obj.rnd.master,'twister');
        randi(2^32-1,obj.nBlk+2,1);
        obj.rnd.intrvl=randi(2^32,obj.nStd,obj.nBlk);
    end
%% BLK TABLE GEN
    function obj=get_blkTable(obj)
        if ~isfield(obj.rnd,'blkRndAll') || isempty(obj.rnd.blkRndAll)
            obj.rnd.blkRndAll=0;
        end
        %if obj.nStd <= 1
            %obj.get_blk_order();
            %% XXX
        if obj.rnd.blkRndAll
            obj.get_blk_table_shuffle_all();
        else
            obj.get_blk_table_shuffle();
        end
    end
    function obj=get_blk_table_shuffle(obj)
        n=obj.nStd*obj.nBlk;
        indStd=zeros(n,1);
        stdXblkTable=zeros(n,1);
        cmpXblkTable=zeros(n,obj.nCmp);
        for i = 1:obj.nBlk
            ind=transpose((1:obj.nStd)+(obj.nStd*(i-1)));
            [indStd(ind),stdXblkTable(ind),cmpXblkTable(ind,:)]=obj.get_blk_table_shuffle_blk_ind(i);
        end
        obj.indStd=indStd;
        obj.indBlk=repelem(transpose(1:obj.nBlk),obj.nStd,1);
        if isfield(obj.methodVars,'stdXunqAll')
            obj.methodVars.stdXblkTable=stdXblkTable;
            obj.methodVars.cmpXblkTable=cmpXblkTable;
        end
    end
    function [indStd,stdXblkTable,cmpXblkTable]=get_blk_table_shuffle_blk_ind(obj,blknum)
        ind=transpose(1:obj.nStd);
        if obj.rnd.bBlk
            rng(obj.rnd.blk(blknum));
            ind=datasample(ind,obj.nStd,'Replace',false);
        end
        indStd=ind;
        if isfield(obj.methodVars,'stdXunqAll')
            stdXblkTable=obj.methodVars.stdXunqAll(ind);
            cmpXblkTable=obj.methodVars.cmpXunqAll(ind,:);
        else
            stdXblkTable=ind;
            cmpXblkTable=ind;
        end
    end
    function obj=get_blk_table_shuffle_all(obj)
        n=obj.nStd*obj.nBlk;
        ind=transpose(1:n);
        if obj.rnd.bBlk
            rng(obj.rnd(1));
            randInd=datasample(ind,n,'Replace',false);
        end
        indStd=repelem(1:obj.nStd,obj.nBlk,1);
        indBlk=repeelm(transpose(1:obj.nBlk),obj.nStd,1);
        obj.indStd=indStd(ind);
        obj.indBlk=cumCount(indStd);
        if isfield(obj.methodVars.stdXunqAll)
            obj.methodVars.stdXblkTable = repelem(obj.methodVars.stdXunqAll,obj.nBlk);
            obj.methodVars.cmpXblkTable = repelem(obj.methodVars.cmpXunqAll,obj.nBlk);

            obj.methodVars.stdXblkTable = obj.methodVars.stdXblkTable(ind);
            obj.methodVars.cmpXblkTable = obj.methodVars.cmpXblkTable(ind,:);
        end
    end
    function obj=get_indTrl(obj)
        rng(obj.rnd.trl,'twister');
        indTrl=1:obj.nTrl;
        obj.indTrl=zeros(obj.nTrlPerBlk,obj.nBlk,obj.nStd);
        if obj.rnd.bTrl
            for i = 1:obj.nStd
                tmp=randsample(indTrl,obj.nTrl);
                obj.indTrl(:,:,i)=reshape(tmp,obj.nTrlPerBlk,obj.nBlk);
            end
        end
        obj.indTrl=permute(obj.indTrl,[1 3 2]);
    end
    function obj=get_intrvl(obj)
        tpb=obj.nTrlPerBlk;

        INDS=distribute(1:obj.nStd,1:obj.nBlk);
        cmpIntrvl=zeros(tpb,obj.nStd,obj.nBlk);
        stdIntrvl=zeros(tpb,obj.nStd,obj.nBlk);
        for i=1:size(INDS,1)
            std=INDS(i,1);
            blk=INDS(i,2);
            [cmpIntrvl(:,std,blk),stdIntrvl(:,std,blk)]=obj.get_intrvl_ind(std,blk);
        end
        obj.methodVars.cmpIntrvl=cmpIntrvl;
        obj.methodVars.stdIntrvl=stdIntrvl;
    end
    function [cmpIntrvl,stdIntrvl]=get_intrvl_ind(obj,std,blk)
        tpb=obj.nTrlPerBlk;
        stdIntrvl = [ones( tpb/2,1); zeros(tpb/2,1)];
        cmpIntrvl = [zeros(tpb/2,1); ones( tpb/2,1)];
        if obj.rnd.bIntrvl
            rng(obj.rnd.intrvl(std,blk),'twister');
            cmpIntrvl= datasample(cmpIntrvl,size(cmpIntrvl,1),'Replace',false);
            stdIntrvl= (~cmpIntrvl);
        end
    end
    function obj=get_cmpXind(obj)
        if ~isfield(obj.rnd,'bCmp') || isempty(obj.rnd.bCmp)
            obj.rnd.bCmp=1;
        end
        cmpXind=repelem(transpose(1:obj.nCmp),obj.nTrlPerBlk/obj.nCmp,1);
        cmpXind=repmat(cmpXind,1,obj.nStd,obj.nBlk);
        if obj.rnd.bCmp
            rng(obj.rnd.cmp,'twister');
            INDS=distribute(1:obj.nStd,1:obj.nBlk);
            for i = 1:size(INDS,1)
                std=INDS(i,1);
                blk=INDS(i,2);
                cmpXind(:,std,blk)=datasample(cmpXind(:,std,blk),size(cmpXind,1),'Replace',false);
            end
        end
        obj.methodVars.cmpXind=cmpXind;
    end
    function obj=get_indTrl_from_data(obj)
       % TODO
    end
    function obj=apply_rnd(obj)
        obj.get_blkTable();
        obj.get_indTrl();
        obj.get_intrvl();
        obj.get_cmpXind();
    end
    function db=gen_DB_flds(obj)
       % TODO
        flds={'prjCode'...
             ,'Alias' ...
             ,'imgDTB'...
             ,'natORflt'...
             ,'imgDim'...
             ,'method'...
             ,'prjInd'...
             ,'subjs'...
             ,'nStd'...
             ,'nCmp'...
             ,'nBlk'...
             ,'nTrl'...
             ,'nTrlPerBlk'...
             ,'nTrlPerLvl'...
             ,'expHost'...
             ,'creationData'...
             ,'description'...
             ,'relPapers'...
             ,'language'...
             };
    end
    function EXP=select_detail(obj)
        EXP=exp_detail(obj);
    end
    function out=rm_subj_pruned(obj,out,rmsubjs)
        out=rm_dim_pruned(obj,out,'subjs',rmsubjs);
    end
    function out=rm_std_pruned(obj,out,rmstds)
        out=rm_dim_pruned(obj,out,'stds',rmstds);
    end
    function out=rm_dim_pruned(obj,out,fld,vals2rm)
        switch fld
            case 'stds'
                dim=2;
                list=obj.get_std_nums;
                vals2rm=obj.auto_std_num(vals2rm);
            case 'subjs'
                dim=3;
                list=obj.(fld);
        end
        if any(ismember(vals2rm,list))
            ind=sort(find(ismember(list,vals2rm)),'descend');
            out=structPrune(out,dim,ind);
        end
    end
end
end
