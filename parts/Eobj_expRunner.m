classdef Eobj_expRunner < handle
properties
    lastPrj
end
methods
    function [I,exitflag]=get_next_run_parms(obj,subj, mode,std,blk)
        I=struct();

        I.subj=obj.auto_subj(subj);
        if ~exist('mode','var') || isempty(mode)
            I.mode=obj.subjStatus.(subj).status;
        else
            I.mode=obj.auto_mode(mode);
        end
        if ~exist('std','var')
            std=[]
        end
        if ~exist('blk','var')
            blk=[]
        end

        [I.ind,I.std,I.blk,I.r1,I.r2,exitflag]=obj.get_next_std_block(I.subj,I.mode,std,blk);
        if exitflag
            return
        end

        I.fname=obj.gen_fname_block(I.subj,I.mode,I.std,I.blk,I.r1,I.r2);
        I.name=obj.gen_name_block(I.subj,I.mode,I.std,I.blk,I.r1,I.r2);

    end
    function [ind,std,blk,r1,r2,exitflag]=get_next_std_block(obj,subj,mode,std,blk)
        exitflag=0;
        % also check if complete
        % TODO handle flags
        %f=obj.get_block_flag(subj,mode,std,blk);
        r1=[];
        r2=[];
        cind=obj.get_std_block_ind_completion(subj,mode,std,blk);
        if all(cind)
            exitflag=1;
            return
        end
        std=obj.auto_std_ind(std);

        k=ismember(obj.indStd,std);
        if exist('std','var') && ~isempty(std)
            cind=(cind==0 & k);
        end
        ind=find(cind,1,'first');
        if isempty(ind)
            exitflag=1;
            return
        end
        std=obj.indStd(ind);

        if ~exist('blk') || isempty(blk)
            blk=obj.indBlk(ind);
        end

        [r1,r2]=obj.get_block_redo(subj,mode,std,blk);
    end
    function [ind,INDS]=get_std_block_ind_completion(obj,subj,mode,std,blk)
        std=obj.auto_std_ind(std);
        if ~exist('std','var') || isempty(std)
            std=1:obj.nStd;
        end
        if ~exist('blk','var') || isempty(blk)
            blk=1:obj.nBlk;
        end
        INDS=distribute(std,blk);
        ind=zeros(size(INDS,1),1);
        for i = 1:size(INDS,1)
            std=INDS(i,1);
            blk=INDS(i,2);
            ind(i)=obj.get_block_val(subj,mode,std,blk);
        end
    end
%%
    function [obj,exitflag]= run_loop(obj,subj, std)
        if ~exist('std','var')
            std=[];
        end
        exitflag=0;
        while ~exitflag
            [obj,exitflag]=obj.run(subj,[],std,[], 1);
        end
    end
    function [obj,exitflag] = run(obj,subj, mode,std,blk, bLoop, bTest)
        exitflag=0;
        if ~exist('bTest','var') || isempty(bTest)
            bTest=0;
        end
        if ~exist('bLoop','var') || isempty(bLoop)
            bLoop=0;
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

        if strcmp(obj.subjStatus.(subj).status,'lock')
            disp(['Subject ' subj  ' is locked']);
            exitflag=1;
            return
        end
        if strcmp(obj.expHost,hostname);
            disp(['Experiment can only run on ' obj.expHots '. Run in test mode instead'])
            exitflag=1;
            return
        end

        [I,exitflag]=obj.get_next_run_parms(subj,mode,std,blk);
        if exitflag
            disp(['Subject complete'])
            return
        end

        % TODO check data parity, no overwriting
        S=obj.load_exp(I.mode,I.std,I.blk);
        Def=obj.load_def('EXP');

        obj.save_prj_path();
        obj.apply_prj_path();

        obj.print_exp(I,S);
        if ~bLoop || bTest==2
            obj.wait_exp();
        end
        [exp,ME,statusCode]=obj.run_exp(S,Def,bTest);

        % ASSIGN IN EXP
        obj.assignin_exp(exp,statusCode,I);
        obj.assignin_pre();

        if ~bTest
            [obj,bSuccess]=obj.save_handler(exp,statusCode,I,0);
        else
            bSuccess=1;
        end

        if bSuccess && bTest < 2
            obj.assignin_post(); 
        end
        obj.handle_exp_error(ME);
        obj.restore_prj_path();
        if (bSuccess && bTest == 0) || bTest==1
            Eobj.plot_exp_results(exp);
        end
    end
    function [obj,bSuccess]=save_handler(obj,exp,statusCode,I, bRetry)
        if ~exist('bRetry','var') || isempty(bRetry)
            bRetry=0;
        end

        bSuccess=0;
        if statusCode~=1
            bSucces=0;
            return
        end

        try
            obj.save_exp_data(exp,statusCode,I);
            obj.update_completion(statusCode,I);
            bSuccess=1;
        catch ME
            if statusCode==1    
                disp(['FIX SAVING ERROR THEN RUN: ' newline '    preEobj.retry_save(exp,statusCode,I)']);
                rethrow(ME);
            else
                disp(ME);
            end 
        end
    end
    function [obj,bSuccess]=retry_save(obj,exp,statusCode,I)
        [obj,bSuccess]=obj.save_handler(exp,statusCode,I,1);
        if bSuccess
            obj.assignin_post();
        end
        obj.restore_prj_path();
        if bSuccess
            Eobj.plot_exp_results(exp);
        end
    end
    function obj=assignin_pre(obj)
        assignin('base','preEobj',obj);
        disp(['Variable ''preEobj'' moved to base workspace']);
    end
    function obj=assignin_exp(obj,exp,statusCode,I)
        assignin('base','exp',exp);
        disp(['Data moved to base workspace as ''exp'''])

        assignin('base','I',I);
        assignin('base','statusCode',statusCode);
        disp(['Variables ''I'', and ''statusCode'' moved to base workspace ']);
    end
    function obj=assignin_post(obj)
        assignin('base','postEobj',obj);
        disp(['Variable ''postEobj'' moved to base workspace']);
    end
    function [exp,ME,statusCode] = run_exp(obj,S,Def,bTest)
        % 2=run
        % 1=complete
        % 0=
        ename=obj.fnames.CODE.EXP;
        if bTest==2
            exp=[];
            ME=[];
            statusCode=1;
        elseif strcmp(ename,'psycho_simple')
            [exp,ME,statusCode]=obj.run_psycho_simple(S,Def,bTest);
        end
    end
    function [exp,ME,statusCode]=run_psycho_simple(obj,S,Def,bTest)
    % 1  complete
    % 0  run
    % -1 error
    % -2 exited
        ME=[];
        exp=[];
        statusCode=0;
        exp=psycho_simple(S,Def,1,bTest);
        exp.run();
        statusCode=exp.returncode;
        ME=exp.ME;
    end
    function [obj,exp]=run_exp_tst(obj,S,Def)
        ename=obj.fnames.CODE.EXP;
        str=[ename  '(S,Def)'];
        disp(str);
        if strcmp(ename,'psycho_simple')
            exp=feval(ename,S,Def);
        end
    end
    function obj=save_exp_data(obj,exp,statusCode,I)
        if statusCode~=1;
            return
        end
        out=obj.check_raw_save(I.fname,I.subj,I.mode,I.std,I.blk);
        if out == 1
            disp('Saving...');
            save(I.fname,'exp');
            disp(['Saved data:' newline '    ' I.fname]);
        end
    end
    function obj=update_completion(obj,statusCode,I)
        % TODO other codes
        subj=obj.auto_subj(I.subj);
        mode=obj.auto_mode(I.mode);
        std=obj.auto_std_str(I.std);
        stdf=obj.auto_std_fld(I.std);
        blk=obj.auto_blk_num(I.blk);

        obj.rawData.(subj).(mode).(stdf).data{blk}=statusCode;
        if statusCode==1
            obj.rawData.(subj).(mode).(stdf).data{blk}=I.name;
            if ~isfield(obj.rawData.(subj).(mode).(stdf),'date')
                 obj.rawData.(subj).(mode).(stdf).date=cell(size(obj.rawData.(subj).(mode).(stdf).data));
            end
            obj.rawData.(subj).(mode).(stdf).date{blk}=date;
            obj.rawData.(subj).(mode).(stdf).blk(blk)=1;
        end
        obj.save();
        disp(['Database updated:' newline '    obj.' subj '.' mode '.' stdf '.' 'data{' num2str(blk) '}' I.name ])
    end
    function obj=handle_exp_error(obj,ME)
        if ~isempty(ME)
            try
                obj.restore_prj_path;
            end
            rethrow(ME);
        end
    end
%% PATH
    function obj=save_prj_path(obj)
        obj.lastPrj=pxCur;
    end
    function obj=restore_prj_path(obj)
        if ~isempty(obj.lastPrj)
            px(obj.lastPrj);
        end
    end
    function obj=apply_prj_path(obj)
        if ispc
            fs='\\';
        else
            fs=filesep;
        end
        dir=sed('s',obj.dir.prjDir,[fs '$'],'');
        px(dir);
    end
%/
    function obj=wait_exp(obj)
        freq = 0.73; 
        sound(sin(freq.*[0:1439]).*cosWindowFlattop([1 1440],720,720,0));
        input('');
    end
    function obj=print_exp(obj,I,S)
    % PRINT DATA TO SCREEN
        ename=obj.fnames.CODE.EXP;
        sfname=obj.get_fname_exp(I.mode,I.std,I.blk);
        dfname=[obj.dir.DEF obj.fnames.DEF.EXP];

        subj=obj.auto_subj(I.subj);
        mode=obj.auto_mode(I.mode);
        std=obj.auto_std_str(I.std);
        stdf=obj.auto_std_fld(I.std);
        blk=obj.auto_blk_num(I.blk);
        nTrial=num2str(S.trlPerRun); % TODO

        sstd=num2str(unique(roundDec(S.stdX,1e-12)));
        scmp=num2str(length(S.cmpXindUnq));

        spc='    ';
        disp([ ...
               spc 'prj name   : ' obj.name newline ...
               spc 'db name    : ' I.name  newline ...
               spc 'def  fname : ' dfname newline ...
               spc 'save fname : ' I.fname  newline ...
               spc 'stim fname : ' sfname  newline ...
               newline ...
               spc 'subj       : ' subj newline ... 
               spc 'mode       : ' I.mode newline...
               spc 'db std     : ' std newline ...
               spc 'stim std   : ' sstd newline ...
               spc 'nCmp       : ' scmp newline ...
               spc 'block      : ' num2str(blk) newline ...
               newline ...
               spc 'host:      : ' hostname newline ...
               spc 'exp handler: ' ename newline ...
               spc 'nTrial:    : ' nTrial newline ...
               'PRESS RETURN' newline ...
        ]);
             
              
    end
%% LATER
    function I=get_IFC_info(obj,I)
        % TODO
        I.indTrl=obj.indTrl(:,std,blk);
        I.trlPerBlk=obj.trlPerBlk;
        [~,I.name,~]=filePartsSane(I.fname);

        if isfield(obj.methodVars,'cmpXunq')
            I=obj.get_IFC_info(I);
        end
        if ~exist('I','var') || isempty(I)
            I=struct();
        end
        I.cmpXblk  =obj.methodVars.cmpXblkTable(ind,:);
        I.stdXblk  =obj.methodVars.stdXblkTable(ind,:);
        I.cmpIntrvl=obj.methodVars.cmpIntrvl(:,std,blk);
        I.stdIntrvl=obj.methodVars.stdIntrvl(:,std,blk);
        I.cmpXind  =obj.methodVars.cmpXind(:,std,blk);
        I.stdXunq  =obj.methodVars.stdXunq;
        I.cmpXunq  =obj.methodVars.cmpXunq;
    end
    function obj=print_exp_complete_status(obj)
        % TODO
    end
    function obj=verify_exp(obj)
        % TODO
        % make sure out data looks like in data
    end
    function part=get_partial_exp(obj,stim,mode,std,blk)
        part=[];
        % TODO
    end
    function [] = px(obj)
        px(obj.dir.prjDir);
    end
    % XXX press key to open def
    % XXX press key to open exp
end
end
