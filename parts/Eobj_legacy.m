classdef Eobj_legacy < handle
methods
    function obj=check_db_legacy(obj)
        chkDirAll(obj.dbDirLegacy,1);
    end
    function E=load_legacy(obj)
        obj.check_db_legacy();
        fnE=[obj.dir.META obj.name '.mat'];
        if exist(fnE,'file')
            load(fnE);
        else
            warning(['legacy db file dos not exist: ' fnE]);
        end
    end
    function E=save_legacy(obj)
        E=obj.convert_to_legacy();
        out=obj.check_data_dir('META');
        if out==0
            return
        end
        fnE=[obj.name '.mat'];
        if exist('fnE','file')
            save(E);
        else
            warning(['legacy db file dos not exist: ' fnE]);
        end
    end
    function obj =get_keyPairs(obj)

        obj.keyPairs={...
            'pass'              ,'pass'      ; ...
            'prjCode'           ,'prjCode'      ; ...
            'expID'             ,'expID'        ; ...
            'creationDate'      ,'creationDate' ; ...
            'relPapers'         ,'relPapers'    ; ...
            'trlPerBlk'         ,'nTrlPerBlk'   ; ...
            'trlPerExp'         ,'nTrl'         ; ...
            'blkPerExp'         ,'nBlk'         ; ...
            'nCmp'              ,'nCmp'         ; ...
            'ncmp'              ,'nCmp'         ; ...
            'bSKIPSYNCTEST'     ,'bSkipSyncTest'; ...
            'pAuthor'           ,'pAuthor'      ; ...
            'expHost'           ,'expHost'      ; ...
            'description'       ,'description'  ; ...
            'Subjs'             ,'subjs'        ; ...
            'SubjData'          ,'rawData'      ; ...
            'fullname'          ,'name'         ; ...
            'burge'             ,'burge'        ; ...
            'imgDTB'            ,'imgDTB'       ; ...
            'natORflt'          ,'natORflt'     ; ...
            'imgDim'            ,'imgDim'       ; ...
            'method'            ,'method'       ; ...
            'prjInd'            ,'prjInd'       ; ...
            'prjHier'           ,'prjHier'      ; ...
            'language'          ,'language'     ; ...
            'nTrlPerLvl'        ,'nTrlPerLvl'   ; ...
        };
    end
    function obj =convert_from_legacy(obj,E)
        obj.nStd=1;
        obj.get_keyPairs();
        c=containers.Map(obj.keyPairs(:,1),obj.keyPairs(:,2));
        flds=fieldnames(E);
        obj.notused=zeros(length(flds),1);
        twoIFC=struct();
        burgest=struct();
        for i = 1:length(flds)
            fld=flds{i};
            if c.isKey(fld)
                val=c(fld);
                obj.(val)=E.(fld);
                continue
            end

            switch fld
            case 'bRandomize'
                obj.rnd.bIntrvl=E.(fld);
                obj.rnd.bBlk=E.(fld);
                obj.rnd.bTrl=E.(fld);
            case 'rndSdEX'
                obj.rnd.master=E.(fld);
            case 'rndSdTO'
                obj.rnd.trl=E.(fld);
            case 'RndSdSTTable'
                obj.rnd.blk=E.(fld);
            case 'bsameRndSdST'
                obj.rnd.bSameBlk=E.(fld);
            case 'defFname'
                obj.fnames.DEF.EXP=['DEF' und E.(fld)];
            case 'StimfName'
                obj.fnames.IN.PRJ_TEST=E.(fld);
            case 'PsyStimfName'
                obj.fnames.IN.EXP=E.(fld);
            case 'experimentfName'
                obj.fnames.CODE.EXP=E.(fld);
            case 'serverDataDir'
                obj.dir.dataDir2=E.(fld);
            case 'localDataDir'
                burgest.localDataDir=E.(fld);
            case 'prjDir'
                obj.dir.prjDir=E.(fld);
            case 'stdXblkTable'
                twoIFC.stdXblkTable=E.(fld);
            case 'bUseFeedback'
                twoIFC.bUseFeedback=E.(fld);
            case 'cmpXblkTable'
                twoIFC.cmpXblTable=E.(fld);
            case 'stdXunqAll'
                twoIFC.stdXunqAll=E.(fld);
            case 'cmpXunqAll'
                twoIFC.cmpXunqAll=E.(fld);
            case 'cmpXSpacing'
                twoIFC.cmpXspacing=E.(fld);
            case 'minDist'
                twoIFC.minDist=E.(fld);
            case 'redo'
                obj.rawData.redo=E.(fld);
            otherwise
                obj.notused(i)=1;
            end
        end
        if isfield(E,'prjInd1') && isfield(E,'prjInd2')
            obj.prjInd=[str2double(E.prjInd1) str2double(E.prjInd2)];
        end
        if strcmp(obj.method,'2IFC')
            obj.methodVars=twoIFC;
            obj.nStd=length(twoIFC.stdXunqAll);
        end
        bflds=fieldnames(burgest);
        for i =1:length(bflds)
            fld=bflds{i};
            obj.burge.(fld)=burgest.(fld);
        end

        SUBJS=fieldnames(E.SubjData);
        for s = 1:length(SUBJS)
            subj=SUBJS{s};
            obj.subjStatus.(subj)=E.SubjData.(subj);
            if isfield(obj.rawData.(subj),'status')
                obj.rawData.(subj)=rmfield(obj.rawData.(subj),'status');
            end
            if isfield(obj.rawData.(subj),'lockmsg')
                obj.rawData.(subj)=rmfield(obj.rawData.(subj),'lockmsg');
            end

            METH=fieldnames(obj.rawData.(subj));
            for m = 1:length(METH)
                meth=METH{m};
                if strcmp(meth,'redo')
                    continue
                end
                if ~isfield(obj.rawData.(subj).(meth),'stdX')
                    continue
                end
                obj.rawData.(subj).(meth)=obj.rawData.(subj).(meth).stdX;
            end
        end

        obj.notused=flds(logical(obj.notused));
        obj.Eflds=flds;

        % XXX
        obj.E = E;
    end
    function E=convert_to_legacy(obj)
        E=struct();
        obj.get_keyPairs();
        c=containers.Map(obj.keyPairs(:,2),obj.keyPairs(:,1));
        flds=fieldnames(obj);
        obj.notused=zeros(length(flds),1);
        if strcmp(obj.method,'2IFC')
            Vars=obj.methodVars;
        else
            Vars=obj.legacy;
        end
        for i = 1:length(flds)
            fld=flds{i};
            if c.isKey(fld)
                val=c(fld);
                E.(val)=obj.(fld);
                continue
            end

            switch fld
            case 'rnd'
                E.bRandomize=obj.rnd.bIntrvl;
                E.bRandomize=obj.rnd.bBlk;
                E.bRandomize=obj.rnd.bTrl;
                E.rndSdEX= obj.rnd.master;
                E.RndSdSTTable= obj.rnd.blk;
                E.bsameRndSdST= obj.rnd.bSameBlk;
            case 'fnames'
                E.defFname= obj.fnames.DEF.EXP;
                E.StimfName= obj.fnames.IN.PRJ_TEST;
                E.PsyStimfName= obj.fnames.IN.PRJ_TEST;
                E.experimentfName= obj.fnames.CODE.EXP;
            case 'dir'
                E.dataDir= obj.dir.RAW;
                E.serverDataDir= obj.dir.dataDir2;
                E.localDataDir= obj.burge.localDataDir;
                E.prjDir= obj.dir.prjDir;
            case 'methodVars'
                E.stdXblkTable= Vars.stdXblkTable;
                E.bUseFeedback= Vars.bUseFeedback;
                E.cmpXblkTable= Vars.cmpXblTable;
                E.stdXunqAll= Vars.stdXunqAll;
                E.cmpXunqAll= Vars.cmpXunqAll;
                E.cmpXSpacing=Vars.cmpXspacing;
                E.minDist=Vars.minDist;
            case 'subjStatus'
                SUBJS=fieldnames(obj.subjStatus);
                for j = 1:length(SUBJS)
                    subj=SUBJS{j};
                    METH=fieldnames(obj.subjStatus.(subj));
                    for k = 1:length(METH)
                        fld=METH{k};
                        E.SubjData.(subj).(fld)=obj.subjStatus.(subj).(fld);
                    end
                end
            otherwise
                obj.notused(i)=1;
            end

        end

        E.redo=E.SubjData.redo;
        E.SubjData=rmfield(E.SubjData,'redo');

        SUBJ=fieldnames(obj.rawData);
        for i = 1:length(SUBJ)
            subj=SUBJ{j};
            METH=fieldnames(obj.rawData.(subj));
            for j = 1:length(METH)
                meth=METH{j};
                if strcmp(meth,'redo')
                    continue
                end
                E.SubjData.(subj).(meth)=[];
                E.SubjData.(subj).(meth).stdX=obj.rawData.(subj).(meth);
            end
        end

        E.burge=obj.burge;

        E.prjHier=obj.prjHier;
        E.expID=obj.expID;
        E.bSKIPSYNCTEST=obj.bSkipSyncTest;
        E.prjInd1=E.prjInd(1);
        E.prjInd2=E.prjInd(2);

        obj.notused=flds(logical(obj.notused));
        fldsNew=fieldnames(E);
        obj.notused;
        setdiff(obj.Eflds,fldsNew)

    end
end
end
