classdef EobjBlk < handle
properties
    E

    alias
    hashes

    dImap
    dPsy
    dBlk
    dTrk

    % trk
    hashAlias
    prjCode
    imgDTB
    natORflt
    imgDim
    method
    prjInd
    pass
    subjs
    pAuthor
    language
    expHost
    Xname
    Xunits
    description

    rnd_bBlk
    rnd_bTrl
    rnd_bIntrvl


    % psy
    psyLib
    expType

    % blk
    Blk
    Sblk
end
methods
    function obj=EobjBlk(alias)
        obj.alias=sed('s',alias,'^D_[a-zA-z]+_','');

        obj.get_def_files();
        obj.load_defs();

        obj.get_hashes();
        obj.load_blk();

        obj.populate_Eobj();
        obj.init_Eobj();
        obj.save_Eobj();

    end
    function get_def_files(obj)
        obj.dImap=which(['D_imap_' obj.alias]);
        obj.dBlk =which(['D_blk_' obj.alias]);
        obj.dTrk =which(['D_trk_' obj.alias]);
        obj.dPsy =which(['D_psy_' obj.alias]);
    end
    function get_hashes(obj)
        obj.hashes=imapCommon.alias2hashes(obj.hashAlias,obj.imgDTB);
    end
    function obj=populate_Eobj(obj)
        base=dbDirs('base');
        loc=dbDirs('loc');

        obj.E=Eobj();
        obj.E.alias=obj.alias;
        obj.E.prjCode=obj.prjCode;
        obj.E.imgDTB=makeUpperCase(obj.imgDTB);
        obj.E.natORflt=makeUpperCase(obj.natORflt);
        obj.E.imgDim=obj.imgDim;
        obj.E.method=obj.method;
        obj.E.prjInd=obj.prjInd;
        obj.E.pass=obj.pass;
        obj.E.subjs=obj.subjs;
        obj.E.pAuthor=obj.pAuthor;
        obj.E.language=obj.language;
        obj.E.expHost=obj.expHost;
        obj.E.Xname=obj.Xname;
        obj.E.Xunits=obj.Xunits;
        obj.E.description=obj.description;


        obj.E.nStd=obj.Sblk.nStd;
        obj.E.nCmp=obj.Sblk.nCmp;
        obj.E.nBlk=obj.Sblk.nBlk;
        obj.E.nTrl=obj.Sblk.nTrl;
        obj.E.nTrlPerBlk=obj.Sblk.nTrlPerBlk;
        obj.E.nTrlPerLvl=obj.Sblk.nTrlPerLvl;
        obj.E.methodVars=obj.Sblk.methodVars;


        % RND
        obj.E.rnd=struct();
        obj.E.rnd.bBlk=obj.rnd_bBlk;
        obj.E.rnd.bTrl=obj.rnd_bTrl;
        obj.E.rnd.bIntrvl=obj.rnd_bIntrvl;

        obj.E.fnames.DEF=struct();
        obj.E.fnames.IN=struct();
        obj.E.fnames.MAP=struct();
        obj.E.fnames.CODE=struct();


        % DEF

        obj.E.fnames.DEF.EXP= strrep(obj.dPsy,loc,'');
        obj.E.fnames.DEF.IMAP=strrep(obj.dImap,loc,'');
        obj.E.fnames.DEF.TRK=strrep(obj.dTrk,loc,'');
        obj.E.fnames.DEF.SEL=strrep(obj.dBlk,loc,'');
        obj.E.creationDate=date();

        obj.E.bBlk=1;
        obj.E.bPtchs=1;
        % IN
        flds={...
             'DMP_PILOT' ...
            ;'DMP_TRAIN' ...
            ;'DMP_TEST' ...
            ;'SEL_PILOT' ...
            ;'SEL_TRAIN' ...
            ;'SEL_TEST' ...
            ;'EXP_PILOT' ...
            ;'EXP_TRAIN' ...
            ;'EXP_TEST' ...
        };
        for  i = 1:length(flds)
            fld=flds{i};
            if startsWith(fld,'DMP')
                % PATCHS

                obj.E.fnames.IN.(fld)=strrep(ptchs.get_dmp_fname_p(obj.imgDTB,obj.hashes.pch),base,'');
            elseif startsWith(fld,'SEL')
                dire=Blk.get_dir(obj.alias);
                obj.E.fnames.IN.(fld)=strrep([dire 'blk'],base,'');
            else
                obj.E.fnames.IN.(fld)='';
            end
        end

        % CODE
        obj.E.fnames.CODE.EXP=obj.psyLib;

        % MAP
        flds=fieldnames(obj.hashes);
        for i = 1:length(flds)
            fld=flds{i};
            FLD=makeUpperCase(fld);
            obj.E.fnames.MAP.(FLD)=obj.hashes.(fld);
        end
    end
    function obj=init_Eobj(obj);
        obj.E.get_name();
        obj.E.init();
    end
    function save_Eobj(obj);
        obj.E.save();
    end
%% LOAD
    function obj=load_blk(obj)
        obj.Blk=Blk.load(obj.alias);
        obj.Sblk=obj.Blk.ret_blk_struct();
    end
    function load_defs(obj)
        obj.load_trk_def();
        obj.load_imap_def();
        obj.load_psy_def();
        obj.load_blk_def();

    end
%% LOAD DEF
    function load_imap_def(obj)
        fname=which(['D_imap_' obj.alias]);
        if isempty(fname)
            return
        end
        run(fname);

        P=EobjBlk.get_imap_parse_opts();

        for i = 1:length(P)
            fld=P{i};
            if exist(fld,'var')
                obj.(fld)=eval([ fld ';']);
            end
        end
    end
    function load_psy_def(obj)
        fname=which(['D_psy_' obj.alias]);
        if isempty(fname)
            return
        end
        run(fname);

        P=EobjBlk.get_psy_parse_opts();

        for i = 1:length(P)
            fld=P{i};
            if exist(fld,'var') && strcmp(fld,'lib')
                obj.psyLib=lib;
            elseif exist(fld,'var')
                obj.(fld)=eval([ fld ';']);
            end
        end
    end
    function load_blk_def(obj)
        fname=which(['D_blk_' obj.alias]);
        if isempty(fname)
            return
        end
        run(fname);

        P=EobjBlk.get_blk_parse_opts();

        for i = 1:length(P)
            fld=P{i};
            if exist(fld,'var')
                obj.(fld)=eval([ fld ';']);
            end
        end
    end
    function load_trk_def(obj)
        fname=which(['D_trk_' obj.alias]);
        if isempty(fname)
            return
        end
        run(fname);

        P=EobjBlk.get_trk_parse_opts();

        for i = 1:length(P)
            fld=P{i};
            if exist(fld,'var')
                obj.(fld)=eval([ fld ';']);
            end
        end
    end

end
methods(Static=true)
    function P=get_psy_parse_opts()
        P={ ...
           'lib',...
           ;'expType'....
          };
    end
    function P=get_imap_parse_opts()
        P={...
          };
    end
    function P=get_blk_parse_opts()
        P={...
          };
    end
    function P=get_trk_parse_opts()

        P={...
            'hashAlias' ...
            ;'prjCode'  ...
            ;'imgDTB'   ...
            ;'natORflt' ...
            ;'imgDim'   ...
            ;'method'   ...
            ;'prjInd'   ...
            ;'pass'     ...
            ;'subjs'    ...
            ;'pAuthor'  ...
            ;'language' ...
            ;'expHost'  ...
            ;'Xname'    ...
            ;'Xunits'   ...
            ;'description' ...
            ;'rnd_bBlk' ...
            ;'rnd_bTrl' ...
            ;'rnd_bIntrvl' ...
        };
    end
end
end
