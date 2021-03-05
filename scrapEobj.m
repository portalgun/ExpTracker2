    function C=cleanup_raw(obj,subj,mode,std,gdSubj,gdStd,gdPss,gdBlk,fromGd,indFld,newInd)
        SS=struct();
        if exist('fromGd','var') && ~isempty(fromGd)
            SS.fromGd=frromGd;
        else
            newInd=[];
        end
        if exist('indFld','var') && ~isempty(indFld)
            SS.indFld=indFld;
        end
        if exist('newInd','var') && ~isempty(newInd)
            SS.newInd=newInd;
        else
            newInd=[];
        end

        chks=zeros(4,1);
        if exist('gdSubj','var') && ~isempty(gdSubj)
            SS.gdSubj=gdSubj
            chks(1)=1;
        end
        if exist('gdStd','var') && ~isempty(gdStd)
            SS.gdStd=gdStd
            chks(2)=1;
        end
        if exist('gdPss','var') && ~isempty(gdPss)
            SS.gdPss=SS.gdPss;
            chks(3)=1;
        end
        if exist('gdBlk','var') && ~isempty(gdBlk)
            SS.gdBlk=gdBlk;
            chks(4)=1;
        end
        if all(chks)
            [SS.SallGd,INDSgd]=obj.load_raw_test_blocks_exemplar(gdSubj,gdStd,gdPss,gdBlk);
        %    SS.gdBlocksAll=[]; % XXX
             SS.gdBlkInd=vertcat(INDSgd{:,4})
        elseif any(chks)
            error('Must specify all Gd params');
        end



        [SS.Sall,SS.fnames,INDS]=obj.load_raw_test_blocks(subj,std);
        SS.blocks=obj.nBlk;
        SS.blkInd=vertcat(INDS{:,4});

        %SS.blocksAll %  all indTrl % XXX

        C=blockCleanup(SS);
        C.main();
        S=C.S

        obj.save_blocks(subj,method,std,S);

    end
