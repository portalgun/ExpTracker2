classdef Eobj_status < handle
methods
    function []=check_meta_lock(obj)
        % TODO
    end
%% COMPLETION
    function [cnt,blkTotal,stdTotal,subjTotal]=get_test_completion(obj)
        [cnt,blkTotal,stdTotal,subjTotal] =obj.get_completion('test');
    end
    function [cnt,blkTotal,stdTotal,subjTotal] =get_train_completion(obj)
        [cnt,blkTotal,stdTotal,subjTotal]=obj.get_completion('train');
    end
    function [cnt,blkTotal,stdTotal,subjTotal]=get_pilot_completion(obj)
        [cnt,blkTotal,stdTotal,subjTotal]=obj.get_completion('pilot');
    end
    function [cnt,blkTotal,stdTotal,subjTotal]=get_completion(obj,mode)
    % all
        [cnt,blkTotal,stdTotal]=get_completion_ind(obj,mode);
        subjTotal=numel(cnt);
        cnt=sum(cnt);
        stdTotal=sum(stdTotal);
        blkTotal=sum(blkTotal);
    end
    function [cnt,blkTotal,stdTotal]=get_completion_ind(obj,mode)
    % 3 subj
        nSubjs=length(obj.subjs);
        stdTotal=zeros(nSubjs,1);
        blkTotal=zeros(nSubjs,1);
        cnt=zeros(nSubjs,1);
        for i = 1:nSubjs
            subj=obj.subjs{i};
            [cnt(i),blkTotal(i),stdTotal(i)]=obj.get_subj_completion(subj,mode);
        end
    end
    function [cnt,blkTotal,stdTotal]=get_subj_completion(obj,subj,mode)
    % 1 subj
        [cnt,blkTotal]=get_subj_completion_ind(obj,subj,mode);
        stdTotal=numel(cnt);
        cnt=sum(cnt);
        blkTotal=sum(blkTotal);
    end

    function [cnt,blkTotal]=get_subj_completion_ind(obj,subj,mode)
    % 5 stds
        stds=obj.get_std_flds();
        cnt=zeros(length(stds),1);
        blkTotal=zeros(length(stds),1);
        for i = 1:length(stds)
            std=stds{i};
            [cnt(i),blkTotal(i)]=obj.get_subj_std_completion(subj,mode,std);
        end
    end

    function [cnt,total]=get_subj_std_completion(obj,subj,mode,std)
    %  1 std completion
        cnt=obj.get_subj_std_completion_ind(subj,mode,std);
        total=numel(cnt);
        cnt=sum(cnt);
    end
    function [cnt]=get_subj_std_completion_ind(obj,subj,mode,std)
    % 5 blocks
        cnt=zeros(obj.nBlk,1);
        for blk = 1:obj.nBlk
            cnt(blk)=get_block_val(obj,subj,mode,std,blk);
        end
    end
%% SUMMARY
    function []=summary(obj)
        % TODO
    end

end
end
