classdef Eobj_mod
    function []=change_subj_status(obj,subj,status)
        validStatus={'lock','train','test','pilot'};
        if ~ismember(status,validStatus)
            error(['Incorrect Status: ' status])
        end
        obj.subjStatus.(subj).status=status;
    end
    function obj=m_mod_nBlocks(obj)
        % TODO
    end
    function obj=m_mod_flag(obj)
        % TODO
    end
    function obj=m_mod_subj_status(obj)
        % TODO
    end
    function obj=m_mod_lock_exp(obj)
        % TODO
    end
    function obj=m_reset_subj(obj)
        % TODO
    end
    function obj=m_add_subj(obj)
        % TODO
    end
    function obj=m_rm_subj(obj)
        % TODO
    end
    function obj=m_mod_cmp_lvls(obj)
        % TODO
    end
    function obj=m_mod_mode_data(obj)
        % TODO
    end
end
