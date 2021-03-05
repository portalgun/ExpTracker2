classdef Eobj_psy < handle
methods
    function plot_block_prog_all(obj,mode,subjs)
        figure(nFn)
        if ~exist('mode','var') || isempty(mode)
            mode='test';
        end
        if ~exist('subjs','var') || isempty(subjs)
            subjs=1:obj.nSubj;
        end
        N=numel(subjs);

        name=strrep(obj.name,'_',' ');
        sgtitle(name);
        for s = subjs
        for i = 1:obj.nStd
            subPlot([N,obj.nStd],s,i);
            obj.plot_block_prog(s,mode,i);
            if i == 1
                rtitl=obj.auto_subj(s);
            else
                rtitl=[];
            end
            if s ==N
                ctitl=obj.auto_std_str(i);
            else
                ctitl=[];
            end
            formatFigure(ctitl,rtitl);
            if s ~= N
                xticks('');
            end
            if i ~= 1
                yticks('');
            end

            xlim([1,obj.nStd]);
            ylim([0 3]);
        end
        end
    end
    function plot_block_prog(obj,subj,mode,std)
        S=cell(obj.nBlk,1);
        for blk = 1:obj.nBlk
            try
                raw=obj.load_block(subj,mode,std,blk);
            catch
                disp(['Cannot load ' obj.auto_subj(subj) obj.auto_mode(mode) obj.auto_std_str(std) num2str(blk)])
                continue
            end
            S{blk}=Eobj.get_2IFC_fit(raw);
        end

        S=cellStructMerge(S); hold on
        %plot(S.sFit,'k','LineWidth',2)
        interv(1:5,S.Tfit-S.TCI,S.Tfit+S.TCI);
        plot(S.Tfit,'k','LineWidth',2)
    end
end
methods(Static=true)
    function S=get_2IFC_fit(exp)
        if isfield(exp,'S')
            exp=exp.S;
        end
        if isfield(exp,'RSP')
            R=exp.RSP.responses;
            if size(R,2) > 1
                R(:, ~any(R,1))=[];
            end
        end
        if isfield(exp,'RcmpChosen')
            RcmpChosen=exp.RcmpChosen;
        elseif isfield(exp,'cmpIntrvl')
            cmpIntrvl=exp.cmpIntrvl;
            RcmpChosen=R==cmpIntrvl;
        end

        stdX=exp.stdX;
        cmpX=exp.cmpX;

        S=struct();
        [S.mFit,S.sFit,S.bFit,S.Tfit,S.PCdta,S.PCfit,S.negLL]=psyfitgengauss(stdX,cmpX,RcmpChosen,[],[],1,1.36,2,0);
        [~,~,~,S.TCI]=psyfitgengaussBootstrap(stdX,cmpX,RcmpChosen,[],[],1,1.36,2,100,95);
        S.TCI=S.TCI(1);

    end
    function []=plot_exp_results(exp)
        if strcmp(exp.expType,'2IFC')
            Eobj.plot_2IFC(exp);
        end
    end
    function []=plot_2IFC(exp)
        if isfield(exp,'S')
            exp=exp.S;
        end
        if isfield(exp,'RSP')
            r=exp.RSP.responses;
            if size(R,2) > 1
                R(:, ~any(R,1))=[];
            end
        end
        if isfield(exp,'RcmpChosen')
            RcmpChosen=exp.RcmpChosen;
        elseif isfield(exp,'cmpIntrvl')
            cmpIntrvl=exp.cmpIntrvl;
            RcmpChosen=R==cmpIntrvl;
        end

        stdX=exp.stdX;
        cmpX=exp.cmpX;

        try
            psyfitgengauss(stdX,cmpX,RcmpChosen,[],[],1,1.36,2,1);

        catch ME
            display(ME.message);
        end
    end
end
end
