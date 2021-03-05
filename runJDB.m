function obj= runJDB()
    subj='JDB';

    rule='rand_all';
    mode=1;
    std={'-11.25', '-7.5','-3.75'};
    Aliases={'DSP1D_1n','DSP2D_1n','DSP1D_2n','DSP2D_2n'};

    EAobj.run_loop_all(subj,Aliases,rule,mode,std);
end
