function obj= runJDB()
    subj='JDB';

    rule='rand_all'
    mode=1;
    std={'-11.25', '-7.5','-3.75'};

    EAobj.run_loop_all(subj,Aliases,rule,mode,std);
end
