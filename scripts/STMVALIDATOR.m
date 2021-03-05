obj1=Eobj.load('DSP1D_1');
obj2=Eobj.load('DSP2D_2');

INDS=distribute(1:5,1:5);
for i =1:length(INDS)
    s=INDS(i,1);
    b=INDS(i,2);

    S1=obj1.load_exp(1,s,b);
    S2=obj2.load_exp(1,s,b);


    for j = 1:S1.trlPerRun
        Ls1=S1.stdIphtL(:,:,j);
        Rs1=S1.stdIphtR(:,:,j);

        Lc1=S1.cmpIphtL(:,:,j);
        Rc1=S1.cmpIphtR(:,:,j);

        Ls2=S2.stdIphtL(:,:,j);
        Rs2=S2.stdIphtR(:,:,j);
        Lc2=S2.cmpIphtL(:,:,j);
        Rc2=S2.cmpIphtR(:,:,j);

        img=[Ls2 Rs2 Lc2 Rc2; Ls1 Rs1 Lc1 Rc1];

        colormap gray;
        imagesc(img);
        waitforbuttonpress
    end
end
