
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a test code to implement Zhu's (2010) method of FMVCM. Note here we use FA and MD.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Analysis Flu Study   %%% 

 names=cell(24, 1);
 names{1}='Genu-';
 names{2}='LeftThal-Ant-'; 
 names{3}='RightPons-Ant-';  
 names{4}='LeftCingulum-'; 
 names{5}='LeftThal-Occ-'; 	
 names{6}='RightPons-Occ-'; 
 names{7}='LeftExtreme-'; 
 names{8}='LeftThal-Par-'; 	
 names{9}='RightPons-Par-';
 names{10}='LeftFornix-';
 names{11}='LeftUncinate-';	
 names{12}='RightThal-Ant-'; 
 names{13}='LeftILF-';		
 names{14}='RightCingulum-';
 names{15}='RightThal-Occ-'; 
 names{16}='LeftPons-Ant-';
 names{17}='RightExtreme-';
 names{18}='RightThal-Par-';
 names{19}='LeftPons-Occ-';
 names{20}='RightFornix-';	
 names{21}='RightUncinate-';
 names{22}='LeftPons-Par-';
 names{23}='RightILF-'; 
 names{24}='Splenium-'; 
 
 Flunames=cell(3, 1);
 Flunames{1}='Fluearly';
 Flunames{2}='Flulate'; 
 Flunames{3}='FluEL';  
 
 dii=1;
 djj=1;
 

addpath 'F:\Study\PDF\UNC\Fibertract\website\code';

dataname=sprintf('F:/Study/PDF/UNC/Fibertract/Hongtu/fibertracts/%s.txt', Flunames{djj});
designdata=load(dataname);

dataname=sprintf('F:/Study/PDF/UNC/Fibertract/Hongtu/fibertracts/tracts/tract%sMD.txt', names{dii});
tractdata=load(dataname);

        %%%  Analysis of lambda1, lambda2, lambda3 together  %%%
        diffusionFiles=cell(2,1);
        dataname=sprintf('F:/Study/PDF/UNC/Fibertract/Hongtu/fibertracts/FA/%s%sFA.txt', Flunames{djj}, names{dii});
        diffusionFiles{1}=load(dataname);
        dataname=sprintf('F:/Study/PDF/UNC/Fibertract/Hongtu/fibertracts/MD/%s%sMD.txt', Flunames{djj}, names{dii});
        diffusionFiles{2}=load(dataname);  
        nofeatures=size(diffusionFiles, 1);
        [NoSetup, arclength, Xdesign, Ydesign]=MVCM_read(tractdata, designdata,diffusionFiles, nofeatures);
        
        n=NoSetup(1);
        L0=NoSetup(2);
        p=NoSetup(3);
        m=NoSetup(4);
        
    %% plot tract
    Mnames=cell(2,1);
    Mnames{1}='FA';
    Mnames{2}='MD';
    plotpath='F:\Study\PDF\UNC\Fibertract\website\figure';
    for mii=1:m
    figure()
    set(gcf,'PaperUnits','centimeters');
    xSize=12;
    ySize=8;
    xLeft=(21-xSize)/2; 
    yTop= (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
    set(gcf,'Position',[24 36 xSize*50 ySize*50]);
        

    for nii=1:n
        plot(arclength,Ydesign(nii,:,mii),'r-','LineWidth', 2);
        hold on;
    end
    hold off;

    xlabel('arc-length');
    ylabel(Mnames{mii});
    set(gca,'XLim',[0 size(arclength,1)-1]);
    
    figurename=sprintf('%s/%s_tract.bmp',plotpath,Mnames{mii});
    saveas(gcf,figurename);
    close()
    end
        %% 2.1 estimate betas using local polynomial kernel smoothing
       
        [mh] = MVCM_lpks_wob(NoSetup,arclength,Xdesign,Ydesign);    
     
        [efitBetas,efitBetas1,InvSigmats,efitYdesign]=MVCM_lpks_wb1(NoSetup,arclength,Xdesign,Ydesign,mh);
    
    %% plot betas
    
        figure()
    set(gcf,'PaperUnits','centimeters');
    xSize=12;
    ySize=8;
    xLeft=(21-xSize)/2; 
    yTop= (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
    set(gcf,'Position',[24 36 xSize*50 ySize*50]);
    
    for mii=1:m
    figure()
    set(gcf,'PaperUnits','centimeters');
    xSize=12;
    ySize=8;
    xLeft=(21-xSize)/2; 
    yTop= (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
    set(gcf,'Position',[24 36 xSize*50 ySize*50]);
    
    plot(arclength,efitBetas(1,:,mii),'-b','LineWidth', 2);
    hold on;
    plot(arclength,efitBetas(2,:,mii),'-r','LineWidth', 2);
    hold on;
    plot(arclength,efitBetas(3,:,mii),'-g','LineWidth', 2);
    hold on;
    plot(arclength,efitBetas(4,:,mii),'-k','LineWidth', 2);
    hold off;
    h=legend('Intercerpt','Gender','Flu','Age',4);
    set(h,'Location','NorthWest');
    set(h,'interpreter','latex');
    set(h,'FontSize',10);
    
    xlabel('arc-length');
    ylabel(Mnames{mii});
    set(gca,'XLim',[0 size(arclength,1)-1]);
    
    figurename=sprintf('%s/%s_beta.bmp',plotpath,Mnames{mii});
    saveas(gcf,figurename);
    close()

    end
    
        %% 2.2 smoothing individual function
        
        ResYdesign=Ydesign-efitYdesign;
        [ResEtas,efitEtas,eSigEta]=MVCM_sif(arclength,ResYdesign);
        
        [mSigEtaEig, mSigEta]=MVCM_eigen(efitEtas);
        
    %% plot eigenvalue

    figure()
    set(gcf,'PaperUnits','centimeters');
    xSize=12;
    ySize=8;
    xLeft=(21-xSize)/2; 
    yTop= (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
    set(gcf,'Position',[24 36 xSize*50 ySize*50]);
    
    plot(mSigEtaEig(1:12,1,1),'-r.','LineWidth', 2);
    hold on;
    plot(mSigEtaEig(1:12,1,2),'-g.','LineWidth', 2);
    hold off; 
    
    set(gca,'XLim',[1 12]);
    set(gca,'YLim',[0 .7]);
    h=legend('FA','MD',2);
    set(h,'Location','NorthEast');
    set(h,'interpreter','latex');
    set(h,'FontSize',14)
    
    figurename=sprintf('%s/eigenvalue.bmp',plotpath);
    saveas(gcf,figurename);
    close()
    
    %% plot eigenvector 

    for mii=1:m
        
        figure()
        set(gcf,'PaperUnits','centimeters');
        xSize=12;
        ySize=8;
        xLeft=(21-xSize)/2; 
        yTop= (30-ySize)/2;
        set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
        set(gcf,'Position',[24 36 xSize*50 ySize*50]);
        
        %subplot(2,3,mii+1);
        plot(arclength,mSigEtaEig(1,2:end,mii),'-r.','LineWidth', 2);
        hold on;
        plot(arclength,mSigEtaEig(2,2:end,mii),'-g.','LineWidth', 2);
        hold on;
        plot(arclength,mSigEtaEig(3,2:end,mii),'-b.','LineWidth', 2);
        hold on;
        h=legend('1st eigenvector','2nd eigenvector','3rd eigenvector',3);
        set(gca,'XLim',[0 L0-1])
        set(gca,'YLim',[-.4 .6])
        set(h,'Location','NorthWest');
        set(h,'interpreter','latex');
        set(h,'FontSize',14);
         xlabel('arc-length');
         ylabel(Mnames{mii});
   
        figurename=sprintf('%s/%s_eigenvector.bmp',plotpath,Mnames{mii});
        saveas(gcf,figurename);
        close()
        
    end
    
        %% 2.3 estimating covariance matrix
        
        [eSigE]=MVCM_ecm(arclength,ResEtas);
        
        %% 2.4 Hypothesis test
        
        n=NoSetup(1); 
        L0=NoSetup(2); 
        p=NoSetup(3); 
        m=NoSetup(4);
        
        Gstats=zeros(1,p-1);
        Lstats=zeros(L0,p-1);
        Gpvals=zeros(1,p-1);
                 
        [ebiasBetas] = MVCM_bias(NoSetup,arclength,Xdesign,Ydesign,InvSigmats,mh);
           
         for pp=2:p
            
            %individual and global statistics calculation
            cdesign=zeros(1,p);
            cdesign(pp)=1;
            Cdesign=kron(eye(m),cdesign);
            B0vector=zeros(m,L0);
            [Gstat,Lstat] = MVCM_ht_stat(NoSetup,arclength,Xdesign,efitBetas,eSigEta,Cdesign,B0vector,ebiasBetas);
            Gstats(1,pp-1)=Gstat;
            Lstats(:,pp-1)=Lstat;
           
                                           
            % Generate random samples and calculate the corresponding statistics and pvalues
            GG=1000;
            [Gpval] = MVCM_bstrp_pvalue3(NoSetup,arclength,Xdesign,Ydesign,efitBetas1,InvSigmats,mh,Cdesign,B0vector,Gstat,GG); 
            Gpvals(1,pp-1)=Gpval;
                        
        end
        
        Lpvals=1-chi2cdf(Lstats,m);
        
    %% plot local p values
    
    figure()
    set(gcf,'PaperUnits','centimeters');
    xSize=12;
    ySize=8;
    xLeft=(21-xSize)/2; 
    yTop= (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
    set(gcf,'Position',[24 36 xSize*50 ySize*50]);
    
    plot(arclength,-log10(Lpvals(:,1)),'-r.','LineWidth', 2,'MarkerSize',15);
    hold on;
    plot(arclength,-log10(Lpvals(:,2)),'-g.','LineWidth', 2,'MarkerSize',15);
    hold on;
    plot(arclength,-log10(Lpvals(:,3)),'-b.','LineWidth', 2,'MarkerSize',15);
    hold off;
    xlabel('arc-length');
    ylabel('-log10(p)');
    set(gca,'XLim',[0 size(arclength,1)-1]);
    h=legend('Gender','Flu','Age',4);
    set(h,'Location','NorthWest');
    set(h,'interpreter','latex');
    set(h,'FontSize',12);
    
    figurename=sprintf('%s/local_p_value.bmp',plotpath);
    saveas(gcf,figurename);
    close()
    
              
        %% 2.5 confidence band
        %GG=1000;
         [Gvalue] = MVCM_cb_Gval(arclength,Xdesign,ResYdesign,InvSigmats,mh,GG);
         alpha=0.05;
         [CBands] = MVCM_CBands(n,alpha,Gvalue,efitBetas,ebiasBetas);

    %% plot confidence bands
    covnames=cell(4,1);
    covnames{1}='inter';
    covnames{2}='gender';
    covnames{3}='flu';
    covnames{4}='age';
    ylabs=cell(4,1);
    ylabs{1}='95% confidence band for intercepts';
    ylabs{2}='95% confidence band for gender';
    ylabs{3}='95% confidence band for flu';
    ylabs{4}='95% confidence band for age';
    
    for mii=1:m
    for pii=1:p
    figure()
    set(gcf,'PaperUnits','centimeters');
    xSize=12;
    ySize=8;
    xLeft=(21-xSize)/2; 
    yTop= (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
    set(gcf,'Position',[24 36 xSize*50 ySize*50]);
    
    %subplot(2,2,2);
    plot(arclength,efitBetas(pii,:,mii),'-b','LineWidth', 2);
    hold on;
    plot(arclength,CBands(2*pii-1,:,mii),'--r','LineWidth', 2);
    hold on;
    plot(arclength,CBands(2*pii,:,mii),'--r','LineWidth', 2);
    xlabel('arc-length');
    ylabel(ylabs{mii})
    set(gca,'XLim',[0 size(arclength,1)-1]);
    
    figurename=sprintf('%s/%s_%s.bmp',plotpath,Mnames{mii},covnames{pii});
    saveas(gcf,figurename);
    close()
    
    end
    end


 


  