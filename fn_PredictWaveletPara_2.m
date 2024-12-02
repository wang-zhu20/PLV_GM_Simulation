%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get predicted wavelet parameter from emperical equation 
%%% Sept 21, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 %% Note: totalEnergy should *0.01
function [outprm]=fn_PredictWaveletPara_2(M,Rrup,Rhyp,Vs30)  %%%% two columns
% M=7;
% Rrup=100;
% Rhyp=100;
% Vs30=270;
outprm = fn_initDBs(3);
%%
load covrho_3dim.mat

%%%%%%%%%%%%%%%%%%  Baker's Prediction Equation          %%%%%%%%%%%%%%%%
%%%% Table 3.2: Coefficients of the prediction equation, p. 59  %%%%%%%%%
%%%% alpha  beta1  beta2  beta3  beta4  beta5  beta6  h  sigma  tau  %%%%


flag=1;

%%%% WZ's regression %%%%%%%%%%%%%%%%%
while true

    epst = mvnrnd(zeros(36,1),cov_total);
    
    coeff=[2.470	0	0	0.00056	-0.0022	0.269	-0.200	1	0.223	0.247
        2.645	0	0	0.00052	-0.0072	0.164	-0.183	1	0.285	0.277
        1.282	-0.1416	0	0	-0.0040	-0.235	0.356	10	0.349	0.252
        0.975	0.0167	0	0	-0.0038	-0.318	0.303	10	0.428	0.322
        -0.445	0.0158	0	0	-0.0007	-0.030	0.053	10	0.059	0.032
        2.262	0	0	0.00060	-0.0023	0.312	-0.224	1	0.255	0.266
        2.527	0	0	0.00061	-0.0084	0.185	-0.296	1	0.323	0.289
        1.125	-0.2284	0	0	-0.0040	-0.183	0.386	10	0.368	0.233
        0.755	-0.1355	0	0	-0.0043	-0.226	0.300	10	0.391	0.260
        -0.703	0.0014	0	0	0.0002	-0.043	0.096	10	0.129	0.056
        -33.007	-3.4944	30.896	0	0	-1.729	-0.970	10	1.076	0.664
        -23.882	-2.2134	23.954	0	0	-1.748	-0.860	10	0.815	0.476
        1.2144	0	0	0	0	0	0	0	0.0994	0];
    alpha=coeff(:,1);
    beta1=coeff(:,2);
    beta2=coeff(:,3);
    beta3=coeff(:,4);
    beta4=coeff(:,5);
    beta5=coeff(:,6);
    beta6=coeff(:,7);
    h=coeff(:,8);
    sigma=coeff(:,9);
    tau=coeff(:,10);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      for jj=1:24
% intra_w(jj)=sqrt(covm(jj,jj));
% inter_w(jj)=sqrt(covi(jj,jj));
%      end
    %%%%%%%%%%%%%%%%%%%  Compute Median Prediction of Wavelet Parameters %%%
    
    Y=alpha+beta1*M + beta2*log(M)+beta3*exp(M)+beta4*(Rhyp-Rrup)+beta5.*log(sqrt(Rrup.^2+h.^2))+beta6*log(Vs30);
    if flag==0
        minorEx = exp(Y(1));
        minorSx = exp(Y(2));
        minorEy = exp(Y(3));
        minorSy = exp(Y(4));
        minorRxy = Y(5);      minorRxy =2*normcdf(minorRxy,0,1)-1;
        majorEx = exp(Y(6));
        majorSx = exp(Y(7));
        majorEy = exp(Y(8));
        majorSy = exp(Y(9));
        majorRxy =Y(10);        majorRxy=2*normcdf(majorRxy,0,1)-1;
        totalEnergy = exp(Y(12));
        majorEa  = exp(Y(11));
    else

        minorEx = exp(Y(1)+epst(1));
        minorSx = exp(Y(2)+epst(2));
        minorEy = exp(Y(3)+epst(3));
        minorSy = exp(Y(4)+epst(4));
        minorRxy = Y(5)+epst(5);      minorRxy =2*normcdf(minorRxy,0,1)-1;
        majorEx = exp(Y(6)+epst(6));
        majorSx = exp(Y(7)+epst(7));
        majorEy = exp(Y(8)+epst(8));
        majorSy = exp(Y(9)+epst(9));
        majorRxy =Y(10)+epst(10);        majorRxy=2*normcdf(majorRxy,0,1)-1;
        totalEnergy = exp(Y(12)+epst(12));
        majorEa  = exp(Y(11)+epst(11));
    end
    %%%% Version 1: random is longmornal
    %minorRnd = exp(Y(13));
    mn13=coeff(13,1);
    sn13=coeff(13,9);
    m13=log(mn13^2/sqrt(sn13^2+mn13^2));
    s13=sqrt(log(sn13^2/mn13^2+1));
    minorRnd=lognrnd(m13,s13);
    
    %% to lxly
    majorVx = majorSx^2;
    majorVy = majorSy^2;
    majorExy = majorRxy*majorSx*majorSy + majorEx*majorEy;
    
    majorElx = log(majorEx^2 / sqrt(majorVx+majorEx^2));
    majorVlx = log(majorVx/majorEx^2 + 1);
    majorSlx = sqrt(majorVlx);
    majorEly = log(majorEy^2 / sqrt(majorVy+majorEy^2));
    majorVly = log(majorVy/majorEy^2 + 1);
    majorSly = sqrt(majorVly);
    
    majorCovlxly = log(majorExy/majorEx/majorEy);
    majorRlxly = majorCovlxly/majorSlx/majorSly;
    
    minorVx = minorSx^2;
    minorVy = minorSy^2;
    minorExy = minorRxy*minorSx*minorSy + minorEx*minorEy;
    
    minorElx = log(minorEx^2 / sqrt(minorVx+minorEx^2));
    minorSlx = sqrt(log(minorVx/minorEx^2 + 1));
    minorVlx = minorSlx^2;
    minorEly = log(minorEy^2 / sqrt(minorVy+minorEy^2));
    minorSly = sqrt(log(minorVy/minorEy^2 + 1));
    minorVly = minorSly^2;
    
    minorCovlxly = log(minorExy/minorEx/minorEy);
    minorRlxly = minorCovlxly/minorSlx/minorSly;
    
    majorLLCov=[majorVlx majorCovlxly; majorCovlxly majorVly];
    [T,err] = cholcov(majorLLCov);
    if err ~= 0
        continue
    end
    outprm(1) = struct( 'M',M,'hdist',Rrup,'vs30',Vs30,...
        'minorElx',minorElx, 'minorSlx',minorSlx, 'minorVlx',minorVlx, 'minorEly',minorEly, 'minorSly',minorSly, 'minorVly',minorVly,'minorRlxly',minorRlxly, ...
        'majorElx',majorElx, 'majorSlx',majorSlx, 'majorVlx',majorVlx, 'majorEly',majorEly, 'majorSly',majorSly, 'majorVly',majorVly,'majorRlxly',majorRlxly, ...
        'totalEnergy',totalEnergy, 'majorEa',majorEa, 'minorRnd',minorRnd);
%%

    Y=alpha+beta1*M + beta2*log(M)+beta3*exp(M)+beta4*(Rhyp-Rrup)+beta5.*log(sqrt(Rrup.^2+h.^2))+beta6*log(Vs30);
    if flag==0
        minorEx = exp(Y(1));
        minorSx = exp(Y(2));
        minorEy = exp(Y(3));
        minorSy = exp(Y(4));
        minorRxy = Y(5);      minorRxy =2*normcdf(minorRxy,0,1)-1;
        majorEx = exp(Y(6));
        majorSx = exp(Y(7));
        majorEy = exp(Y(8));
        majorSy = exp(Y(9));
        majorRxy =Y(10);        majorRxy=2*normcdf(majorRxy,0,1)-1;
        totalEnergy = exp(Y(12));
        majorEa  = exp(Y(11));
    else

        minorEx = exp(Y(1)+epst(1+12));
        minorSx = exp(Y(2)+epst(2+12));
        minorEy = exp(Y(3)+epst(3+12));
        minorSy = exp(Y(4)+epst(4+12));
        minorRxy = Y(5)+epst(5+12);      minorRxy =2*normcdf(minorRxy,0,1)-1;
        majorEx = exp(Y(6)+epst(6+12));
        majorSx = exp(Y(7)+epst(7+12));
        majorEy = exp(Y(8)+epst(8+12));
        majorSy = exp(Y(9)+epst(9+12));
        majorRxy =Y(10)+epst(10+12);        majorRxy=2*normcdf(majorRxy,0,1)-1;
        totalEnergy = exp(Y(12)+epst(12+12));
        majorEa  = exp(Y(11)+epst(11+12));
    end
    %%%% Version 1: random is longmornal
    %minorRnd = exp(Y(13));
    mn13=coeff(13,1);
    sn13=coeff(13,9);
    m13=log(mn13^2/sqrt(sn13^2+mn13^2));
    s13=sqrt(log(sn13^2/mn13^2+1));
    minorRnd=lognrnd(m13,s13);
    
    %% to lxly
    majorVx = majorSx^2;
    majorVy = majorSy^2;
    majorExy = majorRxy*majorSx*majorSy + majorEx*majorEy;
    
    majorElx = log(majorEx^2 / sqrt(majorVx+majorEx^2));
    majorVlx = log(majorVx/majorEx^2 + 1);
    majorSlx = sqrt(majorVlx);
    majorEly = log(majorEy^2 / sqrt(majorVy+majorEy^2));
    majorVly = log(majorVy/majorEy^2 + 1);
    majorSly = sqrt(majorVly);
    
    majorCovlxly = log(majorExy/majorEx/majorEy);
    majorRlxly = majorCovlxly/majorSlx/majorSly;
    
    minorVx = minorSx^2;
    minorVy = minorSy^2;
    minorExy = minorRxy*minorSx*minorSy + minorEx*minorEy;
    
    minorElx = log(minorEx^2 / sqrt(minorVx+minorEx^2));
    minorSlx = sqrt(log(minorVx/minorEx^2 + 1));
    minorVlx = minorSlx^2;
    minorEly = log(minorEy^2 / sqrt(minorVy+minorEy^2));
    minorSly = sqrt(log(minorVy/minorEy^2 + 1));
    minorVly = minorSly^2;
    
    minorCovlxly = log(minorExy/minorEx/minorEy);
    minorRlxly = minorCovlxly/minorSlx/minorSly;
    
    majorLLCov=[majorVlx majorCovlxly; majorCovlxly majorVly];
    [T,err] = cholcov(majorLLCov);
    if err ~= 0
        continue
    end
    outprm(2) = struct( 'M',M,'hdist',Rrup,'vs30',Vs30,...
        'minorElx',minorElx, 'minorSlx',minorSlx, 'minorVlx',minorVlx, 'minorEly',minorEly, 'minorSly',minorSly, 'minorVly',minorVly,'minorRlxly',minorRlxly, ...
        'majorElx',majorElx, 'majorSlx',majorSlx, 'majorVlx',majorVlx, 'majorEly',majorEly, 'majorSly',majorSly, 'majorVly',majorVly,'majorRlxly',majorRlxly, ...
        'totalEnergy',totalEnergy, 'majorEa',majorEa, 'minorRnd',minorRnd);

    %%  for vertical  
    
      Coeff = [-2.342426488	0.68799677	0.227414135	-0.097600939
-1.657789454	0.596876915	0.286321993	-0.148944489
2.205870259	0.035744168	-0.098301376	0.018698363
1.277279461	0.149512253	0.014687628	0.050559948
-0.544966685	0.043923012	-0.01824002	0.021124879
-3.055804437	0.747193353	0.267777643	-0.09716879
-2.269888287	0.586818353	0.363700077	-0.143378833
2.986360677	-0.004213658	-0.109479299	-0.106823394
2.262305805	0.097405858	-0.004620375	-0.122844183
-0.533706594	0.034498168	-0.025133858	0.025108618
-6.904231234	1.251257109	-1.372244069	-0.323508914
-9.847664715	0.700461509	-1.470451386	-0.093534688];




base_x = [1, M, log(sqrt(Rrup^2+36)),log(Vs30)];
Y = Coeff*base_x';

    if flag==0
        minorEx = exp(Y(1));
        minorSx = exp(Y(2));
        minorEy = exp(Y(3));
        minorSy = exp(Y(4));
        minorRxy = Y(5);      minorRxy =2*normcdf(minorRxy,0,1)-1;
        majorEx = exp(Y(6));
        majorSx = exp(Y(7));
        majorEy = exp(Y(8));
        majorSy = exp(Y(9));
        majorRxy =Y(10);        majorRxy=2*normcdf(majorRxy,0,1)-1;
        totalEnergy = exp(Y(11));
        majorEa  = exp(Y(12));
    else
        minorEx = exp(Y(1)+epst(13+12));
        minorSx = exp(Y(2)+epst(14+12));
        minorEy = exp(Y(3)+epst(15+12));
        minorSy = exp(Y(4)+epst(16+12));
        minorRxy = Y(5)+epst(17+12);      minorRxy =2*normcdf(minorRxy,0,1)-1;
        majorEx = exp(Y(6)+epst(18+12));
        majorSx = exp(Y(7)+epst(19+12));
        majorEy = exp(Y(8)+epst(20+12));
        majorSy = exp(Y(9)+epst(21+12));
        majorRxy =Y(10)+epst(22+12);        majorRxy=2*normcdf(majorRxy,0,1)-1;
        totalEnergy = exp(Y(11)+epst(24+12));
        majorEa  = exp(Y(12)+epst(23+12));
    end
    %%%% Version 1: random is longmornal
    %minorRnd = exp(Y(13));
    mn13=coeff(13,1);
    sn13=coeff(13,9);
    m13=log(mn13^2/sqrt(sn13^2+mn13^2));
    s13=sqrt(log(sn13^2/mn13^2+1));
    minorRnd=lognrnd(m13,s13);
    
    %% to lxly
    majorVx = majorSx^2;
    majorVy = majorSy^2;
    majorExy = majorRxy*majorSx*majorSy + majorEx*majorEy;
    
    majorElx = log(majorEx^2 / sqrt(majorVx+majorEx^2));
    majorVlx = log(majorVx/majorEx^2 + 1);
    majorSlx = sqrt(majorVlx);
    majorEly = log(majorEy^2 / sqrt(majorVy+majorEy^2));
    majorVly = log(majorVy/majorEy^2 + 1);
    majorSly = sqrt(majorVly);
    
    majorCovlxly = log(majorExy/majorEx/majorEy);
    majorRlxly = majorCovlxly/majorSlx/majorSly;
    
    minorVx = minorSx^2;
    minorVy = minorSy^2;
    minorExy = minorRxy*minorSx*minorSy + minorEx*minorEy;
    
    minorElx = log(minorEx^2 / sqrt(minorVx+minorEx^2));
    minorSlx = sqrt(log(minorVx/minorEx^2 + 1));
    minorVlx = minorSlx^2;
    minorEly = log(minorEy^2 / sqrt(minorVy+minorEy^2));
    minorSly = sqrt(log(minorVy/minorEy^2 + 1));
    minorVly = minorSly^2;
    
    minorCovlxly = log(minorExy/minorEx/minorEy);
    minorRlxly = minorCovlxly/minorSlx/minorSly;
    
    majorLLCov=[majorVlx majorCovlxly; majorCovlxly majorVly];
    [T,err] = cholcov(majorLLCov);
    if err == 0
        break
    end
end
outprm(3) = struct( 'M',M,'hdist',Rrup,'vs30',Vs30,...
                'minorElx',minorElx, 'minorSlx',minorSlx, 'minorVlx',minorVlx, 'minorEly',minorEly, 'minorSly',minorSly, 'minorVly',minorVly,'minorRlxly',minorRlxly, ...
                'majorElx',majorElx, 'majorSlx',majorSlx, 'majorVlx',majorVlx, 'majorEly',majorEly, 'majorSly',majorSly, 'majorVly',majorVly,'majorRlxly',majorRlxly, ...
                'totalEnergy',totalEnergy, 'majorEa',majorEa, 'minorRnd',minorRnd); 
             
%% save('MedianWaveletParam.mat','MedianWaveletParam');