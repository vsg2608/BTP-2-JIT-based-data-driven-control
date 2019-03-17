
I_0= 0.00258;              %Initial moles of initiator
M_0= 1E4;                  %Initial moles of monomer
T=319.55;
R_lm =0;
Tf=10;
C=MMA_Simulation(I_0, M_0, T, R_lm, Tf);

Data= importdata('Data@10.txt');    %Training Data
X= Data(:,1:2);
Y = Data(:,3);
Xtrain= X(1:6000,:);
Ytrain= Y(1:6000);
[XL,YL,XS,YS,BETA] = plsregress(Xtrain, Ytrain);
Yact = Y(6001);
Ypred = [ones(1,1) X(6001,:)]*BETA;
