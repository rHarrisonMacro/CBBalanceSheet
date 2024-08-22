tidyUpAndSetPath;

%% OPTIONS
varsToPlot =    {'pie';'x';'q';'ir'};
scalingFactors = [100; 100; 1; 400];
varNames = {'Quarterly inflation, %';'Output gap, %';'QE (q)'; 'Policy rate, %'};
resultsFileName = [outputsFolder 'normalisationSim.mat'];

nCols = 2;

plotHorizon = 8;

includeZeroOnMacroVarCharts = true;

%% LOAD RESULTS
load(resultsFileName,'xSimBase','xSimAlt','ExSimBase','ExSimAlt',...
    'rstarSequence','pars','q0','q0alt');

constants =      [0;  0;   0;  400*log(1/pars.beta)]; 

%% BASIC INFO
endogVarMnems = {'pie';'x';'q';'mux';'ir';'muxL';'muxH';'varsigChk'};
nEndogVars = size(endogVarMnems,1);
nVarsToPlot = size(varsToPlot,1);
varInds = lookup_model_index_numbers(endogVarMnems,varsToPlot);
nRows = ceil((nVarsToPlot)/nCols);
T = size(rstarSequence,2);
T = min(T,plotHorizon);

%% SUMMARY CHART
fontSize = 12;
summaryFigHandle = figure;
% subplot(1,4,1);
subplot(2,2,1);
hold on;
plot(1:T,xSimBase(3,1:T),'k','linewidth',2);
plot(1:T,xSimAlt(3,1:T),'linestyle','--','color',0.5*ones(1,3),'linewidth',2);
title('QE stock');
xlim([1 T]);
ylim([0 0.75]);
xlabel('Quarters');
set(gca,'FontSize',fontSize);
% subplot(1,4,2);
subplot(2,2,2);
hold on;
plot(1:T,xSimBase(5,1:T)*400+400*log(1/pars.beta),'k','linewidth',2);
plot(1:T,xSimAlt(5,1:T)*400+400*log(1/pars.beta),'linestyle','--',...
    'color',0.5*ones(1,3),'linewidth',2);
title('Policy rate, %');
xlim([1 T]);
xlabel('Quarters');
set(gca,'FontSize',fontSize);
% subplot(1,4,3);
subplot(2,2,3);
hold on;
plot(1:T,xSimBase(1,1:T)*100,'k','linewidth',2);
plot(1:T,xSimAlt(1,1:T)*100,'linestyle','--','color',0.5*ones(1,3),...
    'linewidth',2);
title('Quarterly inflation, %');
if includeZeroOnMacroVarCharts
    plot(1:T,zeros(1,T),'k','linewidth',0.5);
end
xlim([1 T]);
xlabel('Quarters');
set(gca,'FontSize',fontSize);
% subplot(1,4,4);
subplot(2,2,4);
hold on;
plot(1:T,xSimBase(2,1:T)*100,'k','linewidth',2);
plot(1:T,xSimAlt(2,1:T)*100,'linestyle','--','color',0.5*ones(1,3),...
    'linewidth',2);
if includeZeroOnMacroVarCharts
    plot(1:T,zeros(1,T),'k','linewidth',0.5);
end
xlim([1 T]);
xlabel('Quarters');
title('Output gap, %');
set(gca,'FontSize',fontSize);
summaryLegHandle = legend('With state contingent QE','No state contingency');
legend boxoff;

% figDims = [12 5.75]*1;
figDims = [12 10]*0.8;

summaryFigName = [figureFolder 'FigureA1.pdf'];
adjx = 1;
adjy = 0;

set(summaryFigHandle,'PaperUnits','inches',...
    'PaperSize',(figDims),...
    'PaperPosition',[-adjx -adjy figDims+[adjx*1.5 adjy]]);


summaryLegHandle.Orientation = 'horizontal';
summaryLegHandle.FontSize = fontSize;
oldLeft = summaryLegHandle.Position(1);
oldBottom = summaryLegHandle.Position(2);
legWidth = summaryLegHandle.Position(3);
legHeight = summaryLegHandle.Position(4);
newLeft = 0.5 - legWidth/2;
newBottom = 0;
summaryLegHandle.Position = [newLeft newBottom legWidth legHeight];

print(summaryFigHandle,'-dpdf',summaryFigName);
open(summaryFigName);

%% HEADROOM
zeta = pars.epsilon + pars.epsilonDelta*(1+pars.beta);
xiL = pars.xiL;
xiH = pars.xiH;
ExiMuxCont = -pars.beta/zeta*pars.sigma*...
    (xiL*ExSimBase(7,:)+xiH*ExSimBase(6,:));

ExiMuxHatCont = -pars.beta/zeta*pars.sigma*...
    (xiL*ExSimAlt(7,:)+xiL*ExSimAlt(6,:));

H = ExiMuxCont - ExiMuxHatCont; 

potencyCont = -pars.beta/zeta*pars.sigma*(xiH-xiL)*ExSimBase(6,:);
stabCont = ...
    pars.beta/zeta*pars.sigma*xiL*(ExSimAlt(4,:)-ExSimBase(4,:));

%% RESCALE THE FIGURE ABOVE
debtToGDP = 0.5;
scaling = 100/debtToGDP; 

%% BAR CHART
decompFigHandle = figure;
allBars = [scaling*potencyCont(1:T); scaling*stabCont(1:T)];
barContnsPos = allBars;
barContnsPos(barContnsPos<0)=0;
barContnsNeg = allBars;
barContnsNeg(barContnsNeg>0)=0;
hold on;
barHandle = bar(barContnsPos','stacked');
negBarHandle = bar(barContnsNeg','stacked');
barHandle(1).FaceColor = 0.75*ones(1,3);
barHandle(2).FaceColor = 0.25*ones(1,3);
% Force colours to be aligned
nBars = size(barHandle,2);
for iBar = 1:nBars
    negBarHandle(iBar).FaceColor = barHandle(iBar).FaceColor;
end
xlim([0.5 T+0.5]);
xlabel('Quarters');
lineHandle = plot(1:T,scaling*H(1:T),'k-o','linewidth',2);
legHandle = legend([barHandle,lineHandle],'Potency','Stabilisation',...
    'Effect of state contingency on incentive to build headroom','Location','NorthEast');

legend boxoff;
grid on;
set(gca,'FontSize',fontSize);
legHandle.FontSize = fontSize;

% Add annotations to the chart
fasterUnwind = annotation('arrow',[0.9125 0.9125],[0.575 0.125]);
textboxFast = annotation('textbox',[0.92 0.4 0.08 0.05],...
    'string','Faster QE unwind','LineStyle','none');

slowerUnwind = annotation('arrow',[0.9125 0.9125],[0.625 0.9125]);
textboxSlow = annotation('textbox',[0.92 0.75 0.08 0.05],...
    'string','Slower QE unwind','LineStyle','none');

textboxAxis = annotation('textbox',[0.0625 0.9375 0.2 0.05],...
    'string','% annual GDP','LineStyle','none');


%% SAVE FIGURE
figDims = [12 8]*0.6;

figureName = [figureFolder 'Figure5.pdf'];
adjx = 0*0.5;
adjy = 0;
set(decompFigHandle,'PaperUnits','inches',...
    'PaperSize',(figDims),...
    'PaperPosition',[-adjx -adjy figDims+[adjx*1.5 adjy]]);
print(decompFigHandle,'-dpdf',figureName);

%% CALCULATIONS FOR PAPER
q4Hpercent = mean(scaling*H(4));

million = 1000000;
billion = 1000*million;
GDP2019 = 2255283*million;
meanHpounds = GDP2019*q4Hpercent/100;
meanHpoundsBillions = meanHpounds/billion;
disp(['H in percent of GDP (Q4): ' num2str(q4Hpercent)]);
disp(['H in 2019 GDP £ billions (Q4): ' num2str(meanHpoundsBillions)]);

