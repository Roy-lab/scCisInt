function scCisInt_EP_predict(region_id, ntrees, indir, output_path, verbose)

%% ---- hanndle default input args --------
if nargin < 3 || isempty(indir)
    indir = "Results/";
end

if nargin < 4 || isempty(output_path)
    output_path = "Results/";
end

if nargin < 5 || isempty(verbose)
    verbose = false;
end

%% ----- make output directory -------
if ~exist(output_path)
    mkdir(output_path)
end

%% ----- begin RF procedure --------
if verbose
    sprintf('region_id=%s;treecnt=%s;', region_id, ntrees)
end

if(isstring(ntrees))
    ntrees=str2num(ntrees);
end
treecnt=ntrees;

%% ---- Load data ---------
infile = sprintf('%s/%s_enhancer_data.txt',indir, region_id);
infile1 = sprintf('%s/%s_promoter_data.txt',indir, region_id);
if exist(infile1, 'file')==0
    error('Input file %s does not exist. Please check the file path.', infile1);
end

d=importdata(infile);
enhancers=d.textdata(1,2:end)';
X0=d.data;
p=importdata(infile1);
y0=p.data;

k=5;
validation_set_size=floor(size(X0,1)/k);
all_ynews=[];
all_ys=[];
%pid=1;
%fig=figure;
varimp=[];
forest={};
for fold=1:k

    if (verbose)
        fprintf('Doing CV in %d\n',fold);
    end
    

    %% --- generate train test sets -------
    if fold==k 
        train_range=1:(fold-1)*validation_set_size;
        validation_range=(fold-1)*validation_set_size+1:size(X0,1);
    else
        train_range=[1:(fold-1)*validation_set_size,fold*validation_set_size+1:size(X0,1)];
        validation_range=(fold-1)*validation_set_size+1:fold*validation_set_size;
    end
    
    %% ----- Run RF ------------
    X=X0(train_range,:); %X0(:,1:size(train_data,2)-1);
    Y=y0(train_range,:); %y0(:,size(train_data,2));
    rand('seed',fold);
    B=TreeBagger(treecnt,X,Y,'NVarToSample','all','Method','regression','OOBVarImp','on');

    %% compute Var Importance
    varimp(:,fold)=B.OOBPermutedVarDeltaError';
    forest{fold}=B;

    ytrain=predict(B,X);
    traincc=corrcoef(Y,ytrain);

    if(verbose)
        fprintf('Fold %d train cc=%.3f\n',fold,traincc(1,2));
    end

    Xtest=X0(validation_range,:);
    ynew=predict(B,Xtest);
    Ytest=y0(validation_range,:);
    testcc=corrcoef(Ytest,ynew);
    
    if(verbose)
        fprintf('Fold %d test cc=%.3f\n',fold,testcc(1,2));
    end

    all_ynews=[all_ynews,ynew'];
    all_ys=[all_ys,Ytest'];
end

varimp_mean=mean(varimp')';
[ig,vind]=sort(varimp_mean,'descend');

file_name = sprintf('%s/%s_consensus_EP_predictions.txt', output_path, region_id);
fid=fopen(file_name,'w');
for i=1:size(ig,1)
    fprintf(fid,'%s\t%s\t%f\n',enhancers{vind(i)},region_id,varimp_mean(vind(i)));
end
fclose(fid);
end

