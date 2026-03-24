function [consensus_matrix]=makeConsensusMatrix(dir, U_file, V_file, num_iterations, outdir)
for i=1:num_iterations
    i;
    cell_clust=importdata(sprintf('%s/I_%i/%s', dir, i, U_file));
    gene_clust=importdata(sprintf('%s/I_%i/%s', dir, i, V_file));
    
    if(isstruct(gene_clust))
        gene_clust = gene_clust.data; % Extract data if gene_clust is a struct
    end

    if(isstruct(cell_clust))
        cell_clust = cell_clust.data;   % Extract data if cell_clust is a struct
    end

    if(size(cell_clust, 2) > 1 || size(gene_clust, 2) > 1)
        error("cell and gene cluster must only contain a single cluster assignment per cell.")
    end 
    
    if any(mod(cell_clust(:), 1) ~= 0) || any(mod(gene_clust(:), 1) ~= 0)
        error("cell and gene cluster values must be integers")

    end
 
    if i==1;
        consensus_matrix=zeros(length(cell_clust), length(gene_clust));
    end
    for j=1:length(unique(cell_clust))
        j;
        cellidx=find(cell_clust==j);
        geneidx=find(gene_clust==j);
        rowidx=repmat(cellidx, length(geneidx),1);
        columnidx=reshape(repmat(geneidx,1, length(cellidx))',[],1); 
        idc=sub2ind(size(consensus_matrix), rowidx,columnidx);
        consensus_matrix(idc)=consensus_matrix(idc)+1;        
    end
end

mkdir(outdir)

save(sprintf('%s/consensus_matrix.mat',outdir), 'consensus_matrix');
writematrix(consensus_matrix, sprintf('%s/consensus_matrix.txt', outdir), "Delimiter", 'tab');

[~,sort_gene_idx]=sort(gene_clust);
[~,sort_cell_idx]=sort(cell_clust);

imagesc(consensus_matrix(sort_cell_idx, sort_gene_idx));
saveas(gcf(),sprintf('%s/consensus_matrix_sorted.png', dir));


