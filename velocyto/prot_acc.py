#Gennady Gorin, Pachter Laboratory, Caltech, 9/17/19

#version for velocyto prot_acc.
#adds a suite of tools to do protein acceleration based on feature barcoding data. 
import numpy as np
import matplotlib.pyplot as plt
import csv
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from scipy import sparse
import louvain
import igraph as ig
from .estimation import fit_slope_weighted_offset,  colDeltaCorSqrt
import bezier
from sklearn.manifold import TSNE

def import_prot_data(file_path):
    """
    Imports and parses csv file with feature barcoding information, located at file_path
    Format: rows = ADT (antibody-derived tag) names, columns = cells.  
    Outputs numpy array of gene counts, cell indices, and observed ADT names.
    """
    with open(file_path,'r') as dest_f:
        data_iter = csv.reader(dest_f,
                               delimiter = ',',
                               quotechar = '"')
        data = [data for data in data_iter]
    prot_array = np.asarray(data)
    prot_cells = prot_array[0,1:]
    adt_names = prot_array[1:,0]
    prot_count_array = prot_array[1:,1:]
    return [prot_count_array,prot_cells,adt_names]

def enforce_protein_filter(vlm, mRNA_names,adt_names):
    """
    Forces velocyto filters to retain genes that have protein information, even if they are
    substandard for RNA velocity. Stores ADT names and corresponding mRNA names in the velocyto loom structure vlm.
    """
    chosen_genes = np.asarray([np.where(vlm.ra['Gene'] == mRNA_names[i])[0][0] for i in range(len(mRNA_names))])
    if hasattr(vlm, 'cv_mean_selected'):
        vlm.cv_mean_selected[chosen_genes] = True
    if hasattr(vlm, 'clu_avg_selected'):
        vlm.clu_avg_selected[chosen_genes] = True
    if hasattr(vlm, 'detection_level_selected'):
        vlm.cv_mean_selected[chosen_genes] = True
    vlm.mRNA_names = np.asarray(mRNA_names)
    vlm.adt_names = adt_names
    
def shared_cells_filter(vlm, prot_cells, prot_count_array, first_char, last_char):
    """
    Retains only cells with protein and RNA information. Reports the number of cells in each array.
    Requires the location of the first and last characters of the cell barcode in the vlm structure.
    The cell tages in csv are presumde to be pre-processed.
    Saves protein count array as vlm.P.
    """
    if last_char == 0:
        rna_cells = np.asarray([cell_str[first_char:] for cell_str in vlm.ca['CellID']])
    else:
        rna_cells = np.asarray([cell_str[first_char:last_char] for cell_str in vlm.ca['CellID']])
    shared_cells = [cell_tag for cell_tag in rna_cells if cell_tag in prot_cells]

    print('ADT cell number: '+str(len(prot_cells)))
    print('RNAseq cell number: '+str(len(rna_cells)))
    print('Shared cells: '+str(len(shared_cells)))
    
    shared_cell_prot_ind = [np.where(prot_cells==cell)[0][0] for cell in shared_cells]
    shared_cell_rna_ind = [np.where(rna_cells==cell)[0][0] for cell in shared_cells]
    
    shared_cells_bool = np.zeros(len(rna_cells), dtype=bool)
    shared_cells_bool[shared_cell_rna_ind] = True
    
    for col_att in vlm.ca:
        vlm.ca[col_att] = vlm.ca[col_att][shared_cell_rna_ind]
    vlm.S = vlm.S[:,shared_cell_rna_ind]
    vlm.U = vlm.U[:,shared_cell_rna_ind]
    prot_count_array = prot_count_array[:,shared_cell_prot_ind].astype('float')
    prot_cells = np.asarray(shared_cells) #for idempotence
    
    vlm.P = prot_count_array
    return [prot_count_array, shared_cells, prot_cells]

def fit_pcs(vlm, space_name, pc_name, n_pcs):
    """
    Calculates PCA for top n_pcs components using used-selected space "space_name". Stores as "pc_name". 
    Stores the linear transformation (pc_name_fit) used to generate the PCA space.
    """
    pca = PCA(n_components=n_pcs)
    space = getattr(vlm,space_name)
    pca_fit = pca.fit(space.T)
    setattr(vlm,pc_name,pca.transform(space.T))
    setattr(vlm,pc_name+'_fit',pca_fit)

def fit_tsne(vlm, pc_space_name, ts_name, n_pcs):
    bh_tsne = TSNE()
    setattr(vlm,ts_name,bh_tsne.fit_transform(getattr(vlm,pc_space_name)[:, :n_pcs]))

    
def impute(vlm, P, k=400, impute_in_prot_space=True, size_norm=False, impute_in_pca_space=False, n_pcs = 10):
    """
    Nearest-neighbor imputation function that smoothes using a kNN connection graph. The connection basis may be either
    protein space P or spliced molecule space S. 
    If size-normalization is used, imputation is performed on log2-normalized counts rather than log counts.
    If imputation in PCA space is used, PCA-transformed space is used to construct kNN connection graph.
    """
    P = P.astype('float')
    if impute_in_prot_space or size_norm:
        P_cellsize = P.sum(0)
        P_sz = P / P_cellsize * np.median(P_cellsize)
        P_norm = np.log2(P_sz+1)
        vlm.P_norm = P_norm

    if impute_in_prot_space:
        if impute_in_pca_space:
            if ~hasattr(vlm,prot_pcs):
                fit_pcs(vlm, 'P_norm', 'prot_pcs', n_pcs)
            space = vlm.prot_pcs
        else:
            space = P_norm.T
        

        nn=NearestNeighbors(n_neighbors=k, n_jobs=8, )
        nn.fit(space)
        knn=nn.kneighbors_graph(X=None, mode='distance')
        conn = (knn>0).astype(float)
        vlm.connectivity=conn
        smooth = conn.multiply(1./sparse.csr_matrix.sum(conn, axis=1))
        
        
        if size_norm:
            vlm.Px = sparse.csr_matrix.dot(P_norm, smooth.T)
            vlm.knn_imputation_precomputed(smooth)
        else:
            vlm.Px = sparse.csr_matrix.dot(P, smooth.T)
            vlm.Sx = sparse.csr_matrix.dot(vlm.S, smooth.T)
            vlm.Ux = sparse.csr_matrix.dot(vlm.U, smooth.T)
            
    else:
        vlm.normalize(which="both",size=size_norm)
        if impute_in_pca_space and not hasattr(vlm,'pcs'):
            fit_pcs(vlm, 'S_norm', 'pcs', n_pcs)
     
        vlm.knn_imputation(k=k,size_norm=size_norm,pca_space=impute_in_pca_space)
        if size_norm:
            vlm.Px = sparse.csr_matrix.dot(P_norm, vlm.knn_smoothing_w.T)
        else:
            vlm.Px = sparse.csr_matrix.dot(P, vlm.knn_smoothing_w.T)
            
def identify_clusters(vlm,conn,correct_tags = False, tag_correction_list = [], method_name='ModularityVertexPartition'):
    """
    Cluster identification via the Louvain algorithm. Can be used for cluster discovery. If clusters are manually identified
    (e.g. by visualize_protein_markers()), clusters can be renumbered or combined using the tag correction list.
    Method names are any used in louvain.find_partition method.
    """
    g=ig.Graph.Adjacency(conn.todense().tolist())
    method=getattr(louvain,method_name)
    partition=louvain.find_partition(g,method)
    tag_list = np.zeros(conn.shape[0])
    for x in range(len(partition)):
        tag_list[partition[x]]=int(x)
    if correct_tags:
        cluster_ID = [tag_correction_list[int(X)] for X in tag_list]
    else:
        cluster_ID = [int(X) for X in tag_list]
    

    num_clusters = max(cluster_ID)+1
    
    vlm.cluster_ID = cluster_ID
    vlm.num_clusters = int(num_clusters)
    return [cluster_ID, num_clusters]

def visualize_pcs(vlm, pc_targets, pc_space='prot_pcs', write_labels=True):
    """
    Visualizes principal components of a PC space. pc_targets (list of principal components to plot) is not zero-indexed.
    Presupposes that clusters have already been computed by identify_clusters(), or manually input. 
    write_labels adds a legend with labels, if they have been stored as a vlm attribute. 
    """
    pc_zi = [pc_targets[0]-1, pc_targets[1]-1]
    
 #   print(dir(vlm))

    cluster_ID = vlm.cluster_ID
    num_clusters = vlm.num_clusters
    if hasattr(vlm,'COLORS'):
        color = vlm.COLORS[cluster_ID]
    else: 
        color = 'k'
    
    
    pcs = getattr(vlm, pc_space)
    plt.figure(figsize=(8,6))
    plt.scatter(pcs[:,pc_zi[0]],pcs[:,pc_zi[1]],s=3,c=color,alpha=0.9)
    plt.xlabel('PC'+str(pc_targets[0]))
    plt.ylabel('PC'+str(pc_targets[1]))
    yl = plt.ylim()
    xl = plt.xlim()

    if hasattr(vlm,'COLORS') and hasattr(vlm,'labels') and write_labels:
        for i in range(num_clusters):
            col=np.reshape(vlm.COLORS[i], (-1, 3))
            plt.scatter(100,100,s=0.1,c=col,alpha=1,label=vlm.labels[i])
        plt.ylim(yl)
        plt.xlim(xl)
        plt.legend(markerscale=20)
    plt.axis('off')
    
def visualize_protein_markers(vlm, protein_markers, pc_targets, visualize_clusters=False,
                              colormap='inferno'):
    """
    Visualizes an array of protein markers in protein principal component space. The components to plot are given by
    the list pc_targets. If visualize_clusters is selected, an additional cluster-colored plot is generated.
    Useful for iterative manual procedure to identify clusters based on characteristic markers.
    """
    array_proteins = vlm.adt_names
    pcs = vlm.prot_pcs
    pc_zi = [pc_targets[0]-1, pc_targets[1]-1]
    
    n_addit = int(visualize_clusters)
    
    nrows = int(np.ceil((len(protein_markers)+n_addit)/5))
#     print(nrows)
    
    f, ax = plt.subplots(nrows=nrows,ncols=5,figsize=(12,0.25+2*nrows))
    ax = ax.flatten()


    for j in range(len(protein_markers)):
        prot_name=protein_markers[j]
        ax[j].scatter(pcs[:,pc_zi[0]],pcs[:,pc_zi[1]],s=3,c=np.log(vlm.P[array_proteins == prot_name][0]+1),
                      alpha=0.2,cmap=colormap)
        ax[j].set_title(prot_name)
    
    if visualize_clusters:
        
        if hasattr(vlm,'cluster_ID') and hasattr(vlm,'COLORS'):
            col=vlm.COLORS[vlm.cluster_ID]
        else:
            COLORS = np.rand(np.amax(cluster_ID)+1,3)
            col=COLORS[vlm.cluster_ID]
            
                    
        ax[-1].scatter(pcs[:,pc_zi[0]],pcs[:,pc_zi[1]],s=3,c=col,alpha=0.9)

        
    for k in range(len(ax)):
        ax[k].axis('off')

        
def visualize_phase_portraits(vlm, markers, target='protein', imputed=True, prot_dict=False, plot_fit=False):
    """
    Plots imputed or raw phase portraits for protein or RNA velocity using a panel of genes (markers). If protein
    velocity plots are desired, a dictionary (prot_dict) must be passed into the function to define the relationship 
    between each gene and the protein markers.
    plot_fit plots the linear fit for RNA velocity.
    """
    nrows = int(np.ceil((len(markers))/6))
    
    f, ax = plt.subplots(nrows=nrows,ncols=6,figsize=(16,0.25+2.5*nrows))
    ax = ax.flatten()

    present_gene_list = vlm.ra['Gene']
    present_protein_list = vlm.adt_names

    if hasattr(vlm,'cluster_ID') and hasattr(vlm,'COLORS'):
        col=vlm.COLORS[vlm.cluster_ID]
    else:
        col='k'
                    
    for gene_ind in range(len(markers)):
        gene_name = markers[gene_ind]
        if target == 'mrna':
            if gene_name in present_gene_list:
                if imputed:
                    attr=['Sx','Ux']
                else:
                    attr=['S','U']
                    
                gene_filt = present_gene_list == gene_name
                    
                ax[gene_ind].scatter(getattr(vlm,attr[0])[gene_filt][0],
                                     getattr(vlm,attr[1])[gene_filt][0],s=1,c=col,alpha=0.5)
                ax[gene_ind].set_title(gene_name,fontsize=12)
                ax[gene_ind].set_xlabel('s')
                ax[gene_ind].set_ylabel('u')
                
                if hasattr(vlm,'gammas') and hasattr(vlm,'q') and plot_fit==True:
                    ax_rang = np.asarray([np.amin(getattr(vlm,attr[0])[gene_filt][0]), 
                               np.amax(getattr(vlm,attr[0])[gene_filt][0])])
                    fit_line = vlm.gammas[gene_filt] * ax_rang + vlm.q[gene_filt]
                    ax[gene_ind].plot(ax_rang,fit_line,'r',linewidth=2)
            else:
                ax[gene_ind].set_title('No '+gene_name+'!',fontsize=12)
        if target == 'protein':
            gene_name = markers[gene_ind]
            if gene_name in prot_dict.keys():
                protein_name = prot_dict[gene_name]
            else:
                protein_name = '?'
            
            if gene_name != protein_name:
                titlestr = gene_name+'/'+protein_name
            else:
                titlestr = gene_name
                    
            if gene_name in present_gene_list and protein_name in present_protein_list:

                
                if imputed: 
                    attr=['Px','Sx']
                else:
                    attr=['P','S']
                
                gene_filt = present_gene_list == gene_name
                prot_filt = present_protein_list == protein_name
                
                ax[gene_ind].scatter(getattr(vlm,attr[0])[prot_filt][0],
                                     getattr(vlm,attr[1])[gene_filt][0],s=1,c=col,alpha=0.5)
                ax[gene_ind].set_title(titlestr,fontsize=12)
                ax[gene_ind].set_xlabel('p')
                ax[gene_ind].set_ylabel('s')
                    
            else:
                ax[gene_ind].set_title('No '+titlestr+'!',fontsize=12)
    for ax_ind in range(gene_ind+1,len(ax)):
        ax[ax_ind].axis('off')
    f.tight_layout()

    
def gamma_fit(vlm, x, y, vel_type, genes_used_for_prot_velocity=False, adt_used_for_prot_velocity=False, maxmin_perc=[2, 98]):
    """
    Identifies steady-state genes from imputed data. Specific data to use are defined by names x and y.
    For RNA velocity, use x='Sx', y='Ux'; for protein velocity, use x='Px',y='Sx'. 
    Protein velocity requires explicit panels of gene and ADT names defined in the corresponding variables.
    These panels and their corresponding data are stored as vlm attributes.
    This code is a truncated version of the default workflow in velocyto 0.17.
    """

    Xarr = getattr(vlm,x)
    Yarr = getattr(vlm,y)
        
    if vel_type == 'protein':
        present_gene_list = vlm.ra['Gene']
        present_protein_list = vlm.adt_names
        
        prot_velo_gene_ind = [np.where(present_gene_list==gene)[0][0] for gene in genes_used_for_prot_velocity]
        prot_velo_prot_ind = [np.where(present_protein_list==prot)[0][0] for prot in adt_used_for_prot_velocity]
        
        Xarr = Xarr[prot_velo_prot_ind,:]
        Yarr = Yarr[prot_velo_gene_ind,:]
    
    denom_X = np.percentile(Xarr, 99.9, 1)
    if np.sum(denom_X == 0):
        denom_X[denom_X == 0] = np.maximum(np.max(Xarr[denom_X == 0, :], 1), 0.001)
        
    denom_Y = np.percentile(Yarr, 99.9, 1)
    if np.sum(denom_Y == 0):
        denom_Y[denom_Y == 0] = np.maximum(np.max(Yarr[denom_Y == 0, :], 1), 0.001)
        
    X_maxnorm = Xarr / denom_X[:, None]
    Y_maxnorm = Yarr / denom_Y[:, None]
    
    normsum = X_maxnorm + Y_maxnorm
    down, up = np.percentile(normsum, maxmin_perc, axis=1)
    W = ((normsum <= down[:, None]) | (normsum >= up[:, None])).astype(float)

    [gammas, q, R2] = fit_slope_weighted_offset(Yarr, Xarr, W, return_R2=True, limit_gamma=False)    
    
    if vel_type=='rna':
        vlm.gammas = gammas
        vlm.q = q
    if vel_type=='protein':
        vlm.gammas_p = gammas
        vlm.q_p = q
        
        vlm.Px_prot = Xarr
        vlm.Sx_prot = Yarr
        vlm.genes_used_for_prot_velocity = genes_used_for_prot_velocity
        vlm.adt_used_for_prot_velocity   = adt_used_for_prot_velocity

def extrapolate(vlm, vel_type='rna'):
    """
    Determines direction of movement in spliced RNA or protein space.
    """
    if vel_type == 'rna':
        vlm.delta_S = vlm.Ux - (vlm.Sx * vlm.gammas[:, None] + vlm.q[:,None])
    if vel_type == 'protein':
        vlm.delta_P = vlm.Sx_prot - (vlm.Px_prot * vlm.gammas_p[:, None] + vlm.q_p[:,None])
        
def linear_velocity_projection(vlm, vel_type='rna'):
    """
    Plots linear projection onto principal components. So far only implemented for RNA velocity,
    because the scRNA-seq workflow overwhelmingly uses the spliced RNA embedding.
    """
    if vel_type == 'rna':
        vlm.rna_lin_proj = vlm.pcs_fit.transform(vlm.delta_S.T)
        vlm.embedding = vlm.pcs

def identify_embedding_knn(vlm,embedding_name, embedding_dimensions, n_neighbors=500):
    """
    Identifies nearest neighbors in the embedding's specified dimensions. 
    Both RNA and protein velocity extrapolation hypothesize that each cell is likely to 
    transition to cells neighboring it in the embedding space. 
    """
    vlm.embedding = getattr(vlm,embedding_name)[:,embedding_dimensions]
    nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=8)
    nn.fit(vlm.embedding)
    vlm.embedding_knn = nn.kneighbors_graph(mode="connectivity")        
    
def calculate_embedding_delta(vlm, high_dim_space_name, delta_high_dim_name, delta_name):
    """
    Identifies the net embedding direction, given the assumption of transitioning to neighbors. 
    Largely adapts and truncates the procedure from velocyto 0.17 using a particular set of parameters.
    """
    high_dim = getattr(vlm,high_dim_space_name)
    delta_high_dim = getattr(vlm,delta_high_dim_name)

    psc = 1e-10
    corrcoef=colDeltaCorSqrt(high_dim,np.sqrt(np.abs(delta_high_dim) + psc) * np.sign(delta_high_dim), 
                                            threads=8, psc=psc)

    np.fill_diagonal(corrcoef, 0)

    sigma_corr = 0.05
    transition_prob = np.exp(corrcoef / sigma_corr) * vlm.embedding_knn.A  
    transition_prob /= transition_prob.sum(1)[:, None]


    unitary_vectors = vlm.embedding.T[:, None, :] - vlm.embedding.T[:, :, None]
    with np.errstate(divide='ignore', invalid='ignore'):
        unitary_vectors /= np.linalg.norm(unitary_vectors, ord=2, axis=0)  # divide by L2
        for j in range(unitary_vectors.shape[0]):
            np.fill_diagonal(unitary_vectors[j, ...], 0)  # fix nans

    delta_embedding = (transition_prob * unitary_vectors).sum(2)
    delta_embedding -= (vlm.embedding_knn.A * unitary_vectors).sum(2) / vlm.embedding_knn.sum(1).A.T
    delta_embedding = delta_embedding.T
    
    setattr(vlm,delta_name, delta_embedding)
    
def cluster_specific_plot(vlm, delta_name, draw_cells=False):
    """
    Plots cell-specific velocities for the set of clusters identified in the data. Plots clusters on individual axes
    to facilitate visualization of cluster-specific dynamics.
    """
    nrows = int(np.ceil((vlm.num_clusters)/2))
    
    cluster_labels = vlm.labels
    f, ax = plt.subplots(nrows=nrows,ncols=2,figsize=(16,1+8*nrows))
    ax = ax.flatten()

    for j in range(vlm.num_clusters):
        plt.sca(ax[j])
        subset_true = np.asarray(vlm.cluster_ID) == j
        subset_false = ~ subset_true

        visualize_velocity_projection(vlm, delta_name, subset=subset_false, use_subset=True, use_color=False, draw_cells=draw_cells, alpha=0.2)

        visualize_velocity_projection(vlm, delta_name, subset=subset_true, use_subset=True, use_color=True, draw_cells=draw_cells)
        
        ax[j].set_title(cluster_labels[j],fontsize=12)
        
    for k in range(len(ax)):
        ax[k].axis('off')
        
def visualize_velocity_projection(vlm, delta_name, subset=[], use_subset=False, use_color=True, draw_cells=True, alpha=1):
    """
    Plots cell-specific velocities in the embedding for a particular vector direction, defined by delta_name.
    If draw_cells is chosen, a point is placed at the root of each arrow.
    If use_color is chosen, the arrows/cells are plotted with cluster-specific colors.
    """
    if hasattr(vlm,'COLORS'):
        color = vlm.COLORS[vlm.cluster_ID]
    if not hasattr(vlm,'COLORS') or use_color == False:
        color='k'

    if use_subset == False:
        embedding = vlm.embedding
        vector = getattr(vlm,delta_name)
    else:
        embedding = vlm.embedding[subset,:]
        vector = getattr(vlm,delta_name)[subset,:]
        if len(color) > 1:
            color = color[subset]
    if draw_cells:
        plt.scatter(embedding[:,0],embedding[:,1],s=3,c=color,alpha=alpha)
            
    plt.quiver(embedding[:,0],
               embedding[:,1],
               vector[:,0],
               vector[:,1], color=color,pivot='tail',minlength=0.1,alpha=alpha)
    
    plt.axis('off')

def initialize_grid_embedding(vlm, steps=(20,20,20), n_neighbors=200):

    """
    Intializes the grid for and Gaussian kernel for downstream velocity aggregation. 
    Adapted from velocyto 0.17.
    """
    smooth=0.5
    grs = []
    for dim_i in range(vlm.embedding.shape[1]):
        m, M = np.min(vlm.embedding[:, dim_i]), np.max(vlm.embedding[:, dim_i])
        m = m - 0.025 * np.abs(M - m)
        M = M + 0.025 * np.abs(M - m)
        gr = np.linspace(m, M, steps[dim_i])
        grs.append(gr)
    meshes_tuple = np.meshgrid(*grs)
    gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T

    vlm.XY = gridpoints_coordinates

    from sklearn.neighbors import NearestNeighbors
    from scipy.stats import norm as normal

    
    nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=8)
    nn.fit(vlm.embedding)
    dists, grid_neighbors = nn.kneighbors(gridpoints_coordinates)

    std = np.mean([(g[1] - g[0]) for g in grs])
    # isotropic gaussian kernel
    gaussian_w = normal.pdf(loc=0, scale=smooth * std, x=dists)
    
    vlm.total_p_mass = gaussian_w.sum(1)
    vlm.gaussian_w = gaussian_w
    vlm.grid_neighbors = grid_neighbors
    
def calculate_grid_arrows(vlm, delta_name, UV_suffix, min_mass = 1, uv_multiplier=1):
    """
    Given a pre-calculated grid, calculates the aggregated velocity direction at the grid points. 
    Throws away arrows with an insufficient probability mass
    """
    
    delta_embedding = getattr(vlm,delta_name)

    UV = (delta_embedding[vlm.grid_neighbors] * vlm.gaussian_w[:, :, None]).sum(1) / np.maximum(1, vlm.total_p_mass)[:, None]

    mass_filter = vlm.total_p_mass < min_mass
    UV[mass_filter,:] = 0
    UV = UV * uv_multiplier
    
    if hasattr(vlm,'mass_filter'):
        vlm.mass_filter = mass_filter & vlm.mass_filter
    else:
        vlm.mass_filter = mass_filter
    setattr(vlm,'UV'+UV_suffix,UV)
    
def plot_grid_arrows(vlm, UV_name, plot_cells=False, arr_col='k', color_cells_by_cluster=False, pivot='tail',
                    plot_grid_points = True, arr_scale = 0.6, cell_alpha=0.8):
    """
    Plots a particular precalculated set of grid arrows, denoted by UV_name. 
    If plot_cells is chosen, cells are plotted at their emebdding locations.
    If color_cells_by_cluster is chosen, they are colored by cluster.
    The option "pivot" should be set to "tail" for RNA velocity (forward extrapolation), and
    "head" for protein velocity (backward extrapolation).
    """
    if hasattr(vlm,'cluster_ID') and hasattr(vlm,'COLORS') and color_cells_by_cluster:
        col=vlm.COLORS[vlm.cluster_ID]
    else:
        col='k'
        
        
    if plot_cells:
        plt.scatter(vlm.embedding[:,0],vlm.embedding[:,1],s=30,c=col,alpha=cell_alpha)
        
    _quiver_kwargs = {"angles": 'xy', "scale_units": 'xy', "minlength": 0, "pivot": pivot, "color": arr_col,
                     "headaxislength": 2.75, "headlength": 5, "headwidth": 4.8}
    UV = getattr(vlm,UV_name)
    plt.quiver(vlm.XY[:, 0], vlm.XY[:, 1], UV[:, 0], UV[:, 1], scale=arr_scale, **_quiver_kwargs)
    
    if plot_grid_points:
        plt.scatter(vlm.XY[~vlm.mass_filter, 0], vlm.XY[~vlm.mass_filter, 1],c='k',s=10,zorder=30000)
    plt.axis('off')
    
def plot_bezier(vlm, plot_cells=False, color_cells_by_cluster=False, arr_len_scal=0.25,cell_alpha=0.8):
    """
    Plots aggregated Bezier curves assumed to represent underlying embedding curvature based on RNA and protein velocity grid arrows.
    """        
    if hasattr(vlm,'cluster_ID') and hasattr(vlm,'COLORS') and color_cells_by_cluster:
        col=vlm.COLORS[vlm.cluster_ID]
    else:
        col='k'
    
    if plot_cells:
        plt.scatter(vlm.embedding[:,0],vlm.embedding[:,1],s=30,c=col,alpha=cell_alpha)

        
    n_curves = vlm.XY.shape[0]


    s_vals = np.linspace(0.0, 1.0, 15)



    XYM=vlm.XY
    UVP=vlm.UV_prot
    UVM=vlm.UV_rna
    
    ARR_SCAL = 100;
    for i in range(n_curves):
        nodes = np.asfortranarray([[XYM[i, 0]-UVP[i, 0], XYM[i, 0], XYM[i, 0]+UVM[i, 0]],
                 [XYM[i, 1]-UVP[i, 1], XYM[i, 1], XYM[i, 1]+UVM[i, 1]]])
        curve = bezier.Curve(nodes, degree=1)
        a = curve.evaluate_multi(s_vals)
        arr_len = [a[0,-1]-a[0,-2],a[1,-1]-a[1,-2]]    
        arr_vec = np.asarray(arr_len)
        ARRLENSCAL = np.linalg.norm(arr_vec, ord=2, axis=0)
        if ARRLENSCAL>0.0 and not vlm.mass_filter[i]:
            arr_vec /= ARRLENSCAL
            arr_vec_perp = np.asarray( [arr_vec[1],-arr_vec[0]])
            LEN = arr_len_scal
            WID = LEN*0.7
            BK = arr_len_scal*0.1/0.25
            CC = a[:,-1]
            C_off = LEN*0.5
            CC -= arr_vec*C_off

            patch_x = np.asarray([CC[0], 
                                  CC[0] - WID * arr_vec_perp[0] - arr_vec[0]*BK, 
                                  CC[0] + LEN * arr_vec[0],
                                  CC[0] + WID * arr_vec_perp[0] - arr_vec[0]*BK,
                                  CC[0]])
            patch_y = np.asarray([CC[1], 
                                  CC[1] - WID * arr_vec_perp[1] - arr_vec[1]*BK, 
                                  CC[1] + LEN * arr_vec[1],
                                  CC[1] + WID * arr_vec_perp[1] - arr_vec[1]*BK,
                                  CC[1]])
            plt.fill(patch_x,patch_y,'k',alpha=1)


            plt.plot(a[0],a[1],'k',linewidth=3)
    _=plt.axis('off')
