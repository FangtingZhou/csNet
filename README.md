# csNet: Cell-specific Gene Regulatory Networks Construction via Latent Trajectory

Single-cell rRNA sequencing data provide unprecedented opportunities to decipher complex biological processes such as cell differentiation, embryonic development, disease evolution, and so on. Network construction and trajectory inference are among the most popular topics in single-cell data analysis. Current researches mainly concentrate on addressing the two aspects separately. Although innumerable state-of-the-art methods are proposed through the last few years to construct the overall gene regulatory relationship or cell linkage through pseudotime inference, few consider the problem of combining the two ideas and build the dynamic network through latent trajectory. In this project, we focus on construct cell-specific regulatory network through latent time trajectory. By borrowing information from gene network and pseudotime, and perform inference through rigorous Bayesian statistical modeling, we expect a better performance than individual analysis and a detailed uncertainty measurement than algorithmic methods.

### Model Assumption
We assume the simultaneous equation with factor and varying coefficient model which represents an directed graph. The sequencing data form a matrix where each element is the expression level of gene in cell. Parameters are the pseudotime we want to infer, the effect of gene, the effect of gene- gene ineraction which varies with time, and the random noise. In this setting, non-zero gene-gene coefficients constitute an directed graph, and the graph envolves with time. We model coefficient functions through b-splines which is flexible to capture the functional form of regulatory relationship. Furthermore, we propose to threshold the coefficient to achieve a dynamic network, that if the coefficient is below a threshold, we set it to zero.

However, in the current package, we only consider a fixed network case by assign a spike and slab prior on gene-gene coefficient. Future work is on efficient modeling and calculating graph parameters under the dynamic network case.

### Model Inference
The model is constructed under the Bayesian framework and inference is done with Markov chain Monte Carlo (MCMC). 

### Package Usage
We provide R code to generate simulation data and perform simulation study, and C code to do MCMC updates. In the future work, we hope to visualize the constructed dynamic network and pseudotime trajectory, and provide comparison with other methods.