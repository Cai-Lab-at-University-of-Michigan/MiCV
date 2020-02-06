def custom_plot_multi_gene_trends(gene_trends, genes=None, line_colors=None):
    """ Plot the gene trends: groups of gene trends in the same panel
    :param: gene_trends: Results of the compute_marker_trends function
    """
    
    # Branches and genes
    branches = list(gene_trends.keys())
    if genes is None:
        genes = gene_trends[branches[0]]['trends'].index
    
    if (line_colors is None):
        colors = pd.Series(sns.color_palette('Set2', len(genes)).as_hex(), index=genes)
    elif len(line_colors) < len(branches):
        print("[ERROR] too few colors passed; need " + str(len(genes)) + " unique colors.")
    else:
        colors = pd.Series(sns.color_palette(line_colors[0:len(genes)]), index=genes)
        

    # Set up figure
    fig = plt.figure(figsize=[8, 6])
    for i, branch in enumerate(branches):
        ax = fig.add_subplot(len(branches), 1, i + 1)
        for gene in genes:
            trends = gene_trends[branch]['trends']
            stds = gene_trends[branch]['std']
            ax.plot(trends.columns, trends.loc[gene, :],
                    color=colors[gene], label=gene)
            ax.set_xticks([0, 1])
            ax.fill_between(trends.columns, trends.loc[gene, :] - stds.loc[gene, :],
                            trends.loc[gene, :] + stds.loc[gene, :], alpha=0.1, color=colors[gene])
            ax.set_title(branch)
        # Add legend
        if i == 0:
            ax.legend(ncol=2)
        plt.tight_layout()

    sns.despine()
d.palantir.plot.plot_multi_gene_trends = custom_plot_multi_gene_trends