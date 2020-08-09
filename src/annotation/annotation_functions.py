from Bio import Entrez
Entrez.email = "nigelmic@umich.edu"

def get_gene_data(gene_symbol, organism, is_ENSEMBL=False, find_hsa_ortholog=False):
    """
    A function used to get gene:
    * functional description
    * human orthologs

    Parameters
    ---------- 
    gene_symbol: str
        A gene symbol, i.e. ACTB, Actb
    organism: str
        An organism name, i.e. musmusculus

    Returns
    -------
    dict of string values, with keys: "gene_description",
    "human_orthologues"

    Examples
    --------
    >>> annotations = get_gene_data("Actb")
    >>> print(annotations)
    {"gene_description": "Codes for protein that makes up cellular cytoskeleton",
    "human_orthologues": "ACTB"}
    """
    
    out = {}
    
    # Assign the NCBI standard organism IDs
    if organism == 'musmusclus' or organism == 'mmu':
        organism = '"Mus musculus"[porgn:__txid10090]'
    elif organism == 'homosapien' or organism == 'hsa':
        organism = '"Homo sapiens"[porgn:__txid9606]'
    else:
        return None #unsupported lookup
    
    # Compile NCBI query string
    query = '(' + gene_symbol + '[Gene Name]) AND ' + organism
    
    # Do query
    handle = Entrez.esearch(db='gene', term=query)
    record = Entrez.read(handle)
    handle.close()
    
    try: # the below lines will fail if the lookup failed
        ncbi_gene_id = record['IdList'][0]
    except:
        return None
        
    handle = Entrez.esummary(db="gene", id=ncbi_gene_id)
    record = Entrez.read(handle)
    handle.close()
            
    try: # the below lines will fail if the lookup failed
        out['gene_description'] = record['DocumentSummarySet']['DocumentSummary'][0]['Description']
        out['gene_summary'] = record['DocumentSummarySet']['DocumentSummary'][0]['Summary']
        out['gene_aliases'] = record['DocumentSummarySet']['DocumentSummary'][0]['OtherDesignations'].split( '|' )
        out['nucleotide_uid'] = record['DocumentSummarySet']['DocumentSummary'][0].attributes['uid']
    except:
        return None
        
    if not find_hsa_ortholog:
        return out
    
    # Run an ortholog search
    try:
        query = out['nucleotide_uid'] + '[Gene ID]'
        
        handle = Entrez.esearch(db='HomoloGene', term=query)
        record = Entrez.read(handle)
        handle.close()
        ortho_record = record['IdList'][0]

        handle = Entrez.esummary(db="HomoloGene", id=ortho_record)
        record = Entrez.read(handle)
        handle.close()
        
        out['homology'] = [(x['TaxName'],x['Symbol'],x['Title'],str(x['GeneID'])) for x in record[0]['HomoloGeneDataList']]
    except:
        pass # Replace with a `return out` if below changes
    
    return out