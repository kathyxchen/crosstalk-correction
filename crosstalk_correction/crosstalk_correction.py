"""Donato et al.'s (2013) maximum impact estimation algorithm.
"""
import numpy as np

NEAR_ZERO = 1e-18


def crosstalk_correction(pathway_definitions,
                         random_seed=2015,
                         gene_set=set(),
                         all_genes=True,
                         max_iters=1000):
    """A wrapper function around the maximum impact estimation algorithm.

    Parameters
    -----------
    pathway_definitions : dict(str -> set(str))
        The original pathway definitions.
        A pathway (key) is defined by a set of genes (value).
    random_seed : int (default=2015)
        Sets the numpy random seed
    gene_set : set(str) (default=set())
        Donato et al. (2013) uses this algorithm to remove crosstalk from
        definitions in case-control studies. Here, `gene_set` is equivalent
        to the DE (differentially expressed) genes they refer to in their
        paper. Because crosstalk removal is a preprocessing step
        applicable to pathway analyses in general, we keep the variable
        name nonspecific.
    all_genes : bool (default=True)
        This value is checked if `gene_set` is not empty.
        If False, crosstalk correction is only applied to annotations in
        the `gene_set`.
    max_iters : int (default=1000)
        The maximum number of expectation-maximization steps to take in the
        maximum impact estimation algorithm.

    Returns
    ----------
    dict(str -> tup(set(str), set(str))), where the (str) keys are the
      pathway names.
      tup[0] : crosstalk-correction applied to genes in the pathway
               definition that are also in `gene_set`
      tup[1] :
      - `all_genes` is True.
        Correction is applied to genes outside of `gene_set`.
      - `all_genes` is False.
         The second element in the tuple is all genes
         remaining in the original definition (definition - `gene_set`).
    """
    np.random.seed(seed=random_seed)

    genes_in_pathway_definitions = set.union(*pathway_definitions.values())
    pathway_column_names = index_element_map(pathway_definitions.keys())

    corrected_pathway_defns = {}
    if gene_set:
        gene_set = gene_set & genes_in_pathway_definitions
        if not gene_set and not all_genes:
            print("`gene_set` parameter was {0}, returning original"
                  "pathway definitions".format(gene_set))
            for pathway, definition in pathway_definitions.items():
                corrected_pathway_defns[pathway] = (set(), definition)
            return pathway_definitions

        corrected_pathway_defns = _apply_correction_on_genes(
            gene_set, pathway_column_names, pathway_definitions)

        # crosstalk correction is _only_ applied to `gene_set`
        if not all_genes:
            for pathway, definition in pathway_definitions.items():
                if pathway not in corrected_pathway_defns:
                    corrected_pathway_defns[pathway] = set()
                gene_set_defn = corrected_pathway_defns[pathway]
                remaining_defn = definition - gene_set
                corrected_pathway_defns[pathway] = (
                    gene_set_defn, remaining_defn)
            return corrected_pathway_defns

    remaining_genes = genes_in_pathway_definitions - gene_set
    if not remaining_genes:
        for pathway, definition in corrected_pathway_defns.items():
            corrected_pathway_defns[pathway] = (definition, set())
        return corrected_pathway_defns

    pathway_remaining_defns = _apply_correction_on_genes(
        remaining_genes, pathway_column_names, pathway_definitions)

    for pathway, definitions in pathway_definitions.items():
        if pathway not in corrected_pathway_defns:
            corrected_pathway_defns[pathway] = set()

        if pathway not in pathway_remaining_defns:
            pathway_remaining_defns[pathway] = set()

        corrected_pathway_defns[pathway] = (
            corrected_pathway_defns[pathway],
            pathway_remaining_defns[pathway])
    return corrected_pathway_defns


def maximum_impact_estimation(membership_matrix, max_iters=1000):
    """An expectation maximization technique that produces pathway definitions
    devoid of crosstalk. That is, each gene is mapped to the pathway in
    which it has the greatest predicted impact; this removes any overlap
    between pathway definitions.

    Parameters
    -----------
    membership_matrix : numpy.array(float), shape = [n, k]
        The observed gene-to-pathway membership matrix, where n is the number
        of genes and k is the number of pathways we are interested in.
    max_iters : int (default=1000)
        The maximum number of expectation-maximization steps to take.

    Returns
    -----------
    dict(int -> set(int)), a dictionary mapping a pathway to a set of genes.
    These are the pathway definitions after the maximum impact estimation
    procedure has been applied to remove crosstalk.
      - The keys are ints corresponding to the pathway column indices in the
        membership matrix.
      - The values are sets of ints corresponding to gene row indices in the
        membership matrix.
    """
    # Initialize the probability vector as the sum of each column in the
    # membership matrix normalized by the sum of the entire membership matrix.
    # The probability at some index j in the vector represents the likelihood
    # that a pathway (column) j is defined by the current set of genes (rows)
    # in the membership matrix.
    pr_0 = np.sum(membership_matrix, axis=0) / np.sum(membership_matrix)
    pr_1 = _update_probabilities(pr_0, membership_matrix)
    epsilon = np.linalg.norm(pr_1 - pr_0)/100.

    pr_old = pr_1
    check_for_convergence = epsilon
    count = 0
    while epsilon > NEAR_ZERO and check_for_convergence >= epsilon:
        count += 1
        if count > max_iters:
            print("Reached the maximum number of iterations {0}".format(
                  max_iters))
            break
        pr_new = _update_probabilities(pr_old, membership_matrix)
        check_for_convergence = np.linalg.norm(pr_new - pr_old)
        pr_old = pr_new

    pr_final = pr_old  # renaming for readability

    corrected_pathway_definitions = {}
    n, k = membership_matrix.shape
    for gene_index in range(n):
        gene_membership = membership_matrix[gene_index]
        denominator = np.dot(gene_membership, pr_final)
        # Approximation is used to prevent divide by zero warning.
        # Since we are only looking for the _most_ probable pathway in which a
        # gene contributes its maximum impact, precision is not as important
        # as maintaining the relative differences between each
        # pathway's probability.
        if denominator < NEAR_ZERO:
            denominator = NEAR_ZERO
        # This is equivalent to one row in what Donato et al. (2013) refer
        # to as the underlying (latent) Z matrix.
        conditional_pathway_pr = (np.multiply(gene_membership, pr_final) /
                                  denominator)
        all_pathways_at_max = np.where(
            conditional_pathway_pr == conditional_pathway_pr.max())[0]
        gene_in_pathways = np.where(gene_membership == 1)[0]
        all_pathways_at_max = np.intersect1d(
            all_pathways_at_max, gene_in_pathways)
        pathway_index = np.random.choice(all_pathways_at_max)
        if pathway_index not in corrected_pathway_definitions:
            corrected_pathway_definitions[pathway_index] = set()
        corrected_pathway_definitions[pathway_index].add(gene_index)
    return corrected_pathway_definitions


def initialize_membership_matrix(gene_row_names, pathway_definitions):
    """Create the binary gene-to-pathway membership matrix that
    will be considered in the maximum impact estimation procedure.

    Parameters
    -----------
    gene_row_names : set(str)
        The genes for which we want to assess pathway membership
    pathway_definitions : dict(str -> set(str))
        Pathway definitions, pre-crosstalk-removal.
        A pathway (key) is defined by a set of genes (value).

    Returns
    -----------
    numpy.array, shape = [n, k], the membership matrix
    """
    membership = []
    for pathway, full_definition in pathway_definitions.items():
        pathway_genes = list(full_definition & gene_row_names)
        membership.append(np.in1d(list(gene_row_names), pathway_genes))
    membership = np.array(membership).astype("float").T
    return membership


def index_element_map(arr):
    """Map the indices of the array to the respective elements.

    Parameters
    -----------
    arr : list(a)
        The array to process, of generic type a

    Returns
    -----------
    dict(int -> a), a dictionary corresponding the index to the element
    """
    index_to_element = {}
    for index, element in enumerate(arr):
        index_to_element[index] = element
    return index_to_element


def _apply_correction_on_genes(genes,
                               pathway_column_names,
                               pathway_definitions):
    """Helper function to create the gene-to-pathway
    membership matrix and apply crosstalk correction on that
    matrix. Returns the crosstalk-corrected pathway definitions
    for the input `genes.`
    """
    gene_row_names = index_element_map(genes)
    membership_matrix = initialize_membership_matrix(
        genes, pathway_definitions)
    crosstalk_corrected_index_map = maximum_impact_estimation(
        membership_matrix)

    updated_pathway_definitions = _update_pathway_definitions(
        crosstalk_corrected_index_map,
        gene_row_names, pathway_column_names)
    return updated_pathway_definitions


def _update_pathway_definitions(crosstalk_corrected_index_map,
                                gene_row_names,
                                pathway_column_names):
    """Helper function to convert the mapping of int
    (pathway id -> list of gene ids) to the corresponding pathway
    names and gene identifiers.
    """
    corrected_pathway_definitions = {}
    for pathway_index, gene_indices in crosstalk_corrected_index_map.items():
        pathway = pathway_column_names[pathway_index]
        genes = set([gene_row_names[index] for index in list(gene_indices)])
        corrected_pathway_definitions[pathway] = genes
    return corrected_pathway_definitions


def _update_probabilities(pr, membership_matrix):
    """Updates the probability vector for each iteration of the
    expectation maximum algorithm in maximum impact estimation.

    Parameters
    -----------
    pr : numpy.array(float), shape = [k]
        The current vector of probabilities. An element at index j,
        where j is between 0 and k - 1, corresponds to the probability that,
        given a gene g_i, g_i has the greatest impact in pathway j.
    membership_matrix : numpy.array(float), shape = [n, k]
        The observed gene-to-pathway membership matrix, where n is the number
        of genes and k is the number of pathways we are interested in.

    Returns
    -----------
    numpy.array(float), shape = [k], a vector of updated probabilities
    """
    n, k = membership_matrix.shape
    pathway_col_sums = np.sum(membership_matrix, axis=0)

    weighted_pathway_col_sums = np.multiply(pathway_col_sums, pr)
    sum_of_col_sums = np.sum(weighted_pathway_col_sums)
    try:
        new_pr = weighted_pathway_col_sums / sum_of_col_sums
    except FloatingPointError:
        # In the event that we encounter underflow or overflow issues,
        # apply this approximation.
        cutoff = 1e-150 / k
        log_cutoff = np.log(cutoff)

        weighted_pathway_col_sums = _replace_zeros(
            weighted_pathway_col_sums, cutoff)
        log_weighted_col_sums = np.log(weighted_pathway_col_sums)
        log_weighted_col_sums -= np.max(log_weighted_col_sums)

        below_cutoff = log_weighted_col_sums < log_cutoff
        geq_cutoff = log_weighted_col_sums >= log_cutoff

        print("{1} adjustments made to a vector of length {0}"
              " containing the raw weight values"
              " in a call to 'update_probabilities'".format(
                  k, len(log_weighted_col_sums[below_cutoff])))

        new_pr = np.zeros(k)
        new_pr[below_cutoff] = cutoff
        col_sums_geq_cutoff = log_weighted_col_sums[geq_cutoff]
        new_pr[geq_cutoff] = np.exp(
            col_sums_geq_cutoff) / np.sum(np.exp(sorted(col_sums_geq_cutoff)))

    difference = np.abs(1. - np.sum(new_pr))
    assert difference < 1e-12, "Probabilities sum to {0}.".format(
           np.sum(new_pr))
    return new_pr


def _replace_zeros(arr, default_min_value):
    """Substitute 0s in the list with a near-zero value.

    Parameters
    -----------
    arr : numpy.array(float)
    default_min_value : float
        If the smallest non-zero element in `arr` is greater than the default,
        use the default instead.

    Returns
    -----------
    numpy.array(float)
    """
    min_nonzero_value = min(default_min_value, np.min(arr[arr > 0]))
    closest_to_zero = np.nextafter(min_nonzero_value, min_nonzero_value - 1)
    arr[arr == 0] = closest_to_zero
    return arr
