from itertools import islice
import pandas as pd
import re


def over_slice(sequence, K):
    """
    Splices the string into iterative chunks of K sections. Returns a generator
    For example AGTTGA with K = 3 will become a generator with [AGT, GTT, TTG, TGA]
    """

    itr = iter(sequence)
    result = tuple(islice(itr, K))
    if len(result) == K:
        yield result
    for element in itr:
        result = result[1:] + (element,)
        yield result


def clone_kmer_show(K, fr_clone):
    """
    Gives a list of unique kmer sequences generated from splicing the clone into groups of K
    """
    pd_res = pd.DataFrame(["".join(kmer_sequence) for kmer_sequence in over_slice(fr_clone, K)])
    pd_res = pd_res.drop_duplicates()
    pd_res.rename(columns={0: 'Seq_kmer'}, inplace=True)
    return pd_res['Seq_kmer'].tolist()


def remove_extra_characters(clone):
    """
    Removes extra brackets and new lines
    """
    clone = re.sub('[ ]', "", clone)
    clone = re.sub('\n', '', clone)
    return clone


def clean_and_create_kmer(K, clone):
    """
    Cleans the original clone and then calls the clone_kmer_show method.
    """
    cleaned_clone = remove_extra_characters(clone)
    if len(cleaned_clone) >= K:
        return clone_kmer_show(K, cleaned_clone)
    else:
        raise Exception("Error: Your K length is longer than your clone length")


def create_clone_dictionary_with_k_splice(clone_dictionary, K_value):
    """
    Turns a clone dictionary's sequence into splices of k and assigns it to the same variable.
    """
    k_spliced_dictionary = {}
    for key in clone_dictionary:
        k_spliced_dictionary[key] = clean_and_create_kmer(K_value, clone_dictionary[key])
    return k_spliced_dictionary
