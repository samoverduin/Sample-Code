#!/usr/bin/env python3
"""
Author: Sam Overduin

Implementation of the SSAHA algorithm:
SSAHA: a fast search method for large DNA databases
https://doi.org/10.1101/gr.194201
"""
# import statements
import re


def parse_fasta_generator(input_file):
    """Reads a fasta file and yields label, sequence as output

    Input:
    input_file = str location of fasta file
    Output:
    yields label, sequence as strings
    """
    with open(input_file) as file_obj:
        label = None
        for line in file_obj:
            if not line:
                continue
            if line.startswith('>'):
                if label: #make sure it's assigned
                    yield label, ''.join(seq_lst)
                label = line[1:].strip()
                seq_lst = []
            else:
                seq_lst.append(line.strip().upper())
    yield label, ''.join(seq_lst) #final one


def fasta_to_sequence_lst(input_file):
    """ Converts a fasta file to a list of sequences

    :param input_file: str path to fasta file in relative directory
    :return: list of sequences in fasta file
    """
    seq_lst = []
    for label, seq in parse_fasta_generator(input_file):
        seq_lst.append(seq)
    return seq_lst


def split_to_kmers(seq, k):
    """Takes in a sequence and k to make non-overlapping k-mers from sequence.

    :param seq: string of nucleotides (could be anything I guess)
    :param k: int of how long k-mers should be
    :return: list of k-mer strings
    If there is any more space at the end then the last almost-kmer is dropped.
    """
    kmer_lst = []
    for i in range(0, len(seq), k):
        # This if drops last faKe-mer
        if i+k <= len(seq):
            kmer_lst += [seq[i:i+k]]
    return kmer_lst


def make_hash_table(seq_lst, k, clean_sequences=False):
    """Creates a hash table (dict) containing list of (seq_nr, pos_in_seq)

    :param seq_lst: list of sequences to hash
    :param k: int of size for k-mers
    :param clean_sequences: Bool to replace non-AGCT chars with A

    :return: dict with k-mers as keys and as value a list of tuples:
            (seq_nr, pos_in_seq)
    """
    hash_dct = {}
    for seq_num, seq in enumerate(seq_lst):
        # clean seqs
        if clean_sequences:
            seq = seq.upper()
            nuc_scrubber = re.compile('[^ACGT]')
            if re.match(nuc_scrubber, seq):
                print('Non-AGCT char replaced in input sequence.')
                seq = re.sub(nuc_scrubber, 'A', seq)
        for i in range(0, len(seq), k):
            # This if drops final faKe-mer
            if i+k <= len(seq):
                kmer = seq[i:i+k]
                # setdefault looks for key kmer and if not found inserts []
                hash_dct.setdefault(kmer, []).append((seq_num, i))
    return hash_dct


def query_kmer_hits(query, hash_dct, k):
    """Returns a sorted list of hits in the hash table for a query seq.

    :param query: string of nucleotides
    :param hash_dct: dict with k-mers as keys and as values a lst of
                                                (seq_nr, pos_in_seq)

    :return: sorted list of hits in hash_dct
    """
    assert k == len(next(iter(hash_dct))), "hash_dct and function call have " \
                                           "different k-length"
    match_lst = []
    for i in range(0, len(query)-k+1):
        #print(query[i:i+k])
        # try-except in case the query k-mer is not in hash_dct
        try:
            for index, offset in hash_dct[query[i:i+k]]:
                # add index, shift and offset
                match_lst.append((index, offset - i, offset))
        except KeyError:
            pass
    # Sort list by first and second value of tuple.
    match_lst.sort(key=lambda tup: (tup[0], tup[1]))
    return match_lst


def mode_list(lst):
    return max(set(lst), key=lst.count)


def find_longest_match(match_lst, k):
    """Takes in match_lst as sorted list of matches and returns longest index
        & shift run

    :param match_lst: sorted list of query hits in hash_dct
    :param k: int of k-mer size to use

    :return: index, start offset, end offset, match_start
    """
    index_lst, shift_lst, offset_lst = zip(*match_lst)
    index_shift_lst = list(zip(index_lst, shift_lst))
    best_match = mode_list(index_shift_lst)

    index_match = best_match[0]
    offset_start = offset_lst[index_shift_lst.index(best_match)]
    # find last best_match in list
    offset_end = offset_lst[len(index_shift_lst) - 1
                            - index_shift_lst[::-1].index(best_match)] + k

    # match start in query
    i, shift, offset = match_lst[index_shift_lst.index(best_match)]
    match_start = -(shift - offset)

    #print(match_lst)
    #print(best_match)
    #for i, line in enumerate(match_lst):
    #    print(i, match_lst)

    return index_match, offset_start, offset_end, match_start


def get_similarity_bars(seq1, seq2):
    """Produces a str of '  ' and '|' indicating identity between two sequences

    :param seq1: str
    :param seq2: str (same length as 1!)
    :return: str of '  ' and '|' indicating identity
    """
    assert len(seq1) == len(seq2), "seq1 and seq2 must be same length"
    bars = []
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            bars += '|'
        else:
            bars += ' '
    return ''.join(bars)


def print_match(query, seq_lst, index, offset_start, offset_end, match_start,
                padded_bases = 3, match_info=[]):
    """Function to print query, part of matched sequence with | lines in between.

    :param query: str of match query
    :param seq_lst: lst of sequences that were matched
    :param index: int python index of found match
    :param offset_start: int python index within matched sequence of start match
    :param offset_end: int python index within matched sequence of start match
    :param match_start: int index of match start in query
    :param match_info: default False to not print extra info. from wrapper func
                    lst of parameters: !!!!TO FILL!!!!!

    :return: doesn't.
    """
    print('Longest match found:')
    query_padded = padded_bases*' ' + query + padded_bases*' '
    print(query_padded)
    # Calc new spots to display
    start = offset_start - match_start - padded_bases
    end = offset_start - match_start + len(query) + padded_bases
    match_seq = seq_lst[index][start:end]

    print(get_similarity_bars(query_padded, match_seq))
    print(match_seq)
    #print(seq_lst[index])
    print('In sequence:', index + 1, '\nFrom position:', offset_start + 1)

    if match_info:
        print('!!!TO DO!!')


def find_match_wrap(query, seq_lst, k = 2, print_enabled=True,
                    scrub_seq_lst=True, hash_dct={}, match_lst=[]):
    """Finds longest match in seq_lst (AGCT chars allowed) and returns info

    :param query: str of AGCT query (returns exception if not)
    :param seq_lst: lst of AGCT sequences (other chars become A)
    :param k: int of k-mer size to use (default = 2)
    :param print_enabled: Bool to enable printing of found match
    :param scrub_seq_lst: Bool to enable non-AGCT
    :param hash_dct: OPTIONAL dict with k-mers as keys and as values a lst of
                                                (seq_nr, pos_in_seq)
    :param match_lst: OPTIONAL sorted list of query hits in hash_dct
    :return: hash_dct,
             int of python-based index of matched sequence in seq_lst,
             int of match offset start in matched sequence,
             int of match offset end in matched sequence
    """
    query = query.upper()
    assert set(query) == set('AGCT'), "query not a nucleotide seq of 'ACGT'"

    # make hash table
    if not hash_dct:
        #print('making hash dict')
        hash_dct = make_hash_table(seq_lst, k, clean_sequences=scrub_seq_lst)

    # get sorted matches with hash_dct
    if not match_lst:
        #print('making match_lst')
        match_lst = query_kmer_hits(query, hash_dct, k)

    # get longest match
    index, offset_start, offset_end, match_start = find_longest_match(match_lst, k)

    if print_enabled:
        print_match(query, seq_lst, index, offset_start, offset_end, match_start,
                    padded_bases=3, match_info=[])

    return hash_dct, index, offset_start, offset_end


def reverse_complement(dna_seq):
    """Reverses + complements dna_seq

    dna_seq is uppercase DNA sequence
    :returns: reverse complement of input
    """
    comp_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    reverse_seq = dna_seq[::-1]

    rev_comp = []
    for ch in reverse_seq:  # loop over characters in reversed to complement
        rev_comp.append(comp_dict[ch])  # replace with complements
    return ''.join(rev_comp)


def main():
    query = 'TGCAACAT'
    s1 = 'GTGACGTCACTCTGAGGATCCCCTGGGTGTGG'
    s2 = 'GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT'
    s3 = 'GGATCCCCTGTCCTCTCTGTCACATA'
    seqs = [s1, s2, s3]

    # Find query in seqs:
    find_match_wrap(query, seqs, k=2)

    try:
        print(10*'-')
        ara_seqs = fasta_to_sequence_lst('misc/TAIR10.fasta')
        print('Arabidopsis fasta load done. \nNumber of chromosomes:', len(ara_seqs),
              '\ntotal sequence length:', int(sum([len(seq) for seq in ara_seqs])))
        ara_hash_dct = make_hash_table(ara_seqs, k=15)
        ara_queries = fasta_to_sequence_lst('misc/athal_query.fasta')
        for i, query in enumerate(ara_queries):
            print(10*'-')
            print('Query', i+1)
            find_match_wrap(query, ara_seqs, hash_dct=ara_hash_dct, k=15)

    except IOError:
        print("TAIR10.fasta or atha1_query.fasta files not found, did you download misc folder from repository?")


if __name__ == "__main__":
    main()
