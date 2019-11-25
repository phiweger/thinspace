#!/usr/bin/env python3
"""
This script parses a Centrifuge-generated SAM file to give some stats about the classifications.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""

import argparse
import collections
import gzip
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Summarise classifications from Centrifuge SAM '
                                                 'files')
    parser.add_argument('--centrifuge', type=str, required=True,
                        help='a Centrifuge output file (can be gzipped)')
    parser.add_argument('--tree', type=str, required=True,
                        help='a taxonomy tree file (can be gzipped)')

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()

    tax_id_to_parent, tax_id_to_rank = load_tax_info(args.tree)
    tax_ids_per_read = load_tax_ids_per_read(args.centrifuge)

    missing_tax_ids = set()
    for tax_ids in tax_ids_per_read.values():
        for tax_id in tax_ids:
            if tax_id not in tax_id_to_rank:
                missing_tax_ids.add(tax_id)

    for tax_id in sorted(missing_tax_ids):
        print(tax_id)


def load_tax_info(tree_filename):
    """
    This function reads through the tree file, and returns three dictionaries:
      1) Where the keys are tax IDs and the values are the parent tax IDs. This dictionary allows
         for tracing 'upward' (toward the root) through the tree, starting at any node.
      2) Where the keys are tax IDs and the values are taxonomic ranks. Importantly, this
         dictionary does not include all taxonomic ranks in the tree, but only the standard levels
         (phylum, class, order, etc). When a tax ID has a non-standard rank (e.g. subfamily), that
         ID is given the first standard rank found in its ancestors (e.g. family).
      3) Where the keys are tax IDs and the values are tax IDs for the first ancestor that has a
         standard taxonomic rank.
    """
    tree_data = []
    open_func = get_open_func(tree_filename)
    with open_func(tree_filename, 'rt') as tree_file:
        for line in tree_file:
            parts = line.strip().split('\t')
            tree_data.append([int(parts[0]), int(parts[2]), parts[4].lower()])

    tax_id_to_parent = {}
    for tax_id, parent_id, _ in tree_data:
        tax_id_to_parent[tax_id] = parent_id

    tax_id_to_rank = {0: 'unclassified'}
    for tax_id, parent_id, rank in tree_data:
        tax_id_to_rank[tax_id] = rank
        if tax_id == 1:  # special case for the root node
            assert tax_id == parent_id  # the root is its own parent
            tax_id_to_rank[tax_id] = 'root'

    return tax_id_to_parent, tax_id_to_rank


def load_tax_ids_per_read(centrifuge_filename):
    """
    Returns a dictionary where the key is the read name and the value is a list of all tax IDs that
    read was classified to.
    """
    tax_ids_per_read = collections.OrderedDict()
    with get_open_func(centrifuge_filename)(centrifuge_filename, 'rt') as centrifuge_file:
        for line in centrifuge_file:
            parts = line.strip().split('\t')
            read_name = parts[0]
            if read_name == 'readID':  # header
                continue
            seq_id, tax_id = parts[1], int(parts[2])
            if seq_id == 'no rank':
                tax_id = 1
            if tax_id == 0:
                # A taxID of 0 with a non-unclassified seqID implies something went wrong in the
                # index building.
                assert seq_id == 'unclassified'
            if read_name not in tax_ids_per_read:
                tax_ids_per_read[read_name] = set()
            tax_ids_per_read[read_name].add(tax_id)
    return tax_ids_per_read


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


if __name__ == '__main__':
    main()
