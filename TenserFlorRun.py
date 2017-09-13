#!python

import NucleotideCodeSubstitution
import FastaFile
import GffFile
import MyUtil

import numpy as np

import math

_buckets = []

def read_data(gffFile, fastaFile):
    chromosome_gene_dict, chromosome_gene_list = GffFile.readGff( gffFile )
    chromosome_names, fastas = FastaFile.readFastaFile( fastaFile )
    GffFile.update_sequence_information(fastas, chromosome_gene_dict)

    #get the gene with longest length
    #delete ORF-shirt transcript
    #delete gene without transcript
    #keep only on transcript for each gene
    longest_gene = 0
    for chromosome_name in chromosome_names:
        gene_names_to_delete = []
        for gene_name in chromosome_gene_dict[chromosome_name]:
            if chromosome_gene_dict[chromosome_name][gene_name].length > longest_gene:
                longest_gene = chromosome_gene_dict[chromosome_name][gene_name].length

            # delete ORF-shift transcript begin
            transcript_number = 0
            for transcript in chromosome_gene_dict[chromosome_name][gene_name].transcripts:
                NucleotideCodeSubstitution.checkOrfState(transcript, fastas)
                if not transcript.if_orf_conserved:
                    # print("15 " + transcript.name + transcript.meta_informaiton)
                    chromosome_gene_dict[chromosome_name][gene_name].transcripts=np.delete(chromosome_gene_dict[chromosome_name][gene_name].transcripts, transcript_number, 0)
                    transcript_number = transcript_number - 1
                transcript_number = transcript_number + 1
            #delete ORF-shift transcript end

            ##delete transcript so that the there is only one transcript left for each gene
            while len(chromosome_gene_dict[chromosome_name][gene_name].transcripts)>1:
                chromosome_gene_dict[chromosome_name][gene_name].transcripts = np.delete(chromosome_gene_dict[chromosome_name][gene_name].transcripts, 1, 0)
            #delete those genes contrains no transcript
            if len(chromosome_gene_dict[chromosome_name][gene_name].transcripts) == 0 :
                gene_names_to_delete.append(gene_name)

        for gene_name in gene_names_to_delete:
            del chromosome_gene_dict[chromosome_name][gene_name]

    windowns_size = math.floor(longest_gene * 2.5) + 1

    _buckets.append((windowns_size, windowns_size))
    _buckets.append((longest_gene, longest_gene))
    _buckets.append((windowns_size * 2, windowns_size * 2))
    _buckets.append((windowns_size * 3, windowns_size * 3))
    _buckets.append((windowns_size * 4, windowns_size * 4))
    overLapSize = math.floor(longest_gene * 1.1)+1   ## a bit of larger than one, so that the longest gene could be located in one window

    data_set = [[] for _ in _buckets]
    source_output = open('source.txt', 'w')
    target_output = open('target.txt', 'w')
    target_negative_output = open('target_negative_output.fasta', 'w')
    for chromosome_name in chromosome_names:
        # positive strand
        states = MyUtil.get_genetic_region_states(chromosome_name, "+", chromosome_gene_dict, fastas)
        print(">" + chromosome_name + "+")
        seq = fastas[chromosome_name].seq
        seq_code = NucleotideCodeSubstitution.dna_to_matix3(seq)
        #begin to splice the sequence into windows
        start = 1 - overLapSize
        end = 1
        while end < len(seq):
            start = start + overLapSize
            end = start + windowns_size
            if start > len(seq):
                start = len(seq)
            if end > len(seq):
                end = len(seq)
            overlap_result = MyUtil.overlap_with_certain_gene(start, chromosome_name, "+", chromosome_gene_dict)
            if not None == overlap_result:
                start = chromosome_gene_dict[chromosome_name][overlap_result].end -1
            overlap_result = MyUtil.overlap_with_certain_gene(end, chromosome_name, "+", chromosome_gene_dict)
            if not None == overlap_result:
                end = chromosome_gene_dict[chromosome_name][overlap_result].end + 1
            slice_start = start - 1
            slice_end = end
            source_ids = seq_code[slice_start:slice_end]
            target_ids = states[slice_start:slice_end]
            source_output.writelines(["%s " % item for item in source_ids])
            source_output.write("\n")
            target_output.writelines(["%s " % item for item in target_ids])
            target_output.write("\n")
            for bucket_id, (source_size, target_size) in enumerate(_buckets):
                if (len(source_ids) < source_size) and (len(target_ids) < target_size):
                    data_set[bucket_id].append([source_ids, target_ids])
                    break # stop the for loop

        # negative strand
        states = MyUtil.get_genetic_region_states(chromosome_name, "-", chromosome_gene_dict, fastas)
        print(">" + chromosome_name + "-")
        seq_rev_c = NucleotideCodeSubstitution.getReverseComplementary(seq)
        seq_code_rev_c = NucleotideCodeSubstitution.dna_to_matix3(seq_rev_c)
        seq_code = seq_code_rev_c
        start = 1 - overLapSize
        end = 1
        target_negative_output.write(">" + chromosome_name + "\n")
        target_negative_output.write(''.join(states) + "\n")
        while end < len(seq):
            start = start + overLapSize
            end = start + windowns_size
            if start > len(seq):
                start = len(seq)
            if end > len(seq):
                end = len(seq)
            overlap_result = MyUtil.overlap_with_certain_gene(start, chromosome_name, "-", chromosome_gene_dict)
            if not None == overlap_result:
                start = chromosome_gene_dict[chromosome_name][overlap_result].end - 1
            overlap_result = MyUtil.overlap_with_certain_gene(end, chromosome_name, "+", chromosome_gene_dict)
            if not None == overlap_result:
                end = chromosome_gene_dict[chromosome_name][overlap_result].end + 1
            slice_start = start - 1
            slice_end = end
            source_ids = seq_code[slice_start:slice_end]
            target_ids = states[slice_start:slice_end]
            source_output.writelines(["%s " % item for item in source_ids])
            source_output.write("\n")
            target_output.writelines(["%s " % item for item in target_ids])
            target_output.write("\n")
            for bucket_id, (source_size, target_size) in enumerate(_buckets):
                if (len(source_ids) < source_size) and (len(target_ids) < target_size):
                    data_set[bucket_id].append([source_ids, target_ids])
                    break  # stop the for loop

    source_output.close()
    target_output.close()
    return data_set

print ("begin to run")
read_data("/Users/song/tair11/TAIR10_GFF3_genes_no_UM.gff", "/Users/song/Col.fa")
print ("stop running")
