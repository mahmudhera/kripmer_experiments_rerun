import pandas as pd
import random
from Bio import SeqIO

num_targets = 50
target_length = 150
random_seed = 1
fasta_filename = '../../data/yeast/GCF_000146045.2_R64_genomic.fna'
output_dir = '../../data/yeast/yeast_targets'

if __name__ == "__main__":
    df = pd.read_csv('yeast_gRNAs_multiple.sam', delimiter='\t', header=None)
    random.seed(random_seed)

    contig_names = df[2].tolist()
    start_pos = df[3].tolist()
    list_positions = list( zip(contig_names, start_pos) )

    # test start
    contig_name_test = contig_names[0]
    start_pos_test = start_pos[0]
    for seq_record in SeqIO.parse(fasta_filename, "fasta"):
        if seq_record.id == contig_name_test:
            print( str( seq_record.seq[start_pos_test-10:start_pos_test+30] ) )
    # test end

    random.shuffle(list_positions)
    num_targets_generated = 0
    file_id = 101

    for position in list_positions:
        contig_name = position[0]
        start_pos = position[1]
        if start_pos < 80:
            continue

        for seq_record in SeqIO.parse(fasta_filename, "fasta"):
            if seq_record.id == contig_name:
                fname = output_dir + '/' + f'target{file_id}.fasta'
                f = open(fname, 'w')
                f.write('> ' + f'target{file_id}_' + contig_name + "_" + str(start_pos-75) + '_' + str(start_pos+75) + '\n')
                file_id += 1
                f.write(str(seq_record.seq[start_pos-75:start_pos+75]).upper())
                f.close()

        num_targets_generated+=1
        if num_targets_generated == num_targets:
            break
