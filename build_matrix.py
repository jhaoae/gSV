import numpy as np

def ref_matrix(ref):
    length = len(ref)
    rf_array = np.zeros((4, length)).astype(np.int32)
    for i in range(length):
        ref_base = ref[i]
        if ref_base == 'A' or ref_base == 'a':
            rf_array[0][i] = 1
        elif ref_base == 'T' or ref_base == 't':
            rf_array[1][i] = 1
        elif ref_base == 'C' or ref_base == 'c':
            rf_array[2][i] = 1
        elif ref_base == 'G' or ref_base == 'g':
            rf_array[3][i] = 1
    return rf_array

def read_matrix(cigar, read, chrom, start, end, r_name, min_len, f):
    r_array = np.zeros((4, end-start)).astype(np.int32)
    pos_matrix = 0
    pos_read = 0
    for i in range(len(cigar)):
        cigar_label = cigar[i][0]
        cigar_len = cigar[i][1]
        cigar_seq = read[pos_read:(pos_read+cigar_len)]
        if (cigar_label == 0) or (cigar_label == 7) or (cigar_label == 8):
            #match or mismatch
            r_array[:, pos_matrix:(pos_matrix + cigar_len)]=cigar_matrix(cigar_seq)
            if (cigar_label == 8) and (cigar_len > min_len):
                line = str(chrom)+'\t'+str(start+pos_matrix)+'\t'+r_name+'\t'+str(cigar_len)+'\tMISS\t.\n'
                f.write(line)
            pos_matrix = pos_matrix + cigar_len
            pos_read = pos_read + cigar_len

        elif cigar_label == 1:
            #insertion
            r_array[:, pos_matrix-1:pos_matrix] = r_array[:,pos_matrix-1:pos_matrix]+ins_matrix(cigar_seq)
            if cigar_len > min_len:
                line = str(chrom)+'\t'+str(start+pos_matrix)+'\t'+r_name+'\t'+str(cigar_len)+'\tINS\t'+cigar_seq+'\n'
                f.write(line)
            pos_read = pos_read + cigar_len
            
        elif cigar_label == 2:
            #deletion
            if cigar_len > min_len:
                line = str(chrom)+'\t'+str(start+pos_matrix)+'\t'+r_name+'\t'+str(cigar_len)+'\tDEL\t.\n'
                f.write(line)
            pos_matrix = pos_matrix+cigar_len

        elif cigar_label == 4:
            #softclip
            pos_read = pos_read + cigar_len
            if cigar_len > min_len:
                if i==0:
                    line = str(chrom)+'\t'+str(start+pos_matrix)+'\t'+r_name+'\t'+str(cigar_len)+'\tSoftclip_1\t.\n'
                elif i == (len(cigar)-1):
                    line = str(chrom)+'\t'+str(start+pos_matrix)+'\t'+r_name+'\t'+str(cigar_len)+'\tSoftclip_2\t.\n'
                else:
                    line = str(chrom)+'\t'+str(start+pos_matrix)+'\t'+r_name+'\t'+str(cigar_len)+'\tSoftclip_3\t.\n'
                f.write(line)
        else:
            print('special cigar label:',cigar_label)
    return r_array


def read_matrix_reverse(cigar, read, chrom, start, end, r_name, flag, min_len, f):
    r_array = np.zeros((4, end-start)).astype(np.int32)
    pos_matrix = 0
    pos_read = 0
    if flag == 1:
        reve_read = DNA_reverse(DNA_complement(read))
    reve_read = DNA_reverse(read)
    for i in range(len(cigar)):
        cigar_label = cigar[i][0]
        cigar_len = cigar[i][1]
        cigar_seq = reve_read[pos_read:(pos_read + cigar_len)]
        if (cigar_label == 0) or (cigar_label == 7) or (cigar_label == 8):
            # match or mismatch
            r_array[:, pos_matrix:(pos_matrix + cigar_len)] = cigar_matrix(cigar_seq)
            if (cigar_label == 8) and (cigar_len > min_len):
                line = str(chrom)+'\t'+str(start+pos_matrix)+'\t'+r_name+'\t'+str(cigar_len)+'\tMISS\t.\n'
                f.write(line)
            pos_matrix = pos_matrix + cigar_len
            pos_read = pos_read + cigar_len

        elif cigar_label == 1:
            # insertion
            r_array[:, pos_matrix - 1:pos_matrix] = r_array[:, pos_matrix - 1:pos_matrix] + ins_matrix(cigar_seq)
            
            if cigar_len > min_len:
                line = str(chrom)+'\t'+str(start+pos_matrix)+'\t'+r_name+'\t'+str(cigar_len)+'\tINS\t'+cigar_seq+'\n'
                f.write(line)
            pos_read = pos_read + cigar_len

        elif cigar_label == 2:
            # deletion
            if cigar_len > min_len:
                line = str(chrom)+'\t'+str(start+pos_matrix)+'\t'+r_name+'\t'+str(cigar_len)+'\tDEL\t.\n'
                f.write(line)
            pos_matrix = pos_matrix + cigar_len

        elif cigar_label == 4:
            # softclip
            pos_read = pos_read + cigar_len
            if cigar_len > min_len:
                if i==0:
                    line = str(chrom)+'\t'+str(start+pos_matrix)+'\t'+r_name+'\t'+str(cigar_len)+'\tSoftclip_1\t.\n'
                elif i == (len(cigar)-1):
                    line = str(chrom)+'\t'+str(start+pos_matrix)+'\t'+r_name+'\t'+str(cigar_len)+'\tSoftclip_2\t.\n'
                else:
                    line = str(chrom)+'\t'+str(start+pos_matrix)+'\t'+r_name+'\t'+str(cigar_len)+'\tSoftclip_3\t.\n'
                f.write(line)
        else:
            print('special cigar label:', cigar_label)
    return r_array

def DNA_reverse(sequence):
    return sequence[::-1]

def DNA_complement(sequence):
    trantab = str.maketrans('ATCGatcgRYMKrymkVBHDvbhd','TAGCtagcYRKMyrkmBVDHbvdh')
    string = sequence.translate(trantab)
    return string

def cigar_matrix(read):
    length = len(read)
    ci_array = np.zeros((4, length)).astype(np.int32)
    for i in range(length):
        ref_base = read[i]
        if ref_base == 'A' or ref_base == 'a':
            ci_array[0][i] = 1
        elif ref_base == 'T' or ref_base == 't':
            ci_array[1][i] = 1
        elif ref_base == 'C' or ref_base == 'c':
            ci_array[2][i] = 1
        elif ref_base == 'G' or ref_base == 'g':
            ci_array[3][i] = 1
    return ci_array

def ins_matrix(read):
    ins_array = np.zeros((4, 1)).astype(np.int32)
    for i in range(len(read)):
        ref_base = read[i]
        if ref_base == 'A' or ref_base == 'a':
            ins_array[0][0] = ins_array[0][0]+1
        elif ref_base == 'T' or ref_base == 't':
            ins_array[1][0] = ins_array[1][0]+1
        elif ref_base == 'C' or ref_base == 'c':
            ins_array[2][0] = ins_array[2][0]+1
        elif ref_base == 'G' or ref_base == 'g':
            ins_array[3][0] = ins_array[3][0]+1
    return ins_array
