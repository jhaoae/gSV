import pysam
import Levenshtein
import os

def cal_position(r, ref_pos, length):

    cigar = r.cigar
    pos_rf = 0
    pos_read = 0
    start = 0
    end = r.query_length
    length_now = length
    ssoftclip = 0
    if r.reference_start <= ref_pos:
        for k in range(len(cigar)):
            cigar_label = cigar[k][0]
            cigar_len = cigar[k][1]
            if (cigar_label == 0) or (cigar_label == 7) or (cigar_label == 8):
                # match or mismatch
                pos_rf = pos_rf + cigar_len
                pos_read = pos_read + cigar_len
            elif cigar_label == 1:
                # insertion
                pos_read = pos_read + cigar_len
            elif cigar_label == 2:
                # deletion
                pos_rf = pos_rf + cigar_len
            elif cigar_label == 4:
                # softclip
                pos_read = pos_read + cigar_len
            k += 1
            if r.reference_start + pos_rf >= ref_pos:
                start = pos_read - (r.reference_start + pos_rf - ref_pos) - 1
                break
        if start + length > r.query_length:
            length_now = end-start
        else:
            end = start+length
            length_now = length
        rf_start = ref_pos
    else:
        for k in range(len(cigar)):
            cigar_label = cigar[k][0]
            cigar_len = cigar[k][1]
            if (cigar_label == 0) or (cigar_label == 7) or (cigar_label == 8):
                # match or mismatch
                pos_rf = pos_rf + cigar_len
                pos_read = pos_read + cigar_len
            elif cigar_label == 1:
                # insertion
                pos_read = pos_read + cigar_len
            elif cigar_label == 2:
                # deletion
                pos_rf = pos_rf + cigar_len
            elif cigar_label == 4:
                # softclip
                pos_read = pos_read + cigar_len
                if k==0:
                    ssoftclip = cigar_len
            k += 1
            if r.reference_start + pos_rf >= ref_pos + length :
                end = pos_read - (r.reference_start + pos_rf - ref_pos - length) - 1# + ssoftclip
                break
        rf_start = r.reference_start #- ssoftclip
        length_now = end - start
    return start, end, length_now, rf_start, r.query_length, ssoftclip


def takeend(elem):
    return elem[1]


def sort(reads, ref_pos, length):
    reads_info = list()
    min_mapq = 5
    for r in reads:
        if r.mapping_quality < min_mapq:
            continue
        start, end, length_now, rf_start, read_length, _= cal_position(r, ref_pos, length)
        reads_info.append([r, length_now, rf_start, read_length])
    reads_info.sort(key=lambda x: (x[2], -int(x[1]), -int(x[3])))# reverse=True)
    return reads_info



def seq_overlap(s1,s2):
    length = 100000
    s1_start,s1_end,s1_length_now, rf_start, read_length,_  =  cal_position(s1, s2.reference_start, length)
    if s2.query_length <= s1_length_now:
        s1_end = s1_start + s2.query_length
        s1_seq = s1.query_sequence[s1_start:s1_end]
        s2_seq = s2.query_sequence
    else:
        s1_seq = s1.query_sequence[s1_start:s1_end]
        s2_seq = s2.query_sequence[0:s1_length_now]
    return s1_seq, s2_seq

def clustering(reads, ref_pos, length, min_ratio):
    class1 = list()
    class2 = list()
    if length < 50:
        length = length + 100
    reads_info = sort(reads, ref_pos, 2*length+100)
    i = 0
    j = 0
    num = 0
    ref_pos_new = ref_pos
    
    for read in reads_info:
        num += 1
        r,length_now,_,_ = read
        if i == 0:
            class1.append(r.query_name)
            class1_read = r
            i += 1
        else:
            if j == 0:
                if r.reference_start <= ref_pos:
                    class1_start,class1_end, class1_length,_,_,_ = cal_position(class1_read, ref_pos, 2*length+100)
                    read_start, read_end, read_length,_,_,_ = cal_position(r, ref_pos, 2*length+100)
                    class1_seq = class1_read.query_sequence[class1_start:class1_start+min(read_length,class1_length)]
                    read_seq = r.query_sequence[read_start:read_start+min(read_length,class1_length)]
                else:
                    difference = r.reference_start-class1_read.reference_start
                    if class1_read.query_length - difference >0:
                        class1_start,class1_end, class1_length,_,_,_ = cal_position(class1_read, r.reference_start+10, 2*length+100)
                        read_start, read_end, read_length,_,_,sc_len = cal_position(r, r.reference_start+10, 2*length+100)
                        dis = min(class1_start,read_start+sc_len)
                        class1_seq = class1_read.query_sequence[class1_start-dis:class1_start-dis+min(read_length,class1_length)]
                        read_seq = r.query_sequence[read_start-dis:read_start-dis+min(read_length,class1_length)]
                    else:
                        continue
                sim_1 = Levenshtein.ratio(read_seq, class1_seq)
                if (sim_1 < min_ratio) :
                    class2.append(r.query_name)
                    class2_read = r
                    j += 1
                elif (sim_1 >= min_ratio):
                    class1.append(r.query_name)
                    i += 1
            else:
                if (class1_read.reference_start <= ref_pos) & (class2_read.reference_start <= ref_pos) & (r.reference_start <= ref_pos):
                    class1_start,class1_end, class1_length,rf1,_,_ = cal_position(class1_read, ref_pos, 2*length+100)
                    class2_start,class2_end, class2_length,rf2,_,_ = cal_position(class2_read, ref_pos, 2*length+100)
                    read_start, read_end, read_length,rf3,_,_ = cal_position(r, ref_pos, 2*length+100)
                    length_seq = min([class1_length,class2_length,read_length])
                    class1_seq = class1_read.query_sequence[class1_start:class1_start+length_seq]
                    class2_seq = class2_read.query_sequence[class2_start:class2_start+length_seq]
                    read_seq = r.query_sequence[read_start:read_start+length_seq]
                    sim_1 = Levenshtein.ratio(read_seq, class1_seq)
                    sim_2 = Levenshtein.ratio(read_seq, class2_seq)
                else:
                    ref_sort = sorted([class1_read.reference_start,class2_read.reference_start,r.reference_start])
                    difference1 = r.reference_start-class1_read.reference_start
                    difference2 = r.reference_start-class2_read.reference_start
                    if (class1_read.query_length - difference1 >0) & (class2_read.query_length - difference2 >0):
                        class1_start, class1_end, class1_length,rf1,_,_ = cal_position(class1_read, r.reference_start+10, 2*length+100)
                        class2_start, class2_end, class2_length,rf2,_,_ = cal_position(class2_read, r.reference_start+10, 2*length+100)
                        read_start, read_end, read_length,rf3,_,sc_len = cal_position(r, r.reference_start+10, 2*length+100)
                        dis = min([class1_start,class2_start,read_start+sc_len])
                        length_seq = min([class1_length,class2_length,read_length])
                        class1_seq = class1_read.query_sequence[class1_start-dis:class1_start+length_seq]
                        class2_seq = class2_read.query_sequence[class2_start-dis:class2_start+length_seq]
                        read_seq = r.query_sequence[read_start-dis:read_start+length_seq]
                        sim_1 = Levenshtein.ratio(read_seq, class1_seq)
                        sim_2 = Levenshtein.ratio(read_seq, class2_seq)
                    else:
                        continue
   
                sim_12 = Levenshtein.ratio(class1_seq, class2_seq)

                if (sim_2 < 0.5) & (sim_1 < 0.5):
                    continue
                elif sim_2 > sim_1:
                    class2.append(r.query_name)
                    j += 1
                elif sim_2 < sim_1:
                    class1.append(r.query_name)
                    i += 1


     
    return class1, class2


def bam2fa(Start, End, chro, path, bam, read_set, name):
    outfile = pysam.AlignmentFile(path + 'out.bam','wb',template=bam)
    reads = bam.fetch(chro, Start, End)
    for r in reads:
        if r.query_name in read_set:
            outfile.write(r) 
    outfile.close()
    os.system('samtools bam2fq ' +path + 'out.bam > ' + path + name + '.fastq')
    os.system('rm ' + path + 'out.bam')

def bam2fa_withoutcluster(Start, End, chro, path, bam, name):
    outfile = pysam.AlignmentFile(path + 'out.bam','wb',template=bam)
    reads = bam.fetch(chro, Start, End)
    for r in reads:
        outfile.write(r)
    outfile.close()
    os.system('samtools bam2fq ' +path + 'out.bam > ' + path + name + '.fastq')
    os.system('rm ' + path + 'out.bam')

def assembly(path,depth_mean):
    file_list = os.listdir(path)
    for c in file_list:
        if c.endswith('fastq'):
            name = c.split('/')[-1][:-6].split('-')
            print(name)
            if 'genotype' in c:
                depth = str(int(depth_mean/2))
            else:
                depth = str(int(depth_mean))
 
            print(depth)
            if int(name[2])-int(name[1]) < 100:
                os.system('wtdbg2  -k 22 -p 1  -X '+depth+' -L 1000 -t 16 -i ' + path + c + ' -fo' + path + c.split('/')[-1][:-6] + '.wtdbg')
            else:
                os.system('wtdbg2  -k 2 -p 19 -X '+depth+' -t 16  -i ' + path + c + ' -fo' + path + c.split('/')[-1][:-6] + '.wtdbg')

            os.system('wtpoa-cns -t 16 -i ' + path + c.split('/')[-1][:-6] + '.wtdbg.ctg.lay.gz -fo ' + path + c.split('/')[-1][:-6] + '.fa')
            os.system('rm ' + path + c.split('/')[-1][:-6] + '.wtdbg* ')

            if os.path.getsize(path + c.split('/')[-1][:-6] + '.fa')==0:
                if int(name[2])-int(name[1]) >= 100:
                    os.system('wtdbg2  -k 2 -p 19 -X '+depth+' -R -s -t 16  -i ' + path + c + ' -fo' + path + c.split('/')[-1][:-6] + '.sup.wtdbg')
                    os.system('wtpoa-cns -t 16 -i ' + path + c.split('/')[-1][:-6] + '.sup.wtdbg.ctg.lay.gz -fo ' + path + c.split('/')[-1][:-6] + '.fa')
                    os.system('rm ' + path + c.split('/')[-1][:-6] + '.sup.wtdbg* ')
    os.system('rm ' + path + '*.fastq')
    os.system('find ' + path + ' -name "*" -type f -size 0c | xargs -n 1 rm -rf')




