import pysam
from pysam import VariantFile
from encoding import read_encoding
from build_matrix import ref_matrix
from cluster import clustering, bam2fa, assembly, bam2fa_withoutcluster
from detection import exactmatch, split_fa
from graph_cut_detection import graph_cut_detection
import numpy as np
import os, time
from datetime import datetime
from multiprocessing import Pool
from functools import partial 
from collections import Counter


def p_encode(chrom,bam_path,ref_path,out_path,min_mapq,min_len):

    print('Run Encoding %s (%s)...' % (chrom, os.getpid()))
    start = time.time()

    ref_file = pysam.FastaFile(ref_path)
    rf = ref_file.fetch(chrom)
    rf_array = ref_matrix(rf)
    end = time.time()
    print('Encoding ref %s runs %0.2f seconds.' % (chrom, (end - start)))

    bam_file = pysam.AlignmentFile(bam_path)
    
    sub_num = 200000
    depth_mean, loss_re, support_sum = read_encoding(bam_file, rf_array, out_path, chrom, min_mapq, min_len, sub_num)
    
    print(chrom,depth_mean,loss_re)
    
    if np.isnan(depth_mean):
        print('no reads in '+chrom)
    else:
        data = {'depth_mean':depth_mean,'e0':loss_re}
        file_name = out_path+'/'+chrom+'_graph_cut_input.npy'
        np.save(file_name,data)
    
    bam_file.close()
    end = time.time()
    print('Encoding %s runs %0.2f seconds.' % (chrom, (end - start)))


def graph_cut(chrom,out_path,distance,beta,gamma,amplify):
    
    data = np.load(out_path+'/'+chrom+'_graph_cut_input.npy',allow_pickle=True).item()
    depth_mean = data['depth_mean'].astype(np.double)
    loss_re = data['e0'].astype(np.double)
    print(chrom,depth_mean,loss_re)

    print('Run Graph-Cut %s (%s)...' % (chrom, os.getpid()))
    start = time.time()
    rr,S = graph_cut_detection(loss_re,depth_mean,beta,gamma,amplify)
    end = time.time()
    print('Graph-Cut %s runs %0.2f seconds.' % (chrom, (end - start)))

    del loss_re     
    ### ~ Write graph cut results ~ ###
    f = open(out_path+'/'+chrom+'_rr.vcf', "w+")

    f.write("##fileformat=VCFv4.2\n")
    f.write("##fileDate=" + datetime.today().strftime("%m%d%Y") + "\n")
    f.write("##source=bit_representation_model\n")
    f.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End point of SV\">\n")

    header_vec = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_ids_here']
    f.write("#" + "\t".join(header_vec) + "\n")
        
    i = 0
    rr = np.array(rr)
    rr_merge = []
    for r in rr:
        i = i+1
        body_vec = [chrom, str(int(r[0])+1), "SV" + str(i), "N", ".", ".", "PASS", "END=" + str(int(r[1])+1), ".", "."]
        f.write("\t".join(body_vec) + "\n")
        if not rr_merge or r[0] > rr_merge[-1][1] + distance:
            rr_merge.append(r)
        else:
            rr_merge[-1][1] = max(rr_merge[-1][1], r[1])
    f.close()
        
         
    f_combine = open(out_path+'/'+chrom+'_rr_combine.vcf', "w+")

    f_combine.write("##fileformat=VCFv4.2\n")
    f_combine.write("##fileDate=" + datetime.today().strftime("%m%d%Y") + "\n")
    f_combine.write("##source=bit_representation_model\n")
    f_combine.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End point of SV\">\n")

    header_vec = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_ids_here']
    f_combine.write("#" + "\t".join(header_vec) + "\n")

    i = 1
    for r in rr_merge:
        body_vec = [chrom, str(int(r[0])+1), "SV" + str(i), "N", ".", ".", "PASS", "END=" + str(int(r[1])+1), ".", "."]
        i = i+1
        f_combine.write("\t".join(body_vec) + "\n")
    f_combine.close()
        

def detect(chrom, bam_path, ref_path, out_path, complex_mode, pos_range, L, support_read, cluster_similarity, mempath, memlen, min_SV_len):
    
    print('Run Detection %s (%s)...' % (chrom, os.getpid()))
    start_time = time.time()
    
    ref_file = pysam.FastaFile(ref_path)
    bam_file = pysam.AlignmentFile(bam_path)
        
    preliminary_results = []
    with open(out_path+'/'+chrom+"_preliminary.txt", "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            data = [parts[0], int(parts[1]), parts[2], int(parts[3]), parts[4], parts[5]]
            preliminary_results.append(data)


    path_fa = out_path+'/'+chrom+'/'
    if not os.path.exists(path_fa):
        os.makedirs(path_fa)

    

    vcf_file = open(path_fa+chrom+"_simple_region.vcf", "w+")

    vcf_file.write("##fileformat=VCFv4.2\n")
    vcf_file.write("##source=bit_representation_model\n")
    vcf_file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant described in this record\">\n")
    vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    vcf_file.write("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
    vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

    chromosome = ref_file.references

    for chro in chromosome:
        chro_len = ref_file.get_reference_length(chro)
        vcf_file.write("##contig=<ID="+chro+",length="+str(chro_len)+">\n")

    header_vec = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_ids_here']
    vcf_file.write("#" + "\t".join(header_vec) + "\n")

    

    complex_region = []
    ID = 1
    vcf = VariantFile(out_path+'/'+chrom+'_rr_combine.vcf',"r")
    
    for rec in vcf:
        chrom = rec.chrom
        start = rec.start
        end = rec.stop
        region_len = end - start
        INS_support = []
        DEL_support = []
        INV_support = []
        DUP_support = []
        TRA_support = []
        read_num = 0
        
        if complex_mode:
            complex_region.append([start,end])
        
        else:
            if end - start < 5:
                for read in preliminary_results:
                    if (read[1] < start+pos_range) and (read[1] > start-pos_range):
                        read_num += 1     
                        if read[4] == 'INS':
                            INS_support.append(read)
                        elif read[4] == 'DEL':
                            DEL_support.append(read)
                        elif read[4] == 'INV':
                            INV_support.append(read)
                        elif read[4] == 'DUP':
                            DUP_support.append(read)
                        elif read[4] == 'TRA':
                            TRA_support.append(read)
            else:
                for read in preliminary_results:
                    if (read[1] < end + pos_range) and (read[1] > start-pos_range):
                        read_num += 1
                        if (read[3] > region_len/L):
                            if read[4] == 'INS':
                                INS_support.append(read)
                            elif read[4] == 'DEL':
                                DEL_support.append(read)
                            elif read[4] == 'INV':
                                INV_support.append(read)
                            elif read[4] == 'DUP':
                                DUP_support.append(read)
                            elif read[4] == 'TRA':
                                TRA_support.append(read)
            support = [INS_support, DEL_support, INV_support, DUP_support, TRA_support]
            support_sort = sorted(support, key = len, reverse=True)
            if len(support_sort[0]) >= support_read:
                # one type
                if len(support_sort[1]) < support_read:
                    length = [item[3] for item in support_sort[0]]
                    counter_length = Counter(length)
                    SV_subtype = counter_length.most_common(2)
                    if len(SV_subtype) == 1:
                        SV_subtype.append((0,0))

                    if SV_subtype[0][1] == SV_subtype[1][1]:
                        SV_subtype_candi = [item for item in counter_length.most_common() if item[1] == SV_subtype[0][1]]
                        SV_subtype = sorted(SV_subtype_candi, key=lambda x: x[0], reverse=True)[0:2]

      
                    SV1 = []
                    SV2 = []
                    for item in support_sort[0]:
                        if abs(item[3] - SV_subtype[0][0]) < pos_range/2:
                            SV1.append(item)
                        elif abs(item[3]- SV_subtype[1][0]) < pos_range/2:
                            same_read = [same for same in preliminary_results if (same[2] == item[2] and same[4] == item[4] and abs(same[1]-item[1])<50)]

                            if item[4] == 'INS':
                                sum_length = sum(int(same[3]) for same in same_read)
                            else:
                                same_read = sorted(same_read, key=lambda x: x[1])
                                sum_length = same_read[-1][1]-same_read[0][1]+same_read[-1][3]
                            item[3] = sum_length
                            if len(same_read) > 1 and abs(sum_length-SV_subtype[0][0]) < pos_range/2:
                                SV1.append(item)
                            else:
                                SV2.append(item)
                        
                    SV_list = []
                    if abs(SV_subtype[0][0] - SV_subtype[1][0]) > pos_range/2 and len(SV2) >= support_read:
                        # two SV
                        SV_list = [SV1,SV2]
                    else:
                        # one SV
                        SV_list = [SV1]
                        if len(SV1) < support_read:
                            complex_region.append([start,end])
                            continue
                #two types
                elif len(support_sort[1]) >= support_read: 
                    if len(support_sort[2]) >= support_read:
                        complex_region.append([start,end])
                    SV1_support_name = [item[2] for item in support_sort[0]]
                    SV2_support_name = [item[2] for item in support_sort[1]]
                    intersection = [x for x in SV1_support_name if x in SV2_support_name]
                    if len(intersection) >= 1:
                        complex_region.append([start,end])
                        if intersection == SV2_support_name:
                            SV_list = []
                            SV1 = []
                            SV1_length = [item[3] for item in support_sort[0]]
                            SV1_counter_length = Counter(SV1_length)  
                            SV1_type = SV1_counter_length.most_common(1)
                            SV_subtype = SV1_type
                            for item in support_sort[0]:
                                if abs(item[3] - SV1_type[0][0]) < pos_range/2:
                                    SV1.append(item)
                    
                            if len(SV1) >= support_read:
                                SV_list.append(SV1)
                        else:
                            continue

                    else:
                        SV_list = []
                        SV1 = []
                        SV1_length = [item[3] for item in support_sort[0]]
                        SV1_counter_length = Counter(SV1_length)
                        SV1_type = SV1_counter_length.most_common(1)
                        SV1_pos = [item[1] for item in support_sort[0] if item[3] == SV1_type[0][0]]
            
                        for item in support_sort[0]:
                            if abs(item[3] - SV1_type[0][0]) < pos_range/2 and abs(item[1] - SV1_pos[0]) < pos_range/2:
                                SV1.append(item)
                        SV2 = []
                        SV2_length = [item[3] for item in support_sort[1]]
                        SV2_counter_length = Counter(SV2_length)
                        SV2_type = SV2_counter_length.most_common(1)
                        SV2_pos = [item[1] for item in support_sort[1] if item[3] == SV2_type[0][0]]
                        for item in support_sort[1]:
                            if abs(item[3] - SV2_type[0][0]) < pos_range/2 and abs(item[1] - SV2_pos[0]) < pos_range/2:
                                SV2.append(item)
                        SV_subtype = [SV1_type[0],SV2_type[0]]
                        if len(SV1) < support_read and len(SV2) < support_read:
                            complex_region.append([start,end])
                            continue
                        else:
                            if len(SV1) >= support_read:
                                SV_list.append(SV1)
                            if len(SV2) >= support_read:
                                SV_list.append(SV2)

                else:
                    complex_region.append([start,end])
                    continue


                for i in range(len(SV_list)):
                    SV = SV_list[i]
                    SV_type = SV[0][4]
                    SV_chrom = SV[0][0]
                    SV_support = len(SV)

                    start = [item[1] for item in SV]
                    counter_start = Counter(start)
                    SV_start = counter_start.most_common(1)[0][0]
                    SV_length_candi = []
                    for item in SV:
                        if item[1] == SV_start:
                            SV_length_candi.append(item[3])
                    SV_length = SV_length_candi[len(SV_length_candi)//2]

                    if SV_subtype[i][1] == 1:
                        SV = sorted(SV,key = lambda x:x[3])
                        SV_length = SV[len(SV)//2][3]
                        SV_start = SV[len(SV)//2][1]

                    SV_end = SV_start + SV_length
    
                    GT = '0/1'
        
                    if len(SV_list)==1 and  SV_support > (read_num/2):
                        GT = '1/1'

                    if SV_type == 'DEL':
                        ref_seq = ref_file.fetch(SV_chrom, SV_start, SV_end)
                        ctg_seq = ref_seq[0]
                        i==0
                    else:
                        ref_seq = ref_file.fetch(SV_chrom, SV_start, SV_start+1)
                        if SV_type == 'DUP':
                            ctg_seq = '<DUP>'
                        elif SV_type == 'INV':
                            ctg_seq = '<INV>'
                        elif SV_type == 'TRA':
                            ctg_seq = '<TRA>'
                        elif  SV_type  == 'INS':
                            SV_end = SV_start+1
                            for item in SV:
                                if item[1] == SV_start and abs(item[3]- SV_length) < 5:
                                    ctg_seq = item[5][0:min(item[3],SV_length)]
                                    SV_length = min(item[3],SV_length)
                                    break
                            SV_start = SV_end
                    if abs(SV_length) >= min_SV_len:
                        body_vec = [SV_chrom, str(SV_start), "SV" + str(ID), ref_seq, ctg_seq,".", "PASS", "END=" +str(SV_end)+ ";SVTYPE="+SV_type+";SVLEN="+str(SV_length), "GT", GT]
                        vcf_file.write("\t".join(body_vec) + "\n")
                        ID += 1
            else:
                complex_region.append([start,end])
        
    vcf_file.close()
    
    #complex region
    for i in range(len(complex_region)):
        interval_start = int(complex_region[i][0])
        interval_end = int(complex_region[i][1])
        
        reads = bam_file.fetch(chrom,interval_start, interval_end)
        length = interval_end - interval_start
            
        class1, class2 = clustering(reads, interval_start-100, length, cluster_similarity)
             
        if class2:
            if len(class1)>5 and len(class2)>5:
                if len(class1)/len(class2) >2:
                    name1 = chrom + '-' + str(interval_start) + '-' + str(interval_end) + '-1'
                    bam2fa(interval_start, interval_end, chrom, path_fa, bam_file, class1, name1)
                elif len(class2)/len(class1) >2:
                    name2 = chrom + '-' + str(interval_start) + '-' + str(interval_end) + '-2'
                    bam2fa(interval_start, interval_end, chrom, path_fa, bam_file, class2, name2)
                else:
                    name1 = chrom + '-' + str(interval_start) + '-' + str(interval_end) + '-1-genotype'
                    name2 = chrom + '-' + str(interval_start) + '-' + str(interval_end) + '-2-genotype'
                    bam2fa(interval_start, interval_end, chrom, path_fa, bam_file, class1, name1)
                    bam2fa(interval_start, interval_end, chrom, path_fa, bam_file, class2, name2)
            elif len(class1)>5:
                name1 = chrom + '-' + str(interval_start) + '-' + str(interval_end) + '-1'
                bam2fa(interval_start, interval_end, chrom, path_fa, bam_file, class1, name1)
            elif len(class2)>5:
                name2 = chrom + '-' + str(interval_start) + '-' + str(interval_end) + '-2'
                bam2fa(interval_start, interval_end, chrom, path_fa, bam_file, class2, name2)
        else:
            name1 = chrom + '-' + str(interval_start) + '-' + str(interval_end) + '-1'
            bam2fa(interval_start, interval_end, chrom, path_fa, bam_file, class1, name1)
        
    data = np.load(out_path+'/'+chrom+'_graph_cut_input.npy',allow_pickle=True).item()
    depth_mean = data['depth_mean'].astype(np.double)

    assembly(path_fa,depth_mean)
    
    split_fa(path_fa)
     

    exactmatch(path_fa, ref_file, chrom, mempath, memlen, min_SV_len)
        

    bam_file.close()
    ref_file.close()

    end_time = time.time()
    print('Detection %s runs %0.2f seconds.' % (chrom, (end_time - start_time)))


def summarize(ref_path, out_path, min_SV_len):
    
    ref_file = pysam.FastaFile(ref_path)

    os.system('mkdir '+ out_path+'/workspace')
    os.system('mv '+out_path+'/* '+out_path+'/workspace')
    vcf_path = out_path+'/results_vcf/'
        
    if not os.path.exists(vcf_path):
        os.makedirs(vcf_path)
        
    os.system('cp ' + out_path+'/workspace/*/*.vcf ' + vcf_path)
    os.system('vcf-concat '+ vcf_path +'*complex_SV.vcf > '+ vcf_path+'complex.vcf')
    os.system('vcf-sort -c '+ vcf_path +'complex.vcf > '+ out_path+'/gSV_complexSV.vcf')
    os.system('vcf-concat '+ vcf_path +'*region.vcf > '+ vcf_path+'sum.vcf')
    os.system('vcf-sort -c '+ vcf_path +'sum.vcf > '+ vcf_path+'sum.sort.vcf')
 
    with open(vcf_path+"sum.sort.vcf", "r") as vcf_ori:
        lines = vcf_ori.readlines()

    output = []
    overlap = []
    ID = 1
    for i in range(len(lines)):
        if lines[i].startswith('#'):
            output.append(lines[i])
            continue
        if i in overlap:
            continue

        if len(output) == 0:
            output.append(lines[i])
            
        SV = lines[i].split('\t')
        SV_detail = SV[7].split(';')
        SV_len = int(SV_detail[2].split('=')[1])


        for j in range(i,len(lines)):
            SV_next = lines[j].split('\t')
            SV_next_detail = SV_next[7].split(';')
            SV_next_len = int(SV_next_detail[2].split('=')[1])

            if SV[0] == SV_next[0] and abs(int(SV[1])-int(SV_next[1])) <= 250:
                if (abs(SV_len-SV_next_len)<50):
                    if SV_detail[1] == SV_next_detail[1]:
                        overlap.append(j)
                    else:
                        if SV_detail[1] == 'SVTYPE=DUP' and SV_next_detail[1] == 'SVTYPE=INS':
                            overlap.append(j)
                        elif SV_detail[1] == 'SVTYPE=INS' and SV_next_detail[1] == 'SVTYPE=DUP':
                            SV = SV_next
                            overlap.append(j)
            else:
                break
        if SV_detail[1] == 'SVTYPE=INS':
            seq = open(vcf_path+SV[0]+'_'+SV[1]+'.fa','a')
            seq.write('>'+SV[0]+'_'+SV[1]+'\n'+SV[4]+'\n')
            seq.close()

            ref_start = max(int(SV[1])-max(2000,SV_len),0)
            ref_end = min(int(SV[1])+max(2000,SV_len),ref_file.get_reference_length(SV[0]))
            ref = open(vcf_path+SV[0]+'_'+SV[1]+'_ref.fa','a')
            ref_seq = ref_file.fetch(SV[0], ref_start, ref_end)
            ref.write('>'+SV[0]+'_'+SV[1]+'_ref\n'+ref_seq+'\n')
            ref.close()

            os.system('minimap2 -a -k 9 --secondary=no ' + vcf_path+SV[0]+'_'+SV[1]+'_ref.fa '+vcf_path+SV[0]+'_'+SV[1]+'.fa > '+vcf_path+SV[0]+'_'+SV[1]+'.sam')
            os.system('rm '+ vcf_path+SV[0]+'_'+SV[1]+'*.fa')
            with pysam.AlignmentFile(vcf_path+SV[0]+'_'+SV[1]+'.sam','r') as sam:
                for r in sam:
                    if (r.flag == 0) or (r.flag == 16):
                        maplen = 0
                        for c in r.cigar:
                            if c[0] == 0:
                                maplen = maplen+c[1]
                        if maplen/SV_len > 0.85 and abs(ref_start+r.reference_start-int(SV[1])) < 500:
                            SV[1] = str(ref_start+r.reference_start)
                            SV[4] = '<DUP>'
                            SV[7] = 'END=' +str(int(SV[1])+SV_len)+';SVTYPE=DUP;SVLEN='+str(SV_len)
            os.system('rm '+vcf_path+'*.sam')
        if abs(SV_len) < min_SV_len:
            continue
        
        output.append(SV[0]+'\t'+SV[1]+'\tSV'+str(ID)+'\t'+SV[3]+'\t'+SV[4]+'\t'+SV[5]+'\t'+SV[6]+'\t'+SV[7]+'\t'+SV[8]+'\t'+SV[9])
        ID += 1
    with open(out_path+"/gSV.vcf", "w+") as vcf_final:
        vcf_final.writelines(output)

    vcf_final.close()
    os.system('rm -r '+ vcf_path)




