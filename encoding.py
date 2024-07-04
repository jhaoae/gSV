import numpy as np
import time
from datetime import datetime
from build_matrix import read_matrix, read_matrix_reverse
import pysam
import re
import gc
import scipy.io as scio





def read_encoding(bam_file, rf_array, out_path, chrom, min_mapq, min_len, sub_num):

    chro_len = len(rf_array[0])
    prim, prim_s, supp = bam_iterator(bam_file,min_mapq,chrom,chro_len)
    name = list()
    
    pre_type_info_list = list()
    pre_type_name_list = list()
    s_list = list()
    s_name_list = list()


    ref_len = np.size(rf_array,1)
    loss = np.zeros((1, ref_len)).astype(np.int32)
    depth = np.zeros((1, ref_len)).astype(np.int32)
    support_sum = np.zeros((1, ref_len)).astype(np.int32)
    

    i=0
    currentDateAndTime = datetime.now()
    print(chrom + " The current date and time is", currentDateAndTime)

    
    with open(out_path+'/'+chrom+"_preliminary.txt", "w") as f:

        #primary
        for r in prim:
            i = i + 1
            if i % sub_num == 0:
                print(chrom, ':', i)
            j=0
            k=0
            bf_array = read_matrix(r.cigar, r.query_sequence, chrom, r.reference_start, r.reference_end, r.query_name, min_len, f)
            #name.append(r.query_name)

            support = np.zeros((1,r.reference_end-r.reference_start)).astype(np.int32)
            depth[:,r.reference_start:r.reference_end] += 1
            #print(r.query_name,r.reference_start,r.reference_end)
            difference = abs(bf_array - rf_array[:,r.reference_start:r.reference_end])
            diff_s = sum(difference)
            r_loss = np.where(diff_s>2, diff_s, np.zeros_like(diff_s)) + np.where((diff_s<=2) & (diff_s!=0), np.ones_like(diff_s), np.zeros_like(diff_s))
            support[:,np.where(r_loss!=0)] = 1
            support_sum[:,r.reference_start:r.reference_end] += support
            loss[:,r.reference_start:r.reference_end] += r_loss
    
    
        #supplementary
        for r in prim_s:
            i = i + 1
            #r_name = 'r' + str(len(name)) + '_' + r.query_name.replace('/','_')
            if i%sub_num == 0:
                print(chrom, ':', i)
        
            supp, supp_r, supp_pos, rstart, direct, flag = correspond_supp(supp,r,min_len,f)

            p_start = r.reference_start
            p_end = r.reference_end
            #print(r.query_name, p_start)
            if ((direct == 0) | (flag ==0)):
                j=0
                bf_array = read_matrix(r.cigar, r.query_sequence, chrom, r.reference_start, r.reference_end, r.query_name, min_len, f)
                p_array = bf_array
                idx = 0
                for s in supp_r:
                    s_start = s.reference_start
                    s_end = s.reference_end
                    if s_start > p_start:
                        p1 = p_start
                        p2 = p_end
                        p3 = s_start
                        p4 = s_end
                    else:
                        p1 = s_start
                        p2 = s_end
                        p3 = p_start
                        p4 = p_end
                    if (p3 - p2) > 50000:
                        idx += 1
                        print('larger than 50000')
                        continue
                    combine_array = np.zeros((4, (max(p2, p4) - p1))).astype(np.int32)
                    inv_array = np.zeros((4, (max(p2, p4) - p1))).astype(np.int32)
                    
                    if s.is_reverse == r.is_reverse:
                        sup_array = read_matrix(s.cigar, r.query_sequence, chrom, s_start, s_end, r.query_name, min_len, f)
                        #DEL
                        if (p3-p2 > min_len) and (direct == 0):
                            if (p_start>s_end and rstart>supp_pos[idx]) or (s_start>p_end and supp_pos[idx] > rstart):
                                #print('1',r.query_name,r.reference_start,str(p2),str(p3-p2))
                                line = str(chrom)+'\t'+str(p2)+'\t'+r.query_name+'\t'+str(p3-p2)+'\tDEL\t.\n'
                                f.write(line)
                    else:
                        #INV
                        inv_len = 0
                        sup_array = read_matrix_reverse(s.cigar, r.query_sequence, chrom, s_start, s_end, r.query_name, 0, min_len, f)
                        #print('s is inv')
                        #print(s_start,r.reference_end,r.reference_start,s_end)
                        if s_end - r.reference_end > min_len:
                            if s_start-r.reference_end > 5:
                                inv_array[:,(r.reference_end-p1):(s_start-p1)] += np.fliplr(rf_array[:,r.reference_end:s_start])
                            ss_end = s_start
                            for c in s.cigar:
                                if c[0] == 0 or c[0] == 7 or c[0] == 8:
                                    ss_end += c[1]
                                elif c[0] == 2:
                                    if c[1] > 20:
                                        break
                                    else:
                                        ss_end += c[1]
                            inv_start = r.reference_end
                            inv_len = ss_end - inv_start
                        elif r.reference_start - s_start > min_len:
                            if r.reference_start - s_end > 5:
                                inv_array[:,(s_end-p1):(r.reference_start-p1)] += np.fliplr(rf_array[:,s_end:r.reference_start])
                            ss_start = s_end
                            for c in reversed(s.cigar):
                                if c[0] == 0 or c[0] == 7 or c[0] == 8:
                                    ss_start = ss_start - c[1]
                                elif c[0] == 2:
                                    if c[1] > 20:
                                        break
                                    else:
                                        ss_start =ss_start - c[1]
                            inv_start = ss_start
                            inv_len = r.reference_start - inv_start
                        #print('inv_len',inv_len)
                        if inv_len > min_len:
                            #print('INV1',inv_start,inv_len)
                            line = str(chrom)+'\t'+str(inv_start)+'\t'+r.query_name+'\t'+str(inv_len)+'\tINV\t.\n'
                            f.write(line)

                    if s_start > p_start:
                        combine_array[:, (p1 - p1):(p2 - p1)] = p_array
                        combine_array[:, (p3 - p1):(p4 - p1)] = combine_array[:, (p3 - p1):(p4 - p1)] + sup_array 
                        combine_array += inv_array
                    else:
                        combine_array[:, (p1 - p1):(p2 - p1)] = sup_array
                        combine_array[:, (p3 - p1):(p4 - p1)] = combine_array[:, (p3 - p1):(p4 - p1)] + p_array
                        combine_array += inv_array
                    
                    #DUP
                    if (supp_pos[idx] > rstart) & (r.reference_end > s_start):
                        dup_array = np.zeros((4, (r.reference_end-s_start))).astype(np.int32)
                        if r.reference_end > s_end:
                            dup_array[:, (s_end-s_start):(r.reference_end - s_start)] += rf_array[:,s_end:r.reference_end]
                        if r.reference_start > s_start:
                            dup_array[:,0:(r.reference_start-s_start)] += rf_array[:,s_start:r.reference_start]
                        combine_array[:, (s_start - p1):(r.reference_end - p1)] = combine_array[:, (s_start - p1):(r.reference_end - p1)] + dup_array
                        dup_len = r.reference_end - s.reference_start
                        if dup_len > min_len:
                            line = str(chrom)+'\t'+str(s.reference_start)+'\t'+r.query_name+'\t'+str(dup_len)+'\tDUP\t.\n'
                            f.write(line)
                    
                    elif (supp_pos[idx] < rstart) & (s_end > r.reference_start):
                        dup_array = np.zeros((4, (s_end-r.reference_start))).astype(np.int32)
                        if s_start > r.reference_start:
                            dup_array[:, 0:(s_start - r.reference_start)] += rf_array[:,r.reference_start:s_start]
                        if r.reference_end < s_end:
                            dup_array[:, (r.reference_end - r.reference_start): (s_end - r.reference_start)] += rf_array[:,r.reference_end:s_end]
                        combine_array[:, (r.reference_start - p1):(s_end - p1)] = combine_array[:, (r.reference_start - p1):(s_end - p1)] + dup_array

                        dup_len = s.reference_end - r.reference_start
                        if dup_len > min_len:
                            line = str(chrom)+'\t'+str(r.reference_start)+'\t'+r.query_name+'\t'+str(dup_len)+'\tDUP\t.\n'
                            f.write(line)
                    p_array = combine_array
                    p_start = p1
                    p_end = max(p2, p4)
                    idx += 1

        
            elif flag != 0:
                j=0
                bf_array = read_matrix_reverse(r.cigar, r.query_sequence, chrom, r.reference_start, r.reference_end, r.query_name, 1, min_len, f)
                p_array = bf_array
                idx = 0
                for s in supp_r:
                    #print(idx,len(supp_pos))
                    s_start = s.reference_start
                    s_end = s.reference_end
                    if s_start > p_start:
                        p1 = p_start
                        p2 = p_end
                        p3 = s_start
                        p4 = s_end
                    else:
                        p1 = s_start
                        p2 = s_end
                        p3 = p_start
                        p4 = p_end
                    if (p3 - p2) > 50000:
                        idx += 1
                        continue
                    combine_array = np.zeros((4, (max(p2, p4) - p1))).astype(np.int32)
                    inv_array = np.zeros((4, (max(p2, p4) - p1))).astype(np.int32)
                    if s.is_reverse == r.is_reverse:
                        sup_array = read_matrix_reverse(s.cigar, r.query_sequence, chrom, s_start, s_end, r.query_name, 1, min_len, f)
                        #DEL
                        if p3-p2 > min_len:
                            if (p_start>s_end and rstart>supp_pos[idx]) or (s_start>p_end and supp_pos[idx]>rstart):
                            #print('2',r.query_name,str(p2),str(p3-p2))
                                line = str(chrom)+'\t'+str(p2)+'\t'+r.query_name+'\t'+str(p3-p2)+'\tDEL\t.\n'
                                f.write(line)
                    else:
                        sup_array = read_matrix(s.cigar, s.query_sequence, chrom, s_start, s_end, r.query_name, min_len, f)
                        #INV r is inv
                        inv_len = 0
                        #print('r is inv')
                        if r.reference_end - s_end > min_len:
                            if r.reference_start - s_end > 5:
                                inv_array[:,(s_end-p1):(r.reference_start-p1)] += np.fliplr(rf_array[:,s_end:r.reference_start])
                            rr_end = r.reference_start
                            for c in r.cigar:
                                if c[0] == 0 or c[0] == 7 or c[0] == 8:
                                    rr_end += c[1]
                                elif c[0] == 2:
                                    if c[1] > 20:
                                        break
                                    else:
                                        rr_end += c[1]
                            inv_start = s_end
                            inv_len = rr_end-inv_start

                        elif s_start - r.reference_start > min_len:
                            if s_start - r.reference_end > 5:
                                inv_array[:,(r.reference_end-p1):(s_start-p1)] += np.fliplr(rf_array[:,r.reference_end:s_start])
                            rr_start = r.reference_end
                            for c in reversed(r.cigar):
                                if c[0] == 0 or c[0] == 7 or c[0] == 8:
                                    rr_start = rr_start - c[1]
                                elif c[0] == 2:
                                    if c[1] > 20:
                                        break
                                    else:
                                        rr_start =rr_start - c[1]

                            inv_start = rr_start
                            inv_len = s_start - inv_start
                        if inv_len > min_len:
                            #print('INV2',inv_start,inv_len)
                            line = str(chrom)+'\t'+str(r.reference_start)+'\t'+r.query_name+'\t'+str(inv_len)+'\tINV\t.\n'
                            f.write(line)
                    
                    if s_start > p_start:
                        combine_array[:, (p1 - p1):(p2 - p1)] = p_array
                        combine_array[:, (p3 - p1):(p4 - p1)] = combine_array[:, (p3 - p1):(p4 - p1)] + sup_array 
                        combine_array += inv_array
                    else:
                        combine_array[:, (p1 - p1):(p2 - p1)] = sup_array
                        combine_array[:, (p3 - p1):(p4 - p1)] = combine_array[:, (p3 - p1):(p4 - p1)] + p_array 
                        combine_array += inv_array

                    #DUP
                    if (supp_pos[idx] > rstart) & (r.reference_end > s_start):  
                        dup_array = np.zeros((4, (r.reference_end-s_start))).astype(np.int32)
                        if r.reference_end > s_end:
                            dup_array[:, (s_end-s_start):(r.reference_end - s_start)] += rf_array[:,s_end:r.reference_end]
                        if r.reference_start > s_start:
                            dup_array[:,0:(r.reference_start-s_start)] += rf_array[:,s_start:r.reference_start]
                        combine_array[:, (s_start - p1):(r.reference_end - p1)] = combine_array[:, (s_start - p1):(r.reference_end - p1)] + dup_array
                    
                        dup_len = r.reference_end - s.reference_start
                        if dup_len > min_len:
                            line = str(chrom)+'\t'+str(s.reference_start)+'\t'+r.query_name+'\t'+str(dup_len)+'\tDUP\t.\n'
                            f.write(line)
                    
                    elif (supp_pos[idx] < rstart) & (s_end > r.reference_start):
                        dup_array = np.zeros((4, (s_end-r.reference_start))).astype(np.int32)
                        if s_start > r.reference_start:
                            dup_array[:, 0:(s_start - r.reference_start)] += rf_array[:,r.reference_start:s_start]
                        if r.reference_end < s_end:
                            dup_array[:, (r.reference_end - r.reference_start): (s_end - r.reference_start)] += rf_array[:,r.reference_end:s_end]
                        combine_array[:, (r.reference_start - p1):(s_end - p1)] = combine_array[:, (r.reference_start - p1):(s_end - p1)] + dup_array
                    
                        dup_len = s.reference_end - r.reference_start
                        if dup_len > min_len:
                            line = str(chrom)+'\t'+str(r.reference_start)+'\t'+r.query_name+'\t'+str(dup_len)+'\tDUP\t.\n'
                            f.write(line)
                    
                    p_array = combine_array
                    p_start = p1
                    p_end = max(p2, p4)
                    idx += 1 

            #name.append(r.query_name)
        
            support = np.zeros((1,p_end-p_start)).astype(np.int32)
            depth[:,p_start:p_end] += 1
            difference = abs(p_array - rf_array[:,p_start:p_end])
            diff_s = sum(difference)
            r_loss = np.where(diff_s>2, diff_s, np.zeros_like(diff_s)) + np.where((diff_s<=2) & (diff_s!=0), np.ones_like(diff_s), np.zeros_like(diff_s))
            support[:,np.where(r_loss!=0)] = 1
            support_sum[:,p_start:p_end] += support
            loss[:,p_start:p_end] += r_loss
    
    sub = str(int(i/sub_num+1))

    depth_re = np.where(depth == 0, 1, depth)
    depth_num = (depth !=0)
    depth_mean = depth.sum()/depth_num.sum()
    loss_re = loss/depth_re*depth_mean
    
    return depth_mean, loss_re, support_sum



def bam_iterator(bam_file,min_mapq,chrom,chro_len):
    #alignments = bam_file.fetch(chrom,78240000,78260000)
    alignments = bam_file.fetch(chrom)
    prim = []
    prim_s = []
    supp = []
    for r in alignments:
        if r.mapping_quality <= min_mapq:
            continue
        if ((r.flag == 0) | (r.flag == 16)):
            if r.has_tag('SA'):
                prim_s.append(r)
            else:
                prim.append(r)
        elif ((r.flag == 2048) | (r.flag == 2064)):
            supp.append(r)
    print('Number of primary:', len(prim))
    print('Number of reads with supp:', len(prim_s))
    print('Number of supplementary:', len(supp))
    return prim, prim_s, supp

def correspond_supp(supp,r,min_len,f):
    supp_index = []
    del_inx = []
    supp_r = []
    supp_pos = []
    direct = 0
    flag = 0
    rstart = r.query_alignment_start
    
    for su in supp:
        if su.query_name == r.query_name:
            #print('111')
            if r.reference_name != su.reference_name:
                #TRA!!!!!!!
                sv_len = abs(su.reference_end-su.reference_start)
                if sv_len > min_len:
                    line = su.reference_name+'\t'+str(su.reference_start)+'\t'+r.query_name+'\t'+str(sv_len)+'\tTRA\t.\n'
                    f.write(line)
                del_inx.append(supp.index(su))
                continue

            if ((r.reference_start-su. reference_end) > 50000) or ((su.reference_start-r.reference_end)>50000):
                del_inx.append(supp.index(su))
                continue
            supp_index.append(supp.index(su))
            supp_r.append(su)


            if su.is_reverse != r.is_reverse:
                direct += 1
                if r.is_reverse:
                    #r is reverse
                    rstart = r.infer_read_length()-r.query_alignment_end
                    sstart = su.query_alignment_start
                    if rstart - sstart > 0:
                        if r.reference_start - su.reference_start > 0:
                            flag += 1
                            #r is inv
                    else:
                        if r.reference_start - su.reference_start < 0:
                            flag += 1
                            #r is inv
                else:
                    #r is forward
                    rstart = r.query_alignment_start
                    sstart = su.infer_read_length()-su.query_alignment_end
                    if rstart - sstart > 0:
                        if r.reference_start - su.reference_start < 0:
                            flag += 1
                            #r is inv
                    else:
                        if r.reference_start - su.reference_start > 0:
                            flag += 1
                            #r is inv
                if flag == 0:
                    rstart = r.query_alignment_start
                    supp_pos.append(su.infer_read_length()-su.query_alignment_end)
                else:
                    rstart = r.infer_read_length()-r.query_alignment_end
                    supp_pos.append(su.query_alignment_start)
            else:
                supp_pos.append(su.query_alignment_start)
    sum_index = supp_index + del_inx
    supp = np.delete(supp, sum_index).tolist()
    #print(direct,flag)
    return supp, supp_r, supp_pos, rstart, direct, flag 



def savemat(read_name,read_list,read_info_list,pre_type_name_list,pre_type_info_list,s_name_list,s_list,path,file_name,chrom,sub,file_read,file_readinfo,file_type,file_s):
     #read_list = np.array(read_list,dtype=np.object)
     scio.savemat(path+'results/'+file_name+'/'+file_name+'_'+chrom+'_'+sub+file_read, dict(zip(read_name, read_list)))
     scio.savemat(path+'results/'+file_name+'/'+file_name+'_'+chrom+'_'+sub+file_readinfo, dict(zip(read_name, read_info_list)))
     pre_type_info_list = np.array(pre_type_info_list,dtype=np.object)
     scio.savemat(path+'results/'+file_name+'/'+file_name+'_'+chrom+'_'+sub+file_type, dict(zip(pre_type_name_list,pre_type_info_list)))
     scio.savemat(path+'results/'+file_name+'/'+file_name+'_'+chrom+'_'+sub+file_s, dict(zip(s_name_list,s_list)))

     currentDateAndTime = datetime.now()
     print(chrom +'-'+ sub+" The current date and time is", currentDateAndTime)

