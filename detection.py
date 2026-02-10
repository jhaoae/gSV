import pysam
import numpy as np
import os

def txt2list(file):
    lines = file.split('\n')
    result = []
    len_base = 0
    for line in lines[1:]:
        if line:
            columns = line.split('\t')
            result.append([columns[0].replace(" ",""), int(columns[1]), int(columns[2]), int(columns[3])])
            len_base += int(columns[3])
    return result, len_base


def check_overlap(start, end, SV_start, SV_end):
    start = start - 10
    end = end + 10
    if (start <= SV_start and SV_start <= end) or (start <= SV_end and SV_end <= end) or (SV_start <= start and end <= SV_end):
        return abs(start-SV_start)
    else:
        return -1


def detect_overlap(f):
    f.sort(key=lambda x: x[2])
    to_be_deleted = []
    for i in range(len(f) - 1):    
        if (f[i][5] > f[i + 1][2] and f[i + 1][5] > f[i][5]):
            #partial overlap
            if f[i][2] >= f[i+1][2]:
                to_be_deleted.append(i)
            elif f[i][2] < f[i+1][2] and f[i][3] >= f[i + 1][3]:
                #i longer than i+1
                overlap_len = f[i][5] - f[i + 1][2] + 1
                f[i + 1][1] = f[i + 1][1] + overlap_len
                f[i + 1][2] = f[i][5] + 1
                f[i + 1][3] = f[i + 1][3] - overlap_len
                f[i + 1][4] = f[i + 1][1] + f[i + 1][3] - 1
                f[i + 1][5] = f[i + 1][2] + f[i + 1][3] - 1
            elif f[i][2] < f[i+1][2] and f[i][3] < f[i + 1][3]:
                #i+1 longer than i
                f[i][3] = f[i][3] - (f[i][5] - f[i + 1][2] + 1)
                f[i][4] = f[i][1] + f[i][3] - 1
                f[i][5] = f[i][2] + f[i][3] - 1
            
        elif f[i][5] > f[i + 1][2] and f[i + 1][5] <= f[i][5]:
            #i Cover i+1
            to_be_deleted.append(i)
            f[i + 1] = f[i]
    return f, to_be_deleted

def overlap(f):
    flag = 0
    f.sort(key = lambda x: x[2])
    for i in range(len(f) -1):
        if f[i][5] > f[i+1][2]:
            flag = 1
    if flag == 1:
        return True
    else:
        return False


def combine(f):
    f.sort(key=lambda x: x[2])
    to_be_deleted = []
    diff = 50
    for i in range(len(f) - 1):
        D = f[i+1][1]-f[i][4]
        d = f[i+1][2]-f[i][5]
        if abs(D-d)<10 or (abs(D-d)<diff and abs(D) < diff and abs(d) < diff):
            f[i][3] = f[i][3] + f[i+1][3] + f[i+1][2]- f[i][5]
            f[i][4] = f[i+1][4]
            f[i][5] = f[i+1][5]
            to_be_deleted.append(i)
            f[i + 1] = f[i]
    return f, to_be_deleted


def re_align(f):
    while overlap(f):
        f, f_del = detect_overlap(f)
        f = [row for i, row in enumerate(f) if i not in f_del]
    return f   

def modify_exactmatch(f, distance):
    #distance: connect length
    f_len = len(f)
    f_new = []

    while f_len > 0:
        k = 0
        del_overlap = []
        while f_len > 1:
            for j in range(1, f_len):
                D = f[j][1] - f[0][4]
                d = f[j][2] - f[0][5]
                if (abs(D) < distance and d < distance) or abs(D - d) < distance:
                    k = j
                    break
            if k > 0:
                flag = 0
                if k > 1:
                    for j in range(1, k):
                        if f[j][2] - f[0][5]>-10 and f[k][2] - f[j][5] > -10:
                            flag = 1
                            break
                if flag != 1:
                    f[0][5] = f[k][5]
                    f[0][3] = f[0][5] - f[0][2] + 1
                    f[0][4] = f[0][1] + f[0][3] - 1
                    del f[1:k+1]
                    f_len = len(f)
                index = []
                for i, x in enumerate(f):
                    if i > 0 and x[2] < f[0][5] - 10:
                        index.append(i)
                if len(index) > 0:
                    f = [row for i, row in enumerate(f) if i not in index]
                f_len = len(f)
                k = 0
                if flag == 1:
                    break
            else:
                break
        f_new.append(f[0])
        del f[0]
        f_len = len(f)

    return f_new



def final_results(forward, len_forward, reverse, len_reverse, L, start, end, ref, ctg, genotype, ref_start): 
    
    forward_start = 1000000
    reverse_start = 1000000
    distance = 50

    if len_forward >= len_reverse:
        direction = "forward"
        f = forward
        r = reverse
        read_start = forward_start
    else:
        direction = "reverse"
        f = reverse
        r = forward
        read_start =  reverse_start

    for i in range(len(f)):
        f[i].append(f[i][1] + f[i][3] - 1)
        f[i].append(f[i][2] + f[i][3] - 1)
    for i in range(len(r)):
        r[i].append(r[i][1] + r[i][3] - 1)
        r[i].append(r[i][2] + r[i][3] - 1)
    r =  modify_exactmatch(r,distance)

    for i in range(len(r)):
        r[i][2] = L - (r[i][2] + r[i][3] - 1) + 1
        r[i][5] = r[i][2] + r[i][3] -1

    SV = []
    # Detect INVs
    for i in range(len(r)):
        SV_start = r[i][1]+ref_start
        SV_end = r[i][4]+ref_start
        if check_overlap(start, end, SV_start, SV_end)>0:
            seq = ref.fetch(r[0][0],SV_start,SV_start+1)
            SV.append(['INV',r[0][0], SV_start, SV_end, r[i][3], seq, '<INV>',genotype, r[i][5]+read_start])
    f.extend(r)
    f.sort(key=lambda x: x[2])

    f = modify_exactmatch(f, distance)
    f = re_align(f) 


    max_R = 0
    # Detect other SVs
    for i in range(len(f) - 1):
        max_R = max(max_R, f[i][4])
        D = f[i + 1][1] - f[i][4]
        d = f[i + 1][2] - f[i][5]
        
        D_L = f[i+1][1] - max_R
        D_R = f[i+1][4] - max_R
        if D_L >= 50 and D-d >= 50:
            SV_start_1 = f[i][4] + ref_start
            SV_end_1 = f[i][4]+(D-d+1) + ref_start
            SV_start_2 = f[i+1][1] - (D-d+1) +ref_start
            SV_end_2 = f[i+1][1] + ref_start
            overlap1 = check_overlap(start, end, SV_start_1, SV_end_1)
            overlap2 = check_overlap(start, end, SV_start_2, SV_end_2)
            if (overlap1>0 and overlap2<0) or (overlap1>0 and overlap2>overlap1):
                ref_seq = ref.fetch(f[0][0],SV_start_1,SV_end_1)
                ctg_seq = ref_seq[0]
                SV.append(['DEL', f[0][0], SV_start_1, SV_end_1, -(D-d + 1), ref_seq, ctg_seq, genotype, f[i][5]+read_start])
            elif (overlap2>0 and overlap1<0) or (overlap2>0 and overlap1>overlap2):
                ref_seq = ref.fetch(f[0][0],SV_start_2,SV_end_2)
                ctg_seq = ref_seq[0]#ctg.fetch(ctg.references[0],f[i+1][2]-1,f[i+1][2])
                SV.append(['DEL', f[0][0], SV_start_2, SV_end_2, -(D-d + 1), ref_seq, ctg_seq, genotype, f[i+1][2]-(D-d+1)+read_start])
        elif D_L >= -10 and d >=50 and d-D >=50:
            SV_start = f[i][4] + ref_start
            SV_end = SV_start + 1
            if check_overlap(start, end, SV_start, SV_end)>0:
                ref_seq = ref.fetch(f[0][0],SV_start,SV_start+1)
                ctg_seq = ctg.fetch(ctg.references[0],f[i][5],f[i+1][2])
                SV.append(['INS', f[0][0], SV_start, SV_end, d - D, ref_seq, ctg_seq, genotype, f[i][5]+read_start])
        elif D_L <= -50 and D_R >= -50:
            SV_start = f[i+1][1] + ref_start
            SV_end = max_R + ref_start
            if check_overlap(start, end, SV_start, SV_end)>0:
                seq = ref.fetch(f[0][0],SV_start,SV_start+1)
                SV.append(['DUP', f[0][0], SV_start, SV_end, SV_end-SV_start,seq,'tandemDUP', genotype, f[i][5]+read_start])
        elif D_R < -50:
            # disperedDUP
            SV_start = max_R
            SV_end = SV_start + 1
            disdup_start = f[i+1][1] + ref_start
            disdup_end = SV_start + f[i+1][3]
            if check_overlap(start, end, SV_start, SV_end)>0:
                seq = ref.fetch(f[0][0], SV_start,SV_start+1)
                SV.append(['INS', f[0][0], SV_start, SV_end, SV_end-SV_start,seq,'dispereDUP_'+str(disdup_start)+'_'+str(disdup_end), genotype, f[i+1][2]+read_start])


    
    SV = sorted(SV, key=lambda x: (x[1],x[2]))
    SV_final = []
    for i in range(len(SV)):
        last_index = len(SV_final)-1
        if i != 0 and SV[i][0] == SV_final[last_index][0] and abs(SV[i][2]-SV_final[last_index][2]) <= 250:
            if (abs(SV[i][4]/SV_final[last_index][4])<1.25) and (abs(SV[i][4]/SV_final[last_index][4])>0.8):
                if abs(SV[i][4])>abs(SV_final[last_index][4]):
                    SV_final[last_index] = SV[i]
                else:
                    continue
        SV_final.append(SV[i])
    return SV_final


def exactmatch(path, ref, chrom, mempath, memlen, min_SV_len):
    file_list = os.listdir(path)
    file_list.sort()

    vcf_file = open(path+chrom+"_complex_region.vcf", "w+")

    vcf_file.write("##fileformat=VCFv4.2\n")
    vcf_file.write("##source=bit_representation_model\n")
    vcf_file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant described in this record\">\n")
    vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    vcf_file.write("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
    vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    
    chromosome = ref.references

    for chro in chromosome:
        chro_len = ref.get_reference_length(chro)
        vcf_file.write("##contig=<ID="+chro+",length="+str(chro_len)+">\n")

    header_vec = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_ids_here']
    vcf_file.write("#" + "\t".join(header_vec) + "\n")
 
    
    ID = 1
    SV_list = []
    SV_complex_list = []
    genotype = '1/1'
    for c in file_list:
        if c.endswith('fa'):
            ctg = pysam.FastaFile(path+c)
            L=ctg.get_reference_length(ctg.references[0])            
            SV_name = c.split('/')[-1][:-3]
            chro = SV_name.split('-')[0]
            refL = ref.get_reference_length(chro)
            SV_start = SV_name.split('-')[1]
            SV_end = SV_name.split('-')[2]
            ref_start = max(int(SV_start) - max(3000,int(SV_end)-int(SV_start)),0)
            ref_end = min(int(SV_end) + max(3000,int(SV_end)-int(SV_start)),refL)

            ref_copmem = open(path+SV_start+'_'+str(ref_start)+'_ref.fa','a')
            ref_seq = ref.fetch(chro, ref_start, ref_end)
            ref_copmem.write('>'+chro+'\n'+ref_seq+'\n')
            ref_copmem.close()
            os.system(mempath+'/copmem2 -o '+ path + c.split('/')[-1][:-3] + '_bnew.txt -b -q -l '+ str(memlen)+' '+ path + SV_start+'_'+str(ref_start)+'_ref.fa ' + path + c)
            
            with open(path+c.split('/')[-1][:-3]+'_bnew.txt','r') as file:
                file_content = file.read()
                file_split = file_content.split('> ')[1:]
            start = int(SV_start)
            end = int(SV_end)
            if 'genotype' in c:
                genotype = '0/1'

            for i in range(0,len(file_split),2):
                forward, len_forward = txt2list(file_split[i])
                reverse, len_reverse = txt2list(file_split[i+1])
                if len(forward)==0 and len(reverse)==0:
                    continue
                SV = final_results(forward, len_forward, reverse, len_reverse, L, start, end, ref, ctg, genotype,ref_start)
                if SV:
                    SV_seen = []
                    SV_complex = []
                    SV = sorted(SV, key=lambda x: (x[1],x[2]))
                    for i in range(len(SV)):
                        if abs(SV[i][4])>min_SV_len and abs(SV[i][4])<20000:
                            SV_list.append(SV[i])   
                            if i!=0:
                                for j in range(len(SV_seen)):
                                    if SV_seen[j][1] == SV[i][1]:
                                        #overlap
                                        if SV[i][2]-SV_seen[j][3]<100 or SV[i][3]-SV_seen[j][2]<100 or (SV[i][2]>SV_seen[j][2] and SV[i][2]<SV_seen[j][3]) or (SV[i][3]>SV_seen[j][2] and SV[i][2]<SV_seen[j][2]):
                                            flag = 0 
                                            for k in range(len(SV_complex)):
                                                if  SV_complex[k][1] == SV[i][1]:
                                                    if SV[i][2]-SV_complex[k][2]<10 or SV[i][2]-SV_complex[k][3]<10 or (SV[i][2]>SV_complex[k][2] and SV[i][2]<SV_complex[k][3]) or (SV[i][3]>SV_complex[k][2] and SV[i][2]<SV_complex[k][2]):
                                                        if SV[i][2]<SV_complex[k][2]:
                                                            complex_start =  SV[i][2]
                                                            complex_ref = SV[i][5]
                                                        SV_complex[k][3] =  max(SV[i][3],SV_complex[k][3])
                                                        SV_complex[k][4] =  SV_complex[k][3]-SV_complex[k][2]
                                                        SV_complex[k][8].append([SV[i][0],SV[i][8],SV[i][2],SV[i][3],SV[i][4],SV[i][6]])
                                                        flag = 1
                                                        break
                                            if flag == 0:
                                                complex_end =  max(SV[i][3],SV_seen[j][3])
                                                if SV[i][2]<SV_seen[j][2]:
                                                    complex_start =  SV[i][2]
                                                    complex_ref = SV[i][5]
                                                else:
                                                    complex_start =  SV_seen[j][2]
                                                    complex_ref = SV_seen[j][5]
                                                sv_seen_type = SV_seen[j][0]
                                                sv_type = SV[i][0]
                                                
                                                SV_complex.append(['COMPLEX',SV[i][1], complex_start, complex_end, complex_end - complex_start,complex_ref,'<COM>', genotype,[[sv_seen_type,SV_seen[j][8],SV_seen[j][2],SV_seen[j][3],SV_seen[j][4],SV_seen[j][6]],[sv_type,SV[i][8],SV[i][2],SV[i][3],SV[i][4],SV[i][6]]]])
                                            else:
                                                break
                            SV_seen.append(SV[i])
                    SV_complex_list =  SV_complex_list + SV_complex
            ctg.close()
    SV_sort = sorted(SV_list, key=lambda x: (x[1],x[2]))
    SV_final = []
    for i in range(len(SV_sort)):
        if i != 0 and SV_sort[i][0] == SV_final[len(SV_final)-1][0] and SV_sort[i][1] == SV_final[len(SV_final)-1][1]:
            if abs(SV_sort[i][2]-SV_final[len(SV_final)-1][2]) <= 250:
                if (abs(SV_sort[i][4]/SV_final[len(SV_final)-1][4])<1.25) and (abs(SV_sort[i][4]/SV_final[len(SV_final)-1][4])>0.8):
                    if abs(SV_sort[i][4])>abs(SV_final[len(SV_final)-1][4]):
                        SV_final[len(SV_final)-1] = SV_sort[i]
                        continue
                    else:
                        continue
        
        SV_final.append(SV_sort[i])

    for i in range(len(SV_final)):
        if 'DUP' in SV_final[i][6]:
            body_vec = [SV_final[i][1], str(SV_final[i][2]), "SV" + str(ID),SV_final[i][5], '<'+SV_final[i][0]+'>',".", "PASS", "END=" +str(SV_final[i][3])+ ";SVTYPE="+SV_final[i][6]+";SVLEN="+str(SV_final[i][4]), "GT", SV_final[i][7]]
        else:
            body_vec = [SV_final[i][1], str(SV_final[i][2]), "SV" + str(ID),SV_final[i][5], SV_final[i][6],".", "PASS", "END=" +str(SV_final[i][3])+ ";SVTYPE="+SV_final[i][0]+";SVLEN="+str(SV_final[i][4]), "GT", SV_final[i][7]]
        vcf_file.write("\t".join(body_vec) + "\n")
        ID += 1
    vcf_file.close()
   
    


    vcf_com = open(path+chrom+"_complex_SV.vcf", "w+")

    vcf_com.write("##fileformat=VCFv4.2\n")
    vcf_com.write("##source=bit_representation_model\n")
    vcf_com.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant described in this record\">\n")
    vcf_com.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    vcf_com.write("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
    vcf_com.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

    chromosome = ref.references

    for chro in chromosome:
        chro_len = ref.get_reference_length(chro)
        vcf_com.write("##contig=<ID="+chro+",length="+str(chro_len)+">\n")

    header_vec = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_ids_here']
    vcf_com.write("#" + "\t".join(header_vec) + "\n")
    
    SV_complex_sort = sorted(SV_complex_list, key=lambda x: (x[1],x[2]))
    SV_complex_final = []
    for i in range(len(SV_complex_sort)):
        last_index = len(SV_complex_final)-1
        if i !=0:
            subtype1 = [item[0] for item in SV_complex_sort[i][8]]
            subtype2 = [item[0] for item in SV_complex_final[last_index][8]]
            if SV_complex_sort[i][1] == SV_complex_final[last_index][1] and set(subtype1) == set(subtype2):
                if abs(SV_complex_sort[i][2]-SV_complex_final[last_index][2]) <= 250:
                    if (abs(SV_complex_sort[i][4]/SV_complex_final[last_index][4])<1.25) and (abs(SV_complex_sort[i][4]/SV_complex_final[last_index][4])>0.8):
                        if abs(SV_complex_sort[i][4])>abs(SV_complex_final[last_index][4]):
                            SV_complex_final[last_index] = SV_complex_sort[i]
                            continue
                        else:
                            continue
        SV_complex_final.append(SV_complex_sort[i])
    
    ID = 1 
    for i in range(len(SV_complex_final)):
        detail = sorted(SV_complex_final[i][8], key=lambda x: x[1])
        subtype_list = []
        for item in detail:
            if 'dispereDUP' in item[5]:
                subtype = item[0]+':'+str(item[2])+'-'+str(item[3])+'-'+str(item[4])+'('+item[5]+')'
            else:
                subtype = item[0]+':'+str(item[2])+'-'+str(item[3])+'-'+str(item[4])
            subtype_list.append(subtype)

        subtype = '+'.join(subtype_list)
        body_vec = [SV_complex_final[i][1], str(SV_complex_final[i][2]), "SV" + str(ID),SV_complex_final[i][5], SV_complex_final[i][6],".", "PASS", "END=" +str(SV_complex_final[i][3])+ ";SVTYPE="+SV_complex_final[i][0]+":"+subtype+";SVLEN="+str(SV_complex_final[i][4]), "GT", SV_complex_final[i][7]]
        vcf_com.write("\t".join(body_vec) + "\n")
        ID += 1
    vcf_com.close()
    os.system('rm '+path+'*.fa*')
    

def split_fa(path):
    file_list = os.listdir(path)
    file_list.sort()
    for c in file_list:
        if c.endswith('fa'):
            sequences = {}
            current_sequence = None
            ctg_num = 0

            with open(path+c, 'r') as file:
                for line in file:
                    line = line.strip()

                    if line.startswith('>'):
                        ctg_num += 1
                        sequence_id = line[1:]
                        sequences[sequence_id] = ''
                        current_sequence = sequence_id
                    else:
                        sequences[current_sequence] += line
            if ctg_num == 1:
                continue

            file_num = 1
            
            file_name = c[:-3]
            if 'genotype' not in file_name:
                file_name = file_name + '-genotype'
            
            for sequence_id, sequence in sequences.items(): 
                current_file = path+file_name + '-'+str(file_num)+'.fa'
                file_num += 1
                with open(current_file, 'a') as file:
                    file.write(f'>{sequence_id}\n{sequence}\n')
            os.system('rm '+path + c)
            

