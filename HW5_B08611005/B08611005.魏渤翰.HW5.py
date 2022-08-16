from itertools import groupby
from audioop import reverse
from os import sep
import pandas as pd
from Bio.Seq import Seq
# def Reverse_Complement(base):
#   base = base.replace('A', 'T')
#   base = base.replace('A', 'T')
#   base = base.replace('A', 'T')
#   base = base.replace('A', 'T')


data_1720 = []
data_3079 = []
data = []
dic = {}
position1 = 0
position2 = 0
with open('/Users/weber/Desktop/2022summer_course/Homeowrk/HW5_B08611005/HW5_data/ITGA2B.vcf') as vcf:
    position1 = int(vcf.readline().split()[1])
    position2 = int(vcf.readline().split()[1])
    # print(position1,position2)
vcf.close()
with open('./HW5_data/hw5.sam')as f:
    for line in f.readlines()[94:len(f.readlines())-1]:
        temp = line.split()
        position = int(temp[3])
        index_num = -2
        count = 0
        seq_pos_first = position1 - position
        seq_pos_second = position2 - position
        a = temp[5]
        CIGAR = [''.join(list(g))
                 for k, g in groupby(a, key=lambda x: x.isdigit())]
        if temp[1] == '0' or temp[1] == '2048':
            seq = list(temp[9])
            while count < seq_pos_second + 1:
                index_num += 2
                if CIGAR[index_num+1] == 'I':
                    seq[count] = ''.join(
                        seq[count:count+int(CIGAR[index_num])+1])
                    for r in range(int(CIGAR[index_num])):
                        seq.pop(count+1)
                elif CIGAR[index_num+1] == 'D':
                    for i in range(int(CIGAR[index_num])):
                        seq.insert(count, '.')
                    count += int(CIGAR[index_num])
                else:
                    count += int(CIGAR[index_num])
            # data_1720.append(seq[seq_pos_first])
            # data_3079.append(seq[seq_pos_second])
            if len(seq) > seq_pos_first and len(seq) > seq_pos_second:
                key = seq[seq_pos_first] + "/"+seq[seq_pos_second]
            # print(key)

                if key not in dic:
                    dic[key] = 1
                else:
                    dic[key] += 1
        else:
            CIGAR.reverse()
            # print(CIGAR[0:20], "\n")
            # print(temp[9])
            seq = list(str(Seq(temp[9]).reverse_complement()))
            # print(seq)
            while count < seq_pos_second + 1:
                index_num += 2
                if CIGAR[index_num] == 'I':
                    seq[count] = ''.join(
                        seq[count:count+int(CIGAR[index_num+1])+1])
                    for r in range(int(CIGAR[index_num+1])):
                        seq.pop(count+1)
                elif CIGAR[index_num] == 'D':
                    for i in range(int(CIGAR[index_num+1])):
                        seq.insert(count, '.')
                    count += int(CIGAR[index_num+1])
                else:
                    count += int(CIGAR[index_num+1])
            # print(len(seq),seq_pos_first,seq_pos_second)
            if len(seq) > seq_pos_first and len(seq) > seq_pos_second:
                key = seq[seq_pos_first] + "/"+seq[seq_pos_second]
            # print(key)
                if key not in dic:
                    dic[key] = 1
                else:
                    dic[key] += 1

# print(data_1720,data_3079)
dic
f.close()


sorted_dic = {k: v for k, v in sorted(
    dic.items(), key=lambda item: item[1], reverse=True)}
output = pd.DataFrame(list(sorted_dic.items()), columns=['Haplotype', 'Count'])
output.to_csv("B08611005.魏渤翰.HW5.txt", sep="\t", index=False)
