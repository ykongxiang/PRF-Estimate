from Bio import SeqIO
import os
import itertools
import re

file_path = '/Users/yangyuhao/Desktop/prfect/data/Covid19/gbk/covid19.gb'
file_path = '/Users/yangyuhao/Desktop/prfect/data/FSDB/gbk'
file_path = '/Users/yangyuhao/Desktop/prfect/data/RECODE/gbk'
file_path = '/Users/yangyuhao/Desktop/prfect/data/SEAPHAGES/gbk'
file_path = '/Users/yangyuhao/Desktop/prfect/data/Xu/gbk'
file_path = '/Users/yangyuhao/Desktop/prfect/data/Atkins/gbk'
genbank = SeqIO.read(file_path,'gb')
x =list()
for features in genbank.features:
    if features.type == 'CDS':
        print(features)
        if 'note' in features.qualifiers:
            note_info = features.qualifiers['note']
            x.append(note_info)
flattened_x = list(itertools.chain(*x))
split_x = [s for item in flattened_x for s in item.split(';') if s]
processed_list = [s.lstrip() for s in split_x]

import os
from Bio import SeqIO
import itertools
def find_matches_in_list(strings_list, target_substring = 'ribosomal frameshift'):
    return [string for string in strings_list if target_substring in string]

def extract_frameshift(folder_path):
# 初始化一个空列表来存储所有特征
    all_features = []
    for file_name in os.listdir(folder_path):
        if file_name.endswith('.gb'):
            file_path = os.path.join(folder_path, file_name)
            genbank = SeqIO.read(file_path, 'gb')
            file_features = []
            for features in genbank.features:
                seq_id = genbank.id
                seq = str(genbank.seq)
                strand = genbank.features[0].location.strand
                if features.type == 'CDS':
                    note_list = []
                    if 'note' in features.qualifiers:
                        note_info = features.qualifiers['note']
                        note_list.append(note_info)
            flattened_notes = list(itertools.chain(*note_list))
            split_notes = [s for item in flattened_notes for s in item.split(';') if s]
            processed_notes = [s.lstrip() for s in split_notes if s]

            file_features.append({
                    'id': seq_id,
                    'Sequence' : seq,
                    'Strand' : strand,
                    'note' : processed_notes
                    })
            # 将当前文件的特征列表添加到总列表中
            all_features.extend(file_features)
    return all_features
os.chdir("C:/Users/31598/Desktop/prfect")
Atkins = extract_frameshift('data/Atkins/gbk')
