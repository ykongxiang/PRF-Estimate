import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.preprocessing import OneHotEncoder
file_path ='G:/移码突变/移码数据/all_prf_data.xlsx'
data = pd.read_excel(file_path)

# 假设碱基序列在列 'sequence' 中，目标变量列名为 'PRF'
sequences = data['Sequence'].values
y = data['PRF'].values

# 过滤掉包含未知碱基 'N' 的序列
def is_valid_sequence(seq):
    return all(base in 'ATCG' for base in seq)

filtered_sequences = [seq for seq in sequences if is_valid_sequence(seq)]
filtered_y = [y[i] for i in range(len(sequences)) if is_valid_sequence(sequences[i])]
# 假设碱基序列在列 'Sequence' 中，目标变量列名为 'PRF'
sequences = filtered_sequences
y = filtered_y


# One-hot Encoding
def one_hot_encode(sequence, max_length):
    encoding = {'A': [1, 0, 0, 0], 'T': [0, 1, 0, 0], 'C': [0, 0, 1, 0], 'G': [0, 0, 0, 1]}
    one_hot_sequences = np.zeros((len(sequence), max_length, 4), dtype=int)

    for i, seq in enumerate(sequence):
        for j, base in enumerate(seq):
            if j < max_length:
                one_hot_sequences[i, j] = encoding.get(base, [0, 0, 0, 0])

    return one_hot_sequences


max_length = max(len(seq) for seq in sequences)
X = one_hot_encode(sequences, max_length)
X = X.reshape(len(sequences), -1)  # Flatten for model input

# 划分数据集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y)

# 训练模型
model = RandomForestClassifier()
model.fit(X_train, y_train)

# 评估模型
y_pred = model.predict(X_test)
print(classification_report(y_test, y_pred))
