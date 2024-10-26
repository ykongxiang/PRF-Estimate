import pandas as pd
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from tensorflow.keras.preprocessing.sequence import pad_sequences

# 读取Excel文件
file_path = "C:/Users/31598/Desktop/移码数据/all_blastx.xlsx"
data = pd.read_excel(file_path)

# 提取DNA序列和数值特征
sequences = data['Sequence'].values
frameshift = data['frameshift'].values
numerical_features = data[['pident', 'evalue', 'bitscore', 'gapopen']].values

# 将DNA序列转换为整数编码
def encode_sequence(sequence):
    mapping = {'a': 1, 't': 2, 'c': 3, 'g': 4,'n':5}
    return [mapping[base] for base in sequence]

encoded_sequences = [encode_sequence(seq) for seq in sequences]

# 填充序列
max_length = max(len(seq) for seq in encoded_sequences)
padded_sequences = pad_sequences(encoded_sequences, maxlen=max_length, padding='post')

# 标准化数值特征
scaler = StandardScaler()
numerical_features = scaler.fit_transform(numerical_features)

# 划分训练集和测试集
X_train_seq, X_test_seq, X_train_num, X_test_num, y_train, y_test = train_test_split(
    padded_sequences, numerical_features, frameshift, test_size=0.2, random_state=42)

from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense, Embedding, MultiHeadAttention, LayerNormalization, Dropout, GlobalAveragePooling1D, Concatenate

# Transformer编码器层
class TransformerEncoder(tf.keras.layers.Layer):
    def __init__(self, embed_dim, num_heads, ff_dim, rate=0.1):
        super(TransformerEncoder, self).__init__()
        self.att = MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim)
        self.ffn = tf.keras.Sequential(
            [Dense(ff_dim, activation="relu"), Dense(embed_dim),]
        )
        self.layernorm1 = LayerNormalization(epsilon=1e-6)
        self.layernorm2 = LayerNormalization(epsilon=1e-6)
        self.dropout1 = Dropout(rate)
        self.dropout2 = Dropout(rate)

    def call(self, inputs, training):
        attn_output = self.att(inputs, inputs)
        attn_output = self.dropout1(attn_output, training=training)
        out1 = self.layernorm1(inputs + attn_output)
        ffn_output = self.ffn(out1)
        ffn_output = self.dropout2(ffn_output, training=training)
        return self.layernorm2(out1 + ffn_output)

# 输入层
sequence_input = Input(shape=(max_length,), name='sequence_input')
numerical_input = Input(shape=(X_train_num.shape[1],), name='numerical_input')

# 序列特征处理
embedding_dim = 64
num_heads = 4
ff_dim = 128

x = Embedding(input_dim=5, output_dim=embedding_dim, input_length=max_length)(sequence_input)
x = TransformerEncoder(embed_dim=embedding_dim, num_heads=num_heads, ff_dim=ff_dim)(x)
x = GlobalAveragePooling1D()(x)

# 数值特征处理
y = Dense(64, activation='relu')(numerical_input)
y = Dropout(0.5)(y)
y = Dense(32, activation='relu')(y)

# 合并两个输入的输出
combined = Concatenate()([x, y])
z = Dense(64, activation='relu')(combined)
z = Dropout(0.5)(z)
z = Dense(1, activation='sigmoid')(z)

# 构建模型
model = Model(inputs=[sequence_input, numerical_input], outputs=z)

# 编译模型
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

# 模型结构
model.summary()

# 训练模型
history = model.fit(
    [X_train_seq, X_train_num], y_train,
    validation_data=([X_test_seq, X_test_num], y_test),
    epochs=20,
    batch_size=32
)

# 评估模型
loss, accuracy = model.evaluate([X_test_seq, X_test_num], y_test)
print(f"Test Loss: {loss}")
print(f"Test Accuracy: {accuracy}")
from sklearn.metrics import precision_score, recall_score, f1_score

# 预测
y_pred = (model.predict([X_test_seq, X_test_num]) > 0.5).astype("int32")

# 计算精确率、召回率和F1分数
precision = precision_score(y_test, y_pred)
recall = recall_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred)

print(f"Precision: {precision}")
print(f"Recall: {recall}")
print(f"F1 Score: {f1}")