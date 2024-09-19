import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler

# 读取数据
file_path = "C:/Users/31598/Documents/combined_df.xlsx"
data = pd.read_excel(file_path)

# 选择特征列和目标变量列
X = data[['pident', 'evalue','bitscore','new_distance','Alignment_Score']]
y = data['PRF']  # 二分类标签

# 数据标准化
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 创建逻辑回归模型
#model = LogisticRegression()
model = RandomForestClassifier()
# 创建分层交叉验证对象
stratified_kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# 进行分层交叉验证
cross_val_scores = cross_val_score(model, X_scaled, y, cv=stratified_kfold, scoring='accuracy')

# 拟合模型
model.fit(X_scaled, y)

# 打印结果
print(f"Model Coefficients: {model.coef_}")
print(f"Stratified Cross-validation scores (Accuracy): {cross_val_scores}")
print(f"Mean Stratified Cross-validation score: {np.mean(cross_val_scores)}")
