import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler

# 假设CSV文件名为 'data.csv'，特征列名为 'feature1', 'feature2' 等，目标变量列名为 'target'
file_path = "C:/Users/31598/Documents/combined_df.xlsx"
data = pd.read_excel(file_path)
# 选择特征列和目标变量列
#X = data[['pident', 'evalue','bitscore','distance_score','Alignment_Score']]  # 替换为实际的特征列名
X = data[['pident', 'evalue','bitscore','new_distance','Alignment_Score']]
y = data['PRF']
# 数据标准化
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 创建线性回归模型
model = LinearRegression()

# 拟合模型
model.fit(X_scaled, y)

# 模型系数
coefficients = model.coef_

# 模型评估
r_squared = model.score(X_scaled, y)
cross_val_scores = cross_val_score(model, X_scaled, y, cv=5)

import pandas as pd
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import StratifiedKFold, GridSearchCV, cross_val_score
from sklearn.preprocessing import StandardScaler

# 读取数据
file_path = "C:/Users/31598/Documents/combined_df.xlsx"
data = pd.read_excel(file_path)

# 选择特征列和目标变量列
X = data[['pident', 'evalue','bitscore','new_distance','Alignment_Score']]
y = data['PRF']  # 假设是分类变量

# 数据标准化
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 创建Ridge回归模型
ridge = Ridge()

# 定义参数网格
param_grid = {'alpha': [0.1, 1, 10, 100, 1000]}

# 创建分层交叉验证对象
stratified_kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# 使用网格搜索进行交叉验证，寻找最佳alpha值
grid_search = GridSearchCV(ridge, param_grid, cv=stratified_kfold)
grid_search.fit(X_scaled, y)

# 获取最佳模型
best_ridge = grid_search.best_estimator_

# 模型评估
r_squared = best_ridge.score(X_scaled, y)
cross_val_scores = cross_val_score(best_ridge, X_scaled, y, cv=stratified_kfold)

# 打印结果
print(f"Best Alpha: {grid_search.best_params_['alpha']}")
print(f"Best Ridge Model Coefficients: {best_ridge.coef_}")
print(f"R-squared: {r_squared}")
print(f"Stratified Cross-validation scores: {cross_val_scores}")
print(f"Mean Stratified Cross-validation score: {np.mean(cross_val_scores)}")


# 打印结果
print(f"Coefficients: {coefficients}")
print(f"R-squared: {r_squared}")
print(f"Cross-validation scores: {cross_val_scores}")

# 模型诊断
# 这里可以添加代码来检查残差分布、异方差性等